"""Exercise the saved-pair handoff, branch-specific outgroups and local fallback."""

import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

from awgds.wgd import matched_resampling
from awgds.wgd_linkage import linkage_test, phase_agreement, rng_for
from fixture import CLADES, LABELLED, SPECIES, write


def data(root):
    names = [f"chr{i}" for i in range(1, 7)]
    sequences = {s: {c: [] for c in names} for s in SPECIES}
    annotations = {s: [] for s in SPECIES}
    membership, counts, events = [], [], []
    background = [node for node in CLADES if node != "N2"]
    for i in range(54 + 4 * len(background)):
        family = f"OG{i:07d}"
        focal = i < 54
        node = "N2" if focal else background[(i - 54)//4]
        source = i % 3
        per_species = {s: [] for s in SPECIES}
        a_ids, b_ids = [], []
        for s in SPECIES:
            ingroup = s in CLADES[node]
            ca, cb = names[source], names[source + 3]
            if focal and s == "Taxon_B":
                ca, cb = cb, ca
            for copy, chromosome in (("a", ca), ("b", cb)):
                if copy == "b" and not ingroup:
                    continue
                identifier = f"{family}_{copy}"
                seq = sequences[s][chromosome]
                lo = len(seq) * 75 + 1
                hi = lo + 65
                seq.append("ATG" + "GCT" * 20 + "TAA" + "A" * 9)
                annotations[s].extend((
                    f"{chromosome}\ttest\tgene\t{lo}\t{hi}\t.\t+\t.\tID={identifier}",
                    f"{chromosome}\ttest\tmRNA\t{lo}\t{hi}\t.\t+\t.\tID={identifier}.1;Parent={identifier}",
                    f"{chromosome}\ttest\tCDS\t{lo}\t{hi}\t.\t+\t0\tID=cds-{identifier};Parent={identifier}.1",
                ))
                per_species[s].append(f"{identifier}.1")
                if ingroup:
                    (a_ids if copy == "a" else b_ids).append(f"{s}_{identifier}.1")
        membership.append("\t".join([family, *(", ".join(per_species[s]) for s in SPECIES)]))
        c = [len(per_species[s]) for s in SPECIES]
        counts.append("\t".join([family, *map(str, c), str(sum(c))]))
        events.append("\t".join((family, node, "n0", "1", "Non-Terminal", ", ".join(a_ids), ", ".join(b_ids))))
    for s in SPECIES:
        write(root / "genomes" / f"{s}.softmasked.fna",
              "".join(f">{name}\n{''.join(parts) or 'N'}\n" for name, parts in sequences[s].items()))
        write(root / "gff" / f"{s}.longest.gff3", "##gff-version 3\n" + "\n".join(annotations[s]) + "\n")
    results = root / "orthofinder"
    prefix = "Orthogroup\t" + "\t".join(SPECIES)
    write(results / "Species_Tree/SpeciesTree_rooted_node_labels.txt", LABELLED + "\n")
    write(results / "Orthogroups/Orthogroups.tsv", prefix + "\n" + "\n".join(membership) + "\n")
    write(results / "Orthogroups/Orthogroups.GeneCount.tsv", prefix + "\tTotal\n" + "\n".join(counts) + "\n")
    write(results / "Gene_Duplication_Events/Duplications.tsv",
          "Orthogroup\tSpecies Tree Node\tGene Tree Node\tSupport\tType\tGenes 1\tGenes 2\n"
          + "\n".join(events) + "\n")
    return results


def table(path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


class ValidationTests(unittest.TestCase):
    def setUp(self):
        temp = TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.results = data(self.root)
        self.run_dir = self.root / "run"

    def awgds(self, *args):
        proc = subprocess.run([sys.executable, "-B", "-m", "awgds", *map(str, args), "--threads", "2"],
                              text=True, capture_output=True)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        return proc

    def scan(self, *extras):
        self.awgds(self.root / "genomes", self.root / "gff", "-o", self.run_dir,
                   "--orthofinder-results", self.results, "--no-require-diamond", "--bootstrap", 500,
                   "--wgd-permutations", 99, "--wgd-max-references", 2, *extras)

    def test_integrated_validation_shares_chromosomes_and_events(self):
        self.scan()
        scores = table(self.run_dir / "summary.tsv")
        self.assertEqual(len(scores), 6)
        result, = [x for x in scores if x["recommended"] == "yes"]
        self.assertEqual(result["node"], "N2")
        self.assertEqual(result["wgd_status"], "shared_macrosynteny")
        self.assertEqual(result["wgd_shared_reference_chromosomes"], "3")
        self.assertEqual(result["wgd_selected_references"], "Taxon_C,Taxon_E")
        self.assertEqual(result["wgd_reference_coverage"], "Taxon_C=3/3;Taxon_E=3/3")
        self.assertEqual(result["wgd_positive_target_reference_combinations"], "4")
        self.assertEqual(len(result["wgd_positive_pairs"].split(";")), 4)
        self.assertEqual(result["wgd_local_check"], "skipped")
        self.assertEqual({x["wgd_status"] for x in scores if x["recommended"] == "no"}, {"not_recommended"})
        self.assertEqual(list(self.run_dir.rglob("summary.tsv")), [self.run_dir / "summary.tsv"])
        events = (self.run_dir / "logs/run.log").read_text(encoding="utf-8")
        self.assertTrue(events.startswith("  ###"))
        self.assertIn("AWGDS v", events)
        self.assertIn("[1/6]", events)
        self.assertIn("[6/6]", events)
        self.assertIn(str(self.run_dir / "summary.tsv"), events)
        links = table(self.run_dir / "work/evidence/shared_linkage.tsv")
        self.assertEqual(len(links), 6)
        self.assertEqual({x["phase_consistent_families"] for x in links}, {"18"})
        self.assertEqual({x["swap_target_chromosomes"] for x in links}, {"1"})

    def test_revalidation_reads_saved_pairs_and_falls_back_to_local(self):
        self.scan("--no-wgd-validation")
        pairs = self.run_dir / "work/branches/pairs/N2.pairs.tsv.gz"
        checksum = hashlib.sha256(pairs.read_bytes()).digest()
        write(self.results / "Gene_Duplication_Events/Duplications.tsv", "Unused during revalidation.\n")
        out = self.root / "validation"
        self.awgds("validate", self.run_dir, "-o", out, "--wgd-permutations", 99,
                   "--wgd-min-shared-chromosomes", 4, "--wgd-max-references", 2)
        result, = [x for x in table(out / "summary.tsv") if x["recommended"] == "yes"]
        self.assertEqual(result["wgd_status"], "insufficient_data")
        self.assertEqual(result["wgd_local_check"], "local_signal_detected")
        self.assertEqual(result["wgd_local_target_pairs"], "Taxon_A|Taxon_B")
        local, = table(out / "work/evidence/local_chains.tsv")
        self.assertEqual(local["node"], "N2")
        self.assertEqual(local["ordered_same_copy_events"], "18")
        self.assertEqual(local["copy_segments_swapped"], "1")
        self.assertEqual(hashlib.sha256(pairs.read_bytes()).digest(), checksum)
        self.assertEqual(json.loads((out / "run.json").read_text())["commands"], [])
        original, = [x for x in table(self.run_dir / "summary.tsv") if x["recommended"] == "yes"]
        self.assertEqual(original["wgd_status"], "validation_disabled")
        self.assertEqual(original["adi_raw"], result["adi_raw"])

    def test_whitelist_excludes_unlisted_genomes_and_sequences(self):
        self.scan("--no-wgd-validation")
        info = self.root / "chromosomes.txt"
        lines = [f"{s}.softmasked.fna\tchr{i}\n" for s in SPECIES if s != "Taxon_B" for i in range(1, 7)]
        write(info, "species\tseqid\n" + "".join(lines))
        out = self.root / "whitelist"
        self.awgds("validate", self.run_dir, "-o", out, "--chromosome-info", info, "--wgd-permutations", 99)
        result, = [x for x in table(out / "summary.tsv") if x["recommended"] == "yes"]
        self.assertEqual(result["wgd_status"], "insufficient_data")
        self.assertEqual(result["wgd_right_eligible_targets"], "0")
        inventory = {x["species"]: x for x in table(out / "work/evidence/selection.tsv")}
        self.assertEqual(inventory["Taxon_B"]["reason"], "fewer_than_two_annotated_chromosomes")
        self.assertEqual(inventory["Taxon_B"]["chromosomes_with_genes"], "0")

    def test_legacy_layout_can_be_validated_into_one_summary(self):
        self.scan("--no-wgd-validation")
        saved = self.root / "legacy"
        saved.mkdir()
        for folder, file in (("inputs", "input_files.tsv"), ("trees", "species_tree.nwk")):
            shutil.copy2(self.run_dir / "work" / folder / file, saved / file)
        shutil.copytree(self.run_dir / "work/branches/pairs", saved / "pairs")
        original = json.loads((self.run_dir / "run.json").read_text())
        original.pop("output_layout")
        write(saved / "run.json", json.dumps(original))
        scores = table(self.run_dir / "summary.tsv")
        columns = [k for k in scores[0] if not k.startswith("wgd_")]
        write(saved / "summary.tsv", "\t".join(columns) + "\n" + "\n".join(
            "\t".join(row[k] for k in columns) for row in scores) + "\n")
        out = self.root / "legacy_validated"
        self.awgds("validate", saved, "-o", out, "--wgd-permutations", 99, "--wgd-max-references", 2)
        result, = [x for x in table(out / "summary.tsv") if x["recommended"] == "yes"]
        self.assertEqual(result["wgd_status"], "shared_macrosynteny")
        self.assertEqual(len(table(out / "summary.tsv")), len(scores))

    def test_revalidation_can_use_the_retained_chromosome_info(self):
        info = self.root / "chromosome_input.txt"
        write(info, "species\tseqid\n" + "".join(f"{s}\tchr{i}\n" for s in SPECIES for i in range(1, 7)))
        self.scan("--no-wgd-validation", "--chromosome-info", info)
        info.rename(self.root / "moved_chromosome_input.txt")
        out = self.root / "saved_map_revalidation"
        self.awgds("validate", self.run_dir, "-o", out, "--wgd-permutations", 99, "--wgd-max-references", 2)
        result, = [r for r in table(out / "summary.tsv") if r["recommended"] == "yes"]
        self.assertEqual(result["wgd_status"], "shared_macrosynteny")
        self.assertTrue((self.run_dir / "work/inputs/chromosome_info.txt").is_file())

    def test_each_family_counts_one_even_with_ambiguous_pairs(self):
        families = {f"g{i}" for i in range(40)}
        reference = {f"g{i}": f"chr{i%2}" for i in range(40)}
        target = {f"g{i}": {("a", str(i % 2))} for i in range(40)}
        target["g0"].update({("a", "1"), ("a", "2")})
        hits = linkage_test(families, reference, target, 99, rng_for(42, "test"))
        self.assertAlmostEqual(sum(x.weight for x in hits), 40)
        self.assertEqual(sum(x.families for x in hits), 42)

    def test_shared_families_require_the_same_duplication_and_one_orientation(self):
        left = {"g1": {"e1": True}, "g2": {"e2": True}, "g3": {"e3": True, "e4": False}}
        right = {"g1": {"e5": False}, "g2": {"e2": False}, "g3": {"e3": True, "e4": True}}
        shared, swapped, events = phase_agreement(left, right)
        self.assertEqual(shared, 2)
        self.assertTrue(swapped)
        self.assertEqual(events, [("g2", "e2")])

    def test_family_count_matching_only_uses_overlapping_comparisons(self):
        first = {f"g{i}" for i in range(48)}
        second = {f"g{i}" for i in range(32)}
        reference = {f"g{i}": f"chr{i%2}" for i in range(48)}
        target = {f"g{i}": {("target_chr", str(i % 2))} for i in range(48)}
        a, b = ("N2", first, reference, target), ("N3", second, reference, target)
        contexts = {("G1", "G2"): [a, b], ("G1", "G3"): [a]}
        settings = SimpleNamespace(seed=41, wgd_matched_draws=4, wgd_matched_permutations=99,
                                   wgd_min_link_families=5, wgd_alpha=.05)
        output = []
        matched_resampling(contexts, settings, SimpleNamespace(writerow=output.append))
        self.assertEqual([(x[2], x[3], x[4], x[-1]) for x in output], [("N2", 48, 32, 4), ("N3", 32, 32, 4)])


if __name__ == "__main__":
    unittest.main()
