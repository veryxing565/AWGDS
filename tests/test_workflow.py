"""Regression checks that do not need external executables."""

import csv
import gzip
import json
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from Bio.Seq import Seq

from awgds.dating import load_calibrations, topology
from awgds.io import InputError
from awgds.phylogeny import SpeciesTree, fixed_time_clock, posterior_time_clock, read_newick
from awgds.proteins import Model, Segment, translate
from awgds.statistics import percent_value
from fixture import CLADES, LABELLED, SPECIES, fixture, write


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.results = fixture(self.root)

    def cli(self, outdir, *args):
        return subprocess.run([sys.executable, "-B", "-m", "awgds", str(self.root / "genomes"),
                               str(self.root / "gff"), "-o", str(self.root / outdir),
                               "--orthofinder-results", str(self.results), "--bootstrap", "1200", *args],
                              text=True, capture_output=True)

    def test_end_to_end_pairs_and_filters(self):
        for out, flags in (("default_coords", ()), ("all_coords", ("--no-cross-replicon",))):
            r = self.cli(out, "--no-require-diamond", *flags)
            self.assertEqual(r.returncode, 0, r.stderr)
            target = self.root / out
            with open(target / "work/branches/branch_counts.tsv", newline="") as f:
                counts = {x["node"]: x for x in csv.DictReader(f, delimiter="\t")}
            n = 14 if not flags else 15
            self.assertEqual(int(counts["N1"]["family_species_pairs"]), n * 4)
            self.assertEqual(float(counts["N1"]["mean_families_per_pair"]), n)
            with gzip.open(target / "work/branches/pairs/N1.pairs.tsv.gz", "rt") as f:
                pairs = list(csv.DictReader(f, delimiter="\t"))
            self.assertEqual(len(pairs), n * 4)
            self.assertEqual({p["species"] for p in pairs}, set(CLADES["N1"]))
            self.assertEqual({p["used_in_adi"] for p in pairs}, {"1"})
            with open(target / "summary.tsv", newline="") as f:
                summary = list(csv.DictReader(f, delimiter="\t"))
            self.assertEqual(len(summary), 6)
            self.assertEqual(max(float(s["adi_percent"]) for s in summary), 100.)
            self.assertEqual(min(float(s["adi_percent"]) for s in summary), 0.)
            fence = json.loads((target / "run.json").read_text())["adi"]
            self.assertAlmostEqual(fence["upper_fence"], fence["q3"] + 1.5*fence["iqr"])
            for s in summary:
                self.assertEqual(s["recommended"] == "yes", float(s["adi_raw"]) > fence["upper_fence"])

    def test_reused_results_names_are_strict(self):
        tree = self.results / "Species_Tree/SpeciesTree_rooted_node_labels.txt"
        write(tree, LABELLED.replace("Taxon_A", "Taxon_Z"))
        r = self.cli("bad", "--no-require-diamond")
        self.assertEqual(r.returncode, 1, r.stderr)
        self.assertIn("Species names differ", r.stderr)
        self.assertEqual(json.loads((self.root / "bad/run.json").read_text())["status"], "failed")

    def test_fixed_clocks_with_units_and_an_extra_outgroup(self):
        external = self.root / "FigTree.tre"
        tiny = LABELLED[:-1].replace(":5", ":0.05").replace(":10", ":0.10")
        write(external, "#NEXUS\nBEGIN TREES;\nUTREE 1 = (" + tiny + ":0.10,Extra:0.30);\nEND;\n")
        tree = SpeciesTree(self.root / "fixed.nwk", SPECIES)
        a = fixed_time_clock(tree, self.root / "fixed.nwk", "Myr", "a").durations
        b = fixed_time_clock(tree, external, "100Myr", "b").durations
        for i, j in zip(a[0], b[0]):
            self.assertAlmostEqual(i, j)
        r = self.cli("fixed_clocks", "--no-require-diamond", "--time-tree", str(self.root / "fixed.nwk"),
                     "--time-tree", f"paml={external}", "--time-unit", "Myr", "--time-unit", "100Myr")
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(self.root / "fixed_clocks/summary.tsv", newline="") as f:
            for s in csv.DictReader(f, delimiter="\t"):
                self.assertAlmostEqual(float(s["clock_density_fixed"]), float(s["clock_density_paml"]))
                self.assertEqual(s["density_unit"], "families_per_1000_per_Myr")
        target = self.root / "fixed_clocks"
        manifest = json.loads((target / "run.json").read_text())
        saved = manifest["preserved_clock_inputs"]["paml"]["tree_or_report"]
        self.assertEqual((target / saved).read_bytes(), external.read_bytes())

    def test_orthofinder_input_tree_keeps_branch_lengths(self):
        original = read_newick(LABELLED)
        clean = read_newick(topology(original.root, lengths=True) + ";")
        self.assertEqual([c.branch_length for c in original.find_clades()],
                         [c.branch_length for c in clean.find_clades()])

    def test_paml_numbered_nodes_and_trace_pairing(self):
        tree = SpeciesTree(self.root / "fixed.nwk", SPECIES)
        numbered = ("(((1_Taxon_A,2_Taxon_B)11,(3_Taxon_C,4_Taxon_D)12)10,"
                    "((5_Taxon_E,6_Taxon_F)14,(7_Taxon_G,8_Taxon_H)15)13)9;")
        report, trace = self.root / "mcmctree.out", self.root / "mcmc.txt"
        write(report, f"Species tree for FigTree\n{numbered}\n\n{LABELLED}\n")
        header = "Gen " + " ".join(f"t_n{i}" for i in range(9, 16))
        write(trace, header + "\n1 .50 .40 .10 .11 .40 .11 .11\n"
              "20 .20 .10 .05 .05 .10 .05 .05\n40 .30 .12 .06 .06 .10 .05 .05\n"
              "60 .20 .10 .05 .05 .10 .05 .05\n")
        clock, chains, _ = posterior_time_clock(tree, report, [trace], "100Myr", "test")
        self.assertEqual(clock.durations.shape, (3, 6))
        self.assertAlmostEqual(clock.durations[1, tree.scored.index("N1")], 18.)
        self.assertEqual(len(chains), 1)
        r = self.cli("posterior", "--no-require-diamond", "--mcmc-posterior", "imported", str(report),
                     str(trace), "--time-unit", "100Myr")
        self.assertEqual(r.returncode, 0, r.stderr)
        target = self.root / "posterior"
        manifest = json.loads((target / "run.json").read_text())
        saved = manifest["preserved_clock_inputs"]["imported"]
        self.assertEqual((target / saved["tree_or_report"]).read_bytes(), report.read_bytes())
        with gzip.open(target / saved["traces"][0], "rb") as f:
            self.assertEqual(f.read(), trace.read_bytes())

    def test_calibrations_join_to_mrca(self):
        tree = SpeciesTree(self.root / "fixed.nwk", SPECIES)
        result = load_calibrations(tree, self.root / "divergence.txt")
        self.assertEqual(set(result), {"N0", "N1", "N4"})
        with self.assertRaises(InputError):
            write(self.root / "divergence.txt", "species_a\tspecies_b\tmin_mya\tmax_mya\nTaxon_A\tTaxon_B\t1\t2\n")
            load_calibrations(tree, self.root / "divergence.txt")
        write(self.root / "divergence.txt", "species_a\tspecies_b\tage_mya\nTaxon_A\tTaxon_B\t10\n")
        with_root = load_calibrations(tree, self.root / "divergence.txt", point_width=.25, root_max=50)
        self.assertEqual((with_root["N2"].lower, with_root["N2"].upper), (7.5, 12.5))
        self.assertEqual((with_root["N0"].lower, with_root["N0"].upper), (None, 50))

    def test_replicon_aliases_can_join_scaffolds(self):
        mapping = self.root / "replicons.tsv"
        write(mapping, "species\tseqid\treplicon\nTaxon_A\tchr1\tch1\nTaxon_A\tchr2\tch1\n"
              "Taxon_B\tchr1\tch1\nTaxon_B\tchr2\tunplaced\n")
        r = self.cli("map", "--no-require-diamond", "--replicon-map", str(mapping))
        self.assertEqual(r.returncode, 0, r.stderr)
        with gzip.open(self.root / "map/work/branches/pairs/N1.pairs.tsv.gz", "rt") as f:
            pairs = list(csv.DictReader(f, delimiter="\t"))
        self.assertEqual(pairs, [])

    def test_negative_strand_codon_can_span_cds_segments(self):
        # ATG GAA CTT TAA, split after 4 bases on the coding strand.
        dna = Seq("TTAAAGTT" + "TTTT" + "CCAT")
        s1 = Segment("chr1", 13, 16, "-", 0, "1")
        s2 = Segment("chr1", 1, 8, "-", 2, "1")
        m = Model("t1", "g1", "t1", "p1", frozenset(), (s1, s2))
        self.assertEqual(translate(m, dna, 1), ("MEL", ""))

    def test_constant_percentage_and_above_scale_fence(self):
        self.assertEqual(percent_value(0.75, 0.75, 0.75), 0.)
        self.assertEqual(percent_value(2., 0., 1.), 200.)


if __name__ == "__main__":
    unittest.main()
