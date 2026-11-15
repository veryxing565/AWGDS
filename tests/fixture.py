"""Build a small data set with known branch assignments and two assembly sequences."""

import gzip
from pathlib import Path
import random
import sys


SPECIES = tuple(f"Taxon_{name}" for name in "ABCDEFGH")
LABELLED = ("(((Taxon_A:5,Taxon_B:5)N2:5,(Taxon_C:5,Taxon_D:5)N3:5)N1:10,"
            "((Taxon_E:5,Taxon_F:5)N5:5,(Taxon_G:5,Taxon_H:5)N6:5)N4:10)N0;")
CLADES = {"N1": SPECIES[:4], "N2": SPECIES[:2], "N3": SPECIES[2:4],
          "N4": SPECIES[4:], "N5": SPECIES[4:6], "N6": SPECIES[6:]}
AMINOS = "ACDEFGHIKLMNPQRSTVWY"
CODON = dict(zip(AMINOS, "GCT TGT GAT GAA TTT GGT CAT ATT AAA CTT ATG AAT CCT CAA CGT TCT ACT GTT TGG TAT".split()))


def write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        with gzip.open(path, "wt", encoding="utf-8", newline="\n") as f:
            f.write(data)
    else:
        path.write_text(data, encoding="utf-8", newline="\n")


def fixture(root):
    records, counts, events = [], [], []
    genomes = {s: {"chr1": [], "chr2": []} for s in SPECIES}
    gffs = {s: [] for s in SPECIES}
    rng = random.Random(24601)
    for i in range(56):
        node = ("N1" if i < 15 else "N2" if i < 18 else "N3" if i < 22 else
                "N4" if i < 24 else "N5" if i < 28 else "N6" if i < 31 else None)
        group = f"OG{i:07d}"
        source = tuple(rng.choices(AMINOS, k=110))
        alt = list(source)
        for j in rng.sample(range(len(alt)), 16):
            alt[j] = rng.choice(AMINOS)
        membership = {s: [] for s in SPECIES}
        for s in SPECIES:
            for copy, sequence in (("a", source), ("b", alt)):
                if copy == "b" and (not node or s not in CLADES[node]):
                    continue
                gene = f"{group}_{copy}"
                r = "chr2" if copy == "b" and i != 0 else "chr1"
                aa = list(sequence)
                for j in rng.sample(range(len(aa)), 8):
                    aa[j] = rng.choice(AMINOS)
                coding = "ATG" + "".join(CODON[x] for x in aa) + "TAA"
                start = sum(map(len, genomes[s][r])) + 1
                end = start + len(coding) - 1
                genomes[s][r].append(coding + "A" * 10)
                gffs[s].extend((f"{r}\ttest\tgene\t{start}\t{end}\t.\t+\t.\tID={gene}",
                                f"{r}\ttest\tmRNA\t{start}\t{end}\t.\t+\t.\tID={gene}.1;Parent={gene}",
                                f"{r}\ttest\tCDS\t{start}\t{end}\t.\t+\t0\tID=cds-{gene};Parent={gene}.1"))
                membership[s].append(f"{gene}.1")
        records.append("\t".join([group] + [", ".join(membership[s]) for s in SPECIES]))
        ns = [len(membership[s]) for s in SPECIES]
        counts.append("\t".join([group, *map(str, ns), str(sum(ns))]))
        if node:
            g1 = ", ".join(f"{s}_{group}_a.1" for s in CLADES[node])
            g2 = ", ".join(f"{s}_{group}_b.1" for s in CLADES[node])
            events.append(f"{group}\t{node}\tn0\t1.0\tNon-Terminal\t{g1}\t{g2}")
    for s in SPECIES:
        zipped = ".gz" if s == SPECIES[0] else ""
        write(root / "genomes" / f"{s}.softmasked.fna{zipped}",
              "".join(f">{r} additional information\n{''.join(chunks)}\n" for r, chunks in genomes[s].items()))
        write(root / "gff" / f"{s}.longest.gff3{zipped}", "##gff-version 3\n" + "\n".join(gffs[s]) + "\n")
    results = root / "fixture_orthofinder"
    prefix = "Orthogroup\t" + "\t".join(SPECIES)
    write(results / "Orthogroups/Orthogroups.tsv", prefix + "\n" + "\n".join(records) + "\n")
    write(results / "Orthogroups/Orthogroups.GeneCount.tsv", prefix + "\tTotal\n" + "\n".join(counts) + "\n")
    write(results / "Species_Tree/SpeciesTree_rooted_node_labels.txt", LABELLED + "\n")
    write(results / "Gene_Duplication_Events/Duplications.tsv",
          "Orthogroup\tSpecies Tree Node\tGene Tree Node\tSupport\tType\tGenes 1\tGenes 2\n" + "\n".join(events) + "\n")
    write(root / "fixed.nwk", LABELLED + "\n")
    write(root / "divergence.txt", "species_a\tspecies_b\tmin_mya\tmax_mya\tsource\n"
          "Taxon_A\tTaxon_H\t17\t26\tSynthetic calibration\n"
          "Taxon_A\tTaxon_D\t7\t14\tSynthetic calibration\n"
          "Taxon_E\tTaxon_H\t7\t14\tSynthetic calibration\n")
    return results


if __name__ == "__main__":
    fixture(Path(sys.argv[1]))
