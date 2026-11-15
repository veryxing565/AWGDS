"""Validate OrthoFinder results and connect its gene identifiers to GFF3 models."""

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from Bio.SeqIO.FastaIO import SimpleFastaParser

from .io import InputError, LOG, ids, open_text, read_table, same_names
from .phylogeny import SpeciesTree

TREE_FILE = Path("Species_Tree/SpeciesTree_rooted_node_labels.txt")
DUP_FILE = Path("Gene_Duplication_Events/Duplications.tsv")
COUNT_FILE = Path("Orthogroups/Orthogroups.GeneCount.tsv")
OG_FILE = Path("Orthogroups/Orthogroups.tsv")


def tab_file(directory, relative):
    for path in [directory / relative, directory / (str(relative) + ".gz")]:
        if path.is_file():
            return path
    return None


def find_results(directory):
    if not directory.is_dir():
        raise InputError(f"OrthoFinder directory does not exist: {directory}")
    if tab_file(directory, DUP_FILE):
        return directory
    options = sorted({p.parent.parent for pattern in ("Duplications.tsv", "Duplications.tsv.gz")
                      for p in directory.rglob(pattern) if p.parent.name == "Gene_Duplication_Events"})
    if len(options) == 1:
        return options[0]
    if options:
        raise InputError(f"{directory}: found {len(options)} OrthoFinder runs. Specify one Results_* directory.")
    raise InputError(f"{directory}: cannot find Gene_Duplication_Events/Duplications.tsv. "
                     "Use results from a complete OrthoFinder run.")


@dataclass(frozen=True)
class Group:
    genes: dict
    counts: dict


class OrthoResults:
    def __init__(self, directory, species, read_working_ids=True):
        self.root = find_results(directory)
        path = self.root / TREE_FILE
        if not path.is_file():
            raise InputError(f"Missing {path}; the node labels link duplications to branches.")
        self.tree = SpeciesTree(path, species)
        self.duplications = tab_file(self.root, DUP_FILE)
        self.groups = self.load_groups()
        if not self.groups:
            raise InputError(f"{self.root}: orthogroup table is empty.")
        self.numeric = {}
        self.protein_files = self.index_working_proteins() if read_working_ids else {}

    def index_working_proteins(self):
        wd = self.root / "WorkingDirectory"
        species_ids = wd / "SpeciesIDs.txt"
        sequence_ids = wd / "SequenceIDs.txt"
        if not species_ids.is_file() or not sequence_ids.is_file():
            return {}
        codes = {}
        with open_text(species_ids) as f:
            for line in f:
                if not line.strip() or line.lstrip().startswith("#"):
                    continue
                code, sep, name = line.partition(":")
                if sep and code.strip().isdigit():
                    basename = name.strip().replace("\\", "/").rsplit("/", 1)[-1].removesuffix(".gz")
                    codes[code.strip()] = basename.rsplit(".", 1)[0] if "." in basename else basename
        if set(codes.values()) != self.tree.species:
            LOG.warning("WorkingDirectory species labels differ from the results tree. Skipping its ID map and proteins.")
            return {}
        with open_text(sequence_ids) as f:
            for line in f:
                number, sep, name = line.partition(":")
                number, fields = number.strip(), name.strip().split(maxsplit=1)
                species = codes.get(number.split("_", 1)[0])
                if sep and species and fields:
                    self.numeric[number] = species, fields[0]
        files = {species: tab_file(wd, Path(f"Species{code}.fa")) for code, species in codes.items()}
        return files if all(files.values()) and self.numeric else {}

    def read_working_proteins(self, minimum):
        proteins = {}
        for s, path in sorted(self.protein_files.items()):
            by_id = {}
            with open_text(path) as f:
                for header, sequence in SimpleFastaParser(f):
                    token = header.split()[0] if header.strip() else ""
                    assigned = self.numeric.get(token)
                    if assigned is None or assigned[0] != s:
                        raise InputError(f"{path}: FASTA ID {token!r} cannot be linked to {s} via SequenceIDs.txt. "
                                         "Use --retranslate for the supplied genome/GFF proteins.")
                    gene = assigned[1]
                    protein = sequence.upper().removesuffix("*")
                    if len(protein) < minimum or "*" in protein:
                        continue
                    if gene in by_id and by_id[gene] != protein:
                        raise InputError(f"{path}: {gene!r} has multiple archived protein sequences.")
                    by_id[gene] = protein
            proteins[s] = by_id
        LOG.info("Read %d usable archived proteins from OrthoFinder WorkingDirectory.",
                 sum(map(len, proteins.values())))
        return proteins

    def load_groups(self):
        group_path = tab_file(self.root, OG_FILE)
        counts_path = tab_file(self.root, COUNT_FILE)
        if not group_path:
            # Some result layouts provide the ancestral groups in the root HOG table.
            group_path = tab_file(self.root,
                                  Path(f"Phylogenetic_Hierarchical_Orthogroups/{self.tree.root}.tsv"))
        if not group_path and not counts_path:
            raise InputError(f"{self.root}: requires Orthogroups.tsv or Orthogroups.GeneCount.tsv.")
        taxa = sorted(self.tree.species)

        def read(path, members):
            with open_text(path) as f:
                header = f.readline().rstrip("\r\n").split("\t")
            identifier = "Orthogroup" if "Orthogroup" in header else "OG"
            columns = set(header) - {"Orthogroup", "OG", "HOG", "Gene Tree Parent Clade", "Total"}
            same_names(taxa, columns, "OrthoFinder tree tips", path.name)
            result = {}
            for r in read_table(path, [identifier, *taxa]):
                label = r[identifier].strip()
                if not label:
                    raise InputError(f"{path}: empty orthogroup ID.")
                if members:
                    d = result.setdefault(label, {})
                    for s in taxa:
                        if r[s].strip():
                            d.setdefault(s, set()).update(ids(r[s]))
                else:
                    if label in result:
                        raise InputError(f"{path}: repeated orthogroup {label}.")
                    try:
                        result[label] = {s: int(r[s]) for s in taxa}
                    except ValueError as exc:
                        raise InputError(f"{path}: invalid copy count in {label}.") from exc
                    if min(result[label].values()) < 0:
                        raise InputError(f"{path}: negative copy count in {label}.")
            return result

        # Root HOG rows can represent only fragments of a family. Prefer the standard
        # membership/count tables for the IDs used by Duplications.tsv.
        genes = read(group_path, True) if group_path else {}
        counts = read(counts_path, False) if counts_path else {}
        if group_path and group_path.name.startswith("Orthogroups.tsv") and counts_path:
            if set(counts) != set(genes) or any(
                len(genes[g].get(s, ())) != n for g,c in counts.items() for s,n in c.items()
            ):
                raise InputError("OrthoFinder membership and copy-count tables disagree. "
                                 "Use tables from the same complete result folder.")
        self.by_gene = defaultdict(set)
        groups = {}
        for g in dict.fromkeys([*counts, *genes]):
            members = {s: frozenset(a) for s,a in genes.get(g, {}).items()}
            for s, items in members.items():
                for gene in items:
                    self.by_gene[gene].add(s)
            groups[g] = Group(members, counts[g] if g in counts else
                              {s: len(a) for s,a in members.items()})
        if not genes:
            LOG.warning("No OrthoFinder gene membership table. Same-OG evidence comes from "
                        "the Duplications.tsv row and available name mappings.")
        return groups


def load_gene_map(path, genomes):
    if path is None:
        return {}
    lookup = {}
    for row in read_table(path, ("species", "orthofinder_id", "gff_id")):
        s, a, b = (row[k].strip() for k in ("species", "orthofinder_id", "gff_id"))
        if s not in genomes or not a or not b or (s, a) in lookup:
            raise InputError(f"{path}: invalid or repeated gene map entry {s}/{a}.")
        if genomes[s].locate(b) is None:
            raise InputError(f"{path}: GFF ID {b!r} is missing for {s}.")
        lookup[s, a] = b
    return lookup


@dataclass
class Gene:
    species: str
    identifier: str
    raw: str
    model: object

    @property
    def key(self):
        return self.species, self.identifier


class GeneResolver:
    def __init__(self, results, genomes, aliases, allow_missing=False):
        self.results, self.genomes, self.aliases = results, genomes, aliases
        self.species = sorted(genomes, key=lambda x: (-len(x), x))
        self.cache = {}
        self.allow_missing = allow_missing
        self.unmapped = set()
        self.numeric = results.numeric

    def resolve(self, raw, og, descendants):
        g = self.results.groups.get(og)
        known = g.genes if g is not None else {}
        assigned = self.numeric.get(raw)
        if assigned is None:
            exact_og = [(s, raw) for s in descendants if raw in known.get(s, ())]
            prefixes = [(s, raw[len(s)+1:]) for s in self.species
                        if s in descendants and raw.startswith(s + "_")]
            prefixed_og = [(s, gene) for s, gene in prefixes if gene in known.get(s, ())]
            family_matches = set(exact_og + prefixed_og)
            if len(family_matches) == 1:
                assigned = next(iter(family_matches))
            elif family_matches:
                raise InputError(f"OrthoFinder identifier {raw!r} is ambiguous within {og}. "
                                 "Use unambiguous species-prefixed gene names.")
            else:
                for s, gene in prefixes:
                    if ((s, gene) in self.aliases or (s, raw) in self.aliases or
                            self.genomes[s].locate(gene) is not None):
                        assigned = s, gene
                        break
                if assigned is None:
                    matches = [(s, raw) for s in self.results.by_gene.get(raw, ()) if s in descendants]
                    if not matches:
                        matches = [(s, raw) for s in descendants
                                   if (s, raw) in self.aliases or self.genomes[s].locate(raw) is not None]
                    if len(matches) > 1:
                        raise InputError(f"OrthoFinder identifier {raw!r} is shared by multiple species. "
                                         "Use Duplications.tsv with species prefixes.")
                    if matches:
                        assigned = matches[0]
                    elif prefixes:
                        assigned = prefixes[0]
        if assigned is None or assigned[0] not in descendants:
            raise InputError(f"Cannot identify the species for OrthoFinder gene {raw!r} in {og}. "
                             "Check the input species names and use a species-prefixed duplication table.")
        species, identifier = assigned
        if known and identifier not in known.get(species, ()):
            prefixed = species + "_" + identifier
            if prefixed in known.get(species, ()):
                identifier = prefixed
        key = species, identifier
        if key in self.cache:
            return self.cache[key]
        lookup = self.aliases.get((species, raw), self.aliases.get(key, identifier))
        model = locate_annotation(self.genomes[species], lookup)
        if model is None:
            self.unmapped.add(key)
            if not self.allow_missing:
                raise InputError(f"OrthoFinder gene {raw!r} ({species}/{identifier}) has no GFF CDS match. "
                                 "Check that the assembly/annotation versions match, or use --gene-id-map.")
        gene = Gene(species, identifier, raw, model)
        self.cache[key] = gene
        return gene

    def in_group(self, gene, og):
        group = self.results.groups.get(og)
        if group is None or not group.genes:
            return True   # stated in the OrthoFinder duplication row
        names = group.genes.get(gene.species, ())
        return gene.identifier in names or gene.raw in names


def protein_sequence(archives, species, identifier, model=None, alternate=None):
    if archives is not None:
        sequences = archives[species]
        prefixed = identifier.removeprefix(species + "_")
        return sequences.get(identifier) or sequences.get(alternate) or sequences.get(prefixed, "")
    return model.translation if model is not None else ""


def locate_annotation(proteome, identifier, allow_ambiguous_skip=False):
    try:
        model = proteome.locate(identifier)
        unprefixed = identifier.removeprefix(proteome.species + "_")
        if model is None and unprefixed != identifier:
            model = proteome.locate(unprefixed)
        return model
    except InputError:
        if not allow_ambiguous_skip:
            raise
        return None

