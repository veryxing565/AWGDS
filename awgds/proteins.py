"""Translate GFF3 coding sequences and retain their identifiers and coordinates."""

from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path
from urllib.parse import unquote

from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq

from .io import InputError, LOG, open_text, tsv_writer


@dataclass(frozen=True)
class Segment:
    seqid: str
    start: int
    end: int
    strand: str
    phase: int | None
    table: str


@dataclass
class Model:
    target: str
    gene: str
    transcript: str
    protein_id: str
    aliases: frozenset
    segments: tuple
    translation: str = ""

    @cached_property
    def seqid(self):
        return self.segments[0].seqid

    @cached_property
    def start(self):
        return min(s.start for s in self.segments)

    @cached_property
    def end(self):
        return max(s.end for s in self.segments)

    @cached_property
    def strand(self):
        return self.segments[0].strand


def attributes(text, path, number):
    data = {}
    for item in text.split(";"):
        if not item.strip() or item == ".":
            continue
        field, equals, value = item.partition("=")
        if not equals or not field:
            raise InputError(f"{path}:{number}: expected GFF3 key=value attributes, not GTF.")
        values = value.split(",") if field in {"Parent", "Alias"} else [value]
        data[field] = [unquote(part) for part in values if part]
    return data


def first(fields, *keys):
    return next((fields[k][0] for k in keys if fields.get(k)), "")


def valid_protein_id(value):
    return bool(value) and not any(c.isspace() or c in ",;()[]:'\"=" for c in value)


def parse_gff(path, default_table):
    features, coding = {}, []
    with open_text(path) as stream:
        for number, line in enumerate(stream, 1):
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) != 9:
                raise InputError(f"{path}:{number}: expected 9 GFF3 columns.")
            seqid, _, kind, start, end, _, strand, phase, raw = parts
            if kind in {"exon", "region", "match", "match_part"}:
                continue
            fields = attributes(raw, path, number)
            feature_id = first(fields, "ID")
            if kind == "CDS":
                try:
                    begin, stop = int(start), int(end)
                    offset = int(phase) if phase != "." else None
                    if begin < 1 or stop < begin or strand not in {"+", "-"} or offset not in {0, 1, 2, None}:
                        raise ValueError
                except ValueError as exc:
                    raise InputError(f"{path}:{number}: invalid CDS coordinates, strand or phase.") from exc
                coding.append((Segment(unquote(seqid), begin, stop, strand, offset,
                                       first(fields, "transl_table") or str(default_table)), fields))
            elif feature_id:
                if feature_id in features and features[feature_id] != (kind, fields):
                    # Some discontinuous non-coding features share an ID. Only the attributes
                    # used for locating a CDS parent must agree.
                    old_kind, old_fields = features[feature_id]
                    if old_kind != kind or old_fields.get("Parent") != fields.get("Parent"):
                        raise InputError(f"{path}:{number}: conflicting parent information for {feature_id!r}.")
                else:
                    features[feature_id] = (kind, fields)
    if not coding:
        raise InputError(f"{path}: no CDS features. A gene-only annotation cannot produce proteins.")

    grouped = defaultdict(list)
    for piece, fields in coding:
        pid = first(fields, "protein_id")
        tid = first(fields, "transcript_id")
        targets = fields.get("Parent") or [tid or pid or first(fields, "ID")]
        for target in targets:
            if not target:
                raise InputError(f"{path}: a CDS is missing Parent, transcript_id, protein_id and ID.")
            kind = features.get(target, ("", {}))[0]
            grouped[(target, (tid or pid) if kind in {"gene", "pseudogene"} else "")].append((piece, fields))

    models = []
    for (target, isoform), chunks in grouped.items():
        kind, feature = features.get(target, ("", {}))
        parent_fields = features.get(first(feature, "Parent"), ("", {}))[1]
        pid_set = {first(f, "protein_id") for _, f in chunks} - {""}
        tid_set = {first(f, "transcript_id") for _, f in chunks} - {""}
        if len(pid_set) > 1 or len(tid_set) > 1:
            raise InputError(f"{path}: inconsistent protein_id/transcript_id for {target!r}.")
        pid = next(iter(pid_set), first(feature, "protein_id"))
        transcript = (next(iter(tid_set), "") or first(feature, "transcript_id")
                      or isoform or target)
        gene = (target if kind in {"gene", "pseudogene"} else first(feature, "Parent") or
                first(feature, "gene_id") or first(chunks[0][1], "gene_id") or target)
        gene_fields = features.get(gene, ("", {}))[1]
        all_aliases = {target, isoform, transcript, pid, gene}
        for source in (feature, parent_fields, gene_fields, *(parts for _, parts in chunks)):
            for key in ("ID", "Alias", "protein_id", "transcript_id", "gene_id", "locus_tag"):
                all_aliases.update(source.get(key, ()))
        # Permit NCBI prefixes when they resolve uniquely, while keeping exact aliases first.
        for alias in tuple(all_aliases):
            for prefix in ("rna-", "cds-", "gene-"):
                if alias.startswith(prefix):
                    all_aliases.add(alias[len(prefix):])
        parts = tuple(sorted({p for p, _ in chunks}, key=lambda x: (x.start, x.end),
                             reverse=chunks[0][0].strand == "-"))
        models.append(Model(target, gene, transcript, pid, frozenset(all_aliases - {""}), parts))
    return models


def translate(model, chromosome, min_aa):
    parts = model.segments
    if len({(p.seqid, p.strand, p.table) for p in parts}) != 1:
        return "", "mixed_sequences_strands_or_codes"
    if any(p.phase is None for p in parts):
        return "", "missing_phase"
    pieces, coding_bases = [], 0
    previous = None
    for index, segment in enumerate(parts):
        if segment.end > len(chromosome):
            raise InputError(f"{segment.seqid}:{segment.start}-{segment.end}: CDS lies outside FASTA sequence.")
        if previous and ((segment.start <= previous.end) if model.strand == "+"
                         else (segment.end >= previous.start)):
            return "", "overlapping_cds_parts"
        if index and segment.phase != (-coding_bases) % 3:
            return "", "inconsistent_phase"
        piece = chromosome[segment.start-1:segment.end]
        if model.strand == "-":
            piece = piece.reverse_complement()
        if not index:
            piece = piece[segment.phase:]
        pieces.append(str(piece))
        coding_bases += len(piece)
        previous = segment
    seq = "".join(pieces).upper().replace("U", "T")
    seq = seq[:len(seq) - len(seq) % 3]
    try:
        protein = str(Seq(seq).translate(table=int(parts[0].table))).removesuffix("*")
    except (CodonTable.TranslationError, ValueError, KeyError):
        return "", "invalid_codon_or_code"
    if "*" in protein:
        return "", "internal_stop"
    if len(protein) < min_aa:
        return "", "short_protein"
    return protein, ""


class Proteome:
    def __init__(self, species, genome, annotation, code, min_aa, sequences=True):
        self.species = species
        self.genome = genome
        self.gff = annotation
        self.models = parse_gff(annotation, code)
        self.aliases = defaultdict(list)
        self.preferred = defaultdict(list)
        self.by_seqid = defaultdict(list)
        self.replicons = {}
        self.rejected = Counter()
        for model in self.models:
            if len({(part.seqid, part.strand) for part in model.segments}) == 1:
                self.by_seqid[model.seqid].append(model)
            else:
                self.rejected["mixed_sequences_or_strands"] += 1
                continue
            for name in model.aliases:
                self.aliases[name].append(model)
            for name in {model.transcript, model.protein_id, model.target} - {""}:
                self.preferred[name].append(model)
        def check_record(seqid, length, sequence=None):
            if seqid in self.replicons:
                raise InputError(f"{genome}: duplicate FASTA sequence ID {seqid!r}.")
            self.replicons[seqid] = length
            for model in self.by_seqid.get(seqid, ()):
                if model.end > length:
                    raise InputError(f"{annotation}: CDS {model.target!r} extends beyond "
                                     f"{seqid!r} in {genome.name}.")
                if sequence is not None:
                    model.translation, failure = translate(model, sequence, min_aa)
                    if failure:
                        self.rejected[failure] += 1

        index = Path(str(genome) + ".fai")
        use_index = (not sequences and genome.suffix.lower() != ".gz" and index.is_file() and
                     index.stat().st_mtime >= genome.stat().st_mtime)
        if use_index:
            with open_text(index) as f:
                for number, line in enumerate(f, 1):
                    fields = line.rstrip("\r\n").split("\t")
                    try:
                        if len(fields) != 5 or not fields[0] or int(fields[1]) < 0:
                            raise ValueError
                        length = int(fields[1])
                    except ValueError as exc:
                        raise InputError(f"{index}:{number}: malformed FASTA index.") from exc
                    check_record(fields[0], length)
        else:
            with open_text(genome) as stream:
                for record in SeqIO.parse(stream, "fasta"):
                    check_record(record.id, len(record), record.seq if sequences else None)
        if not self.replicons:
            raise InputError(f"{genome}: no FASTA sequences or index records.")
        missing = sorted(self.by_seqid.keys() - self.replicons.keys())
        if missing:
            raise InputError(f"{annotation}: GFF CDS seqid not present in {genome.name}: "
                             f"{', '.join(missing[:8])}. Match GFF column 1 to the FASTA first token.")
        self.primary = {}
        for m in self.models:
            if m.translation and (m.gene not in self.primary or
                                  (-len(m.translation), m.transcript, m.target) <
                                  (-len(self.primary[m.gene].translation),
                                   self.primary[m.gene].transcript, self.primary[m.gene].target)):
                self.primary[m.gene] = m
        if sequences and not self.primary:
            raise InputError(f"{species}: no usable translated proteins. Check CDS phase/genetic code "
                             f"({dict(self.rejected)}).")

    def locate(self, gene_id):
        matches = self.preferred.get(gene_id) or self.aliases.get(gene_id, [])
        if not matches:
            return None
        if len(matches) == 1:
            return matches[0]
        # A reused protein may be named by its gene. Pick the documented longest isoform
        # only if all matches belong to the same locus, never across loci.
        if len({m.gene for m in matches}) == 1:
            gene = matches[0].gene
            chosen = self.primary.get(gene)
            if chosen is not None and chosen in matches:
                return chosen
            if len({(m.seqid, m.start, m.end, m.strand) for m in matches}) == 1:
                return matches[0]
        raise InputError(f"{self.species}: GFF identifier {gene_id!r} matches multiple CDS models. "
                         "Use unique transcript/protein IDs or --gene-id-map.")

    def write_fasta(self, path, id_kind, prefix_ids=False):
        selected = {}
        with path.open("w", encoding="utf-8", newline="\n") as fasta:
            for model in sorted(self.primary.values(), key=lambda m: (m.seqid, m.start, m.transcript)):
                requested = {"gene": model.gene, "protein": model.protein_id or model.transcript,
                             "transcript": model.transcript}[id_kind]
                if not valid_protein_id(requested):
                    raise InputError(f"{self.gff}: protein ID {requested!r} is unsafe for a gene tree.")
                if prefix_ids:
                    requested = f"{self.species}_{requested}"
                if requested in selected:
                    raise InputError(f"{self.gff}: duplicate output protein ID {requested!r}. "
                                     "Choose another --protein-id scheme.")
                selected[requested] = model
                fasta.write(f">{requested}\n{model.translation}\n")
        return selected


def extract(inputs, folder, table, min_aa, id_kind, sequences=True, write_proteins=True, prefix_ids=False):
    if write_proteins:
        folder.mkdir(parents=True, exist_ok=True)
    genomes = {}
    with (folder.parent / "extraction.tsv").open("w", encoding="utf-8", newline="") as stream:
        out = tsv_writer(stream, ("species", "genome", "annotation", "cds_models",
                                  "primary_proteins", "rejected_models"))
        for i, (species, (genome, annotation)) in enumerate(inputs.items(), 1):
            LOG.info("Reading %s CDS annotations (%d/%d): %s",
                     "and translating" if sequences else "without translation", i, len(inputs), species)
            p = Proteome(species, genome, annotation, table, min_aa, sequences)
            if write_proteins:
                p.write_fasta(folder / f"{species}.faa", id_kind, prefix_ids)
            genomes[species] = p
            summary = ",".join(f"{k}={v}" for k, v in sorted(p.rejected.items()))
            out.writerow((species, genome, annotation, len(p.models), len(p.primary), summary))
            if summary:
                LOG.warning("%s: %s (see extraction.tsv)", species, summary)
    return genomes

