"""DIAMOND support for intragenomic pairs and optional chromosome aliases."""

from collections import defaultdict
from dataclasses import dataclass
from itertools import product
import math

from .io import InputError, LOG, executable, file_key, open_text, read_table, tsv_writer
from .orthology import protein_sequence


@dataclass(frozen=True)
class Hit:
    evalue: float
    bitscore: float
    identity: float | None = None
    coverage_1: float | None = None
    coverage_2: float | None = None

    @property
    def priority(self):
        return self.bitscore, -self.evalue, self.identity if self.identity is not None else -1


HIT_FIELDS = ("species", "gene1", "gene2", "evalue", "bitscore", "pident", "qcov", "scov")


def pair_key(x, y):
    return (x, y) if x <= y else (y, x)


def read_replicon_map(path, genomes):
    if path is None:
        return {}
    result = defaultdict(dict)
    for r in read_table(path, ("species", "seqid", "replicon")):
        s, seqid, chromosome = (r[k].strip() for k in ("species", "seqid", "replicon"))
        if s not in genomes or seqid not in genomes[s].replicons or seqid in result[s]:
            raise InputError(f"{path}: unknown or duplicate sequence {s}/{seqid}.")
        result[s][seqid] = "" if chromosome.casefold() in {"", ".", "na", "unplaced"} else chromosome
    return dict(result)


def read_chromosome_info(path, genomes):
    """Read an explicit chromosome whitelist for every input genome."""
    result = {name: {} for name in genomes}
    with open_text(path) as stream:
        found_row = False
        for line_no, line in enumerate(stream, 1):
            fields = line.split()
            if not fields or fields[0].startswith("#"):
                continue
            if not found_row and fields[:2] in (["species", "seqid"], ["genome", "seqid"]):
                found_row = True
                if len(fields) not in (2, 3) or len(fields) == 3 and fields[2] != "chromosome":
                    raise InputError(f"{path}:{line_no}: expected 'species seqid [chromosome]'.")
                continue
            found_row = True
            if len(fields) not in (2, 3):
                raise InputError(f"{path}:{line_no}: expected 'species seqid [chromosome]'.")
            filename, seqid, *alias = fields
            species = file_key(filename.replace("\\", "/"))
            if species not in genomes or seqid not in genomes[species].replicons:
                raise InputError(f"{path}:{line_no}: unknown genome or FASTA sequence: {filename}/{seqid}.")
            if seqid in result[species]:
                raise InputError(f"{path}:{line_no}: repeated sequence: {species}/{seqid}.")
            chromosome = alias[0] if alias else seqid
            result[species][seqid] = "" if chromosome.casefold() in {".", "na", "unplaced"} else chromosome
    if not any(any(x.values()) for x in result.values()):
        raise InputError(f"{path}: no chromosomes were listed.")
    return result


def sequence_replicon(species, seqid, mapping):
    return mapping[species].get(seqid, "") if species in mapping else seqid


def replicon(gene, mapping):
    if gene.model is None:
        return ""
    return sequence_replicon(gene.species, gene.model.seqid, mapping)


def pass_hit(hit, opts):
    if hit.evalue > opts.evalue:
        return False
    if opts.min_identity and (hit.identity is None or hit.identity < opts.min_identity):
        return False
    if opts.min_coverage and (hit.coverage_1 is None or hit.coverage_2 is None or
                              min(hit.coverage_1, hit.coverage_2) < opts.min_coverage):
        return False
    return True


def parse_hit(vals, path):
    ev, bit, identity, qcov, scov = vals
    try:
        evalue, bitscore = float(ev), float(bit)
        extras = [float(x) if x not in {"", ".", "NA"} else None
                  for x in (identity, qcov, scov)]
        if not all(math.isfinite(x) and x >= 0 for x in (evalue, bitscore)) or any(
            x is not None and (not math.isfinite(x) or not 0 <= x <= 100) for x in extras
        ):
            raise ValueError
    except ValueError as exc:
        raise InputError(f"{path}: invalid DIAMOND numeric fields.") from exc
    return Hit(evalue, bitscore, *extras)


def trusted_hits(path, species, allowed, opts):
    collected = defaultdict(dict)
    supplied, matched = 0, 0
    for r in read_table(path, HIT_FIELDS[:5]):
        s = r["species"]
        if s not in species:
            raise InputError(f"{path}: DIAMOND table has unknown species {s!r}.")
        supplied += 1
        a, b = r["gene1"].strip(), r["gene2"].strip()
        key = pair_key(a, b)
        if key not in allowed.get(s, ()):
            continue
        matched += 1
        hit = parse_hit((r["evalue"], r["bitscore"], r.get("pident", ""),
                         r.get("qcov", ""), r.get("scov", "")), path)
        if key[0] != a:
            hit = Hit(hit.evalue, hit.bitscore, hit.identity, hit.coverage_2, hit.coverage_1)
        if pass_hit(hit, opts) and (key not in collected[s] or collected[s][key].priority < hit.priority):
            collected[s][key] = hit
    LOG.info("Reused DIAMOND hits: %d/%d supplied rows match potential pairs.", matched, supplied)
    return collected


def run_diamond(needed, allowed, folder, runner, opts, archives):
    """Search only proteins in potential A/B events, with both orientations available."""
    collected = {}
    folder.mkdir(parents=True, exist_ok=True)
    if not needed:
        return collected
    binary = executable(opts.diamond)
    for index, (s, genes) in enumerate(sorted(needed.items()), 1):
        if len(genes) < 2:
            continue
        LOG.info("Checking protein similarity for %s (%d/%d; %d candidates)",
                 s, index, len(needed), len(genes))
        sequences_by_id = {name: protein_sequence(archives, s, name, gene.model, gene.raw)
                           for name, gene in genes.items()}
        no_protein = [name for name, seq in sequences_by_id.items() if not seq]
        if no_protein:
            LOG.warning("%s: %d of %d candidate proteins are unavailable or unusable and cannot pass "
                        "DIAMOND validation; first IDs: %s", s, len(no_protein), len(genes),
                        ", ".join(no_protein[:5]))
        work = folder / s
        work.mkdir(exist_ok=True)
        sequences = work / "candidate_proteins.faa"
        names = {}
        with sequences.open("w", encoding="utf-8", newline="\n") as fasta:
            for n, (gene, seq) in enumerate(sorted(sequences_by_id.items())):
                if seq:
                    key = f"p{n}"
                    names[key] = gene
                    fasta.write(f">{key}\n{seq}\n")
        if len(names) < 2:
            continue
        with (work / "sequence_ids.tsv").open("w", encoding="utf-8", newline="") as f:
            tsv_writer(f, ("diamond_id","orthofinder_id")).writerows(names.items())
        raw = work / "raw.tsv"
        runner.run(f"diamond_{s}_db", [binary, "makedb", "--in", sequences, "--db",
                                            "candidates", "--threads", opts.threads], work)
        runner.run(f"diamond_{s}_blastp", [
            binary, "blastp", "--query", sequences, "--db", "candidates",
            "--out", raw, "--outfmt", 6, "qseqid", "sseqid", "evalue", "bitscore", "pident",
            "qcovhsp", "scovhsp", "--threads", opts.threads,
            "--evalue", opts.evalue, "--max-target-seqs", opts.max_targets, "--sensitive",
        ], work)
        data = {}
        with raw.open("r", encoding="utf-8") as f:
            for line in f:
                a, b, *rest = line.rstrip("\r\n").split("\t")
                x, y = names[a], names[b]
                key = pair_key(x, y)
                if key not in allowed[s]:
                    continue
                hit = parse_hit(rest, raw)
                if key[0] != x:
                    hit = Hit(hit.evalue, hit.bitscore, hit.identity, hit.coverage_2, hit.coverage_1)
                if pass_hit(hit, opts) and (key not in data or data[key].priority < hit.priority):
                    data[key] = hit
        collected[s] = data
        LOG.info("%s: %d unique passing DIAMOND pairs", s, len(data))
    return collected


def write_hits(path, hits):
    with open_text(path, "w") as f:
        writer = tsv_writer(f, HIT_FIELDS)
        for s, data in sorted(hits.items()):
            for (a,b), hit in sorted(data.items()):
                writer.writerow((s,a,b,hit.evalue,hit.bitscore,hit.identity,hit.coverage_1,hit.coverage_2))


def potential_lookup(gene_list, opts):
    by_species = defaultdict(list)
    for gene in gene_list:
        if gene.model is not None:
            by_species[gene.species].append(gene)
        elif opts.allow_unmapped_genes and not opts.cross_replicon:
            by_species[gene.species].append(gene)
    return by_species


def compatible_pairs(side_a, side_b, mapping, opts):
    a = potential_lookup(side_a, opts)
    b = potential_lookup(side_b, opts)
    for s in sorted(a.keys() & b.keys()):
        if opts.cross_replicon:
            known_a = [(x, replicon(x,mapping)) for x in a[s]]
            known_b = [(y, replicon(y,mapping)) for y in b[s]]
            for (x,rx),(y,ry) in product(known_a,known_b):
                if rx and ry and rx != ry and x.identifier != y.identifier and x.model.gene != y.model.gene:
                    yield x,y,rx,ry
        else:
            for x,y in product(a[s],b[s]):
                if x.identifier != y.identifier and (
                        x.model is None or y.model is None or x.model.gene != y.model.gene):
                    yield x,y,replicon(x,mapping),replicon(y,mapping)

