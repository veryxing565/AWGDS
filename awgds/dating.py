"""Calibrate a fixed OrthoFinder topology from a protein alignment with MCMCtree."""

from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
import math
import os
import shutil

import numpy as np
from Bio import AlignIO

from .io import LOG, InputError, executable, open_text, read_table, same_names, tsv_writer
from .orthology import locate_annotation, protein_sequence
from .phylogeny import parse_tree_file, posterior_time_clock, project_tree


@dataclass(frozen=True)
class Calibration:
    lower: float | None
    upper: float
    row: dict
    interval_type: str


def load_calibrations(tree, path, point_width=0.2, root_max=None):
    fields = ("species_a", "species_b")
    data = {}
    for row in read_table(path, fields):
        a, b = row["species_a"].strip(), row["species_b"].strip()
        if a == b or a not in tree.species or b not in tree.species:
            raise InputError(f"{path}: expected two different OrthoFinder tip names, found {a!r}, {b!r}.")
        node = min((n for n in tree.branches if {a, b} <= tree.descendants[n]),
                   key=lambda n: len(tree.descendants[n]))
        try:
            low, high = (row.get(key) or "" for key in ("min_mya", "max_mya"))
            if low.strip() or high.strip():
                lo, hi = float(low), float(high)
                interval = "supplied_range"
            else:
                age = float(row.get("age_mya") or "")
                lo, hi = age * (1 - point_width), age * (1 + point_width)
                interval = "range_from_point_age"
            if not (math.isfinite(lo) and math.isfinite(hi) and 0 < lo < hi):
                raise ValueError
        except ValueError as exc:
            raise InputError(f"{path}: expected 0 < min_mya < max_mya or a positive age_mya for {a}/{b}.") from exc
        if node in data:
            prior = data[node]
            if (lo, hi) != (prior.lower, prior.upper):
                raise InputError(f"{path}: conflicting ranges for {node}. Supply one range per MRCA.")
            continue
        data[node] = Calibration(lo, hi, row, interval)
    if not data:
        raise InputError(f"{path}: the divergence-time table has no calibration rows.")
    if tree.root not in data:
        root = tree.branches[tree.root]
        if root_max is None:
            raise InputError(f"{path}: supply a root-age range, for example the MRCA of "
                             f"{root.left[0]} and {root.right[0]}, or add --root-max-mya.")
        data[tree.root] = Calibration(None, root_max, {"species_a": root.left[0], "species_b": root.right[0],
                                      "source": "--root-max-mya"}, "supplied_upper_bound")
    elif root_max is not None:
        raise InputError("Use a root range in --divergence-times or --root-max-mya, with one source for the root prior.")
    for node, bounds in data.items():
        for ancestor, older in data.items():
            if (bounds.lower is not None and tree.descendants[node] < tree.descendants[ancestor]
                    and bounds.lower >= older.upper):
                raise InputError(f"{path}: the minimum age of {node} is at least the maximum age of {ancestor}.")
    return data


def write_calibrations(path, data):
    with open_text(path, "w") as f:
        table = tsv_writer(f, ("node", "species_a", "species_b", "age_mya", "min_mya", "max_mya",
                               "interval_type", "source"))
        for node, c in data.items():
            table.writerow((node, c.row["species_a"], c.row["species_b"], c.row.get("age_mya", ""),
                            c.lower, c.upper, c.interval_type, c.row.get("source", "")))


def topology(clade, labels=None, lengths=False):
    if clade.is_terminal():
        node = clade.name
    else:
        label = "" if labels is None else labels.get(clade, "")
        node = "(" + ",".join(topology(c, labels, lengths) for c in clade.clades) + ")" + label
    if lengths and clade.branch_length is not None:
        if not math.isfinite(clade.branch_length) or clade.branch_length < 0:
            raise InputError("A supplied species tree has an invalid branch length.")
        node += f":{clade.branch_length:.16g}"
    return node


def read_alignment(path, species):
    with open_text(path) as f:
        first = f.readline()
        f.seek(0)
        style = "fasta" if first.lstrip().startswith(">") else "phylip-relaxed"
        try:
            aln = AlignIO.read(f, style)
        except ValueError as exc:
            raise InputError(f"{path}: expected one aligned protein FASTA or relaxed PHYLIP: {exc}") from exc
    same_names(species, [r.id for r in aln], "OrthoFinder tree tips", path.name)
    if len(aln) != len(species):
        raise InputError(f"{path}: duplicate species in protein alignment.")
    if aln.get_alignment_length() < 2:
        raise InputError(f"{path}: the protein alignment has fewer than two columns.")
    return {r.id: str(r.seq).upper() for r in aln}


def write_fasta(path, sequences):
    with open_text(path, "w") as f:
        for key, sequence in sequences.items():
            f.write(f">{key}\n{sequence}\n")


def single_copy_alignment(results, genomes, aliases, directory, runner, opts, archives):
    species = sorted(genomes)
    models = {}
    for og, group in results.groups.items():
        if any(group.counts.get(s) != 1 or len(group.genes.get(s, ())) != 1 for s in species):
            continue
        seqs = {}
        for s in species:
            identifier = next(iter(group.genes[s]))
            model = None
            if archives is None:
                model = locate_annotation(genomes[s], aliases.get((s, identifier), identifier))
            seq = protein_sequence(archives, s, identifier, model)
            if not seq:
                break
            seqs[s] = seq
        if len(seqs) == len(species):
            models[og] = seqs
    if not models:
        raise InputError("No complete single-copy orthogroups have usable protein sequences. "
                         "Supply an aligned protein FASTA with --dating-alignment.")
    groups = list(models)
    if opts.dating_max_ogs and len(groups) > opts.dating_max_ogs:
        indices = np.random.default_rng(opts.seed).choice(len(groups), opts.dating_max_ogs, replace=False)
        groups = [groups[int(i)] for i in sorted(indices)]
    mafft = executable(opts.mafft)
    work = directory / "alignments"
    work.mkdir()
    LOG.info("Aligning %d/%d usable complete single-copy families with MAFFT.", len(groups), len(models))

    def align(number, og):
        name = f"og{number+1:05d}"
        inp = work / f"{name}.faa"
        out = work / f"{name}.aln.faa"
        write_fasta(inp, models[og])
        runner.run(f"mafft_{name}", [mafft, "--auto", "--thread", "1", inp], work,
                   output=out, quiet=True)
        return out

    with ThreadPoolExecutor(max_workers=min(opts.threads, len(groups))) as pool:
        futures = [pool.submit(align, i, og) for i, og in enumerate(groups)]
        for number, future in enumerate(as_completed(futures), 1):
            future.result()
            if number == 1 or number % 25 == 0 or number == len(groups):
                LOG.info("MAFFT alignments: %d/%d", number, len(groups))

    blocks = {s: [] for s in species}
    with open_text(directory / "alignment_partitions.tsv", "w") as f:
        table = tsv_writer(f, ("orthogroup", "start_1based", "end_1based", "alignment"))
        stop = 0
        for og, future in zip(groups, futures):
            path = future.result()
            alignment = read_alignment(path, species)
            start = stop + 1
            stop += len(alignment[species[0]])
            table.writerow((og, start, stop, path.relative_to(directory)))
            for s in species:
                blocks[s].append(alignment[s])
    LOG.info("Protein supermatrix: %d species, %d genes, %s columns.", len(species), len(groups), f"{stop:,}")
    return {s: "".join(blocks[s]) for s in species}


def write_settings(path, clock, seed, root_max, opts):
    settings = {
        "seed": seed, "seqfile": "dummy.phy", "treefile": "calibrated.nwk", "outfile": "mcmctree.out",
        "mcmcfile": "mcmc.txt", "ndata": 1, "seqtype": 0, "usedata": "2 in.BV", "clock": clock,
        "RootAge": f"<{root_max/100:.12g}", "BDparas": "1 1 0.1 M", "rgene_gamma": "2 10 1",
        "sigma2_gamma": "2 4 1", "finetune": "1: 0.1 0.1 0.1 0.01 0.5 0.1", "print": 1,
        "burnin": opts.mcmc_burnin, "sampfreq": opts.mcmc_sample_every, "nsample": opts.mcmc_samples,
    }
    with open_text(path, "w") as f:
        for key, value in settings.items():
            f.write(f"{key} = {value}\n")


def split_rhat(chains):
    if len(chains) < 2:
        return None
    n = min(len(chain) for chain in chains) // 2
    if n < 2:
        return None
    sample = np.asarray([block for chain in chains for block in (chain[:n], chain[-n:])])
    within = sample.var(axis=1, ddof=1).mean(axis=0)
    between = n * sample.mean(axis=1).var(axis=0, ddof=1)
    variance = ((n-1)*within + between) / n
    return np.sqrt(np.divide(variance, within, out=np.full_like(within, np.nan), where=within > 0))


def record_diagnostics(path, name, chains, nodes):
    rhat = split_rhat(chains)
    with open_text(path, "w") as f:
        table = tsv_writer(f, ("node", "clock", "chains", "draws", "split_rhat"))
        for i, node in enumerate(nodes):
            score = None if rhat is None else float(rhat[i])
            table.writerow((node, name, len(chains), sum(map(len, chains)),
                            f"{score:.6g}" if score is not None and math.isfinite(score) else "NA"))
    if rhat is None or not np.isfinite(rhat).all() or np.max(rhat) > 1.05:
        LOG.warning("%s: some split R-hat values are unavailable or >1.05. Inspect %s before interpreting dates.",
                    name, path)


def link_or_copy(source, destination):
    try:
        os.link(source, destination)
    except OSError:
        shutil.copyfile(source, destination)


def run_dating(results, genomes, aliases, runner, directory, opts, bounds, archives):
    directory.mkdir()
    tree = results.tree
    write_calibrations(directory / "calibrations.tsv", bounds)
    iqtree, mcmctree = executable(opts.iqtree), executable(opts.mcmctree)
    alignment = (read_alignment(opts.dating_alignment, tree.species) if opts.dating_alignment else
                 single_copy_alignment(results, genomes, aliases, directory, runner, opts, archives))
    write_fasta(directory / "proteins.aln.faa", alignment)
    with open_text(directory / "fixed.nwk", "w") as out:
        out.write(topology(tree.tree.root) + ";\n")
    runner.run("iqtree_mcmctree_hessian", [iqtree, "-s", "proteins.aln.faa", "-st", "AA", "-m",
                opts.dating_model, "-te", "fixed.nwk", "--dating", "mcmctree", "-T", opts.threads,
                "--prefix", "hessian"], directory)
    hessian = directory / "hessian.mcmctree.hessian"
    rooted = directory / "hessian.rooted.nwk"
    dummy = directory / "hessian.dummy.phy"
    if not all(x.is_file() and x.stat().st_size for x in (hessian, rooted, dummy)):
        raise InputError("IQ-TREE did not write the Hessian, rooted Newick and dummy PHYLIP. "
                         "Use an IQ-TREE 3 build with --dating mcmctree.")
    rooted_tree = parse_tree_file(rooted)
    clades = project_tree(tree, rooted_tree)
    labels = {clades[n]: f"'B({c.lower/100:.12g},{c.upper/100:.12g},0.025,0.025)'"
              for n, c in bounds.items() if c.lower is not None}
    calibrated = directory / "calibrated.nwk"
    with open_text(calibrated, "w") as f:
        f.write(f"{len(tree.species)} 1\n{topology(rooted_tree.root, labels)};\n")
    numbers = {"independent": [2], "autocorrelated": [3], "both": [2, 3]}[opts.dating_clock]
    random = np.random.default_rng(opts.seed + 1)
    clocks = []
    for number in numbers:
        for chain in (1, 2):
            work = directory / f"clock{number}_chain{chain}"
            work.mkdir()
            for source, name in [(hessian, "in.BV"), (calibrated, "calibrated.nwk"), (dummy, "dummy.phy")]:
                link_or_copy(source, work / name)
            write_settings(work / "mcmctree.ctl", number, int(random.integers(1, 2147483647)),
                           bounds[tree.root].upper, opts)
            runner.run(f"mcmctree_clock{number}_chain{chain}", [mcmctree, "mcmctree.ctl"], work)
        name = {2: "independent", 3: "autocorrelated"}[number]
        chain1, chain2 = (directory / f"clock{number}_chain{x}" for x in (1, 2))
        clock, draws, nodes = posterior_time_clock(tree, chain1 / "mcmctree.out",
                        [chain1 / "mcmc.txt", chain2 / "mcmc.txt"], "100Myr", name)
        record_diagnostics(directory / f"{name}_diagnostics.tsv", name, draws, nodes)
        clocks.append(clock)
    return clocks
