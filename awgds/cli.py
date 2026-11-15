"""Ancient Whole-Genome Duplication Search command-line interface."""

import argparse
from collections import Counter
import gzip
import json
import logging
import math
import os
from pathlib import Path
import re
import shutil
import sys
from tempfile import TemporaryDirectory
import time

from Bio.Data import CodonTable
from threadpoolctl import threadpool_limits

from . import __version__
from .dating import load_calibrations, record_diagnostics, run_dating, topology
from .io import (InputError, LOG, Runner, as_json, executable, genomes_and_annotations,
                 open_text, read_table, same_names, tsv_writer)
from .orthology import GeneResolver, OrthoResults, load_gene_map
from .pairs import collect_needed, extract_pairs, prepare_events
from .phylogeny import SpeciesTree, UNIT_TO_MYR, fixed_time_clock, parse_tree_file, posterior_time_clock, tip_names
from .proteins import extract
from .report import RunPaths, read_scores, write_summary
from .similarity import read_chromosome_info, read_replicon_map, run_diamond, trusted_hits, write_hits
from .statistics import assess
from .wgd import run_validation


LOGO = """  ###   #   #   ####  ####    ####
 #   #  #   #  #      #   #  #
 #####  # # #  #  ##  #   #   ###
 #   #  ## ##  #   #  #   #      #
 #   #  #   #   ###   ####   ####"""


def integer(minimum):
    def parse(raw):
        try:
            number = int(raw)
            if number < minimum:
                raise ValueError
            return number
        except ValueError:
            raise argparse.ArgumentTypeError(f"expected an integer >= {minimum}") from None
    return parse


def real(minimum, maximum=None, strict=False):
    def parse(raw):
        try:
            number = float(raw)
            if not math.isfinite(number) or number < minimum or strict and number == minimum:
                raise ValueError
            if maximum is not None and number > maximum:
                raise ValueError
            return number
        except ValueError:
            constraint = f"{'>' if strict else '>='} {minimum}"
            if maximum is not None:
                constraint += f" and <= {maximum}"
            raise argparse.ArgumentTypeError(f"expected a finite number {constraint}") from None
    return parse


def default_threads():
    available = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count() or 1
    return min(available, 8)


def chromosome_arguments(p, standalone=False):
    group = p.add_argument_group("Chromosome input")
    options = group.add_mutually_exclusive_group()
    options.add_argument("--chromosome-info", type=Path, metavar="TXT",
                         help="chromosome whitelist: species (or genome basename), seqid, optional chromosome label")
    options.add_argument("--replicon-map", type=Path, metavar="TSV",
                         help="compatible species/seqid/replicon table; unlisted species use FASTA seqids")
    if standalone:
        options.add_argument("--assume-chromosomes", action="store_true",
                             help="use every FASTA seqid as a chromosome instead of reusing the original whitelist/map")


def wgd_arguments(p):
    wgd = p.add_argument_group("WGD validation of ADI recommendations")
    wgd.add_argument("--wgd-permutations", type=integer(1), default=1999, metavar="N",
                     help="permutations per target/reference comparison (default: 1999)")
    wgd.add_argument("--wgd-alpha", type=real(0, 1, True), default=.05, metavar="P",
                     help="strict cutoff for within-comparison max-adjusted P (default: 0.05)")
    wgd.add_argument("--wgd-min-anchor-families", type=integer(2), default=20, metavar="N",
                     help="minimum anchored orthogroups per target/reference test (default: 20)")
    wgd.add_argument("--wgd-min-link-families", type=real(0, strict=True), default=5., metavar="N",
                     help="minimum sum of orthogroup weights per chromosome link (default: 5)")
    wgd.add_argument("--wgd-min-shared-events", type=integer(1), default=3, metavar="N",
                     help="minimum common, consistently oriented events across daughter sides (default: 3)")
    wgd.add_argument("--wgd-min-shared-chromosomes", type=integer(2), default=2, metavar="N",
                     help="minimum reference chromosomes with matching links across daughter sides (default: 2)")
    wgd.add_argument("--wgd-min-placed-fraction", type=real(0, 1), default=.5, metavar="FRAC",
                     help="minimum fraction of coding loci on assigned chromosomes (default: 0.5)")
    wgd.add_argument("--wgd-max-references", type=integer(0), default=4, metavar="N",
                     help="maximum outgroups per branch (default: 4, 0 means all eligible)")
    wgd.add_argument("--wgd-max-outgroup-depth", type=integer(0), default=0, metavar="N",
                     help="number of successive sister clades to consider (default: 0, all)")
    wgd.add_argument("--wgd-exclude-outgroup", action="append", default=[], metavar="SPECIES",
                     help="exclude a particular species from use as an outgroup; may be repeated")
    wgd.add_argument("--wgd-local", choices=("fallback", "always", "off"), default="fallback",
                     help="when to check shared ordered local chains (default: fallback)")
    wgd.add_argument("--wgd-local-gap", type=integer(0), default=25, metavar="N",
                     help="maximum intervening annotated coding loci in a local chain step (default: 25)")
    wgd.add_argument("--wgd-local-min-chain", type=integer(2), default=5, metavar="N",
                     help="minimum distinct families in each local chain (default: 5)")
    wgd.add_argument("--wgd-matched-draws", type=integer(0), default=0, metavar="N",
                     help="optional family-count-matched resamples for cases shared between branches (default: 0)")
    wgd.add_argument("--wgd-matched-permutations", type=integer(1), default=499, metavar="N",
                     help="permutations per family-count-matched resample (default: 499)")


def parser():
    p = argparse.ArgumentParser(
        prog="awgds", description="Ancient Whole-Genome Duplication Search. Rank branches by duplication enrichment.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Input naming
  Genome: Species_A.softmasked.fna    GFF3: Species_A.longest.gff3.gz
  First dot ends the case-sensitive species key: Species_A. Keys use letters,
  digits, underscores and hyphens, and start with a letter or digit. One genome
  and one GFF3 per key. Subdirectories are searched. FASTA accepts .fna, .fa,
  .fasta, .fas and optional .gz. GFF3 accepts .gff, .gff3 and optional .gz.
  GFF3 column 1 must match the first token of the corresponding FASTA record.
  Reused OrthoFinder tip names must match the complete set of species keys.
  Time-tree tips must match these keys. Extra time-tree tips are pruned by clade.
  Use assemblies and annotations from the same releases as the OrthoFinder proteins.

Examples
  awgds genomes annotations -o awgds_run --threads 24
  awgds genomes annotations -o reuse --orthofinder-results Results_Nov09
  awgds genomes annotations -o dated --orthofinder-results Results_Nov09
        --time-tree ir=clock2/FigTree.tre --time-tree ar=clock3/FigTree.tre --time-unit 100Myr
  awgds genomes annotations -o calibrated --divergence-times divergence_times.txt
  awgds genomes annotations -o scan --chromosome-info chromosomes.txt
  awgds validate previous_run -o revalidation --chromosome-info chromosomes.txt

Open summary.tsv for ADI ranks and WGD evidence together. Gene pairs are
in work/branches/pairs/. See work/trees/ for retained trees, logs/run.log
for the run log and README.md for formats. If no chromosome file is
supplied, all input FASTA sequences are assumed to be chromosomes.
""")
    p.add_argument("genomes", type=Path, metavar="GENOMES", help="genome FASTA directory")
    p.add_argument("annotations", type=Path, metavar="GFFS", help="matching CDS GFF3 directory")
    p.add_argument("-o", "--outdir", type=Path, required=True, metavar="DIR", help="new or empty output directory")
    p.add_argument("-t", "--threads", type=integer(1), default=default_threads(),
                   help="maximum threads for a parallel stage (default: up to 8)")
    p.add_argument("--version", action="version", version=f"AWGDS {__version__}")

    inputs = p.add_argument_group("OrthoFinder and protein extraction")
    inputs.add_argument("--orthofinder-results", type=Path, metavar="DIR", help="reuse a complete Results_* folder")
    inputs.add_argument("--species-tree", type=Path, metavar="TREE",
                        help="supply a rooted species topology when running a fresh OrthoFinder analysis")
    inputs.add_argument("--orthofinder-arg", action="append", default=[], metavar="TOKEN",
                        help="extra OrthoFinder token, repeat as --orthofinder-arg=-FLAG --orthofinder-arg=VALUE")
    inputs.add_argument("--orthofinder", default="orthofinder", metavar="EXE", help="OrthoFinder executable")
    inputs.add_argument("--diamond", default="diamond", metavar="EXE", help="DIAMOND executable")
    inputs.add_argument("--protein-id", choices=("transcript", "protein", "gene"), default="transcript",
                        help="FASTA ID for fresh OrthoFinder proteins (default: transcript)")
    inputs.add_argument("--prefix-protein-ids", action="store_true",
                        help="prefix each fresh protein ID with the species key to disambiguate shared IDs")
    inputs.add_argument("--genetic-code", type=int, choices=sorted(CodonTable.unambiguous_dna_by_id),
                        default=1, metavar="NCBI_TABLE", help="default translation table (default: 1)")
    inputs.add_argument("--retranslate", action="store_true",
                        help="use genome/GFF translations instead of archived OrthoFinder proteins on reuse")
    inputs.add_argument("--min-protein-aa", type=integer(1), default=20, metavar="AA",
                        help="minimum usable translation length (default: 20)")
    inputs.add_argument("--gene-id-map", type=Path, metavar="TSV", help="alias table: species, orthofinder_id, gff_id")
    inputs.add_argument("--allow-unmapped-genes", action="store_true",
                        help="permit missing GFF gene IDs; pairs still must pass all enabled filters")

    flt = p.add_argument_group("Candidate filters and family panel")
    flt.add_argument("--same-og", action=argparse.BooleanOptionalAction, default=True,
                     help="require agreement with the event's OrthoFinder orthogroup (default: on)")
    flt.add_argument("--cross-replicon", action=argparse.BooleanOptionalAction, default=True,
                     help="require two different chromosomes or scaffolds (default: on)")
    flt.add_argument("--require-diamond", action=argparse.BooleanOptionalAction, default=True,
                     help="require protein similarity under the DIAMOND settings (default: on)")
    flt.add_argument("--cross-daughters", action=argparse.BooleanOptionalAction, default=True,
                     help="export only events with a qualifying pair on each daughter side (default: on)")
    flt.add_argument("--min-dup-support", type=real(0, 1), default=0., metavar="FRAC",
                     help="minimum OrthoFinder duplication support, 0 to 1 (default: 0)")
    flt.add_argument("--family-occupancy", type=real(0, 1, strict=True), default=1., metavar="FRAC",
                     help="fraction of sampled species required in an ADI orthogroup (default: 1)")
    flt.add_argument("--diamond-hits", type=Path, metavar="TSV", help="reuse diamond_hits.tsv.gz from a compatible search")
    flt.add_argument("--evalue", type=real(0, strict=True), default=1e-5, metavar="E",
                     help="maximum DIAMOND E-value (default: 1e-5)")
    flt.add_argument("--min-identity", type=real(0, 100), default=0., metavar="PCT",
                     help="minimum DIAMOND percent identity (default: 0)")
    flt.add_argument("--min-coverage", type=real(0, 100), default=0., metavar="PCT",
                     help="minimum DIAMOND coverage of each protein, in percent (default: 0)")
    flt.add_argument("--max-targets", type=integer(0), default=0, metavar="N",
                     help="DIAMOND target cap per query (default: 0, all passing)")
    flt.add_argument("--rbh", action="store_true", help="keep only reciprocal best eligible protein pairs within species")
    flt.add_argument("--keep-search", action="store_true", help="retain native DIAMOND inputs, hits, databases and ID maps")

    dating = p.add_argument_group("Branch lengths and dating")
    dating.set_defaults(time_tree=[], time_unit=[], mcmc_posterior=[])
    dating.add_argument("--time-tree", action="append", metavar="[NAME=]TREE",
                        help="dated Newick, FigTree/NEXUS or MCMCtree report; repeat for multiple trees")
    dating.add_argument("--mcmc-posterior", action="append", nargs=3, metavar=("NAME", "REPORT", "TRACE[,TRACE]"),
                        help="MCMCtree report and one or more raw mcmc.txt traces; repeat for multiple clocks")
    dating.add_argument("--time-unit", action="append", choices=tuple(UNIT_TO_MYR),
                        help="input units: one shared value or one per time tree, then per posterior bundle")
    dating.add_argument("--divergence-times", type=Path, metavar="TXT",
                        help="date MRCA pairs with MCMCtree; TSV columns: species_a, species_b, min_mya/max_mya or age_mya")
    dating.add_argument("--point-age-width", type=real(0, strict=True), default=0.2, metavar="FRAC",
                        help="fraction below and above age_mya for a soft range (default: 0.2; less than 1)")
    dating.add_argument("--root-max-mya", type=real(0, strict=True), metavar="MYR",
                        help="root age upper bound if the calibration table supplies no root range")
    dating.add_argument("--dating-alignment", type=Path, metavar="FILE",
                        help="aligned protein FASTA or relaxed PHYLIP instead of MAFFT single-copy genes")
    dating.add_argument("--dating-max-ogs", type=integer(0), default=0, metavar="N",
                        help="sample at most N complete single-copy groups for dating (default: 0, all)")
    dating.add_argument("--dating-model", default="LG+F+G4", metavar="MODEL", help="IQ-TREE protein model (default: LG+F+G4)")
    dating.add_argument("--dating-clock", choices=("both", "independent", "autocorrelated"), default="both",
                        help="MCMCtree rate model when calibrating (default: both, with two chains per model)")
    dating.add_argument("--mcmc-burnin", type=integer(0), default=100000, metavar="N", help="burn-in per chain (default: 100000)")
    dating.add_argument("--mcmc-samples", type=integer(4), default=18000, metavar="N", help="posterior draws per chain (default: 18000)")
    dating.add_argument("--mcmc-sample-every", type=integer(1), default=20, metavar="N", help="MCMC sampling interval (default: 20)")
    dating.add_argument("--iqtree", default="iqtree3", metavar="EXE", help="IQ-TREE 3 executable with --dating mcmctree")
    dating.add_argument("--mafft", default="mafft", metavar="EXE", help="MAFFT executable for single-copy protein alignments")
    dating.add_argument("--mcmctree", default="mcmctree", metavar="EXE", help="PAML MCMCtree executable")

    stats = p.add_argument_group("ADI calculations")
    stats.add_argument("--bootstrap", type=integer(20), default=12000, metavar="N", help="family bootstrap replicates (default: 12000)")
    stats.add_argument("--seed", type=integer(0), default=410611, metavar="N", help="random seed (default: 410611)")
    p.add_argument("--wgd-validation", action=argparse.BooleanOptionalAction, default=True,
                   help="validate ADI-recommended branches with chromosome evidence (default: on)")
    chromosome_arguments(p)
    wgd_arguments(p)
    return p


def validation_parser():
    p = argparse.ArgumentParser(prog="awgds validate",
        description="Validate recommended branches from a saved AWGDS run without repeating the candidate search.")
    p.add_argument("run", type=Path, metavar="RUN", help="full AWGDS output with summary.tsv and saved branch pairs")
    p.add_argument("-o", "--outdir", type=Path, required=True, metavar="DIR", help="new or empty validation directory")
    p.add_argument("--genomes", type=Path, metavar="DIR", help="override the saved genome paths")
    p.add_argument("--annotations", type=Path, metavar="DIR", help="override the saved GFF paths")
    p.add_argument("--orthofinder-results", type=Path, metavar="DIR", help="override the saved OrthoFinder location")
    p.add_argument("--gene-id-map", type=Path, metavar="TSV", help="override the saved GFF ID map")
    p.add_argument("--seed", type=integer(0), default=None, metavar="N", help="random seed (default: original run's seed)")
    p.add_argument("-t", "--threads", type=integer(1), default=default_threads(), help="parallel comparisons (default: up to 8)")
    chromosome_arguments(p, standalone=True)
    wgd_arguments(p)
    return p


def checked_file(value):
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise InputError(f"Input file does not exist: {path}")
    return path


def external_clocks(opts):
    configs = []
    for value in opts.time_tree:
        label, sep, filename = value.partition("=")
        if not sep or not re.fullmatch(r"[A-Za-z0-9_-]+", label):
            label, filename = None, value
        source = checked_file(filename)
        if label is None:
            base = source.name.removesuffix(".gz").rsplit(".", 1)[0]
            if base in {"FigTree", "mcmctree"}:
                base = source.parent.name
            label = re.sub(r"[^A-Za-z0-9_-]", "_", base)
        configs.append((label, source, []))
    for label, report, traces in opts.mcmc_posterior:
        paths = [checked_file(x) for x in traces.split(",")]
        if len(paths) != len(set(paths)):
            raise InputError(f"{label}: the same MCMCtree trace was given more than once.")
        configs.append((label, checked_file(report), paths))
    labels = [label for label, _, _ in configs]
    if any(not re.fullmatch(r"[A-Za-z0-9_-]+", label) for label in labels) or len(set(labels)) != len(labels):
        raise InputError("Clock labels must be unique and use letters, digits, underscores or hyphens. "
                         "Use NAME=PATH with --time-tree to set a label.")
    has_external = bool(configs)
    if opts.divergence_times and has_external:
        raise InputError("Choose either --divergence-times or supplied time trees/posteriors for this run.")
    count = len(opts.time_unit)
    if has_external and count not in {1, len(configs)}:
        raise InputError("Supply one shared --time-unit or one --time-unit per external tree/posterior.")
    if not has_external and count:
        raise InputError("--time-unit applies to supplied time trees/posteriors. Divergence times are in Myr.")
    return [(*value, opts.time_unit[0 if count == 1 else i]) for i, value in enumerate(configs)]


def check_wgd_settings(opts):
    if opts.wgd_alpha >= 1:
        raise InputError("--wgd-alpha must be less than 1.")
    if 1 / (opts.wgd_permutations + 1) >= opts.wgd_alpha:
        raise InputError("Increase --wgd-permutations so the minimum empirical P is below --wgd-alpha.")
    if opts.wgd_matched_draws and 1 / (opts.wgd_matched_permutations + 1) >= opts.wgd_alpha:
        raise InputError("Increase --wgd-matched-permutations for this --wgd-alpha.")


def check_output(path, inputs):
    for folder in inputs:
        if folder is not None and path.is_relative_to(folder):
            raise InputError(f"The output directory must be outside input directories: {folder}")
    if path.exists() and (not path.is_dir() or any(path.iterdir())):
        raise InputError(f"Output directory must be new or empty: {path}")


def arguments(opts):
    for name in ("genomes", "annotations", "outdir", "orthofinder_results"):
        value = getattr(opts, name)
        if value is not None:
            setattr(opts, name, value.expanduser().resolve())
    for name in ("species_tree", "gene_id_map", "replicon_map", "chromosome_info", "diamond_hits",
                 "divergence_times", "dating_alignment"):
        value = getattr(opts, name)
        if value is not None:
            setattr(opts, name, checked_file(value))
    if opts.species_tree and opts.orthofinder_results:
        raise InputError("--species-tree is for fresh OrthoFinder runs. The reused run already has a species tree.")
    if opts.orthofinder_arg and opts.orthofinder_results:
        raise InputError("--orthofinder-arg applies to fresh OrthoFinder runs.")
    if opts.prefix_protein_ids and opts.orthofinder_results:
        raise InputError("--prefix-protein-ids applies to fresh OrthoFinder runs.")
    if set(opts.orthofinder_arg) & {"-f", "-o", "-s", "-t", "-a", "-S", "-b", "-d", "--assign", "--core"}:
        raise InputError("Genome, result, topology, thread and search arguments are already set by AWGDS.")
    if (opts.dating_alignment or opts.root_max_mya) and not opts.divergence_times:
        raise InputError("--dating-alignment and --root-max-mya require --divergence-times.")
    if opts.point_age_width >= 1:
        raise InputError("--point-age-width must be less than 1.")
    if (opts.rbh or opts.diamond_hits) and not opts.require_diamond:
        raise InputError("--rbh and --diamond-hits require the DIAMOND filter to be enabled.")
    check_output(opts.outdir, (opts.genomes, opts.annotations, opts.orthofinder_results))
    if opts.wgd_validation:
        check_wgd_settings(opts)
    return external_clocks(opts)


def resolve_chromosomes(opts, genomes):
    if opts.chromosome_info:
        mapping = read_chromosome_info(opts.chromosome_info, genomes)
        LOG.info("Using the explicit chromosome whitelist: %s", opts.chromosome_info)
        return mapping, {"mode": "whitelist", "file": opts.chromosome_info}
    if opts.replicon_map:
        return read_replicon_map(opts.replicon_map, genomes), {"mode": "replicon_map", "file": opts.replicon_map}
    if getattr(opts, "wgd_validation", True):
        LOG.info("Chromosome validation assumes every supplied FASTA seqid is a chromosome.")
    return {}, {"mode": "all_input_sequences_assumed_chromosomes"}


def read_clocks(tree, configs, paths, manifest):
    clocks = []
    for name, source, traces, unit in configs:
        LOG.info("Loading clock %s from %s (input unit: %s)", name, source, unit)
        archive = paths.trees / "imported" / f"clock_{name}"
        archive.mkdir(parents=True, exist_ok=True)
        saved_tree = archive / f"tree_{source.name}"
        shutil.copy2(source, saved_tree)
        saved_traces = []
        if traces:
            clock, draws, nodes = posterior_time_clock(tree, source, traces, unit, name)
            for i, trace in enumerate(traces, 1):
                saved = archive / f"trace_{i:02d}.txt.gz"
                if trace.suffix.lower() == ".gz":
                    shutil.copy2(trace, saved)
                else:
                    with trace.open("rb") as src, gzip.open(saved, "wb", compresslevel=5) as dst:
                        shutil.copyfileobj(src, dst)
                saved_traces.append(saved.relative_to(paths.root).as_posix())
            record_diagnostics(archive / "diagnostics.tsv", name, draws, nodes)
        else:
            clock = fixed_time_clock(tree, source, unit, name)
        manifest.setdefault("preserved_clock_inputs", {})[name] = {
            "tree_or_report": saved_tree.relative_to(paths.root).as_posix(), "input_unit": unit,
            "traces": saved_traces,
        }
        clocks.append(clock)
    return clocks


def write_input_files(path, inputs):
    with open_text(path, "w") as f:
        tsv_writer(f, ("species", "genome", "annotation")).writerows((s, *files) for s, files in inputs.items())


def preserve_mapping_files(opts, paths, manifest):
    for attr in ("chromosome_info", "replicon_map", "gene_id_map"):
        source = getattr(opts, attr)
        if source:
            suffix = ".txt" if attr == "chromosome_info" else ".tsv"
            if source.suffix.lower() == ".gz":
                suffix += ".gz"
            saved = paths.inputs / (attr + suffix)
            shutil.copy2(source, saved)
            manifest.setdefault("preserved_inputs", {})[attr] = saved.relative_to(paths.root).as_posix()


def write_metadata(paths, tree, metrics, stages, families):
    with open_text(paths.branches / "branch_counts.tsv", "w") as f:
        output = tsv_writer(f, ("node", "parent", "scored", "left_species", "right_species", "species_pairs",
                               "supported_families", "family_species_pairs", "mean_families_per_pair",
                               "input_events", "low_support", "genes_outside_og", "no_within_species_pair",
                               "no_daughter_coverage", "potential_events", "potential_adi_events", "passing_events",
                               "passing_gene_pairs", "adi_events", "adi_gene_pairs"))
        for b in [*tree.scored, tree.root]:
            data = stages[b]
            parent, left, right, pairs, groups, observations, mean = metrics[b]
            n = b != tree.root
            output.writerow((b, parent, int(n), left, right, pairs,
                             groups if n else "NA", observations if n else "NA", f"{mean:.12g}" if n else "NA",
                             *[data[k] for k in ("input_events", "low_support", "genes_outside_og",
                                "no_within_species_pair", "no_daughter_coverage", "potential_events",
                                "potential_adi_events", "passing_events", "passing_gene_pairs",
                                "adi_events", "adi_gene_pairs")]))
    with open_text(paths.adi / "families.tsv", "w") as f:
        tsv_writer(f, ("orthogroup",)).writerows((og,) for og in families)
    with open_text(paths.branches / "branch_taxa.tsv", "w") as f:
        out = tsv_writer(f, ("node", "parent", "left_species", "right_species"))
        for b, br in tree.branches.items():
            out.writerow((b, br.parent, ",".join(br.left), ",".join(br.right)))


def step(paths, manifest, phase, index, total, message):
    manifest["phase"] = phase
    as_json(paths.manifest, manifest)
    LOG.info("[%d/%d] %s", index, total, message)


def run_wgd(pair_dir, paths, nodes, tree, groups, genomes, aliases, mapping, opts, source, manifest):
    if nodes:
        results = run_validation(pair_dir, paths.evidence, nodes, tree, groups, genomes, aliases, mapping, opts)
    else:
        LOG.info("No branches met the ADI recommendation threshold")
        results = {}
    manifest["wgd"] = {"examined_branches": len(results), "status_counts": dict(Counter(r["status"] for r in results.values())),
                       "chromosome_source": source}
    return results


def workflow(opts, configs, paths, runner, manifest):
    step(paths, manifest, "inputs", 1, 6, "Check input files and analysis settings")
    inputs = genomes_and_annotations(opts.genomes, opts.annotations)
    LOG.info("Matched %d genome/GFF pairs", len(inputs))
    paths.inputs.mkdir(parents=True)
    paths.trees.mkdir(parents=True)
    write_input_files(paths.input_files, inputs)
    preserve_mapping_files(opts, paths, manifest)
    results, bounds, clocks = None, None, None
    ortho, diamond = None, None
    if opts.orthofinder_results:
        results = OrthoResults(opts.orthofinder_results, inputs)
        LOG.info("Reusing OrthoFinder results from %s", results.root)
        if opts.divergence_times:
            bounds = load_calibrations(results.tree, opts.divergence_times, opts.point_age_width, opts.root_max_mya)
        else:
            clocks = read_clocks(results.tree, configs, paths, manifest) or [results.tree.substitution_clock()]
    else:
        ortho = executable(opts.orthofinder)
        diamond = executable(opts.diamond)
        if opts.species_tree:
            supplied = parse_tree_file(opts.species_tree)
            same_names(inputs, tip_names(supplied), "genome/GFF keys", opts.species_tree.name)
            if any(len(c.clades) != 2 for c in supplied.get_nonterminals()):
                raise InputError(f"{opts.species_tree}: supply a rooted, bifurcating species tree.")
            (paths.trees / "orthofinder_input_tree.nwk").write_text(topology(supplied.root, lengths=True) + ";\n",
                                                               encoding="utf-8")
    if opts.require_diamond and not opts.diamond_hits:
        diamond = diamond or executable(opts.diamond)
    if opts.divergence_times:
        executable(opts.iqtree)
        executable(opts.mcmctree)
        if not opts.dating_alignment:
            executable(opts.mafft)
    needs_sequence = (results is None or opts.require_diamond and not opts.diamond_hits or
                      opts.divergence_times and not opts.dating_alignment)
    use_archives = (needs_sequence and results is not None and results.protein_files and not opts.retranslate)
    needs_translation = needs_sequence and not use_archives
    step(paths, manifest, "orthology", 2, 6, "Prepare annotations and orthogroups")
    proteins = paths.inputs / "proteins"
    genomes = extract(inputs, proteins, opts.genetic_code, opts.min_protein_aa, opts.protein_id,
                      sequences=bool(needs_translation), write_proteins=results is None,
                      prefix_ids=opts.prefix_protein_ids)
    if results is None:
        target = paths.orthofinder
        command = [ortho, "-f", proteins, "-o", target, "-S", "diamond", "-t", opts.threads,
                   "-a", opts.threads]
        if opts.species_tree:
            command.extend(["-s", paths.trees / "orthofinder_input_tree.nwk"])
        command.extend(opts.orthofinder_arg)
        try:
            runner.run("orthofinder", command, paths.work, extra_path=[Path(diamond).parent])
        except InputError as exc:
            log_path = runner.log_dir / "orthofinder.log"
            if not opts.species_tree and log_path.is_file():
                with open_text(log_path) as f:
                    if any("No orthogroup is single-copy in enough species" in line for line in f):
                        raise InputError(f"OrthoFinder could not infer a species tree from its single-copy genes. "
                                         f"Supply --species-tree or an alternative OrthoFinder tree method. "
                                         f"See {log_path}.") from exc
            raise
        results = OrthoResults(target, inputs)
        if opts.divergence_times:
            bounds = load_calibrations(results.tree, opts.divergence_times, opts.point_age_width, opts.root_max_mya)
        else:
            clocks = read_clocks(results.tree, configs, paths, manifest) or [results.tree.substitution_clock()]

    tree = results.tree
    archives = results.read_working_proteins(opts.min_protein_aa) if use_archives else None
    manifest["protein_source"] = ("orthofinder_working_directory" if use_archives else "genome_cds_translation"
                                  if needs_translation else "none_required")
    manifest["orthofinder_directory"] = results.root
    if results.root.is_relative_to(paths.root):
        manifest["orthofinder_relative_path"] = results.root.relative_to(paths.root).as_posix()
    manifest["species"] = sorted(tree.species)
    shutil.copy2(tree.path, paths.species_tree)
    aliases = load_gene_map(opts.gene_id_map, genomes)
    replicons, chromosome_source = resolve_chromosomes(opts, genomes)
    step(paths, manifest, "branch_pairs", 3, 6, "Map duplication events and export branch gene pairs")
    resolver = GeneResolver(results, genomes, aliases, opts.allow_unmapped_genes)
    families, events, stages, ignored = prepare_events(results, resolver, opts)
    hits = {}
    if opts.require_diamond:
        paths.similarity.mkdir(parents=True)
        needed, allowed = collect_needed(events, replicons, opts)
        if opts.diamond_hits:
            hits = trusted_hits(opts.diamond_hits, inputs, allowed, opts)
        elif opts.keep_search:
            hits = run_diamond(needed, allowed, paths.similarity / "raw", runner, opts, archives)
        else:
            with TemporaryDirectory(prefix=".diamond-search-", dir=paths.similarity) as scratch:
                hits = run_diamond(needed, allowed, Path(scratch), runner, opts, archives)
        write_hits(paths.similarity / "diamond_hits.tsv.gz", hits)
        if not any(hits.values()):
            LOG.warning("No DIAMOND pairs passed. Inspect the inputs and filters if you expected candidates.")
    paths.adi.mkdir(parents=True)
    matrix, metrics = extract_pairs(events, families, tree, hits, replicons, paths.pairs,
                                    paths.adi / "family_species_pairs.tsv.gz", opts, stages)
    write_metadata(paths, tree, metrics, stages, families)
    if resolver.unmapped:
        with open_text(paths.inputs / "unmapped_genes.tsv", "w") as f:
            tsv_writer(f, ("species", "orthofinder_id")).writerows(sorted(resolver.unmapped))
        LOG.warning("%d genes could not be mapped to CDS models", len(resolver.unmapped))
    step(paths, manifest, "adi", 4, 6, "Prepare branch lengths and calculate ADI")
    if bounds is not None:
        clocks = run_dating(results, genomes, aliases, runner, paths.trees / "mcmctree", opts, bounds, archives)
    LOG.info("Computing ADI on %d internal branches using %d clock(s)", len(tree.scored), len(clocks))
    with threadpool_limits(limits=opts.threads, user_api="blas"):
        threshold, columns, scores = assess(matrix, tree, clocks, opts.bootstrap, opts.seed, paths.adi)
    manifest["terminal_duplication_records"] = ignored["terminal_events"]
    manifest["clocks"] = [{"name": c.name, "length_unit": c.unit, "source": c.source,
                           "joint_posterior_draws": len(c.durations) if c.full_posterior else 0} for c in clocks]
    manifest["adi"] = threshold
    manifest["chromosomes"] = chromosome_source
    validation = None
    if opts.wgd_validation:
        write_summary(paths.summary, columns, scores)
        step(paths, manifest, "wgd_validation", 5, 6, "Validate recommended branches")
        validation = run_wgd(paths.pairs, paths, threshold["recommended"], tree, results.groups,
                             genomes, aliases, replicons, opts, chromosome_source, manifest)
    else:
        step(paths, manifest, "wgd_disabled", 5, 6, "WGD validation disabled for this run")
    step(paths, manifest, "reporting", 6, 6, "Write the unified branch summary")
    write_summary(paths.summary, columns, scores, validation, awaiting="validation_disabled")


def validation_arguments(opts):
    opts.run = opts.run.expanduser().resolve()
    opts.outdir = opts.outdir.expanduser().resolve()
    check_wgd_settings(opts)
    try:
        with open_text(opts.run / "run.json") as f:
            previous = json.load(f)
    except (OSError, ValueError) as exc:
        raise InputError(f"{opts.run}: a valid AWGDS run.json is required.") from exc
    if previous.get("mode") == "validate":
        raise InputError(f"Use the original full run ({previous.get('source_run')}) as the validation input.")
    if previous.get("status") != "completed" and previous.get("phase") not in {"wgd_validation", "reporting"}:
        raise InputError(f"{opts.run}: the ADI phase must finish before validation.")
    old_config = previous.get("configuration", {})
    source_paths = RunPaths(opts.run, legacy=previous.get("output_layout") != 2)

    def saved_input(key):
        relative = previous.get("preserved_inputs", {}).get(key)
        if relative and (opts.run / relative).is_file():
            return opts.run / relative
        return old_config.get(key)

    if bool(opts.genomes) != bool(opts.annotations):
        raise InputError("Use --genomes and --annotations together to override the saved inputs.")
    opts.seed = opts.seed if opts.seed is not None else old_config.get("seed", 410611)
    if not (opts.chromosome_info or opts.replicon_map or opts.assume_chromosomes):
        opts.chromosome_info = saved_input("chromosome_info")
        opts.replicon_map = saved_input("replicon_map") if not opts.chromosome_info else None
    opts.gene_id_map = opts.gene_id_map or saved_input("gene_id_map")
    for name in ("chromosome_info", "replicon_map", "gene_id_map"):
        raw = getattr(opts, name)
        if raw:
            setattr(opts, name, checked_file(raw))
    relative_ortho = previous.get("orthofinder_relative_path")
    retained_ortho = opts.run / relative_ortho if relative_ortho else None
    if retained_ortho and not retained_ortho.is_dir():
        retained_ortho = None
    ortho = (opts.orthofinder_results or retained_ortho or previous.get("orthofinder_directory") or
             old_config.get("orthofinder_results"))
    opts.orthofinder_results = Path(ortho).expanduser().resolve() if ortho else None
    for name in ("genomes", "annotations"):
        path = getattr(opts, name)
        if path:
            setattr(opts, name, path.expanduser().resolve())
    check_output(opts.outdir, (source_paths.pairs, opts.genomes, opts.annotations, opts.orthofinder_results))
    return source_paths, previous


def validate_existing(opts, source, previous, paths, manifest):
    step(paths, manifest, "inputs", 1, 3, "Load the saved branch scores and validation inputs")
    config = previous.get("configuration", {})
    columns, scores = read_scores(source.summary)
    nodes = [x["node"] for x in scores if x["recommended"] == "yes"]
    input_files = {}
    for item in read_table(source.input_files, ("species", "genome", "annotation")):
        s = item["species"]
        if s in input_files:
            raise InputError(f"{source.input_files}: duplicate species {s!r}.")
        input_files[s] = (item["genome"], item["annotation"])
    tree = SpeciesTree(source.species_tree, input_files)
    for node in nodes:
        if node not in tree.scored:
            raise InputError(f"{node!r} in the saved summary does not label a scored branch.")
        checked_file(source.pairs / f"{node}.pairs.tsv.gz")
    genomes, groups, aliases, mapping = {}, {}, {}, {}
    chromosome_source = {"mode": "no_recommended_branches"}
    manifest["source_run"] = source.root
    manifest["adi"] = previous.get("adi") or {"branches": len(scores), "recommended": nodes}
    if nodes:
        if opts.orthofinder_results is None:
            raise InputError("Use --orthofinder-results to find the orthogroup memberships.")
        if opts.genomes:
            inputs = genomes_and_annotations(opts.genomes, opts.annotations)
            same_names(input_files, inputs, "source run", "replacement genome/GFF directories")
        else:
            try:
                inputs = {s: tuple(checked_file(path) for path in paths) for s, paths in input_files.items()}
            except InputError as exc:
                raise InputError(f"{exc} Supply --genomes and --annotations if the inputs have moved.") from exc
        results = OrthoResults(opts.orthofinder_results, inputs, read_working_ids=False)
        if results.tree.descendants != tree.descendants:
            raise InputError("The OrthoFinder species-tree labels or clades differ from the saved run.")
        check_output_inputs = {p.parent for pair in inputs.values() for p in pair}
        for folder in check_output_inputs:
            if opts.outdir.is_relative_to(folder):
                raise InputError(f"The output directory must be outside input directories: {folder}")
        LOG.info("Reusing the ranked branches and candidate pairs in %s", source.root)
        groups = results.groups
        manifest["orthofinder_directory"] = results.root
        paths.inputs.mkdir(parents=True)
        write_input_files(paths.input_files, inputs)
        preserve_mapping_files(opts, paths, manifest)
        genomes = extract(inputs, paths.inputs / "proteins", config.get("genetic_code", 1),
                          config.get("min_protein_aa", 20), config.get("protein_id", "transcript"),
                          sequences=False, write_proteins=False)
        aliases = load_gene_map(opts.gene_id_map, genomes)
        mapping, chromosome_source = resolve_chromosomes(opts, genomes)
    step(paths, manifest, "wgd_validation", 2, 3, "Validate recommended branches")
    validation = run_wgd(source.pairs, paths, nodes, tree, groups, genomes, aliases, mapping, opts,
                         chromosome_source, manifest)
    step(paths, manifest, "reporting", 3, 3, "Write the unified branch summary")
    write_summary(paths.summary, columns, scores, validation)


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    validation_only = bool(argv and argv[0] == "validate")
    opts = validation_parser().parse_args(argv[1:]) if validation_only else parser().parse_args(argv)
    handlers = []
    manifest = None
    started = time.monotonic()
    try:
        configs = validation_arguments(opts) if validation_only else arguments(opts)
        paths = RunPaths(opts.outdir)
        paths.logs.mkdir(parents=True)
        log_file = paths.logs / "run.log"
        fmt = logging.Formatter("%(asctime)s %(levelname)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
        handlers = [logging.StreamHandler(), logging.FileHandler(log_file, encoding="utf-8")]
        LOG.setLevel(logging.INFO)
        LOG.propagate = False
        for handler in handlers:
            handler.stream.write(f"{LOGO}\nAWGDS v{__version__} | Ancient Whole-Genome Duplication Search\n\n")
            handler.flush()
            handler.setFormatter(fmt)
            LOG.addHandler(handler)
        runner = None if validation_only else Runner(paths.logs / "tools")
        manifest = {"software": "AWGDS", "version": __version__, "status": "running",
                    "output_layout": 2, "mode": "validate" if validation_only else "scan",
                    "configuration": vars(opts), "commands": runner.commands if runner else []}
        if validation_only:
            manifest["source_run"] = opts.run
        as_json(paths.manifest, manifest)
        LOG.info("Mode: %s | Threads: %d | Output directory: %s",
                 "validate saved results" if validation_only else "full scan", opts.threads, paths.root)
        if validation_only:
            validate_existing(opts, *configs, paths, manifest)
        else:
            workflow(opts, configs, paths, runner, manifest)
        manifest["status"] = "completed"
        manifest["phase"] = "completed"
        manifest["elapsed_seconds"] = round(time.monotonic() - started, 3)
        as_json(paths.manifest, manifest)
        adi = manifest["adi"]
        LOG.info("Completed in %.1f s | %d ranked branches | %d recommended | %d validated",
                 manifest["elapsed_seconds"], adi["branches"], len(adi["recommended"]),
                 manifest.get("wgd", {}).get("examined_branches", 0))
        LOG.info("Final summary  %s", paths.summary)
        if validation_only:
            LOG.info("Source pairs   %s", configs[0].pairs)
        else:
            LOG.info("Branch pairs   %s", paths.pairs)
            LOG.info("Species trees  %s", paths.trees)
        if paths.evidence.is_dir():
            LOG.info("WGD evidence   %s", paths.evidence)
        LOG.info("Run metadata   %s", paths.manifest)
        LOG.info("Run log        %s", log_file)
        return 0
    except (Exception, KeyboardInterrupt) as exc:
        error = "Interrupted by user." if isinstance(exc, KeyboardInterrupt) else str(exc)
        if handlers:
            if isinstance(exc, (InputError, OSError, KeyboardInterrupt)):
                LOG.error("%s", error)
            else:
                LOG.exception("Analysis failed: %s", error)
        else:
            print(f"awgds: error: {error}", file=sys.stderr)
        if manifest is not None:
            manifest.update(status="interrupted" if isinstance(exc, KeyboardInterrupt) else "failed", error=error)
            try:
                as_json(paths.manifest, manifest)
            except OSError:
                LOG.error("Could not update the run manifest.")
        return 130 if isinstance(exc, KeyboardInterrupt) else 1
    finally:
        for handler in handlers:
            LOG.removeHandler(handler)
            handler.close()
