"""Select outgroups and validate the recommended branches using saved A/B pairs."""

from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from contextlib import ExitStack
from dataclasses import dataclass
from itertools import product

from .io import InputError, LOG, open_text, read_table, tsv_writer
from .orthology import locate_annotation
from .similarity import sequence_replicon
from .statistics import fmt
from .wgd_linkage import family_phases, linkage_test, phase_agreement, rng_for
from .wgd_local import LOCAL_FIELDS, write_local_matches


PAIR_FIELDS = ("node", "orthogroup", "gene_tree_node", "species", "species_side",
               "seqid_a", "gff_gene_a", "seqid_b", "gff_gene_b")
SELECTION_FIELDS = ("node", "role", "species", "daughter_side", "sister_depth", "selected", "reason",
                    "chromosomes_with_genes", "coding_loci_on_chromosomes", "fraction_placed",
                    "qualifying_families", "testable_left_targets", "testable_right_targets")
TEST_FIELDS = ("node", "reference", "target", "daughter_side", "status", "anchored_families",
               "reference_chromosomes", "target_chromosome_pairs", "significant_links")
LINK_FIELDS = ("node", "reference", "target", "daughter_side", "reference_chromosome",
               "target_chromosome_1", "target_chromosome_2", "distinct_families", "family_weight",
               "expected_family_weight", "max_adjusted_p", "significant")
SHARED_FIELDS = ("node", "reference", "reference_chromosome", "left_target", "left_chromosome_1",
                 "left_chromosome_2", "right_target", "right_chromosome_1", "right_chromosome_2",
                 "shared_event_families", "phase_consistent_families", "swap_target_chromosomes",
                 "common_events", "passes")
REFERENCE_FIELDS = ("node", "reference", "sister_depth", "left_targets_tested", "right_targets_tested",
                    "testable_reference_chromosomes", "shared_reference_chromosomes", "shared_chromosome_fraction")
MATCHED_FIELDS = ("target", "reference", "node", "available_families", "matched_families",
                  "resamples", "permutations_per_resample", "resamples_with_a_link")


@dataclass(frozen=True)
class NodePair:
    species: str
    family: str
    event: str
    a: tuple
    b: tuple
    chromosome_a: str
    chromosome_b: str


@dataclass
class TargetComparison:
    species: str
    side: str
    families: set
    reference_chromosomes: set
    target_pairs: set
    status: str


@dataclass
class Reference:
    species: str
    depth: int
    anchors: dict
    comparisons: list

    @property
    def sides(self):
        return Counter(c.side for c in self.comparisons if c.status == "eligible")


class Chromosomes:
    def __init__(self, genomes, mapping):
        self.labels = {}
        self.loci = {}
        self.stats = {}
        for species, genome in genomes.items():
            labels = {s: sequence_replicon(species, s, mapping) for s in genome.replicons}
            self.labels[species] = {s: c for s, c in labels.items() if c}
            counts = {s: len({x.gene for x in models}) for s, models in genome.by_seqid.items()}
            self.loci[species] = counts
            placed = sum(n for s, n in counts.items() if s in self.labels[species])
            distinct = {c for s, c in self.labels[species].items() if counts.get(s, 0)}
            total = sum(counts.values())
            self.stats[species] = (len(distinct), placed, placed / total if total else 0.0)

    def annotation_stats(self, species):
        return self.stats[species]

    def usable(self, species, opts):
        count, _, placed = self.annotation_stats(species)
        if count < 2:
            return "fewer_than_two_annotated_chromosomes"
        return "low_chromosome_coverage" if placed < opts.wgd_min_placed_fraction else "eligible"


def read_node_pairs(path, branch, genomes, chromosomes):
    for row in read_table(path, PAIR_FIELDS):
        s = row["species"]
        if row["node"] != branch.name or s not in branch.left and s not in branch.right:
            raise InputError(f"{path}: pair from {s!r} does not belong to branch {branch.name}.")
        expected_side = "left" if s in branch.left else "right"
        if expected_side != row["species_side"]:
            raise InputError(f"{path}: {s!r} is assigned to the wrong side of branch {branch.name}.")
        seq_a, seq_b = row["seqid_a"], row["seqid_b"]
        if not seq_a or not seq_b or not row["gff_gene_a"] or not row["gff_gene_b"]:
            continue
        for seqid in (seq_a, seq_b):
            if seqid not in genomes[s].replicons:
                raise InputError(f"{path}: {s}/{seqid} is absent from the input FASTA. "
                                 "Use the same genome and GFF releases as the source run.")
        a, b = chromosomes.labels[s].get(seq_a), chromosomes.labels[s].get(seq_b)
        if a and b:
            yield NodePair(s, row["orthogroup"], row["gene_tree_node"],
                           (seq_a, row["gff_gene_a"]), (seq_b, row["gff_gene_b"]), a, b)


def gather_targets(node_pairs, allowed):
    targets = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
    for row in node_pairs:
        a, b = row.chromosome_a, row.chromosome_b
        if row.species in allowed and a != b:
            pair = (a, b) if a < b else (b, a)
            targets[row.species][row.family][pair].add((row.event, a < b))
    return targets


def outgroup_candidates(tree, node):
    """Walk successive sister clades outward from the focal branch."""
    depth, current = 0, node
    while current in tree.parent:
        depth += 1
        parent = tree.parent[current]
        for s in sorted(tree.descendants[parent] - tree.descendants[current]):
            yield depth, s
        current = parent


def single_copy_anchors(species, groups, families, genomes, aliases, chromosomes, cache):
    cached = cache[species]
    for family in families:
        if family in cached:
            continue
        group = groups[family]
        genes = group.genes.get(species, ())
        seqid = ""
        if group.counts.get(species) == 1 and len(genes) == 1:
            identifier = next(iter(genes))
            model = locate_annotation(genomes[species], aliases.get((species, identifier), identifier),
                                      allow_ambiguous_skip=True)
            if model is not None:
                seqid = model.seqid
        cached[family] = chromosomes.labels[species].get(seqid, "")
    return {family: cached[family] for family in families if cached[family]}


def compare_target(species, side, target, anchors, minimum):
    families = target.keys() & anchors.keys()
    ref_chr = {anchors[g] for g in families}
    dest = {pair for g in families for pair in target[g]}
    if len(families) < minimum:
        status = "too_few_shared_families"
    elif len(ref_chr) < 2 or len(dest) < 2:
        status = "too_little_chromosome_variation"
    else:
        status = "eligible"
    return TargetComparison(species, side, families, ref_chr, dest, status)


def select_references(references, maximum):
    if not maximum:
        maximum = len(references)
    # Prefer references that provide coverage in both daughter lineages. Within
    # that set, take one from each successive sister clade before a second there.
    two_sided, partial = defaultdict(list), defaultdict(list)
    for ref in references:
        sides = ref.sides
        if min(sides["left"], sides["right"]) > 0:
            two_sided[ref.depth].append(ref)
        elif sum(sides.values()):
            partial[ref.depth].append(ref)

    def quality(ref):
        counts = ref.sides
        return (-min(counts["left"], counts["right"]), -sum(counts.values()),
                -len(ref.anchors), ref.species)

    ordered = []
    for tier in (two_sided, partial):
        batches = [sorted(tier[depth], key=quality) for depth in sorted(tier)]
        for n in range(max(map(len, batches), default=0)):
            ordered.extend(refs[n] for refs in batches if n < len(refs))
    return ordered[:maximum]


def plan_branch(tree, branch, targets, groups, genomes, aliases, chromosomes, cache, opts, selection):
    all_families = set().union(*(t.keys() for t in targets.values()))
    unknown = all_families - groups.keys()
    if unknown:
        raise InputError(f"{branch.name}: {min(unknown)} in the saved pairs is absent from the orthogroup table.")
    reference_plans = []
    for depth, species in outgroup_candidates(tree, branch.name):
        reason = chromosomes.usable(species, opts)
        if species in opts.wgd_exclude_outgroup:
            reason = "user_excluded"
        if opts.wgd_max_outgroup_depth and depth > opts.wgd_max_outgroup_depth:
            reason = "outside_outgroup_depth_limit"
        if reason != "eligible":
            selection.writerow((branch.name, "reference", species, "", depth, "no", reason,
                                *chr_stats(chromosomes, species), 0, 0, 0))
            continue
        anchors = single_copy_anchors(species, groups, all_families, genomes, aliases, chromosomes, cache)
        comparisons = [compare_target(s, side, targets[s], anchors, opts.wgd_min_anchor_families)
                       for side, taxa in (("left", branch.left), ("right", branch.right)) for s in taxa
                       if s in targets]
        reference_plans.append(Reference(species, depth, anchors, comparisons))
    selected = select_references(reference_plans, opts.wgd_max_references)
    selected_names = {x.species for x in selected}
    for ref in reference_plans:
        counts = ref.sides
        kept = ref.species in selected_names
        reason = "selected" if kept else "reference_limit" if counts else "no_testable_targets"
        selection.writerow((branch.name, "reference", ref.species, "", ref.depth, "yes" if kept else "no", reason,
                            *chr_stats(chromosomes, ref.species), len(ref.anchors), counts["left"], counts["right"]))
    return selected


def chr_stats(chromosomes, species):
    unique, n_loci, fraction = chromosomes.annotation_stats(species)
    return unique, n_loci, fmt(fraction)


def test_comparison(node, ref, c, target, opts):
    if c.status != "eligible":
        return []
    rng = rng_for(opts.seed, "wgd", node, ref.species, c.species)
    return linkage_test(c.families, ref.anchors, target, opts.wgd_permutations, rng)


def shared_links(node, ref, left, right, opts, writer):
    passed = set()
    for chromosome in sorted(left.keys() & right.keys()):
        best = None
        for (sp1, hit1, e1), (sp2, hit2, e2) in product(left[chromosome], right[chromosome]):
            shared, swap, events = phase_agreement(e1, e2)
            if not shared:
                continue
            # Keep the best comparison for this reference chromosome. Nearby
            # species sharing the same chromosome evidence only contribute once.
            key = (len(events), shared, -max(hit1.p_max, hit2.p_max))
            if best is None or best[0] < key:
                best = (key, sp1, hit1, sp2, hit2, shared, swap, events)
        if best is None:
            continue
        _, sp1, hit1, sp2, hit2, shared, swap, events = best
        supported = len(events) >= opts.wgd_min_shared_events
        if supported:
            passed.add(chromosome)
        writer.writerow((node, ref, chromosome, sp1, *hit1.target_pair, sp2, *hit2.target_pair,
                         shared, len(events), int(swap), ";".join(f"{g}:{e}" for g, e in events), int(supported)))
    return passed


def validate_branch(node, references, targets, opts, writers, pool, match_contexts):
    tested, positive, coverage = 0, [], []
    per_reference = []
    for rank, ref in enumerate(references):
        comparable = set(), set()
        by_side = defaultdict(list), defaultdict(list)
        comparisons = [(ref, c) for c in ref.comparisons]
        LOG.info("WGD %s: testing %s against %d eligible target(s), sister depth %d", node, ref.species,
                 sum(c.status == "eligible" for _, c in comparisons), ref.depth)
        futures = [pool.submit(test_comparison, node, r, c, targets[c.species], opts) for r, c in comparisons]
        for (_, c), future in zip(comparisons, futures):
            results = future.result()
            eligible = c.status == "eligible"
            index = int(c.side == "right")
            if eligible:
                sizes = Counter(ref.anchors[og] for og in c.families)
                comparable[index].update(chrom for chrom, n in sizes.items() if n >= opts.wgd_min_link_families)
                tested += 1
                if opts.wgd_matched_draws:
                    match_contexts[c.species, ref.species].append(
                        (node, c.families, ref.anchors, targets[c.species]))
            n_pass = 0
            for hit in results:
                passes = hit.passes(opts.wgd_min_link_families, opts.wgd_alpha)
                writers["links"].writerow((node, ref.species, c.species, c.side, hit.reference_chromosome,
                    *hit.target_pair, hit.families, fmt(hit.weight), fmt(hit.expected), fmt(hit.p_max), int(passes)))
                if passes:
                    n_pass += 1
                    phases = family_phases(targets[c.species], ref.anchors, hit)
                    by_side[index][hit.reference_chromosome].append((c.species, hit, phases))
            if n_pass:
                positive.append((c.species, ref.species))
            writers["tests"].writerow((node, ref.species, c.species, c.side,
                "tested" if eligible else c.status, len(c.families), len(c.reference_chromosomes),
                len(c.target_pairs), n_pass))
        found = shared_links(node, ref.species, *by_side, opts, writers["shared"])
        denominator = comparable[0] & comparable[1]
        counts = ref.sides
        ratio = fmt(len(found) / len(denominator) if denominator else None)
        writers["references"].writerow((node, ref.species, ref.depth, counts["left"], counts["right"],
            ",".join(sorted(denominator)), ",".join(sorted(found)), ratio))
        per_reference.append(f"{ref.species}={len(found)}/{len(denominator)}" if denominator else f"{ref.species}=NA")
        if denominator:
            coverage.append((len(found), len(denominator), -rank, ref.species))
    comparable = [x for x in coverage if x[1] >= opts.wgd_min_shared_chromosomes]
    best = max(comparable or coverage, key=lambda x: (x[0], x[0] / x[1], x[2]), default=None)
    status = "insufficient_data"
    if comparable:
        status = "shared_macrosynteny" if best[0] >= opts.wgd_min_shared_chromosomes else "no_broad_linkage_detected"
    n_shared, n_testable, _, best_ref = best or (0, 0, 0, "")
    return {
        "status": status, "tested_target_reference_combinations": tested,
        "positive_target_reference_combinations": len(positive), "positive_pairs": format_pairs(positive),
        "reference_coverage": ";".join(per_reference), "best_reference": best_ref,
        "shared_reference_chromosomes": n_shared, "testable_reference_chromosomes": n_testable,
        "shared_chromosome_fraction": fmt(n_shared / n_testable if n_testable else None),
    }


def format_pairs(pairs):
    return ";".join(f"{a}|{b}" for a, b in pairs)


def matched_resampling(contexts, opts, writer):
    common = [(combo, inputs) for combo, inputs in sorted(contexts.items()) if len(inputs) > 1]
    for i, ((target, reference), inputs) in enumerate(common, 1):
        LOG.info("WGD matched families (%d/%d): %s vs %s", i, len(common), target, reference)
        size = min(len(x[1]) for x in inputs)
        for node, families, ref, candidate in inputs:
            generator = rng_for(opts.seed, "wgd-matched", node, target, reference)
            population = sorted(families)
            with_link = 0
            for _ in range(opts.wgd_matched_draws):
                sample = generator.choice(population, size, replace=False).tolist()
                associations = linkage_test(sample, ref, candidate, opts.wgd_matched_permutations, generator)
                with_link += any(x.passes(opts.wgd_min_link_families, opts.wgd_alpha) for x in associations)
            writer.writerow((target, reference, node, len(families), size, opts.wgd_matched_draws,
                             opts.wgd_matched_permutations, with_link))


def run_validation(pair_dir, output, nodes, tree, groups, genomes, aliases, mapping, opts):
    invalid = set(nodes) - set(tree.scored)
    if invalid:
        raise InputError(f"The saved summary recommends unknown or unscored nodes: {', '.join(sorted(invalid))}.")
    unknown = set(opts.wgd_exclude_outgroup) - tree.species
    if unknown:
        raise InputError(f"Unknown excluded outgroup species: {', '.join(sorted(unknown))}.")
    output.mkdir(parents=True, exist_ok=True)
    chromosomes = Chromosomes(genomes, mapping)
    reports = {}
    cache, match_contexts = defaultdict(dict), defaultdict(list)
    with ExitStack() as stack:
        tables = {
            "selection": ("selection.tsv", SELECTION_FIELDS),
            "tests": ("target_reference.tsv", TEST_FIELDS), "links": ("chromosome_links.tsv.gz", LINK_FIELDS),
            "shared": ("shared_linkage.tsv", SHARED_FIELDS), "matched": ("matched_families.tsv", MATCHED_FIELDS),
            "references": ("reference_summary.tsv", REFERENCE_FIELDS),
        }
        writers = {k: tsv_writer(stack.enter_context(open_text(output / path, "w")), fields)
                   for k, (path, fields) in tables.items() if k != "matched" or opts.wgd_matched_draws}
        with open_text(output / "chromosomes.tsv", "w") as stream:
            writer = tsv_writer(stream, ("species", "seqid", "chromosome", "coding_loci"))
            for species, labels in sorted(chromosomes.labels.items()):
                writer.writerows((species, seqid, chrom, chromosomes.loci[species].get(seqid, 0))
                                 for seqid, chrom in sorted(labels.items()))
        pool = stack.enter_context(ThreadPoolExecutor(max_workers=opts.threads))
        for i, node in enumerate(nodes, 1):
            LOG.info("Validating recommended branch %s (%d/%d)", node, i, len(nodes))
            branch = tree.branches[node]
            path = pair_dir / f"{node}.pairs.tsv.gz"
            node_pairs = lambda: read_node_pairs(path, branch, genomes, chromosomes)
            quality = {s: chromosomes.usable(s, opts) for s in (*branch.left, *branch.right)}
            usable = {s for s, reason in quality.items() if reason == "eligible"}
            targets = gather_targets(node_pairs(), usable)
            side_counts = Counter()
            for side, taxa in (("left", branch.left), ("right", branch.right)):
                for species in taxa:
                    n_families = len(targets.get(species, ()))
                    reason = quality[species]
                    if reason == "eligible" and n_families < opts.wgd_min_anchor_families:
                        reason = "too_few_cross_chromosome_families"
                    keep = reason == "eligible"
                    side_counts[side] += keep
                    if not keep:
                        targets.pop(species, None)
                    writers["selection"].writerow((node, "target", species, side, "", "yes" if keep else "no",
                        reason, *chr_stats(chromosomes, species), n_families, "", ""))
            references = plan_branch(tree, branch, targets, groups, genomes, aliases, chromosomes, cache, opts,
                                     writers["selection"])
            report = validate_branch(node, references, targets, opts, writers, pool, match_contexts)
            status = report["status"]
            local_run = opts.wgd_local == "always" or opts.wgd_local == "fallback" and status != "shared_macrosynteny"
            report.update(local_check="skipped", local_matches="NA", local_target_pairs="NA")
            if local_run:
                if "local" not in writers:
                    stream = stack.enter_context(open_text(output / "local_chains.tsv", "w"))
                    writers["local"] = tsv_writer(stream, LOCAL_FIELDS)
                local_check, matches = write_local_matches(branch, node_pairs(), genomes, chromosomes,
                                                           writers["local"], opts)
                report.update(local_check=local_check, local_matches=len(matches), local_target_pairs=format_pairs(matches))
            report.update(left_eligible_targets=side_counts["left"], right_eligible_targets=side_counts["right"],
                          selected_references=",".join(r.species for r in references),
                          two_sided_references=sum(len(r.sides) == 2 for r in references))
            reports[node] = report
            LOG.info("WGD %s: %s; %d tested target/reference combinations; local ordered matches: %s",
                     node, status, report["tested_target_reference_combinations"], report["local_matches"])
        if opts.wgd_matched_draws:
            matched_resampling(match_contexts, opts, writers["matched"])
    return reports
