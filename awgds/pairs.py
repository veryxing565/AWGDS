"""Assign gene pairs to OrthoFinder events and count supported family/species pairs."""

from collections import Counter, defaultdict
from contextlib import ExitStack
from dataclasses import dataclass
from itertools import product
import math

import numpy as np

from .io import LOG, InputError, ids, open_text, read_table, tsv_writer
from .similarity import compatible_pairs, pair_key

DUP_FIELDS = ("Orthogroup", "Species Tree Node", "Gene Tree Node", "Support", "Genes 1", "Genes 2")
PAIR_FIELDS = (
    "node", "parent", "orthogroup", "gene_tree_node", "of_support", "of_type", "species",
    "species_side", "used_in_adi", "gene_a", "gene_b", "of_id_a", "of_id_b", "gff_gene_a",
    "gff_gene_b", "gff_transcript_a", "gff_transcript_b", "seqid_a", "start_a", "end_a",
    "strand_a", "seqid_b", "start_b", "end_b", "strand_b", "replicon_a", "replicon_b",
    "diamond_evalue", "diamond_bitscore", "diamond_identity_pct", "diamond_coverage_a_pct",
    "diamond_coverage_b_pct",
)


@dataclass
class Event:
    node: str
    group: str
    split: str
    support: float
    kind: str
    a: list
    b: list


def prepare_events(results, resolver, opts):
    tree = results.tree
    minimum = max(1, math.ceil(len(tree.species)*opts.family_occupancy - 1e-10))
    families = [k for k,g in results.groups.items() if
                sum(n>0 for n in g.counts.values()) >= minimum]
    if not families:
        raise InputError(f"No orthogroups are present in {minimum} of {len(tree.species)} species. "
                         "Check the input data or lower --family-occupancy.")
    LOG.info("Using %d orthogroups present in at least %d/%d species for ADI.",
             len(families), minimum, len(tree.species))
    allowed = set(families)
    by_branch = {b: [] for b in tree.branches}
    ignored, stages, seen = Counter(), defaultdict(Counter), set()
    for row_n, r in enumerate(read_table(results.duplications, DUP_FIELDS), 1):
        if row_n % 100000 == 0:
            LOG.info("Read %s OrthoFinder duplication records", f"{row_n:,}")
        b = r["Species Tree Node"].strip()
        if b not in by_branch:
            if b not in tree.species:
                raise InputError(f"{results.duplications}: unknown species-tree node {b!r}.")
            ignored["terminal_events"] += 1
            continue
        g, s = r["Orthogroup"].strip(), r["Gene Tree Node"].strip()
        if not g or not s:
            raise InputError(f"{results.duplications}: empty family or gene-tree event ID.")
        if g not in results.groups:
            raise InputError(f"{results.duplications}: orthogroup {g!r} is absent from the group tables.")
        if (g, s) in seen:
            raise InputError(f"{results.duplications}: duplicate event key {g}/{s}.")
        seen.add((g, s))
        stage = stages[b]
        stage["input_events"] += 1
        try:
            support = float(r["Support"])
            if not math.isfinite(support) or not 0. <= support <= 1.:
                raise ValueError
        except ValueError as exc:
            raise InputError(f"{results.duplications}: invalid support value for {g}/{s}.") from exc
        if support < opts.min_dup_support:
            stage["low_support"] += 1
            continue
        desc = tree.descendants[b]
        a, c = ([resolver.resolve(x, g, desc) for x in ids(r[k])] for k in ("Genes 1", "Genes 2"))
        if set(gene.key for gene in a) & set(gene.key for gene in c):
            raise InputError(f"OrthoFinder event {g}/{s} assigns the same gene to both copies.")
        a, c = (list({gene.key: gene for gene in side}.values()) for side in (a,c))
        if opts.same_og:
            stage["genes_outside_og"] += sum(not resolver.in_group(x, g) for x in a+c)
            a, c = ([x for x in side if resolver.in_group(x, g)] for side in (a,c))
        both = {x.species for x in a} & {x.species for x in c}
        if not both:
            stage["no_within_species_pair"] += 1
            continue
        br = tree.branches[b]
        if opts.cross_daughters and not (both.intersection(br.left) and both.intersection(br.right)):
            stage["no_daughter_coverage"] += 1
            continue
        stage["potential_events"] += 1
        if g in allowed:
            stage["potential_adi_events"] += 1
        by_branch[b].append(Event(b, g, s, support, r.get("Type", ""), a, c))
    return families, by_branch, stages, ignored


def collect_needed(events, replicons, opts):
    needed = defaultdict(dict)
    allowed = defaultdict(set)
    for group in events.values():
        for event in group:
            for x,y,_,_ in compatible_pairs(event.a, event.b, replicons, opts):
                needed[x.species][x.identifier] = x
                needed[y.species][y.identifier] = y
                allowed[x.species].add(pair_key(x.identifier, y.identifier))
    return needed, allowed


def locate(model):
    if model is None:
        return ("", "", "", "")
    return model.seqid, model.start, model.end, model.strand


def hit_fields(hits, x, y):
    label = pair_key(x.identifier,y.identifier)
    hit = hits.get(x.species, {}).get(label)
    if hit is None:
        return "", "", "", "", ""
    if label[0] == x.identifier:
        q,s = hit.coverage_1,hit.coverage_2
    else:
        q,s = hit.coverage_2,hit.coverage_1
    return hit.evalue, hit.bitscore, hit.identity, q, s


def best_partners(hits):
    # Hits have already been limited to eligible pairs in the input event set.
    best = {}
    for species, pairs in hits.items():
        partners = {}
        for (a, b), hit in pairs.items():
            for query, partner in ((a, b), (b, a)):
                previous = partners.get(query)
                choice = (-hit.bitscore, hit.evalue, partner)
                if previous is None or choice < previous:
                    partners[query] = choice
        best[species] = partners
    return best


def extract_pairs(events, families, tree, hits, replicons, folder, observations_file, opts, stages):
    folder.mkdir(parents=True, exist_ok=True)
    family_index = {g:i for i,g in enumerate(families)}
    branch_index = {b:i for i,b in enumerate(tree.scored)}
    contributes = np.zeros((len(families), len(tree.scored)), dtype=np.float64)
    metric_rows = {}
    with ExitStack() as stack:
        event_rows = tsv_writer(stack.enter_context(open_text(folder.parent/"events.tsv.gz", "w")),
            ("node","parent","orthogroup","gene_tree_node","support","type","left_species_with_pairs",
             "right_species_with_pairs","n_gene_pairs","eligible_species_pairs","used_in_adi"))
        observations = tsv_writer(stack.enter_context(open_text(observations_file, "w")),
            ("node","orthogroup","left_species","right_species"))
        best = best_partners(hits) if opts.rbh else None
        for step,b in enumerate([*tree.scored,tree.root], 1):
            br = tree.branches[b]
            pairs = list(product(br.left, br.right))
            bit = {p: 1 << i for i,p in enumerate(pairs)}
            support_masks = defaultdict(int)
            fpath = folder/f"{b}.pairs.tsv.gz"
            with open_text(fpath, "w") as f:
                writer = tsv_writer(f, PAIR_FIELDS)
                for ei, event in enumerate(events[b], 1):
                    possible = compatible_pairs(event.a, event.b, replicons, opts)
                    selected = []
                    for x,y,rx,ry in possible:
                        key = pair_key(x.identifier,y.identifier)
                        if (opts.require_diamond or opts.rbh) and key not in hits.get(x.species, {}):
                            continue
                        if best is not None:
                            d = best[x.species]
                            if d[x.identifier][-1] != y.identifier or d[y.identifier][-1] != x.identifier:
                                continue
                        selected.append((x,y,rx,ry))
                    both = {x.species for x,_,_,_ in selected}
                    left = tuple(s for s in br.left if s in both)
                    right = tuple(s for s in br.right if s in both)
                    if not selected or opts.cross_daughters and not (left and right):
                        continue
                    used = b in branch_index and event.group in family_index and bool(left and right)
                    stage = stages[b]
                    stage["passing_events"] += 1
                    stage["passing_gene_pairs"] += len(selected)
                    if used:
                        stage["adi_events"] += 1
                        stage["adi_gene_pairs"] += len(selected)
                        support_masks[event.group] |= sum(bit[p] for p in product(left,right))
                    # Prefer reporting all descendants of the qualifying event. Mark its ADI status.
                    for x,y,rx,ry in selected:
                        writer.writerow((b, br.parent, event.group, event.split, event.support, event.kind,
                                         x.species, "left" if x.species in br.left else "right", int(used),
                                         x.identifier,y.identifier,x.raw,y.raw,
                                         x.model.gene if x.model else "",y.model.gene if y.model else "",
                                         x.model.transcript if x.model else "",y.model.transcript if y.model else "",
                                         *locate(x.model), *locate(y.model), rx,ry,*hit_fields(hits,x,y)))
                    event_rows.writerow((b,br.parent,event.group,event.split,event.support,event.kind,
                                         ",".join(left),",".join(right),len(selected),len(left)*len(right),int(used)))
                    if ei % 2000 == 0:
                        LOG.info("Branch %s: processed %d/%d duplication events", b, ei, len(events[b]))
            n_records = sum(mask.bit_count() for mask in support_masks.values())
            mean = n_records / len(pairs) if b in branch_index else 0.
            if b in branch_index:
                for g, mask in sorted(support_masks.items()):
                    contributes[family_index[g],branch_index[b]] = mask.bit_count()/len(pairs)
                    while mask:
                        low = mask & -mask
                        u,v = pairs[low.bit_length()-1]
                        observations.writerow((b,g,u,v))
                        mask ^= low
            metric_rows[b] = (br.parent, len(br.left),len(br.right),len(pairs),len(support_masks),n_records,mean)
            LOG.info("Branch %s (%d/%d): %d passing pairs; %.3f ADI families per species pair",
                     b,step,len(tree.branches),stages[b]["passing_gene_pairs"],mean)
    return contributes, metric_rows

