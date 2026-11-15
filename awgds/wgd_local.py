"""Look for ordered local A/B anchor chains shared across daughter lineages."""

from bisect import bisect_left
from collections import Counter, defaultdict
from dataclasses import dataclass

from .io import LOG, InputError


LOCAL_FIELDS = ("node", "left_target", "right_target", "left_sequences", "right_sequences",
                "left_coding_locus_ranks", "right_coding_locus_ranks", "left_chain_families", "right_chain_families",
                "shared_events", "ordered_same_copy_events", "copy_segments_swapped", "order_reversed", "events")


@dataclass(frozen=True, order=True)
class Anchor:
    x: int
    y: int
    family: str
    event: str
    first_is_a: bool

    @property
    def key(self):
        return self.family, self.event


@dataclass(frozen=True)
class Chain:
    sequences: tuple
    direction: int
    anchors: tuple

    @property
    def events(self):
        return {a.key: (i, a.first_is_a) for i, a in enumerate(self.anchors)}

    def intervals(self):
        x = self.anchors[0].x, self.anchors[-1].x
        y = sorted((self.anchors[0].y, self.anchors[-1].y))
        return x, tuple(y)


def gene_ranks(genome, labels):
    ranks = {}
    for seqid in labels:
        loci = {}
        for model in genome.by_seqid.get(seqid, ()):
            old = loci.get(model.gene)
            value = model.start, model.end
            loci[model.gene] = (min(old[0], value[0]), max(old[1], value[1])) if old else value
        for i, gene in enumerate(sorted(loci, key=lambda g: (*loci[g], g))):
            ranks[seqid, gene] = i + 1
    return ranks


def collect_anchors(node_pairs, genomes, chromosomes):
    all_ranks = {}
    anchors = defaultdict(lambda: defaultdict(set))
    for pair in node_pairs:
        s = pair.species
        if s not in all_ranks:
            all_ranks[s] = gene_ranks(genomes[s], chromosomes.labels[s])
        ranks = all_ranks[s]
        if pair.a not in ranks or pair.b not in ranks:
            raise InputError(f"{s}: a saved pair has a locus absent from the current GFF. "
                             "Use the same annotation as in the source run.")
        a, b = (pair.a[0], ranks[pair.a]), (pair.b[0], ranks[pair.b])
        if a != b:
            first, second = sorted((a, b))
            anchors[s][first[0], second[0]].add(Anchor(first[1], second[1], pair.family, pair.event, first == a))
    return anchors


def chain_on_pair(sequences, values, gap, minimum):
    if len(values) < minimum or len({a.family for a in values}) < minimum:
        return
    xs = [a.x for a in values]
    for direction in (1, -1):
        lengths, previous, families = [], [], []
        for i, a in enumerate(values):
            best, parent = 1, -1
            # Each step can skip at most `gap` intervening annotated loci.
            start = bisect_left(xs, a.x - gap - 1, 0, i)
            for j in range(start, i):
                b = values[j]
                if lengths[j] + 1 <= best:
                    continue
                if a.x > b.x and 0 < direction * (a.y - b.y) <= gap + 1 and a.family not in families[j]:
                    best, parent = lengths[j] + 1, j
            prefix = families[parent] if parent >= 0 else frozenset()
            lengths.append(best)
            previous.append(parent)
            families.append(prefix | {a.family})
        longer_paths = set()
        for i in sorted(range(len(values)), key=lambda j: (-lengths[j], j)):
            if lengths[i] < minimum:
                break
            if i in longer_paths:
                continue
            path = []
            while i >= 0:
                path.append(i)
                i = previous[i]
            chain = Chain(sequences, direction, tuple(values[j] for j in reversed(path)))
            (x1, x2), (y1, y2) = chain.intervals()
            if sequences[0] != sequences[1] or max(x1, y1) > min(x2, y2):
                longer_paths.update(path)
                yield chain


def increasing_events(events):
    """Return the longest increasing subset of ordered (right rank, key) pairs."""
    tails, endpoints, prev = [], [], []
    for i, (rank, key) in enumerate(events):
        k = bisect_left(tails, rank)
        prev.append(endpoints[k - 1] if k else -1)
        if k == len(tails):
            tails.append(rank)
            endpoints.append(i)
        else:
            tails[k] = rank
            endpoints[k] = i
    result = []
    i = endpoints[-1] if endpoints else -1
    while i >= 0:
        result.append(events[i][1])
        i = prev[i]
    return tuple(reversed(result))


def ordered_matches(left, right, minimum):
    lookup = defaultdict(list)
    right_events = [chain.events for chain in right]
    for i, items in enumerate(right_events):
        for key in items:
            lookup[key].append(i)
    best = None
    for a in left:
        source = a.events
        matches = Counter(i for key in source for i in lookup.get(key, ()))
        for i in sorted(matches):
            if matches[i] < minimum or best and matches[i] < best[0][0]:
                continue
            other = right_events[i]
            shared = sorted(source.keys() & other.keys(), key=lambda k: source[k][0])
            for flip in (False, True):
                ranks = [(other[k][0], k) for k in shared if (source[k][1] != other[k][1]) == flip]
                if len(ranks) < minimum:
                    continue
                for reverse in (False, True):
                    ordered = increasing_events([(-rank if reverse else rank, k) for rank, k in ranks])
                    score = len(ordered), len(shared)
                    if len(ordered) >= minimum and (best is None or best[0] < score):
                        best = (score, a, right[i], flip, reverse, ordered)
    return best


def print_intervals(chain):
    return ";".join(f"{start}-{end}" for start, end in chain.intervals())


def write_local_matches(branch, node_pairs, genomes, chromosomes, writer, opts):
    LOG.info("WGD %s: looking for shared ordered anchors on both daughter sides", branch.name)
    anchors = collect_anchors(node_pairs, genomes, chromosomes)
    chains = {s: [c for pair, hits in sorted(by_pair.items())
                  for c in chain_on_pair(pair, sorted(hits), opts.wgd_local_gap, opts.wgd_local_min_chain)]
              for s, by_pair in sorted(anchors.items())}
    matches = []
    for a in branch.left:
        if not chains.get(a):
            continue
        for b in branch.right:
            if not chains.get(b):
                continue
            match = ordered_matches(chains[a], chains[b], opts.wgd_min_shared_events)
            if match is None:
                continue
            (ordered, shared), left, right, flip, reverse, events = match
            writer.writerow((branch.name, a, b, ";".join(left.sequences), ";".join(right.sequences),
                print_intervals(left), print_intervals(right), len(left.anchors), len(right.anchors),
                shared, ordered, int(flip), int(reverse), ";".join(f"{g}:{e}" for g, e in events)))
            matches.append((a, b))
    has_both_sides = all(any(s in anchors for s in side) for side in (branch.left, branch.right))
    status = "local_signal_detected" if matches else "no_local_signal_detected" if has_both_sides else "insufficient_pairs"
    return status, matches
