"""Chromosome-pair permutation tests using one total unit per orthogroup."""

from collections import defaultdict
from dataclasses import dataclass
from hashlib import blake2b

import numpy as np


def rng_for(seed, *labels):
    key = "\0".join([str(seed), *map(str, labels)]).encode()
    return np.random.default_rng(int.from_bytes(blake2b(key, digest_size=16).digest(), "little"))


@dataclass(frozen=True)
class Association:
    reference_chromosome: str
    target_pair: tuple
    families: int
    weight: float
    expected: float
    p_max: float

    def passes(self, minimum, alpha):
        return self.weight + 1e-10 >= minimum and self.p_max < alpha


def count_cells(indices, weights, size, observed_ids):
    if size <= max(4096, 6 * len(indices)):
        counts = np.bincount(indices, weights, minlength=size)
        return counts[observed_ids], counts.max()
    cells, membership = np.unique(indices, return_inverse=True)
    counts = np.bincount(membership, weights)
    present = np.searchsorted(cells, observed_ids)
    clipped = np.minimum(present, len(cells) - 1)
    values = np.where(cells[clipped] == observed_ids, counts[clipped], 0.0)
    return values, counts.max()


def linkage_test(families, reference, target, permutations, rng):
    """Test all observed cells. Each alternative pair receives 1/k of its family's weight."""
    if not families:
        return []
    source = sorted({reference[g] for g in families})
    destinations = sorted({pair for g in families for pair in target[g]})
    if len(source) < 2 or len(destinations) < 2:
        return []
    src_index = {c: i for i, c in enumerate(source)}
    dst_index = {p: i for i, p in enumerate(destinations)}
    families = sorted(families)
    labels = np.array([src_index[reference[g]] for g in families], dtype=np.int64)
    n_choices = np.array([len(target[g]) for g in families], dtype=np.int64)
    family_index = np.repeat(np.arange(len(families)), n_choices)
    pair_index = np.array([dst_index[p] for g in families for p in sorted(target[g])], dtype=np.int64)
    weights = 1.0 / n_choices[family_index]

    width = len(destinations)
    size = len(source) * width
    indices = labels[family_index] * width + pair_index
    cells, membership, distinct = np.unique(indices, return_inverse=True, return_counts=True)
    observed = np.bincount(membership, weights)
    totals = np.zeros(len(cells), dtype=float)
    maxima = np.empty(permutations, dtype=float)
    for i in range(permutations):
        shuffled = rng.permutation(labels)[family_index] * width + pair_index
        values, maxima[i] = count_cells(shuffled, weights, size, cells)
        totals += values
    maxima.sort()
    exceedances = permutations - np.searchsorted(maxima, observed - 1e-10, side="left")
    adjusted = (1 + exceedances) / (1 + permutations)
    return [Association(source[int(cell)//width], destinations[int(cell)%width], int(n),
                        float(score), float(total/permutations), float(p))
            for cell, n, score, total, p in zip(cells, distinct, observed, totals, adjusted)]


def family_phases(target, reference, association):
    """Collect identifiable A/B orientations for the supporting events in a cell."""
    result = {}
    for family, options in target.items():
        if reference.get(family) != association.reference_chromosome:
            continue
        events = options.get(association.target_pair, ())
        phases = defaultdict(set)
        for event, orientation in events:
            phases[event].add(orientation)
        identifiable = {event: next(iter(x)) for event, x in phases.items() if len(x) == 1}
        if identifiable:
            result[family] = identifiable
    return result


def phase_agreement(left, right):
    consistent = [[], []]
    shared = 0
    for family in sorted(left.keys() & right.keys()):
        a, b = left[family], right[family]
        events = sorted(a.keys() & b.keys())
        if not events:
            continue
        shared += 1
        directions = {a[e] != b[e] for e in events}
        if len(directions) == 1:
            consistent[int(next(iter(directions)))].append((family, events[0]))
    swapped = len(consistent[1]) > len(consistent[0])
    return shared, swapped, consistent[int(swapped)]
