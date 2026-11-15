"""Compute family bootstrap intervals and the handoff's Abnormal Duplication Index."""

from dataclasses import dataclass

import numpy as np

from .io import InputError, LOG, open_text, tsv_writer


def bootstrap_counts(matrix, replicates, rng):
    n, n_branches = matrix.shape
    draws = np.empty((replicates, n_branches))
    probabilities = np.full(n, 1/n)
    for i in range(0, replicates, 512):
        stop = min(i+512, replicates)
        weights = rng.multinomial(n, probabilities, size=stop-i)
        draws[i:stop] = weights.astype(float) @ matrix
        if i == 0 or stop == replicates or stop % 4096 == 0:
            LOG.info("Family bootstrap: %s/%s replicates", f"{stop:,}", f"{replicates:,}")
    return draws


def q_linear(a, p, axis=None):
    return np.quantile(a, p, axis=axis, method="linear")


def safe_ratio(value, background):
    return value/background if background > 0 else None


def shared_floor(ratios):
    return min(ratios) if ratios and all(x is not None for x in ratios) else None


@dataclass
class Assessment:
    rank: int
    name: str
    adi: float | None
    percent: float | None
    global_ratio: float | None
    local_ratio: float | None
    recommended: bool | None


def assess(matrix, tree, clocks, bootstrap, seed, folder):
    if not clocks or len({c.unit for c in clocks}) != 1:
        raise InputError("At least one clock and one common density unit are required.")
    rng = np.random.default_rng(seed)
    counts = bootstrap_counts(matrix, bootstrap, rng)
    n_families, n_branches = matrix.shape
    names = tree.scored
    b_index = {b:i for i,b in enumerate(names)}
    direct, two = tree.neighbor_sets()
    densities = {}
    for c in clocks:
        if c.durations.ndim != 2 or c.durations.shape[1] != n_branches or \
                not np.isfinite(c.durations).all() or (c.durations <= 0).any():
            raise InputError(f"{c.name}: invalid branch-duration array.")
        if c.full_posterior:
            drawn = c.durations[rng.integers(len(c.durations), size=bootstrap)]
        else:
            drawn = c.durations
        density = (1000/n_families) * counts/drawn
        intervals = q_linear(density, [.025,.5,.975], axis=0)
        durations = q_linear(c.durations,[.025,.5,.975],axis=0)
        densities[c.name] = (intervals, durations)
        LOG.info("%s: calculated branch density intervals from %s.",
                 c.name, "family and joint time draws" if c.full_posterior else "family draws")

    components = {}
    clock_columns = ("node","clock","density_unit","uncertainty","branch_length_median",
                     "branch_length_low95","branch_length_high95","density_median",
                     "density_low95","density_high95","global_background","neighbor_background",
                     "global_ratio","neighbor_ratio","direct_neighbors","two_step_neighbors_used")
    unit = f"families_per_1000_per_{clocks[0].unit}"
    with open_text(folder/"clock_densities.tsv", "w") as f:
        per_clock = tsv_writer(f, clock_columns)
        for b in names:
            i = b_index[b]
            comp_global, comp_local = [], []
            for c in clocks:
                density, time = densities[c.name]
                lo, mid, up = map(float, density[:,i])
                others = np.delete(density[2], i)
                g = float(q_linear(others,.75))
                b_direct = sorted(direct[b], key=b_index.get)
                b_two = sorted(two[b], key=b_index.get) if len(b_direct) == 1 else []
                loc = (float(np.median([density[2,b_index[x]] for x in b_direct])) if b_direct else 0.)
                if b_two:
                    loc = max(loc,float(np.median([density[2,b_index[x]] for x in b_two])))
                rg, rl = safe_ratio(lo,g), safe_ratio(lo,loc)
                comp_global.append(rg)
                comp_local.append(rl)
                floor, median, ceiling = map(float, time[:,i])
                per_clock.writerow((b,c.name,unit,
                                    "families_and_time_posterior" if c.full_posterior else "families_only",
                                    fmt(median),fmt(floor),fmt(ceiling),fmt(mid),fmt(lo),fmt(up),
                                    fmt(g),fmt(loc),fmt(rg),fmt(rl),",".join(b_direct),",".join(b_two)))
            gg, ll = shared_floor(comp_global),shared_floor(comp_local)
            vals = [r for r in (gg,ll) if r is not None]
            components[b] = max(vals) if vals else None, gg, ll

    finite = [vals[0] for vals in components.values() if vals[0] is not None]
    if finite:
        low, high = min(finite), max(finite)
        q1, q3 = map(float, q_linear(finite, [.25, .75]))
    else:
        low = high = q1 = q3 = None
    fence = q3 + 1.5*(q3-q1) if finite else None
    if fence is None:
        LOG.warning("No branch has a finite ADI. No recommendation cutoff can be estimated.")
    ranking = sorted(names, key=lambda b: (components[b][0] is None,
                                            -(components[b][0] or 0), b_index[b]))
    rows, previous, rank = [], object(), None
    for i,b in enumerate(ranking, 1):
        adi,g,l = components[b]
        if adi != previous:
            rank = i
            previous = adi
        percent = percent_value(adi,low,high)
        rows.append(Assessment(rank,b,adi,percent,g,l,
                               adi > fence if adi is not None and fence is not None else None))
    columns = ("rank", "node", "parent", "density_unit", *(f"clock_density_{c.name}" for c in clocks),
               "adi_raw", "adi_percent", "recommended", "global_ratio", "neighbor_ratio")
    scores = []
    for r in rows:
        flag = {True: "yes", False: "no", None: "NA"}[r.recommended]
        medians = [fmt(densities[c.name][0][1,b_index[r.name]]) for c in clocks]
        scores.append(dict(zip(columns, (r.rank if r.adi is not None else "NA", r.name, tree.parent[r.name], unit,
            *medians, fmt(r.adi), fmt(r.percent), flag, fmt(r.global_ratio), fmt(r.local_ratio)))))

    threshold_pct = percent_value(fence,low,high)
    brief = {
        "branches": n_branches, "branches_with_finite_adi": len(finite), "family_count": n_families,
        "clock_names": [c.name for c in clocks], "bootstrap_replicates": bootstrap, "seed": seed,
        "adi_min": low, "adi_max": high, "q1": q1, "q3": q3,
        "iqr": q3-q1 if q1 is not None else None, "upper_fence": fence,
        "upper_fence_percent": threshold_pct,
        "recommended": [r.name for r in rows if r.recommended],
        "rule": "Recommend when unrounded ADI > Q3 + 1.5 * (Q3 - Q1), linear quantiles.",
    }
    LOG.info("ADI upper fence: %s. Recommended branches: %s", fmt(fence),
             ", ".join(brief["recommended"]) or "none")
    return brief, columns, scores


def percent_value(score, low, high):
    if score is None or low is None:
        return None
    if low == high:
        return 0. if score == low else None
    return 100*(score-low)/(high-low)


def fmt(x):
    return "NA" if x is None else format(float(x), ".12g")

