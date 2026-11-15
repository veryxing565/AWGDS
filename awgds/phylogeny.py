"""Read rooted species trees and map external dating results by clade."""

from dataclasses import dataclass
from io import StringIO
import math
import re

import numpy as np
from Bio import Phylo

from .io import InputError, LOG, SPECIES_NAME, open_text, same_names

UNIT_TO_MYR = {"years": 1e-6, "kyr": 1e-3, "Myr": 1.0, "100Myr": 100.0, "Gyr": 1000.0}


@dataclass(frozen=True)
class Branch:
    name: str
    parent: str
    left: tuple
    right: tuple


@dataclass
class Clock:
    name: str
    durations: np.ndarray
    unit: str
    source: str

    @property
    def full_posterior(self):
        return len(self.durations) > 1


def read_newick(raw):
    try:
        return Phylo.read(StringIO(raw), "newick")
    except (ValueError, TypeError, IndexError) as exc:
        raise InputError(f"Invalid Newick tree: {exc}") from exc


def match_text(path):
    with open_text(path) as f:
        return f.read()


def newick_strings(s):
    # PAML writes ordinary Newick, NEXUS/UTREE lines, or several trees in a report.
    return re.findall(r"\([^;]*;", s, re.DOTALL)


def parse_tree_file(path):
    raw = match_text(path).strip()
    if raw.upper().startswith("#NEXUS"):
        if re.search(r"\bTRANSLATE\b", raw, flags=re.I):
            try:
                return Phylo.read(StringIO(raw), "nexus")
            except (ValueError, TypeError, IndexError) as exc:
                raise InputError(f"{path}: cannot parse the NEXUS tree: {exc}") from exc
        lines = re.findall(r"^\s*(?:UTREE|TREE)\s+[^=\n]+\s*=\s*(.*?;)", raw,
                           re.MULTILINE | re.IGNORECASE | re.DOTALL)
        if len(lines) != 1:
            raise InputError(f"{path}: expected one NEXUS tree. Supply separate --time-tree files.")
        raw = lines[0]
    elif "Species tree for FigTree" in raw:
        return paml_report(path)[1]
    else:
        raw = re.sub(r"^\s*\d+\s+\d+\s*(?:\r?\n)", "", raw, count=1)
    return read_newick(raw)


def node_name(clade):
    if clade.name is not None:
        return str(clade.name).strip("'\"")
    if clade.confidence is not None and float(clade.confidence).is_integer():
        return str(int(clade.confidence))
    return ""


def childsets(tree, taxa=None):
    selected = None if taxa is None else set(taxa)
    sets = {}
    for c in tree.find_clades(order="postorder"):
        if c.is_terminal():
            sets[c] = frozenset([c.name]) if selected is None or c.name in selected else frozenset()
        else:
            sets[c] = frozenset().union(*(sets[x] for x in c.clades))
    return sets


def tip_names(tree):
    names = [c.name for c in tree.get_terminals()]
    if any(not name for name in names) or len(names) != len(set(names)):
        raise InputError("The tree contains unnamed or duplicate tip labels.")
    return set(names)


class SpeciesTree:
    def __init__(self, filename, species):
        self.tree = parse_tree_file(filename)
        self.path = filename
        self.species = set(species)
        same_names(self.species, tip_names(self.tree), "genome/GFF keys", "OrthoFinder tree tips")
        desc = childsets(self.tree)
        self.names = {c: node_name(c) for c in desc}
        if len(set(self.names.values())) != len(self.names):
            raise InputError(f"{filename}: tree nodes must all have unique labels.")
        for c in self.tree.get_nonterminals():
            if len(c.clades) != 2:
                raise InputError(f"{filename}: {node_name(c) or 'root'} is not bifurcating.")
            if not SPECIES_NAME.fullmatch(self.names[c]):
                raise InputError(f"{filename}: missing or unsafe internal node label {self.names[c]!r}. "
                                 "Use OrthoFinder's SpeciesTree_rooted_node_labels.txt.")
        self.root = self.names[self.tree.root]
        self.parent = {self.names[c]: self.names[p] for p in desc for c in p.clades}
        self.clades = {self.names[c]: c for c in desc}
        self.descendants = {self.names[c]: tips for c, tips in desc.items()}
        self.branches = {
            self.names[c]: Branch(self.names[c], self.parent.get(self.names[c], ""),
                                  tuple(sorted(desc[c.clades[0]])), tuple(sorted(desc[c.clades[1]])))
            for c in self.tree.get_nonterminals()
        }
        self.scored = sorted((b for b in self.branches if b != self.root), key=natural)
        if len(self.scored) < 2:
            raise InputError("At least two non-root internal branches are required.")

    def neighbor_sets(self):
        direct = {b: set() for b in self.scored}
        endpoints = {}
        for b in self.scored:
            endpoints.setdefault(self.parent[b], []).append(b)
            endpoints.setdefault(b, []).append(b)
        for connected in endpoints.values():
            for b in connected:
                direct[b].update(x for x in connected if x != b)
        two = {b: (n | set().union(*(direct[i] for i in n))) - {b} for b, n in direct.items()}
        return direct, two

    def substitution_clock(self):
        lengths = [self.clades[b].branch_length for b in self.scored]
        positive_lengths(self.scored, lengths)
        return Clock("orthofinder", np.asarray([lengths], dtype=float),
                     "orthofinder_branch_length", str(self.path))


def natural(name):
    return tuple((1, int(part)) if part.isdigit() else (0, part.casefold())
                 for part in re.split(r"(\d+)", name))


def positive_lengths(names, lengths):
    problem = [b for b, t in zip(names, lengths)
               if t is None or not math.isfinite(t) or t <= 0]
    if problem:
        raise InputError(f"Missing or non-positive lengths on internal branches: {', '.join(problem[:12])}. "
                         "Supply a rooted compatible time tree with positive lengths.")


def project_tree(ref, supplied):
    names = tip_names(supplied)
    if not ref.species <= names:
        missing = sorted(ref.species - names)
        raise InputError("Time tree is missing species: " + ", ".join(missing[:12]) +
                         ". Use the exact genome/GFF species keys.")
    desc = childsets(supplied, ref.species)
    retained = {s: c for c, s in desc.items()
                if s and (c.is_terminal() or sum(bool(desc[x]) for x in c.clades) >= 2)}
    expected = set(ref.descendants.values())
    if set(retained) != expected:
        missing = sorted((x for x in expected if x not in retained), key=lambda x: (len(x), sorted(x)))
        example = ", ".join(sorted(missing[0])[:8]) if missing else "additional internal clades"
        raise InputError("The rooted time-tree topology differs from the OrthoFinder topology "
                         f"after removing extra species (example clade: {example}).")
    return {b: retained[s] for b, s in ref.descendants.items()}


def normalize_numbered_tips(tree, ref=None, force=False):
    recognized = set() if ref is None else ref.species
    leaves = tree.get_terminals()
    stripped = [re.sub(r"^\d+_", "", leaf.name or "", count=1) for leaf in leaves]
    gain = len(set(stripped) & recognized) > len({leaf.name for leaf in leaves} & recognized)
    if force or gain:
        for leaf, name in zip(leaves, stripped):
            leaf.name = name
    return tree


def paml_report(path):
    raw = match_text(path)
    if "Species tree for FigTree" not in raw:
        raise InputError(f"{path}: expected an MCMCtree report containing 'Species tree for FigTree'.")
    region = raw.split("Species tree for FigTree", 1)[1]
    region = re.split(r"(?m)^Posterior mean", region, maxsplit=1)[0]
    numbered, timed = None, None
    for piece in newick_strings(region):
        tree = read_newick(piece)
        terms = tree.get_terminals()
        if (all(c.branch_length is not None for c in tree.find_clades() if c != tree.root)
                and timed is None):
            timed = tree
        if (terms and all(re.match(r"^\d+_", x.name or "") for x in terms) and all(
                node_name(x).isdigit() for x in tree.get_nonterminals()) and numbered is None):
            numbered = tree
    if timed is None:
        raise InputError(f"{path}: MCMCtree report contains no tree with posterior time branch lengths.")
    return numbered, timed


def fixed_time_clock(ref, path, unit, name):
    tree = normalize_numbered_tips(parse_tree_file(path), ref)
    projection = project_tree(ref, tree)
    # Ignore extra outgroups when checking distances, and use projected branch durations.
    root_paths = {}
    projected_root = projection[ref.root]
    for b in ref.species:
        path_to_tip = tree.get_path(projection[b])
        root_paths[b] = path_to_tip[path_to_tip.index(projected_root) + 1:] \
            if projected_root != tree.root else path_to_tip
    for c in set().union(*root_paths.values()):
        if c.branch_length is None or not math.isfinite(c.branch_length) or c.branch_length < 0:
            raise InputError(f"{path}: all projected paths in a time tree need finite nonnegative lengths.")
    depth = tree.depths()
    heights = [math.fsum(c.branch_length for c in path) for path in root_paths.values()]
    if max(heights) <= 0 or max(heights) - min(heights) > max(1e-7, max(heights)*1e-4):
        raise InputError(f"{path}: projected time tree is not ultrametric. Check time units and tree format.")
    lengths = [(depth[projection[b]] - depth[projection[ref.parent[b]]]) * UNIT_TO_MYR[unit]
               for b in ref.scored]
    positive_lengths(ref.scored, lengths)
    return Clock(name, np.asarray([lengths]), "Myr", str(path))


def posterior_draws(path, wanted):
    with open_text(path) as f:
        header = f.readline().split()
        if len(header) != len(set(header)) or not set(wanted) <= set(header):
            missing = sorted(set(wanted) - set(header))
            raise InputError(f"{path}: MCMCtree trace lacks age columns: {', '.join(missing[:12])}")
        cols = [header.index(key) for key in wanted]
        # Some versions write an extra initial state at generation 1, followed by
        # the first sampling interval. A regular every-generation trace starts 1, 2.
        has_generation = "Gen" in header
        if has_generation:
            cols.insert(0, header.index("Gen"))
        try:
            draws = np.loadtxt(f, usecols=cols, ndmin=2)
        except ValueError as exc:
            raise InputError(f"{path}: invalid MCMCtree trace: {exc}") from exc
    if has_generation and len(draws) > 2:
        g = draws[:, 0]
        if g[0] == 1 and g[1] > 1 and g[2] - g[1] == g[1]:
            draws = draws[1:]
    if has_generation:
        draws = draws[:, 1:]
    if len(draws) < 2 or not np.isfinite(draws).all():
        raise InputError(f"{path}: expected at least two finite joint posterior age draws.")
    return draws


def posterior_time_clock(ref, report, traces, unit, name):
    numbered, _ = paml_report(report)
    if numbered is None:
        raise InputError(f"{report}: cannot find the numbered tree for linking t_n columns to clades.")
    projection = project_tree(ref, normalize_numbered_tips(numbered, force=True))
    names = set(ref.scored) | {ref.parent[b] for b in ref.scored}
    ordered = sorted(names, key=natural)
    labels = [f"t_n{node_name(projection[b])}" for b in ordered]
    if any(not label[3:].isdigit() for label in labels):
        raise InputError(f"{report}: numbered tree has missing or invalid MCMCtree node IDs.")
    index = {b: i for i, b in enumerate(ordered)}
    draws = [posterior_draws(path, labels) for path in traces]
    left = [index[ref.parent[b]] for b in ref.scored]
    right = [index[b] for b in ref.scored]
    intervals = []
    for i, (chain, path) in enumerate(zip(draws, traces)):
        dt = chain[:, left] - chain[:, right]
        if (dt < 0).any():
            raise InputError(f"{path}: an ancestral age is younger than its daughter. "
                             f"Check that this trace came from the numbered tree in {report}.")
        valid = (dt > 0).all(axis=1)
        if not valid.all():
            n = np.count_nonzero(valid)
            if n < 2 or n < .99*len(chain):
                raise InputError(f"{path}: too many posterior draws have equal parent and daughter ages. "
                                 "Check the dating model or export ages at higher precision.")
            LOG.warning("%s: omitted %d/%d joint draws with equal parent and daughter ages in the trace.",
                        path, len(chain) - n, len(chain))
            draws[i] = chain[valid]
        intervals.append(dt[valid])
    durations = np.concatenate(intervals)*UNIT_TO_MYR[unit]
    if not np.isfinite(durations).all() or (durations <= 0).any():
        raise InputError(f"{report}: scaled branch durations must be finite and positive.")
    return Clock(name, durations, "Myr", f"{report}; " + "; ".join(map(str, traces))), draws, ordered

