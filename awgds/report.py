"""Run layout and the unified, one-row-per-branch result table."""

import csv

from .io import InputError, open_text, tsv_writer


class RunPaths:
    def __init__(self, root, legacy=False):
        self.root = root
        self.summary = root / "summary.tsv"
        self.manifest = root / "run.json"
        self.logs = root / "logs"
        self.work = root / "work"
        self.inputs = root if legacy else self.work / "inputs"
        self.orthofinder = self.work / "orthofinder"
        self.trees = root if legacy else self.work / "trees"
        self.similarity = self.work / "similarity"
        self.branches = root if legacy else self.work / "branches"
        self.adi = self.work / "adi"
        self.evidence = self.work / "evidence"
        self.pairs = self.branches / "pairs"
        self.species_tree = self.trees / "species_tree.nwk"
        self.input_files = self.inputs / "input_files.tsv"


WGD_FIELDS = (
    "status", "selected_references", "tested_target_reference_combinations",
    "positive_target_reference_combinations", "reference_coverage", "best_reference",
    "shared_reference_chromosomes", "testable_reference_chromosomes", "shared_chromosome_fraction",
    "two_sided_references", "left_eligible_targets", "right_eligible_targets",
    "local_check", "local_matches", "positive_pairs", "local_target_pairs",
)


def read_scores(path):
    """Read ADI columns from either legacy or unified output."""
    with open_text(path) as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        fields = reader.fieldnames or []
        if not {"node", "recommended"} <= set(fields) or len(fields) != len(set(fields)):
            raise InputError(f"{path}: expected unique columns including node and recommended.")
        columns = tuple(k for k in fields if not k.startswith("wgd_"))
        seen, rows = set(), []
        for row in reader:
            if None in row or any(row[k] is None for k in fields) or not row["node"] or row["node"] in seen:
                raise InputError(f"{path}: missing, duplicate or malformed branch at line {reader.line_num}.")
            seen.add(row["node"])
            rows.append({k: row[k] for k in columns})
    return columns, rows


def write_summary(path, columns, rows, validation=None, awaiting="pending"):
    """Publish the ADI checkpoint or final merged summary atomically."""
    validated = {} if validation is None else validation
    temp = path.with_name(f".{path.name}.tmp")
    with open_text(temp, "w") as stream:
        out = tsv_writer(stream, (*columns, *(f"wgd_{key}" for key in WGD_FIELDS)))
        for row in rows:
            evidence = validated.get(row["node"])
            if evidence is None:
                status = {"no": "not_recommended", "yes": awaiting}.get(row["recommended"], "not_scored")
                evidence = {"status": status}
            out.writerow((*[row[k] for k in columns], *[evidence.get(k, "NA") for k in WGD_FIELDS]))
    temp.replace(path)
