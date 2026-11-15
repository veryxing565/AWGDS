"""Shared file, input validation and external process utilities."""

import csv
import gzip
import json
import logging
import os
import re
import shlex
import shutil
import subprocess
import signal
import time
from pathlib import Path

LOG = logging.getLogger("awgds")
SPECIES_NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9_-]*\Z")


class InputError(ValueError):
    """A configuration or input-data error suitable for the command line."""


def open_text(path, mode="r"):
    if Path(path).suffix.lower() == ".gz":
        return gzip.open(path, mode + "t", encoding="utf-8-sig" if mode == "r" else "utf-8",
                         newline="")
    return open(path, mode, encoding="utf-8-sig" if mode == "r" else "utf-8", newline="")


def tsv_writer(stream, fields):
    writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
    writer.writerow(fields)
    return writer


def read_table(path, required):
    with open_text(path) as stream:
        reader = csv.DictReader((line for line in stream if line.strip() and not line.startswith("#")),
                                delimiter="\t")
        if not reader.fieldnames or not set(required) <= set(reader.fieldnames):
            raise InputError(f"{path}: expected tab-separated columns: {', '.join(required)}")
        for row in reader:
            if None in row or any(row[k] is None for k in required):
                raise InputError(f"{path}: malformed row near line {reader.line_num}")
            yield row


def file_key(path):
    return Path(path).name.split(".", 1)[0]


def genomes_and_annotations(genome_dir, gff_dir):
    def scan(folder, suffixes):
        if not folder.is_dir():
            raise InputError(f"Input directory does not exist: {folder}")
        found, folded = {}, {}
        for path in sorted(folder.rglob("*")):
            name = path.name.lower().removesuffix(".gz")
            if not name.endswith(suffixes) or not path.is_file():
                continue
            key = file_key(path)
            if not SPECIES_NAME.fullmatch(key):
                raise InputError(f"Invalid species key {key!r} in {path}. See awgds --help.")
            if key.casefold() in folded:
                raise InputError(f"Duplicate or case-conflicting key {key!r}: "
                                 f"{folded[key.casefold()]} and {path}")
            folded[key.casefold()] = path
            found[key] = path.resolve()
        if not found:
            raise InputError(f"No input files ending in {', '.join(suffixes)} (optionally .gz) in {folder}")
        return found

    genomes = scan(genome_dir, (".fna", ".fa", ".fasta", ".fas"))
    gffs = scan(gff_dir, (".gff", ".gff3"))
    same_names(genomes, gffs, "genome filenames", "GFF filenames")
    if len(genomes) < 4:
        raise InputError("At least four species are needed to compare two non-root internal branches.")
    return {key: (genomes[key], gffs[key]) for key in sorted(genomes)}


def same_names(expected, actual, expected_name, actual_name):
    expected, actual = set(expected), set(actual)
    if expected != actual:
        def short(items):
            ordered = sorted(items)
            return ", ".join(ordered[:12]) + (" ..." if len(ordered) > 12 else "") if ordered else "none"
        raise InputError(f"Species names differ between {expected_name} and {actual_name}. "
                         f"Missing in {actual_name}: {short(expected - actual)}. "
                         f"Unexpected in {actual_name}: {short(actual - expected)}. "
                         "Names are case-sensitive; use the filename prefix before the first dot.")


def as_json(path, data):
    temp = path.with_name(path.name + ".tmp")
    temp.write_text(json.dumps(data, indent=2, ensure_ascii=True, allow_nan=False,
                               default=str) + "\n", encoding="utf-8")
    temp.replace(path)


def executable(command):
    resolved = shutil.which(os.path.expanduser(str(command)))
    if not resolved:
        raise InputError(f"Executable not found: {command}. Install it or supply its path (see --help).")
    return os.path.abspath(resolved)


class Runner:
    """Keep complete logs on disk and report progress for long-running commands."""

    def __init__(self, log_dir):
        self.log_dir = log_dir
        self.commands = []

    def run(self, name, argv, cwd, output=None, quiet=False, extra_path=()):
        argv = [str(x) for x in argv]
        self.log_dir.mkdir(parents=True, exist_ok=True)
        log_file = self.log_dir / f"{name}.log"
        env = dict(os.environ)
        binary_dir = str(Path(argv[0]).parent)
        env["PATH"] = os.pathsep.join([binary_dir, *(str(p) for p in extra_path), env.get("PATH", "")])
        env["PYTHONUNBUFFERED"] = "1"
        record = {"stage": name, "argv": argv, "cwd": str(cwd), "log": str(log_file),
                  "stdout": str(output or log_file)}
        self.commands.append(record)
        if not quiet:
            LOG.info("Starting %s (log: %s)", name, log_file.name)
        started = last = time.monotonic()
        with log_file.open("wb") as diagnostics:
            diagnostics.write(("$ " + shlex.join(argv) + "\n").encode("utf-8"))
            diagnostics.flush()
            target = output.open("wb") if output is not None else diagnostics
            try:
                posix = os.name == "posix"
                with subprocess.Popen(argv, cwd=cwd, stdout=target, stderr=diagnostics,
                                      env=env, start_new_session=posix) as proc:
                    try:
                        while True:
                            try:
                                code = proc.wait(timeout=5)
                                break
                            except subprocess.TimeoutExpired:
                                now = time.monotonic()
                                if now - last >= 30:
                                    LOG.info("%s is running (%.0f s; %s)", name, now-started,
                                             log_file.name)
                                    last = now
                    except BaseException:
                        if posix:
                            try:
                                os.killpg(proc.pid, signal.SIGTERM)
                            except ProcessLookupError:
                                pass
                        else:
                            proc.terminate()
                        try:
                            proc.wait(timeout=5)
                        except subprocess.TimeoutExpired:
                            if posix:
                                os.killpg(proc.pid, signal.SIGKILL)
                            else:
                                proc.kill()
                            proc.wait()
                        raise
            finally:
                if output is not None:
                    target.close()
        record["exit_code"] = code
        if code:
            with log_file.open("rb") as stream:
                stream.seek(max(0, log_file.stat().st_size - 4096))
                tail = stream.read().decode("utf-8", "replace").strip().splitlines()[-6:]
            raise InputError(f"{name} exited with status {code}; see {log_file}\n" + "\n".join(tail))
        if not quiet:
            LOG.info("Finished %s in %.1f s", name, time.monotonic()-started)


def ids(cell):
    return [part.strip() for part in cell.split(",") if part.strip()]

