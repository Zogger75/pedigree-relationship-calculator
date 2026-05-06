#!/usr/bin/env python3
"""
pedigree_relationship_calculator.py

A standalone, portfolio-friendly coefficient of parentage / pedigree relationship
calculator. This version reads a plain text pedigree file instead of querying a
CIERA database.

The goal is to demonstrate:
    - text-file parsing
    - recursive pedigree traversal
    - method-code handling for different germplasm derivation types
    - pairwise relationship matrix generation
    - clear, documented Python suitable for a GitHub portfolio

Input file format
-----------------
The input file is pipe-delimited by default:

    gid|name|dam|sire|method_code|include

Columns:
    gid         Unique germplasm identifier. Can be numeric or text.
    name        Display name for the germplasm.
    dam         Female parent / source parent / group parent. Blank or 0 = unknown.
    sire        Male parent / secondary parent. Blank or 0 = unknown.
    method_code Basic derivation method. Supported examples:
                  FND = founder / no known parents
                  GEN = generative cross; relationship comes from both parents
                  DER = derivative/selection; treated as genetically equivalent
                        to the source parent unless two parents are supplied
                  MAN = management/increase/release; no genetic change, so it
                        unwraps to its source parent
                  SEL = selection; treated like DER
    include     Optional yes/no flag. Lines marked yes are used in the final matrix.

Usage
-----
    python pedigree_relationship_calculator.py example_pedigree.txt

Optional:
    python pedigree_relationship_calculator.py example_pedigree.txt --output cop_matrix.csv
    python pedigree_relationship_calculator.py example_pedigree.txt --debug-pair LINE_A LINE_B

Notes
-----
This is a generic implementation intended for demonstration and portfolio use.
It is not tied to private CIERA data or database tables.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

UNKNOWN_PARENT_VALUES = {"", "0", "NA", "N/A", "NONE", "NULL", "."}


@dataclass(frozen=True)
class GermplasmRecord:
    """One row of pedigree information.

    Attributes:
        gid: Unique germplasm identifier.
        name: Human-readable display name.
        dam: Female/source parent identifier, or None if unknown.
        sire: Male/second parent identifier, or None if unknown.
        method_code: Basic derivation method code such as FND, GEN, DER, or MAN.
        include: Whether this line should be included in the final COP matrix.
    """

    gid: str
    name: str
    dam: Optional[str]
    sire: Optional[str]
    method_code: str
    include: bool = False


class PedigreeError(ValueError):
    """Raised when the pedigree input file contains invalid or inconsistent data."""


def normalize_parent(value: str) -> Optional[str]:
    """Convert blank/unknown parent values to None."""

    cleaned = str(value).strip()
    return None if cleaned.upper() in UNKNOWN_PARENT_VALUES else cleaned


def parse_bool(value: str) -> bool:
    """Parse a human-friendly boolean value from the include column."""

    return str(value).strip().upper() in {"Y", "YES", "TRUE", "T", "1", "INCLUDE"}


def read_pedigree_file(input_path: Path, delimiter: str = "|") -> Dict[str, GermplasmRecord]:
    """Read a pipe-delimited pedigree text file into a dictionary.

    Args:
        input_path: Path to the pedigree text file.
        delimiter: Field separator. Defaults to pipe because line names often
            contain spaces, slashes, or commas.

    Returns:
        Dictionary keyed by germplasm id.

    Raises:
        PedigreeError: If required columns are missing or duplicate GIDs exist.
    """

    required_columns = {"gid", "name", "dam", "sire", "method_code"}
    records: Dict[str, GermplasmRecord] = {}
    input_path = Path(input_path)
    with input_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None:
            raise PedigreeError("The pedigree file appears to be empty.")

        # Normalize headers so the file is forgiving about case and spaces.
        reader.fieldnames = [header.strip().lower() for header in reader.fieldnames]
        missing = required_columns - set(reader.fieldnames)
        if missing:
            raise PedigreeError(f"Missing required columns: {', '.join(sorted(missing))}")

        for row_number, row in enumerate(reader, start=2):
            gid = row["gid"].strip()
            if not gid:
                raise PedigreeError(f"Row {row_number}: gid cannot be blank.")
            if gid in records:
                raise PedigreeError(f"Row {row_number}: duplicate gid '{gid}'.")

            method_code = row["method_code"].strip().upper() or "GEN"
            include = parse_bool(row.get("include", "no"))

            records[gid] = GermplasmRecord(
                gid=gid,
                name=row["name"].strip() or gid,
                dam=normalize_parent(row["dam"]),
                sire=normalize_parent(row["sire"]),
                method_code=method_code,
                include=include,
            )

    validate_parent_references(records)
    return records

def choose_pedigree_file() -> str:
    """
    Open a file picker so the user can select the pedigree input file.

    Returns:
        str: Full path to the selected pedigree file.

    Raises:
        SystemExit: If the user cancels the file picker.
    """
    try:
        import tkinter as tk
        from tkinter import filedialog
    except ImportError as exc:
        raise SystemExit(
            "Tkinter is not available. Please provide the pedigree file path manually."
        ) from exc

    root = tk.Tk()
    root.withdraw()

    file_path = filedialog.askopenfilename(
        title="Select pedigree input file",
        filetypes=[
            ("Text files", "*.txt"),
            ("CSV files", "*.csv"),
            ("All files", "*.*"),
        ],
    )

    root.destroy()

    if not file_path:
        raise SystemExit("No pedigree file selected. Analysis cancelled.")

    return file_path
def validate_parent_references(records: Dict[str, GermplasmRecord]) -> None:
    """Ensure every non-blank parent points to another row in the file."""

    missing: List[Tuple[str, str]] = []
    for record in records.values():
        for parent in (record.dam, record.sire):
            if parent is not None and parent not in records:
                missing.append((record.gid, parent))

    if missing:
        formatted = "; ".join(f"{child} references missing parent {parent}" for child, parent in missing)
        raise PedigreeError(formatted)


def genetic_parents(gid: str, records: Dict[str, GermplasmRecord]) -> Tuple[Optional[str], Optional[str]]:
    """Return the genetically contributing parents for a germplasm.

    This is where the simple method-code logic is handled.

    FND:
        Founder. No known parents.
    GEN:
        Generative cross. Both DAM and SIRE contribute when present.
    DER / SEL:
        Derivative or selected line. If one source parent is supplied, the line is
        treated as genetically equivalent to that source. If two parents are
        supplied, it behaves like a normal cross.
    MAN:
        Management, increase, or release step. It represents no genetic change,
        so it unwraps to the source parent. The DAM/source parent is preferred,
        then SIRE if DAM is blank.
    """

    record = records[gid]
    method = record.method_code

    if method == "FND":
        return None, None

    if method == "MAN":
        source = record.dam or record.sire
        return source, source if source else None

    if method in {"DER", "SEL"}:
        if record.dam and not record.sire:
            return record.dam, record.dam
        if record.sire and not record.dam:
            return record.sire, record.sire

    return record.dam, record.sire


def resolve_genetic_identity(
    gid: str,
    records: Dict[str, GermplasmRecord],
    visited: Optional[set[str]] = None,
) -> str:
    """Resolve method codes that represent no new recombination.

    A single-parent DER, SEL, or MAN row is treated as genetically equivalent to
    its source parent. This mirrors the practical interpretation that an increase,
    release, or simple selection step did not create a new cross.

    Args:
        gid: Germplasm id to resolve.
        records: Pedigree dictionary.
        visited: Used internally to detect circular references.

    Returns:
        The GID that should be used for genetic relationship calculations.
    """

    if visited is None:
        visited = set()
    if gid in visited or gid not in records:
        return gid

    visited.add(gid)
    record = records[gid]

    if record.method_code in {"DER", "SEL", "MAN"}:
        source = record.dam or record.sire
        has_only_one_parent = bool(source) and not (record.dam and record.sire)
        if has_only_one_parent:
            return resolve_genetic_identity(source, records, visited)

    return gid



def generation_depth(
    gid: str,
    records: Dict[str, GermplasmRecord],
    memo: Optional[Dict[str, int]] = None,
) -> int:
    """Return a simple generation depth for choosing the recursive direction.

    Founders have depth 0. A child is one generation deeper than its deepest
    known genetic parent. The relationship function uses this to recurse from
    the more-derived line toward the older line, which produces more intuitive
    pairwise results when both lines have known parents.
    """

    if memo is None:
        memo = {}

    gid = resolve_genetic_identity(gid, records)
    if gid in memo:
        return memo[gid]
    if gid not in records:
        return 0

    parents = [resolve_genetic_identity(parent, records) for parent in genetic_parents(gid, records) if parent]
    if not parents:
        memo[gid] = 0
    else:
        memo[gid] = 1 + max(generation_depth(parent, records, memo) for parent in parents)
    return memo[gid]

def relationship(
    gid_a: str,
    gid_b: str,
    records: Dict[str, GermplasmRecord],
    memo: Optional[Dict[Tuple[str, str], float]] = None,
    active_path: Optional[set[Tuple[str, str]]] = None,
) -> float:
    """Recursively calculate the pedigree relationship between two germplasm IDs.

    The calculation follows a simple recursive pedigree rule:

        relationship(child, other) = average relationship of child's genetic
        parent(s) to the other line.

    Single-parent DER, SEL, and MAN rows are first resolved to their genetic
    source. For example, an increase of a line is treated as the same genetic
    identity as the original line.

    For the same line, the diagonal is 1.0 in this portfolio version. That keeps
    the matrix intuitive for recruiters and avoids adding a full inbreeding
    coefficient model. The off-diagonal values are the primary focus.

    Args:
        gid_a: First germplasm id.
        gid_b: Second germplasm id.
        records: Pedigree dictionary.
        memo: Cache for repeated recursive calls.
        active_path: Used to detect cycles in malformed pedigrees.

    Returns:
        Relationship coefficient between 0.0 and 1.0.
    """

    if memo is None:
        memo = {}
    if active_path is None:
        active_path = set()

    gid_a = resolve_genetic_identity(gid_a, records)
    gid_b = resolve_genetic_identity(gid_b, records)

    if gid_a not in records or gid_b not in records:
        return 0.0

    if gid_a == gid_b:
        return 1.0

    key = tuple(sorted((gid_a, gid_b)))
    if key in memo:
        return memo[key]
    if key in active_path:
        raise PedigreeError(f"Cycle detected while evaluating {gid_a} and {gid_b}.")

    active_path.add(key)

    parents_a = tuple(resolve_genetic_identity(parent, records) for parent in genetic_parents(gid_a, records) if parent)
    parents_b = tuple(resolve_genetic_identity(parent, records) for parent in genetic_parents(gid_b, records) if parent)

    if parents_a and parents_b:
        depth_a = generation_depth(gid_a, records)
        depth_b = generation_depth(gid_b, records)
        if depth_a >= depth_b:
            value = sum(relationship(parent, gid_b, records, memo, active_path) for parent in parents_a) / len(parents_a)
        else:
            value = sum(relationship(gid_a, parent, records, memo, active_path) for parent in parents_b) / len(parents_b)
    elif parents_a:
        value = sum(relationship(parent, gid_b, records, memo, active_path) for parent in parents_a) / len(parents_a)
    elif parents_b:
        value = sum(relationship(gid_a, parent, records, memo, active_path) for parent in parents_b) / len(parents_b)
    else:
        value = 0.0

    active_path.remove(key)
    memo[key] = value
    return value


def included_gids(records: Dict[str, GermplasmRecord]) -> List[str]:
    """Return GIDs marked for inclusion, or all non-founder rows if none are marked."""

    marked = [gid for gid, record in records.items() if record.include]
    if marked:
        return marked
    return [gid for gid, record in records.items() if record.method_code != "FND"]


def build_relationship_matrix(gids: Iterable[str], records: Dict[str, GermplasmRecord]) -> List[List[str]]:
    """Build a printable pairwise relationship matrix."""

    gid_list = list(gids)
    memo: Dict[Tuple[str, str], float] = {}
    header = ["Line"] + [records[gid].name for gid in gid_list]
    matrix = [header]

    for gid_a in gid_list:
        row = [records[gid_a].name]
        for gid_b in gid_list:
            row.append(f"{relationship(gid_a, gid_b, records, memo):.6f}")
        matrix.append(row)

    return matrix


def write_matrix_csv(matrix: List[List[str]], output_path: Path) -> None:
    """Write the relationship matrix to a CSV file."""
    output_path = Path(output_path)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerows(matrix)


def print_matrix(matrix: List[List[str]]) -> None:
    """Print an aligned matrix to the console for quick review."""

    widths = [max(len(row[col]) for row in matrix) for col in range(len(matrix[0]))]
    for row in matrix:
        print("  ".join(value.rjust(widths[idx]) for idx, value in enumerate(row)))


def explain_pair(gid_a: str, gid_b: str, records: Dict[str, GermplasmRecord]) -> None:
    """Print a small debug explanation for one relationship pair."""

    value = relationship(gid_a, gid_b, records)
    print(f"\nDebug pair: {records[gid_a].name} × {records[gid_b].name}")
    print(f"Relationship coefficient: {value:.6f}")
    print(f"{gid_a} genetic identity: {resolve_genetic_identity(gid_a, records)}")
    print(f"{gid_b} genetic identity: {resolve_genetic_identity(gid_b, records)}")
    print(f"{gid_a} parents: {genetic_parents(gid_a, records)}")
    print(f"{gid_b} parents: {genetic_parents(gid_b, records)}")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description="Calculate a generic coefficient of parentage / pedigree relationship matrix from a text file."
    )
    parser.add_argument(
    "pedigree_file",
    nargs="?",
    help="Path to the pedigree input text file. If omitted, a file picker will open."
    )   
    parser.add_argument("--output", type=Path, default=Path("cop_matrix.csv"), help="Output CSV path.")
    parser.add_argument("--delimiter", default="|", help="Input delimiter. Default: '|'.")
    parser.add_argument("--debug-pair", nargs=2, metavar=("GID_A", "GID_B"), help="Print details for one pair.")
    return parser.parse_args()


def main() -> None:
    """Command-line entry point."""

    args = parse_args()
    if not args.pedigree_file:
        args.pedigree_file = choose_pedigree_file()
    records = read_pedigree_file(args.pedigree_file, delimiter=args.delimiter)
    gids = included_gids(records)

    if not gids:
        raise PedigreeError("No lines were selected for the matrix. Set include=yes for at least one row.")

    matrix = build_relationship_matrix(gids, records)
    write_matrix_csv(matrix, args.output)
    print_matrix(matrix)
    print(f"\nSaved matrix to: {args.output.resolve()}")

    if args.debug_pair:
        explain_pair(args.debug_pair[0], args.debug_pair[1], records)


if __name__ == "__main__":
    main()
