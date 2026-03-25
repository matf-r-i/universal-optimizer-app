#!/usr/bin/env python3

from __future__ import annotations

import csv
import re
from pathlib import Path


ROOT_DIR = Path(__file__).resolve().parents[2]
RESULTS_DIR = ROOT_DIR / "comparison" / "knapsack" / "results"
OUTPUT_CSV = ROOT_DIR / "visualization" / "knapsack" / "results_summary.csv"


PATTERNS = {
    "representation": re.compile(r"^Best solution representation:\s*(.+)$"),
    "code": re.compile(r"^Best solution code:\s*(.+)$"),
    "objective": re.compile(r"^Best solution objective:\s*([+-]?\d+(?:\.\d+)?)$"),
    "fitness": re.compile(r"^Best solution fitness:\s*([+-]?\d+(?:\.\d+)?)$"),
    "feasible": re.compile(r"^Best solution feasible:\s*(True|False)$"),
    "iterations": re.compile(r"^Number of iterations:\s*(\d+)$"),
    "evaluations": re.compile(r"^Number of evaluations:\s*(\d+)$"),
}


def parse_result_file(file_path: Path, method: str) -> dict[str, object]:
    text = file_path.read_text(encoding="utf-8")
    lines = text.splitlines()

    data: dict[str, object] = {
        "instance": file_path.stem.replace(f"_{method}", ""),
        "method": method,
        "representation": "",
        "code": "",
        "objective": None,
        "fitness": None,
        "feasible": None,
        "iterations": None,
        "evaluations": None,
    }

    for line in lines:
        for key, pattern in PATTERNS.items():
            match = pattern.match(line.strip())
            if not match:
                continue

            value = match.group(1)
            if key in {"objective", "fitness"}:
                data[key] = float(value)
            elif key in {"iterations", "evaluations"}:
                data[key] = int(value)
            elif key == "feasible":
                data[key] = value == "True"
            else:
                data[key] = value

    missing = [
        key for key in ("objective", "fitness", "feasible", "iterations", "evaluations")
        if data[key] is None
    ]
    if missing:
        raise ValueError(f"Missing fields {missing} in file: {file_path}")

    return data


def collect_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []

    for method in ("ga", "vns"):
        method_dir = RESULTS_DIR / method
        if not method_dir.exists():
            raise FileNotFoundError(f"Missing results directory: {method_dir}")

        for file_path in sorted(method_dir.glob("knapsack_*.txt")):
            rows.append(parse_result_file(file_path, method))

    rows.sort(key=lambda row: (str(row["instance"]), str(row["method"])))
    return rows


def write_csv(rows: list[dict[str, object]]) -> None:
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "instance",
        "method",
        "representation",
        "code",
        "objective",
        "fitness",
        "feasible",
        "iterations",
        "evaluations",
    ]

    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def print_summary(rows: list[dict[str, object]]) -> None:
    print(f"Saved summary CSV to: {OUTPUT_CSV}")
    print()
    print("Parsed results:")
    for row in rows:
        print(
            f"{row['instance']} | {row['method']} | "
            f"objective={row['objective']} | "
            f"fitness={row['fitness']} | "
            f"feasible={row['feasible']} | "
            f"iterations={row['iterations']} | "
            f"evaluations={row['evaluations']}"
        )


def main() -> None:
    rows = collect_rows()
    write_csv(rows)
    print_summary(rows)


if __name__ == "__main__":
    main()
