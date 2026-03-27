#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

DATA_DIR="$ROOT_DIR/opt/single_objective/comb/knapsack_problem/data"

RESULTS_DIR="$ROOT_DIR/comparison/knapsack/results"
GA_DIR="$RESULTS_DIR/ga"
VNS_DIR="$RESULTS_DIR/vns"

mkdir -p "$GA_DIR"
mkdir -p "$VNS_DIR"

SEED=123
EVALUATIONS_MAX=1000

for input_file in "$DATA_DIR"/knapsack_*.txt; do
    base_name="$(basename "$input_file" .txt)"

    echo "Running GA for $base_name"

    poetry run python -m opt.single_objective.comb.knapsack_problem.solver \
        --input-file "$input_file" \
        --method ga \
        --seed "$SEED" \
        --evaluations-max "$EVALUATIONS_MAX" \
        > "$GA_DIR/${base_name}_ga.txt"


    echo "Running VNS for $base_name"

    poetry run python -m opt.single_objective.comb.knapsack_problem.solver \
        --input-file "$input_file" \
        --method vns \
        --seed "$SEED" \
        --evaluations-max "$EVALUATIONS_MAX" \
        > "$VNS_DIR/${base_name}_vns.txt"

done

echo "Done."
