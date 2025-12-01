#!/usr/bin/env bash
# Always fail if any command fails and treat unset vars as errors
set -euo pipefail

usage() {
    cat <<EOF
Usage: $(basename "$0") -i INPUT_FASTA [-c COVERAGE] [-p PERCENTS] [-h]

Arguments:
  -i, --input       Path to input FASTA file (required)
  -c, --coverage    Minimum coverage for mmseqs (--cov). Accepts decimal (0.9) or percent (90). Default: 0.9
  -p, --identity    Comma-separated or space-separated identity percentages to cluster at (default: "1 0.95 0.90 0.80 0.70 0.60 0.50 0.40 0.30")
  -h, --help        Show this help and exit

Examples:
  $(basename "$0") -i input.fasta
  $(basename "$0") --input input.fasta --coverage 90 --identity "0.95,0.90,0.80"
EOF
}

# Parse arguments
input=""
coverage="0.9"
identity="0.95 0.90 0.80 0.70 0.60 0.50 0.40 0.30"

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)
            shift
            input="${1:-}"
            ;;
        -c|--coverage)
            shift
            coverage="${1:-}"
            ;;
        -p|--identity)
            shift
            identity="${1:-}"
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            usage
            exit 2
            ;;
    esac
    shift || true
done

# Basic validations
if [[ -z "${input}" ]]; then
    echo "ERROR: Input file is required." >&2
    usage
    exit 2
fi

if ! command -v mmseqs >/dev/null 2>&1; then
    echo "ERROR: mmseqs not found in PATH. Please install mmseqs and ensure it is available." >&2
    exit 3
fi

if [[ ! -s "${input}" ]]; then
    echo "ERROR: Input file '${input}' does not exist or is empty." >&2
    exit 4
fi

# Normalize identity  list: accept comma or space-separated values
# Convert commas to spaces
identity="${identity//,/ }"
# Trim extra spaces
identity="$(echo "${identity}" | xargs)"

# Validate coverage and convert percent to decimal if needed
# Accept both 0.9 and 90
if ! awk -v c="${coverage}" 'BEGIN{ if (c+0!=c+0) exit 1; exit 0 }' >/dev/null 2>&1; then
    echo "ERROR: Coverage must be numeric (like 0.9 or 90)." >&2
    exit 5
fi

# convert coverage to a number using awk to support floats
coverage_num=$(awk -v c="$coverage" 'BEGIN { if (c > 1) { c = c / 100 } printf "%.6f", c }')

# Ensure coverage between 0 and 1 (exclusive of 0)
if awk -v c="$coverage_num" 'BEGIN { if (c <= 0 || c > 1) exit 1; exit 0 }' >/dev/null 2>&1; then
    :
else
    echo "ERROR: Coverage must be > 0 and <= 1 after normalizing (e.g., 0.9 or 90)." >&2
    exit 6
fi

# Ensure identity values are numeric and between 0 and 1, convert if needed
validated_identities=()
for id in ${identity}; do
    if ! awk -v i="$id" 'BEGIN{ if (i+0!=i+0) exit 1; exit 0 }' >/dev/null 2>&1; then
        echo "ERROR: Identity values must be numeric (like 1, 0.95, or 90)." >&2
        exit 5
    fi
    # Convert percent to decimal if needed
    if (( $(awk -v i="$id" 'BEGIN { print (i > 1) }') )); then
        id_decimal=$(awk -v i="$id" 'BEGIN { printf "%.3f", i / 100 }')
    else
        id_decimal=$(awk -v i="$id" 'BEGIN { printf "%.3f", i }')
    fi
    # Validate range
    if awk -v i="$id_decimal" 'BEGIN { if (i < 0 || i > 1) exit 1; exit 0 }' >/dev/null 2>&1; then
        validated_identities+=("$id_decimal")
    else
        echo "ERROR: Identity values must be between 0 and 1 after normalizing (e.g., 1, 0.95, or 90)." >&2
        exit 6
    fi
done

echo "Input: ${input}"
echo "Coverage: ${coverage_num}"
echo "Identity: ${validated_identities[*]}"

# Cleanup function and trap
tmpdirs_to_clean=()
cleanup() {
    for dir in "${tmpdirs_to_clean[@]:-}"; do
        if [[ -d "$dir" ]]; then
            rm -rf "$dir"
        fi
    done
}
trap cleanup EXIT

# Clustering loop
current_input="${input}"
for min_seq_id in ${validated_identities[*]}; do
    echo "Clustering at ${min_seq_id} sequence identity, coverage ${coverage_num}"
    output="mmseqs_${min_seq_id}_clusters"
    # Create a unique temporary directory for mmseqs workfiles
    tmpdir=$(mktemp -d -t mmseqs_tmp_"${min_seq_id}"_XXXX)
    tmpdirs_to_clean+=("${tmpdir}")
    # Convert to float in 0-1 range for mmseqs
    # Run mmseqs easy-cluster
    # Set the k-mer size to 3 for min_seq_id higher than 0.8
    if (( $(awk -v i="${min_seq_id}" 'BEGIN { print (i > 0.8) }') )); then
        command=(mmseqs easy-cluster "${current_input}" --dbtype 2 -k 3 --min-seq-id "${min_seq_id}" -c "${coverage_num}" --cov-mode 5 "${output}" "${tmpdir}")
    else
        command=(mmseqs easy-cluster "${current_input}" --dbtype 2 -k 4 --min-seq-id "${min_seq_id}" -c "${coverage_num}" --cov-mode 5 "${output}" "${tmpdir}")
    fi
    echo command: "${command[@]}"
    if ! eval "${command[@]}"; then
        echo "ERROR: mmseqs easy-cluster failed for ${min_seq_id} identity." >&2
        exit 7
    fi

    # Expected representative fasta created by mmseqs
    rep_seq="${output}_rep_seq.fasta"
    if [[ ! -s "${rep_seq}" ]]; then
        echo "ERROR: Expected representative sequence file '${rep_seq}' not found after clustering ${min_seq_id} identity." >&2
        exit 8
    fi

    # Use rep_seq as next iteration input
    current_input="${rep_seq}"
done

echo "Clustering completed successfully."

