#!/usr/bin/env bash

set -euo pipefail

repo_root="$(cd "$(dirname "$0")/.." && pwd -P)"
runner="$repo_root/scripts/run_tm_extraction.sh"
fixture_root="$(mktemp -d "${TMPDIR:-/tmp}/tm-runner-test.XXXXXX")"
fixture_root="$(cd "$fixture_root" && pwd -P)"
trap 'rm -rf "$fixture_root"' EXIT INT TERM

patent_dir="$fixture_root/Company A/WO2023000001"
mkdir -p "$patent_dir/output/WO2023000001"
printf '{"patent_id":"WO2023000001","pages":[]}' > \
    "$patent_dir/output/WO2023000001/WO2023000001.json"
printf 'sequence listing fixture\n' > "$patent_dir/WO2023000001_sequence.txt"
printf '%s\n' "$patent_dir" > "$fixture_root/jobs.txt"

output="$($runner --job-list "$fixture_root/jobs.txt" --workers 1 --dry-run)"

printf '%s\n' "$output" | grep -F "JSON     : $patent_dir/output/WO2023000001/WO2023000001.json" >/dev/null
printf '%s\n' "$output" | grep -F "Sequence : $patent_dir/WO2023000001_sequence.txt" >/dev/null
printf '%s\n' "$output" | grep -F "Dry run  : 1" >/dev/null

if [[ -d "$patent_dir/output" && -n "$(find "$patent_dir/output" -maxdepth 1 -type f -print -quit)" ]]; then
    echo "dry run unexpectedly wrote output files" >&2
    exit 1
fi

echo "run_tm_extraction.sh dry-run test passed"
