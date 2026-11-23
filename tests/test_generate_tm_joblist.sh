#!/usr/bin/env bash

set -euo pipefail

repo_root="$(cd "$(dirname "$0")/.." && pwd -P)"
generator="$repo_root/scripts/generate_tm_joblist.sh"
fixture_root="$(mktemp -d "${TMPDIR:-/tmp}/tm-joblist-test.XXXXXX")"
fixture_root="$(cd "$fixture_root" && pwd -P)"
trap 'rm -rf "$fixture_root"' EXIT INT TERM

mkdir -p "$fixture_root/Company A/WO2023000001"
mkdir -p "$fixture_root/Company A/WO2023000001/output/WO2023000001"
mkdir -p "$fixture_root/Company_B/WO2023000002"
mkdir -p "$fixture_root/Company_B/US2023000003"
mkdir -p "$fixture_root/WO2023999999"

output_file="$fixture_root/jobs.txt"
"$generator" "$fixture_root" "$output_file" >/dev/null

expected_file="$fixture_root/expected.txt"
printf '%s\n' \
    "$fixture_root/Company A/WO2023000001" \
    "$fixture_root/Company_B/WO2023000002" > "$expected_file"

diff -u "$expected_file" "$output_file"
echo "generate_tm_joblist.sh test passed"
