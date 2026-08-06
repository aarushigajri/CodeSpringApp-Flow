#!/usr/bin/env bash

set -Eeuo pipefail

repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
test_root="$(mktemp -d /tmp/codespring-fetchngs-smoke.XXXXXX)"
test_home="$test_root/home"
input_file="$test_root/accessions.txt"

cleanup() {
  [[ "$test_root" == /tmp/codespring-fetchngs-smoke.* ]] && rm -rf -- "$test_root"
}
trap cleanup EXIT

mkdir -p "$test_home"
printf 'SRR14593545\nERR1160846\n' > "$input_file"

HOME="$test_home" "$repo_root/bin/codespringflow" fetchngs \
  --input "$input_file" \
  --name smoke_fastq \
  --dry-run

fastq_run="$test_home/csl_results/fetchngs/smoke_fastq"
test -d "$fastq_run/results"
test -d "$test_home/.codespringflow/cache/singularity"
test -d "$test_home/.codespringflow/work/fetchngs/smoke_fastq"
test ! -e "$fastq_run/job_id.txt"
test -f "$fastq_run/input/accessions.csv"
cmp "$input_file" "$fastq_run/input/accessions.csv"
grep -Fq "input: '$fastq_run/input/accessions.csv'" "$fastq_run/params.yml"
bash -n "$fastq_run/run.sbatch"
grep -Fq "download_method: 'sratools'" "$fastq_run/params.yml"
grep -Fq "skip_fastq_download: false" "$fastq_run/params.yml"
grep -Fq $'fetchngs_version\t1.12.0' "$fastq_run/run_manifest.tsv"
grep -Fq $'nextflow_version\t24.04.4' "$fastq_run/run_manifest.tsv"

HOME="$test_home" "$repo_root/bin/codespringflow" fetchngs \
  --input "$input_file" \
  --name smoke_metadata \
  --metadata-only \
  --dry-run

metadata_run="$test_home/csl_results/fetchngs/smoke_metadata"
grep -Fq "skip_fastq_download: true" "$metadata_run/params.yml"
grep -Fq $'metadata_only\ttrue' "$metadata_run/run_manifest.tsv"

if HOME="$test_home" "$repo_root/bin/codespringflow" fetchngs --input "$input_file" --name smoke_fastq --dry-run >/dev/null 2>&1; then
  echo "Duplicate FetchNGS run names were not rejected." >&2
  exit 1
fi

if HOME="$test_home" "$repo_root/bin/codespringflow" fetchngs --input "$input_file" --name '../unsafe' --dry-run >/dev/null 2>&1; then
  echo "Unsafe FetchNGS run names were not rejected." >&2
  exit 1
fi

test "$(HOME="$test_home" "$repo_root/bin/codespringflow" root)" = "$test_home/csl_results/fetchngs"
test "$(HOME="$test_home" "$repo_root/bin/codespringflow" runtime-root)" = "$test_home/.codespringflow"

fake_bin="$test_root/fake-bin"
mkdir -p "$fake_bin"
printf '#!/usr/bin/env bash\nprintf "Submitted batch job 424242\\n"\n' > "$fake_bin/sbatch"
chmod 0755 "$fake_bin/sbatch"

PATH="$fake_bin:$PATH" HOME="$test_home" NEXTFLOW_BIN=/bin/true \
  "$repo_root/bin/codespringflow" fetchngs \
  --input "$input_file" \
  --name smoke_submit

submitted_run="$test_home/csl_results/fetchngs/smoke_submit"
test "$(<"$submitted_run/job_id.txt")" = "424242"
test "$(wc -l < "$submitted_run/job_history.txt")" -eq 1

PATH="$fake_bin:$PATH" HOME="$test_home" NEXTFLOW_BIN=/bin/true \
  "$repo_root/bin/codespringflow" resume smoke_submit
test "$(wc -l < "$submitted_run/job_history.txt")" -eq 2

echo "CodeSpringApp FetchNGS bundle smoke tests passed."
