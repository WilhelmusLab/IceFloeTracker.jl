#!/usr/bin/env bash
# Benchmark the `tracking` Snakemake rule (wall/CPU time, max RSS) for
# before/after comparison of tracker parallelization work.
#
# Usage: ./benchmark_tracking.sh <target> <label> [configfile]
#   target     snakemake target, e.g.
#              track/beaufort_sea-large.250m.2019-03-23.2019-03-23.LopezAcosta2019Tiling.tracked.csv
#   label      free-form tag for this run, e.g. "before" or "after"
#   configfile optional --configfile path (default: configs/large-regions/config.yaml)
#
# Env vars:
#   JULIA_NUM_THREADS  threads Julia will use (default: 10)
#   REPEATS            number of trials to run (default: 1)
set -euo pipefail

TARGET="${1:?usage: $0 <target> <label> [configfile]}"
LABEL="${2:?usage: $0 <target> <label> [configfile]}"
CONFIGFILE="${3:-configs/large-regions/config.yaml}"

export JULIA_NUM_THREADS="${JULIA_NUM_THREADS:-10}"
REPEATS="${REPEATS:-1}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR="${SCRIPT_DIR}/results/benchmarks"
LOG_DIR="${RESULTS_DIR}/logs"
SUMMARY_CSV="${RESULTS_DIR}/tracking_benchmark.csv"

mkdir -p "$LOG_DIR"

if [ ! -f "$SUMMARY_CSV" ]; then
	echo "timestamp,label,branch,commit,julia_num_threads,wall_s,user_s,sys_s,max_rss_mb,exit_code" >"$SUMMARY_CSV"
fi

BRANCH="$(git -C "$SCRIPT_DIR" rev-parse --abbrev-ref HEAD 2>/dev/null || echo unknown)"
COMMIT="$(git -C "$SCRIPT_DIR" rev-parse --short HEAD 2>/dev/null || echo unknown)"

for ((i = 1; i <= REPEATS; i++)); do
	TIMESTAMP="$(date -u +%Y%m%dT%H%M%SZ)"
	TIME_LOG="${LOG_DIR}/${TIMESTAMP}-${LABEL}-run${i}.time.log"

	echo "== Run ${i}/${REPEATS} (label=${LABEL}, threads=${JULIA_NUM_THREADS}) =="

	set +e
	(
		cd "$SCRIPT_DIR"
		/usr/bin/time -v -o "$TIME_LOG" \
			snakemake -R tracking "$TARGET" --configfile "$CONFIGFILE"
	)
	EXIT_CODE=$?
	set -e

	# Parse the fields we need from GNU time's verbose output.
	WALL_RAW="$(grep 'Elapsed (wall clock) time' "$TIME_LOG" | awk -F': ' '{print $2}')"
	USER_S="$(grep 'User time (seconds)' "$TIME_LOG" | awk -F': ' '{print $2}')"
	SYS_S="$(grep 'System time (seconds)' "$TIME_LOG" | awk -F': ' '{print $2}')"
	MAX_RSS_KB="$(grep 'Maximum resident set size' "$TIME_LOG" | awk -F': ' '{print $2}')"

	# Convert wall clock from [h:]mm:ss(.ss) to seconds.
	WALL_S="$(awk -F: '{
		if (NF == 3) print $1*3600 + $2*60 + $3
		else if (NF == 2) print $1*60 + $2
		else print $1
	}' <<<"$WALL_RAW")"
	MAX_RSS_MB="$(awk -v kb="$MAX_RSS_KB" 'BEGIN { printf "%.1f", kb / 1024 }')"

	echo "${TIMESTAMP},${LABEL},${BRANCH},${COMMIT},${JULIA_NUM_THREADS},${WALL_S},${USER_S},${SYS_S},${MAX_RSS_MB},${EXIT_CODE}" >>"$SUMMARY_CSV"

	echo "-> wall=${WALL_S}s user=${USER_S}s sys=${SYS_S}s max_rss=${MAX_RSS_MB}MB exit=${EXIT_CODE} (log: $TIME_LOG)"
done

echo "Summary appended to $SUMMARY_CSV"
