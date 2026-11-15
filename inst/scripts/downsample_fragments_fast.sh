#!/usr/bin/env bash

# Allow explicit `sh this_script ...`: re-enter Bash before Bash-only syntax.
if [ "${DOWNSAMPLE_FRAGMENTS_BASH_REEXEC:-0}" != "1" ]; then
  export DOWNSAMPLE_FRAGMENTS_BASH_REEXEC=1
  exec /usr/bin/env bash "$0" "$@"
fi
unset DOWNSAMPLE_FRAGMENTS_BASH_REEXEC

set -Eeuo pipefail

# Stream a coordinate-sorted 10x fragments file once, retain each fragment
# record independently, write BGZF, and build a tabix BED index.
#
# Packaged copy of scRT_learning/.../downsample_fragments_fast.sh. Changes:
# input/output are required (no hard-coded default paths), bgzip/tabix/cc
# default to the executables found on PATH, and the temporary directory
# honours $TMPDIR. The sampling algorithm is unchanged.
# Call it from R with scRTtools::downsample_fragments().
#
# Usage:
#   bash downsample_fragments_fast.sh [options] input.tsv.gz output.tsv.gz
#
# Options:
#   -p FLOAT  fragment-record keep probability (default: 0.33)
#   -s INT    random seed (default: 20260821)
#   -t INT    bgzip threads for input/output (default: 16)
#   -l INT    bgzip compression level, 0-9 (default: 1; fastest practical)
#   -f        replace an existing output and index
#   -h        show help

KEEP_PROB="0.33"
RANDOM_SEED="20260821"
N_THREADS="16"
COMPRESSION_LEVEL="1"
REPLACE_OUTPUT="0"

usage() {
  sed -n '/^# Usage:/,/^#   -h/p' "$0" | sed 's/^# \{0,1\}//'
}

while getopts ":p:s:t:l:fh" opt; do
  case "$opt" in
    p) KEEP_PROB="$OPTARG" ;;
    s) RANDOM_SEED="$OPTARG" ;;
    t) N_THREADS="$OPTARG" ;;
    l) COMPRESSION_LEVEL="$OPTARG" ;;
    f) REPLACE_OUTPUT="1" ;;
    h) usage; exit 0 ;;
    :) echo "Missing value for -$OPTARG" >&2; usage >&2; exit 2 ;;
    \?) echo "Unknown option: -$OPTARG" >&2; usage >&2; exit 2 ;;
  esac
done
shift $((OPTIND - 1))

if (( $# != 2 )); then
  echo "Exactly two positional arguments (input, output) are required." >&2
  usage >&2
  exit 2
fi
INPUT_PATH="$1"
OUTPUT_PATH="$2"
if (( $# > 2 )); then
  echo "Too many positional arguments." >&2
  usage >&2
  exit 2
fi

BGZIP_BIN="${BGZIP_BIN:-$(command -v bgzip || true)}"
TABIX_BIN="${TABIX_BIN:-$(command -v tabix || true)}"
CC_BIN="${CC_BIN:-$(command -v cc || true)}"

[[ -r "$INPUT_PATH" ]] || { echo "Input is not readable: $INPUT_PATH" >&2; exit 1; }
[[ "$INPUT_PATH" != "$OUTPUT_PATH" ]] || { echo "Input and output must differ." >&2; exit 1; }
[[ -x "$BGZIP_BIN" ]] || { echo "bgzip is not executable: $BGZIP_BIN" >&2; exit 1; }
[[ -x "$TABIX_BIN" ]] || { echo "tabix is not executable: $TABIX_BIN" >&2; exit 1; }
[[ -x "$CC_BIN" ]] || { echo "C compiler is not executable: $CC_BIN" >&2; exit 1; }
[[ "$RANDOM_SEED" =~ ^[0-9]+$ ]] || { echo "Seed must be a non-negative integer." >&2; exit 2; }
[[ "$N_THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "Threads must be a positive integer." >&2; exit 2; }
[[ "$COMPRESSION_LEVEL" =~ ^[0-9]$ ]] || { echo "Compression level must be 0-9." >&2; exit 2; }
awk -v p="$KEEP_PROB" 'BEGIN { exit !(p > 0 && p <= 1) }' || {
  echo "Keep probability must be in (0, 1]." >&2
  exit 2
}

OUTPUT_DIR=$(dirname -- "$OUTPUT_PATH")
OUTPUT_BASE=$(basename -- "$OUTPUT_PATH")
[[ -d "$OUTPUT_DIR" && -w "$OUTPUT_DIR" ]] || {
  echo "Output directory is not writable: $OUTPUT_DIR" >&2
  exit 1
}

if [[ "$REPLACE_OUTPUT" != "1" && ( -e "$OUTPUT_PATH" || -e "${OUTPUT_PATH}.tbi" ) ]]; then
  echo "Output or index already exists; use -f to replace it: $OUTPUT_PATH" >&2
  exit 1
fi

RUN_TMPDIR=$(mktemp -d "${TMPDIR:-/tmp}/downsample_fragments.XXXXXX")
PARTIAL_PATH="$OUTPUT_DIR/.${OUTPUT_BASE%.gz}.partial.$$.gz"
PARTIAL_INDEX="${PARTIAL_PATH}.tbi"
SAMPLER_STATS="$RUN_TMPDIR/sampler.stats"
SAMPLER_BIN="$RUN_TMPDIR/fragment_sampler"

cleanup() {
  status=$?
  trap - EXIT INT TERM
  rm -f -- "$PARTIAL_PATH" "$PARTIAL_INDEX"
  rm -rf -- "$RUN_TMPDIR"
  exit "$status"
}
trap cleanup EXIT INT TERM

"$CC_BIN" -O3 -std=c99 -D_FILE_OFFSET_BITS=64 -x c -o "$SAMPLER_BIN" - <<'C_CODE'
#define _POSIX_C_SOURCE 200809L
#include <errno.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

static uint64_t rng_state;

static inline uint64_t splitmix64(void) {
    uint64_t z = (rng_state += UINT64_C(0x9e3779b97f4a7c15));
    z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
    return z ^ (z >> 31);
}

static inline double uniform01(void) {
    return (splitmix64() >> 11) * 0x1.0p-53;
}

static double wall_seconds(void) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec + (double)tv.tv_usec / 1e6;
}

int main(int argc, char **argv) {
    if (argc != 3) {
        fprintf(stderr, "fragment_sampler requires PROBABILITY and SEED\n");
        return 2;
    }

    char *end = NULL;
    errno = 0;
    double probability = strtod(argv[1], &end);
    if (errno || !end || *end != '\0' || probability <= 0.0 || probability > 1.0) {
        fprintf(stderr, "Invalid probability: %s\n", argv[1]);
        return 2;
    }

    errno = 0;
    rng_state = strtoull(argv[2], &end, 10);
    if (errno || !end || *end != '\0') {
        fprintf(stderr, "Invalid seed: %s\n", argv[2]);
        return 2;
    }

    const size_t io_buffer_size = 8U * 1024U * 1024U;
    char *input_buffer = malloc(io_buffer_size);
    char *output_buffer = malloc(io_buffer_size);
    if (!input_buffer || !output_buffer) {
        fprintf(stderr, "Unable to allocate I/O buffers\n");
        return 1;
    }
    setvbuf(stdin, input_buffer, _IOFBF, io_buffer_size);
    setvbuf(stdout, output_buffer, _IOFBF, io_buffer_size);

    char *line = NULL;
    size_t capacity = 0;
    ssize_t line_length;
    uint64_t records = 0, kept = 0, headers = 0, bytes = 0;
    uint64_t next_report = UINT64_C(25000000);
    const double start_time = wall_seconds();

    while ((line_length = getline(&line, &capacity, stdin)) != -1) {
        bytes += (uint64_t)line_length;
        if (line_length > 0 && line[0] == '#') {
            if (fwrite(line, 1, (size_t)line_length, stdout) != (size_t)line_length) return 1;
            headers++;
            continue;
        }

        records++;
        if (probability >= 1.0 || uniform01() < probability) {
            if (fwrite(line, 1, (size_t)line_length, stdout) != (size_t)line_length) return 1;
            kept++;
        }

        if (records == next_report) {
            double elapsed = wall_seconds() - start_time;
            fprintf(stderr,
                    "progress records=%" PRIu64 " kept=%" PRIu64
                    " input_GiB=%.2f rate=%.2f_Mrecords/s elapsed=%.1f_min\n",
                    records, kept, bytes / 1073741824.0,
                    records / elapsed / 1e6, elapsed / 60.0);
            next_report += UINT64_C(25000000);
        }
    }

    if (ferror(stdin) || fflush(stdout) != 0) {
        fprintf(stderr, "I/O error while streaming fragments\n");
        return 1;
    }

    double elapsed = wall_seconds() - start_time;
    fprintf(stderr,
            "complete records_in=%" PRIu64 " records_kept=%" PRIu64
            " observed_keep_probability=%.8f headers=%" PRIu64
            " input_GiB=%.2f elapsed=%.1f_min\n",
            records, kept, records ? (double)kept / records : 0.0,
            headers, bytes / 1073741824.0, elapsed / 60.0);
    free(line);
    return 0;
}
C_CODE

echo "[$(date '+%F %T')] Streaming fragment downsampling"
echo "input=$INPUT_PATH"
echo "output=$OUTPUT_PATH"
echo "keep_probability=$KEEP_PROB seed=$RANDOM_SEED bgzip_threads=$N_THREADS compression_level=$COMPRESSION_LEVEL"
SECONDS=0

LC_ALL=C "$BGZIP_BIN" -@ "$N_THREADS" -dc "$INPUT_PATH" \
  | "$SAMPLER_BIN" "$KEEP_PROB" "$RANDOM_SEED" 2> "$SAMPLER_STATS" \
  | LC_ALL=C "$BGZIP_BIN" -@ "$N_THREADS" -l "$COMPRESSION_LEVEL" -c \
  > "$PARTIAL_PATH"
cat "$SAMPLER_STATS" >&2

echo "[$(date '+%F %T')] Building tabix index"
LC_ALL=C "$TABIX_BIN" -f -p bed "$PARTIAL_PATH"

{
  printf 'input=%s\noutput=%s\nkeep_probability=%s\nseed=%s\nbgzip_threads=%s\ncompression_level=%s\nelapsed_seconds=%s\n' \
    "$INPUT_PATH" "$OUTPUT_PATH" "$KEEP_PROB" "$RANDOM_SEED" "$N_THREADS" "$COMPRESSION_LEVEL" "$SECONDS"
  cat "$SAMPLER_STATS"
} > "${OUTPUT_PATH}.downsample_info.txt.partial.$$"

mv -f -- "$PARTIAL_INDEX" "${OUTPUT_PATH}.tbi"
mv -f -- "$PARTIAL_PATH" "$OUTPUT_PATH"
mv -f -- "${OUTPUT_PATH}.downsample_info.txt.partial.$$" "${OUTPUT_PATH}.downsample_info.txt"

CONTIGS=$(LC_ALL=C "$TABIX_BIN" -l "$OUTPUT_PATH")
[[ -n "$CONTIGS" ]] || { echo "Output index contains no contigs." >&2; exit 1; }

echo "[$(date '+%F %T')] Done in ${SECONDS}s"
ls -lh "$OUTPUT_PATH" "${OUTPUT_PATH}.tbi" "${OUTPUT_PATH}.downsample_info.txt"
