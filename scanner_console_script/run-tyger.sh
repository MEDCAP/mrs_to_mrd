#!/usr/bin/env bash
#
# Run the conversion -> shift -> recon pipeline on Tyger, over every scan directory under a
# target directory. One tar per scan goes up, one _recon.mrd2 comes back.
#
#   run-tyger.sh /data/kidney_2026_09
#   run-tyger.sh /data/kidney_2026_09 -- -urea 2.3 -pyr_s 9.7 -lac_m 21.8
#   run-tyger.sh /data/kidney_2026_09 --method alloc --method roll
#
# Everything after `--` is passed to mrd2recon.py. Those flags are the experiment's chemistry,
# not a default, so set them with the study. With none given, the peak set committed in
# tyger_deploy/recon_codespec.yml is used unchanged.
#
# Why the codespecs are rendered rather than passed through: `tyger run exec` overrides only
# --replicas, --node-pool, --buffer, --cluster, --timeout and --tag. There is no flag for a
# codespec's args or env, and mrd2recon.py reads its peaks from argv only - it consults no
# environment variable. So a per-study peak set can reach the container only by being in the YAML
# handed to -f, which is what render_spec builds between the ARGS_BEGIN/ARGS_END markers.
#
# Plots stay local: there is no plot codespec, and the recon image ships no matplotlib.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SPEC_DIR="${REPO_ROOT}/tyger_deploy"
PYTHON="${PYTHON:-python}"

die() { echo "Error: $*" >&2; exit 1; }

usage() {
  cat >&2 <<'EOF'
Usage: run-tyger.sh <target_directory> [--method METHOD]... [--no-plot] [-- recon flags...]

  <target_directory>  holds one subdirectory per scan; each is tarred and sent up
  --method METHOD     shift method: shift, regrid, contiguous, roll, alloc.
                      Repeat to sweep; each method costs a full recon. Default: regrid
  --no-plot           skip the local mrdplot.py figures
  --dry-run           print the rendered codespecs and exit, without touching the cluster
  -- recon flags      passed to mrd2recon.py, e.g. -pyr_s 9.7 -lac_m 21.8
EOF
  exit 1
}

# ---------- argument handling --------------------------------------------

TARGET_DIR=""
METHODS=()
RECON_ARGS=()
PLOT=1
DRY_RUN=0

while [ "$#" -gt 0 ]; do
  case "$1" in
    --method)   [ "$#" -ge 2 ] || die "--method needs a value"; METHODS+=("$2"); shift 2 ;;
    --no-plot)  PLOT=0; shift ;;
    --dry-run)  DRY_RUN=1; shift ;;
    -h|--help)  usage ;;
    --)         shift; RECON_ARGS=("$@"); break ;;
    -*)         die "unknown option $1 (recon flags go after --)" ;;
    *)          [ -z "$TARGET_DIR" ] || die "only one target directory"; TARGET_DIR="$1"; shift ;;
  esac
done

[ "${#METHODS[@]}" -gt 0 ] || METHODS=(regrid)
if [ "$DRY_RUN" -eq 0 ]; then
  [ -n "$TARGET_DIR" ] || usage
  [ -d "$TARGET_DIR" ] || die "directory '$TARGET_DIR' does not exist"
  TARGET_DIR="$(cd "$TARGET_DIR" && pwd)"
fi
OUTPUT_DIR="${TARGET_DIR}/png"

for spec in convert_codespec.yml shift_codespec.yml recon_codespec.yml; do
  [ -f "${SPEC_DIR}/${spec}" ] || die "missing codespec ${SPEC_DIR}/${spec}"
done
[ "$DRY_RUN" -eq 1 ] || command -v tyger >/dev/null || die "tyger is not on PATH"

# ---------- codespec rendering -------------------------------------------

# Emit flags as codespec args: one YAML list item per token, every value quoted.
#
# Both halves matter. split_peak_args (mrd2recon.py) reads a peak's value from argv[i+1], so a
# flag and its value folded into one item ('- -pyr_s 9.7') arrives as a single argv element with
# no value after it, falls through to argparse, and dies as "unrecognized arguments". And Tyger's
# args are []string, so an unquoted 9.7 is a YAML float that fails validation before the job is
# even created - as does an unquoted -0.4, which YAML would not read as a string either way.
yaml_args() {
  local indent="      "
  while [ "$#" -gt 0 ]; do
    case "$1" in
      -i|--input|-o|--output|-f|--folder)
        die "$1 is set by the codespec's buffer wiring; it cannot be overridden" ;;
      -*) ;;
      *)  die "expected a flag, got '$1'" ;;
    esac
    [ "$#" -ge 2 ] || die "$1 has no value"
    printf '%s- %s\n%s- "%s"\n' "$indent" "$1" "$indent" "$2"
    shift 2
  done
}

# Every recon flag takes a number: a peak offset, or one of -lb/-df/-dw/-dph/--fidpad. Checking
# it here is what turns a typo into a one-second local failure rather than a scheduled pod that
# exits 2 - split_peak_args only treats a token as a peak if a float follows it, so '-pyr_s x'
# silently becomes an argparse error instead of a peak.
require_numeric_values() {
  # the flag shape, the reserved flags and the pairing are yaml_args' checks; reuse them so a
  # misplaced token is reported as what it is rather than as a non-numeric value
  yaml_args "$@" > /dev/null
  while [ "$#" -gt 0 ]; do
    case "$2" in
      ''|*[!0-9eE.+-]*) die "value for $1 must be a number, got '$2'" ;;
    esac
    shift 2
  done
}

# Replace the ARGS_BEGIN/ARGS_END region of a codespec, leaving the rest of the YAML alone.
# With no flags the committed block is kept, so the rendered spec and the committed one agree.
render_spec() {
  local template="$1"; shift
  if [ "$#" -eq 0 ]; then
    cat "$template"
    return
  fi
  # through a file, not -v: awk rejects a newline inside a -v assignment
  local block_file
  block_file="$(mktemp)"
  yaml_args "$@" > "$block_file"
  awk -v blockfile="$block_file" '
    /# ARGS_BEGIN/ {
      print
      while ((getline line < blockfile) > 0) print line
      close(blockfile)
      skip = 1
      next
    }
    /# ARGS_END/ { skip = 0 }
    !skip        { print }
  ' "$template"
  rm -f "$block_file"
}

# ---------- run ----------------------------------------------------------

if [ "${#RECON_ARGS[@]}" -gt 0 ]; then
  echo "Recon flags: ${RECON_ARGS[*]}" >&2
else
  echo "Recon flags: the peak set committed in recon_codespec.yml" >&2
fi

# render once, up front: a bad flag should fail before anything is tarred or uploaded
WORK_DIR="$(mktemp -d)"
trap 'rm -rf "$WORK_DIR"' EXIT

# bash 3.2 is what macOS ships, and there "${ARR[@]}" on an empty array trips set -u
recon_flags=(${RECON_ARGS[@]+"${RECON_ARGS[@]}"})
if [ "${#recon_flags[@]}" -gt 0 ]; then
  require_numeric_values "${recon_flags[@]}"
  render_spec "${SPEC_DIR}/recon_codespec.yml" "${recon_flags[@]}" > "${WORK_DIR}/recon.yml"
else
  render_spec "${SPEC_DIR}/recon_codespec.yml" > "${WORK_DIR}/recon.yml"
fi

for m in "${METHODS[@]}"; do
  case "$m" in
    shift|regrid|contiguous|roll|alloc) ;;
    *) die "unknown shift method '$m' (METHODS in mrd2shift.py)" ;;
  esac
  render_spec "${SPEC_DIR}/shift_codespec.yml" --method "$m" > "${WORK_DIR}/shift_${m}.yml"
done

if [ "$DRY_RUN" -eq 1 ]; then
  for rendered in "${WORK_DIR}"/*.yml; do
    echo "===== $(basename "$rendered") ====="
    cat "$rendered"
  done
  exit 0
fi

if [ "$PLOT" -eq 1 ]; then
  mkdir -p "$OUTPUT_DIR"
fi

# One scan directory, end to end. Every stage is checked explicitly rather than left to set -e,
# because errexit is disabled inside a function called from an `if`, and a partial output is
# removed so a later stage cannot read a truncated stream and a rerun is not fooled by a file
# that exists but is half written.
run_scan() {
  local base_name="$1"
  local raw_path="${TARGET_DIR}/${base_name}_raw.mrd2"
  local failed=0

  echo "Converting ${base_name}" >&2
  if ! tar cf - -C "$TARGET_DIR" "$base_name" \
       | tyger run exec -f "${SPEC_DIR}/convert_codespec.yml" > "$raw_path"; then
    rm -f "$raw_path"
    echo "  conversion failed" >&2
    return 1
  fi

  local m shift_path recon_path
  for m in "${METHODS[@]}"; do
    shift_path="${TARGET_DIR}/${base_name}_${m}.mrd2"
    recon_path="${TARGET_DIR}/${base_name}_${m}_recon.mrd2"

    echo "Shifting ${base_name} (${m})" >&2
    if ! tyger run exec -f "${WORK_DIR}/shift_${m}.yml" < "$raw_path" > "$shift_path"; then
      rm -f "$shift_path"
      echo "  shift failed for method ${m}" >&2
      failed=1
      continue
    fi

    echo "Reconstructing ${base_name} (${m})" >&2
    if ! tyger run exec -f "${WORK_DIR}/recon.yml" < "$shift_path" > "$recon_path"; then
      rm -f "$recon_path"
      echo "  recon failed for method ${m}" >&2
      failed=1
      continue
    fi

    # a missing figure is not a reason to call the reconstruction failed
    if [ "$PLOT" -eq 1 ]; then
      MPLBACKEND=Agg "$PYTHON" "${REPO_ROOT}/mrdplot.py" -i "$recon_path" -s "$OUTPUT_DIR" \
        || echo "  plotting failed for method ${m}; ${recon_path} is still good" >&2
    fi
  done

  return "$failed"
}

shopt -s nullglob
found_subdirs=0
FAILED=()

for dir in "$TARGET_DIR"/*/; do
  dir="${dir%/}"
  [ -d "$dir" ] || continue
  [ "$dir" = "$OUTPUT_DIR" ] && continue

  found_subdirs=1
  base_name="$(basename "$dir")"

  # one bad scan directory should not cost the rest of an overnight batch
  if run_scan "$base_name"; then
    echo "Finished ${base_name}" >&2
  else
    echo "Skipping ${base_name}: a pipeline stage failed" >&2
    FAILED+=("$base_name")
  fi
done

if [ "$found_subdirs" -eq 0 ]; then
  die "no subdirectories found in '$TARGET_DIR'"
fi

if [ "${#FAILED[@]}" -gt 0 ]; then
  echo "Done, with ${#FAILED[@]} scan(s) skipped: ${FAILED[*]}" >&2
  exit 1
fi
echo "Done. All subdirectories processed." >&2
