#!/usr/bin/env bash
set -euo pipefail

EXE=""
MERCURYDPM_BUILD_DIR="${MERCURYDPM_BUILD_DIR:-}"

log() {
  echo "[run.sh]" "$*" >&2
}

usage() {
  cat >&2 <<'EOF'
Usage:
  script.sh [--exec=/path/to/ChannelTransport] [--build-dir=/path/to/mercurydpm/build]

Options:
  -e=, --exec=        Explicit executable path to run or
  -b=, --build-dir=   MercuryDPM build dir (root). Script will look for:
                      <build-dir>/Drivers/PreCICE/ChannelTransport
EOF
}

# go through CLI
for i in "$@"; do
  case "$i" in
    -e=*|--exec=*)
      EXE="${i#*=}"
      shift
      ;;
    -b=*|--build-dir=*)
      MERCURYDPM_BUILD_DIR="${i#*=}"
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      ;;
  esac
done

# Helper to run and exit
run_and_exit() {
  local cmd="$1"
  log "Using executable: ${cmd}"
  "${cmd}"
  exit 0
}

# 1) Explicit executable
if [[ -n "${EXE}" ]]; then
  log "Explicit --exec provided: ${EXE}"
  if [[ -x "${EXE}" ]]; then
    run_and_exit "${EXE}"
  else
    log "Explicit executable was not found or is not an executable: ${EXE}"
  fi
fi

EXE="ChannelTransport"

# 2) MercuryDPM build dir
if [[ -n "${MERCURYDPM_BUILD_DIR}" ]]; then
  # location in the build dir:
  # /path/to/mercurydpm/build/Drivers/PreCICE/ChannelTransport
  build_candidate="${MERCURYDPM_BUILD_DIR%/}/Drivers/PreCICE/${EXE}"
  log "Checking MERCURYDPM_BUILD_DIR: ${MERCURYDPM_BUILD_DIR}"
  log "Looking for build candidate: ${build_candidate}"
  if [[ -x "${build_candidate}" ]]; then
    run_and_exit "${build_candidate}"
  fi
fi

# 3) Global PATH
log "Checking global PATH for: ${EXE}"
if command -v "${EXE}" >/dev/null 2>&1; then
  path_candidate="$(command -v "${EXE}")"
  log "Found in PATH: ${path_candidate}"
  run_and_exit "${EXE}"
fi

# 4) Local directory
log "Checking local directory: ./${EXE}"
if [[ -x "./${EXE}" ]]; then
  run_and_exit "./${EXE}"
fi

# Not found
cat >&2 <<EOF
Unable to find the executable ${EXE}.
Searched (in order):
  1) --exec=... (explicit path)
  2) MERCURYDPM_BUILD_DIR (or --build-dir=...): <build>/Drivers/PreCICE/${EXE}
  3) PATH: command -v ${EXE}
  4) Local: ./${EXE}

Hints:
  - Specify explicitly: --exec=/path/to/${EXE}
  - Or export MERCURYDPM_BUILD_DIR=/path/to/mercurydpm/build
  - Or make it discoverable: export PATH=...:\$PATH
  - MercuryDPM must be built with preCICE coupling:
      -D MercuryDPM_PreCICE_COUPLING="ON"
EOF
exit 1
