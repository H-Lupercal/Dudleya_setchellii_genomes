#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPOSITORY_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
export PYTHONDONTWRITEBYTECODE=1
export PYTHONPYCACHEPREFIX="${REPOSITORY_ROOT}/supplementary_analysis/work/pycache"
export MPLCONFIGDIR="${REPOSITORY_ROOT}/supplementary_analysis/work/matplotlib"
export PYTHONPATH="${REPOSITORY_ROOT}/supplementary_analysis/pipeline/src:${REPOSITORY_ROOT}/canonical_publication/pipeline/src${PYTHONPATH:+:${PYTHONPATH}}"
export PATH="${REPOSITORY_ROOT}/.tools/bioconda-env/bin:${PATH}"
exec python "${REPOSITORY_ROOT}/supplementary_analysis/pipeline/scripts/run_pipeline.py" "$@"
