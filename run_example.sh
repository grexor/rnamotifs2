#!/usr/bin/env bash
#
# Run a single RNAmotifs2 analysis end to end.
#
#   ./run_example.sh                 # runs comps/paper.bh
#   ./run_example.sh my_comparison   # runs comps/my_comparison
#
# One-time environment setup:
#   micromamba create -n rnamotifs2 -c conda-forge -c bioconda python=3.12 \
#       numpy scipy matplotlib pandas pysam psutil beautifulsoup4 requests pip
#   micromamba run -n rnamotifs2 pip install fisher pybio
#   micromamba run -n rnamotifs2 pip install -e .   # from the repo root
#
set -euo pipefail

COMPS="${1:-paper.bh}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="${RNAMOTIFS2_ENV:-rnamotifs2}"

# comps/ is resolved relative to the working directory (or $RNAMOTIFS2_COMPS)
cd "$REPO"
exec micromamba run -n "$ENV_NAME" rnamotifs2 run -comps "$COMPS"
