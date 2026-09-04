#!/usr/bin/env bash
#
# Run a single RNAmotifs2 analysis end to end.
#
#   ./run_example.sh                 # runs comps/paper.bh
#   ./run_example.sh my_comparison   # runs comps/my_comparison
#
# Requires the `rnamotifs2` micromamba environment:
#   micromamba create -n rnamotifs2 -c conda-forge -c bioconda python=3.12 \
#       numpy scipy matplotlib pandas pysam psutil beautifulsoup4 requests pip
#   micromamba run -n rnamotifs2 pip install fisher pyliftover /home/gregor/pybio
#
set -euo pipefail

COMPS="${1:-paper.bh}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="${RNAMOTIFS2_ENV:-rnamotifs2}"

# rnamotifs2 is imported as a top-level package living at $REPO, so its parent
# must be on PYTHONPATH; the bin/ helper scripts must be on PATH because the
# pipeline shells out to them (rnamotifs2.motif, rnamotifs2.draw, ...).
export PYTHONPATH="$(dirname "$REPO")${PYTHONPATH:+:$PYTHONPATH}"
export PATH="$REPO/bin:$PATH"

cd "$REPO"
exec micromamba run -n "$ENV_NAME" python bin/rnamotifs2 -comps "$COMPS"
