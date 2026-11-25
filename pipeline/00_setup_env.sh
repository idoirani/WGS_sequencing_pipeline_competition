#!/usr/bin/env bash
# Create the analysis environment from environment.yml, which pins the exact
# versions that produced the published results.  All tools are native arm64 on
# Apple Silicon -- no Rosetta, unlike the comparative-genomics pipeline in this
# project (geNomad/Roary/Pharokka), which needed osx-64 emulation.
set -euo pipefail
cd "$(dirname "$0")/.."
ENV_NAME="${1:-phage-evo}"
if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo "environment '$ENV_NAME' already exists; updating"
  conda env update -n "$ENV_NAME" -f environment.yml --prune
else
  conda env create -n "$ENV_NAME" -f environment.yml
fi
echo "activate with: conda activate $ENV_NAME"
