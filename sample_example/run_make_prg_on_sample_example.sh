#!/usr/bin/env bash
set -eu

# configs
make_prg_URL="https://github.com/leoisl/make_prg/releases/download/v0.3.0/make_prg_0.3.0"
make_prg_executable="./make_prg_0.3.0"

wget "${make_prg_URL}" -O "${make_prg_executable}"
chmod +x "${make_prg_executable}"

echo "Building PRGs from MSAs..."
"${make_prg_executable}" from_msa --input msas/ --output_prefix msas_output/sample

echo "Updating PRGs with denovo paths..."
"${make_prg_executable}" update --update_DS msas_output/sample.update_DS --denovo_paths denovo_paths/denovo_paths.txt --output_prefix msas_updated/updated_sample

echo "All done!"
