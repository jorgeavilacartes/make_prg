#!/usr/bin/env bash
set -eu

# configs
make_prg_URL="https://github.com/leoisl/make_prg/releases/download/v0.2.0_prototype/make_prg_0.2.0_prototype"
make_prg_executable="./make_prg_0.2.0_prototype"
make_prg_md5sum_file="./make_prg_0.2.0_prototype.md5sum.txt"

if md5sum -c "${make_prg_md5sum_file}"; then
    # The MD5 sum match
    echo "${make_prg_executable} has correct MD5 sum, proceeding..."
else
    # The MD5 sum didn't match
    echo "${make_prg_executable} does not exist or does not have correct MD5 sum, downloading..."
    wget "${make_prg_URL}" -O "${make_prg_executable}"
    chmod +x "${make_prg_executable}"
fi

echo "Building PRGs from MSAs..."
"${make_prg_executable}" from_msa --input msas/ --output_prefix msas_output/sample

echo "Updating PRGs with denovo paths..."
"${make_prg_executable}" update --update_DS msas_output/sample.update_DS --denovo_paths denovo_paths/denovo_paths.txt --output_prefix msas_updated/updated_sample

echo "All done!"
