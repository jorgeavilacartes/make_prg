#!/usr/bin/env bash
set -eu

# configs
version="1.0.0"
make_prg_URL="docker://leandroishilima/make_prg:${version}"

if [ ! -f "make_prg.sif" ]; then
  echo "make_prg image not found, pulling it..."
  singularity pull --name make_prg.sif "${make_prg_URL}"
fi

echo "Building PRGs from MSAs..."
from_msa_command_line="singularity exec make_prg.sif make_prg from_msa --input msas/ --output_prefix msas_output/sample"
echo "Running ${from_msa_command_line}"
${from_msa_command_line}

echo "Updating PRGs with denovo paths..."
update_command_line="singularity exec make_prg.sif make_prg update --update_DS msas_output/sample.update_DS.zip --denovo_paths denovo_paths/denovo_paths.txt --output_prefix msas_updated/updated_sample"
echo "Running ${update_command_line}"
${update_command_line}

echo "All done!"
