#!/usr/bin/env bash

# 1. build container (see Dockerfile)
# 2. run from project root: scripts/build_precompiled_binary/build_precompiled_binary.sh
version="0.2.0_prototype"
sudo docker run --rm -v "$(pwd)":/make_prg make_prg_precompiled_binary_builder:0.0.1 /bin/bash -c "cd make_prg && pyinstaller \
  --distpath precompiled_binary/dist --workpath precompiled_binary/build --noconfirm --clean --onefile \
  --specpath precompiled_binary/spec --name make_prg_${version} --hidden-import=\"sklearn.utils._weight_vector\" \
  make_prg/__main__.py"
