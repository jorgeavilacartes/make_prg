#!/usr/bin/env bash
# run from project root: scripts/build_precompiled_binary/build_precompiled_binary.sh

set -eu
CURRENT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SCRIPTS_DIR="$(dirname "${CURRENT_DIR}")"
MAKE_PRG_DIR="$(dirname "${SCRIPTS_DIR}")"
PORTABLE_EXECUTABLE_BUILD_DIR="${MAKE_PRG_DIR}/precompiled_binary"

cd "$MAKE_PRG_DIR"

if [ -d "${PORTABLE_EXECUTABLE_BUILD_DIR}" ]; then
  echo "Please remove ${PORTABLE_EXECUTABLE_BUILD_DIR} before proceeding."
  exit 1
fi

version="0.2.0_prototype"
sudo docker run --rm \
  -v "$(pwd)":/src \
  -e PYTHONHASHSEED=42 \
  leandroishilima/make_prg_precompiled_binary_builder:0.0.1 \
  --noconfirm \
  --onefile \
  --log-level DEBUG \
  --clean \
  --distpath /src/precompiled_binary/dist \
  --workpath /src/precompiled_binary/build \
  --specpath /src/precompiled_binary/spec \
  --name make_prg_${version} \
  --hidden-import=\"sklearn.utils._weight_vector\" \
  make_prg/__main__.py
