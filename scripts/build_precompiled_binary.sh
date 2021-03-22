#!/usr/bin/env bash

version="0.2.0_prototype"

# To run: scripts/build_precompiled_binary.sh
# Precompiled binary will be on the dist/ folder
pyinstaller -F --hidden-import="sklearn.utils._weight_vector" make_prg/__main__.py
mv dist/__main__ "dist/make_prg_${version}"
chmod +x "dist/make_prg_${version}"