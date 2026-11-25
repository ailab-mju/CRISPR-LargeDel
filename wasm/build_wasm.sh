#!/bin/bash
# Build this core as the website's WebAssembly module: out/kmer_align.js + out/kmer_align.wasm,
# the two files the site's docs/ folder loads (docs/worker.js calls run_analysis()).
# Same sources as the native build; wasm_entry.cpp replaces kmer_align.cpp (the CLI main).
# Flags as /home/seb/web/build.sh, which built the live site's core. Needs emsdk (3.1.18).
set -e
cd "$(dirname "$0")"
command -v emcc >/dev/null || source "$HOME/emsdk/emsdk_env.sh" >/dev/null 2>&1
mkdir -p out
# no -DUSE_HTSLIB: the browser build writes SAM, which docs/viewer.js reads. samout.cpp still sorts.
SRC="../align.cpp ../align_wrap.cpp ../arg_parse.cpp ../func.cpp ../post.cpp ../samout.cpp ../report.cpp wasm_entry.cpp"
emcc -O3 -std=gnu++14 -I.. $SRC -o out/kmer_align.js \
  -s MODULARIZE=1 \
  -s EXPORT_NAME=CrisprLargeDel \
  -s ENVIRONMENT=web,worker \
  -s INVOKE_RUN=0 \
  -s EXIT_RUNTIME=0 \
  -s ALLOW_MEMORY_GROWTH=1 \
  -s INITIAL_MEMORY=64MB \
  -s TOTAL_STACK=8MB \
  -s FORCE_FILESYSTEM=1 \
  -s EXPORTED_RUNTIME_METHODS='["ccall","FS"]' \
  -s EXPORTED_FUNCTIONS='["_run_analysis","_malloc","_free"]' \
  -s ASSERTIONS=0
echo "built: $(pwd)/out/kmer_align.js + kmer_align.wasm"
