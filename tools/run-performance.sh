#!/usr/bin/env bash
set -euo pipefail
root="${1:?benchmark root}"
output="${2:?new output directory}"
[[ ! -e "$output" ]] || { echo 'Refusing to overwrite evidence'; exit 1; }
mkdir -p "$output"
trap 'echo $? > "$output/exit-status.txt"; date -u +%FT%TZ > "$output/finished-at.txt"' EXIT
date -u +%FT%TZ > "$output/started-at.txt"
echo $$ > "$output/driver.pid"
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
export CUDAVERSE_CUBLAS_PATH=/lib/x86_64-linux-gnu/libcublas.so.12
export CUDAVERSE_CUSOLVER_PATH=/lib/x86_64-linux-gnu/libcusolver.so.11
unset NVIDIA_TF32_OVERRIDE
dependencies=/home/exouser/cudaverse-validation-20261010/lib
export R_LIBS_USER="$root/library:$dependencies"
cd "$root/source"
sha256sum tools/performance*.R tools/compare-performance.R tools/run-performance.sh > "$output/scripts.sha256"
nvidia-smi -q -x > "$output/gpu-at-start.xml"
for case in legacy_1000 legacy_10000 legacy_50000 sparse_50000 pbmc_10k_v3; do
  for variant in baseline candidate; do
    if [[ "$variant" == baseline ]]; then
      export R_LIBS_USER="$root/baseline-library:$root/library:$dependencies"
    else
      export R_LIBS_USER="$root/library:$dependencies"
    fi
    for mode in timing profile; do
      echo "$case $variant $mode" > "$output/stage.txt"
      Rscript --vanilla tools/performance-worker.R "$root/data" "$case" "$variant" "$mode" \
        "$output/$case-$variant-$mode.json" > "$output/$case-$variant-$mode.log" 2>&1
    done
  done
done
nvidia-smi -q -x > "$output/gpu-at-finish.xml"
Rscript --vanilla tools/compare-performance.R "$output" > "$output/comparison.log" 2>&1
echo complete > "$output/stage.txt"
