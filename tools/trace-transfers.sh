#!/usr/bin/env bash
set -euo pipefail
root="${1:?benchmark root}"
output="${2:?new trace directory}"
[[ ! -e "$output" ]] || { echo 'Refusing to overwrite evidence'; exit 1; }
mkdir -p "$output"
trap 'echo $? > "$output/exit-status.txt"' EXIT
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
export CUDAVERSE_CUBLAS_PATH=/lib/x86_64-linux-gnu/libcublas.so.12
export CUDAVERSE_CUSOLVER_PATH=/lib/x86_64-linux-gnu/libcusolver.so.11
unset NVIDIA_TF32_OVERRIDE
dependencies=/home/exouser/cudaverse-validation-20261010/lib
candidate_library="${CUDACELLR_CANDIDATE_LIBRARY:-$root/library}"
cd "${CUDACELLR_BENCHMARK_SOURCE:-$root/source}"
nsys --version > "$output/nsys-version.txt"
sha256sum tools/performance-worker.R tools/trace-transfers.sh > "$output/scripts.sha256"
for case in sparse_50000 pbmc_10k_v3; do
  for variant in baseline candidate; do
    if [[ "$variant" == baseline ]]; then
      export R_LIBS_USER="$root/baseline-library:$root/library:$dependencies"
    else
      export R_LIBS_USER="$candidate_library:$root/library:$dependencies"
    fi
    name="$output/$case-$variant"
    echo "$case $variant" > "$output/stage.txt"
    nsys profile --trace=cuda --sample=none --cpuctxsw=none --force-overwrite=false \
      --output="$name" Rscript --vanilla tools/performance-worker.R \
      "$root/data" "$case" "$variant" profile "$name.json" > "$name.log" 2>&1
    nsys stats --report gpumemsizesum,gpumemtimesum,cudaapisum,gpukernsum \
      --format csv "$name.nsys-rep" > "$name-stats.csv" 2> "$name-stats.log"
  done
done
echo complete > "$output/stage.txt"
