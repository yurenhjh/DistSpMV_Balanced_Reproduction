#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p results
LIST=matrices/matrix_list.csv
NP_LIST=(2)
THREADS_LIST=(1)
NITER=5
LOWER_BOUND=0.8
CSV=results/bench_results.csv
echo "matrix,np,threads,mode,total_time,compute_time,comm_time,gflops,log" > "$CSV"

get_val(){
  local key="$1"; local log="$2"
  local v; v=$(grep -E "${key}" "$log" | head -n1 | sed -E 's/.*: *([0-9.+-eE]+).*/\1/') || true
  if [ -z "$v" ]; then echo NA; else echo "$v"; fi
}

while IFS=, read -r collection name; do
  mtx="matrices/${name}.mtx"
  if [ ! -f "$mtx" ]; then echo "Missing $mtx, skipping"; continue; fi
  for np in "${NP_LIST[@]}"; do
    for threads in "${THREADS_LIST[@]}"; do
      export OMP_NUM_THREADS="$threads"
      # baseline (naive)
      log="results/${name}_naive_np${np}_t${threads}.log"
      echo "Running naive: $name np=$np threads=$threads"
      mpirun -np "$np" ./dist_spmv_naive "$mtx" "$NITER" < /dev/null 2>&1 | tee "$log"
      total=$(get_val "Total time" "$log")
      comp=$(get_val "Compute time" "$log")
      comm=$(get_val "Comm time" "$log")
      gflops=$(get_val "GFlops" "$log")
      echo "${name},${np},${threads},naive,${total},${comp},${comm},${gflops},${log}" >> "$CSV"

      # METIS reorder + naive
      log="results/${name}_metis_np${np}_t${threads}.log"
      echo "Running METIS+naive: $name np=$np threads=$threads"
      mpirun -np "$np" ./dist_spmv_naive "$mtx" "$NITER" METIS < /dev/null 2>&1 | tee "$log"
      total=$(get_val "Total time" "$log")
      comp=$(get_val "Compute time" "$log")
      comm=$(get_val "Comm time" "$log")
      gflops=$(get_val "GFlops" "$log")
      echo "${name},${np},${threads},metis,${total},${comp},${comm},${gflops},${log}" >> "$CSV"

      # balanced
      log="results/${name}_balanced_np${np}_t${threads}.log"
      echo "Running balanced: $name np=$np threads=$threads"
      mpirun -np "$np" ./dist_spmv_balanced "$mtx" "$NITER" "$LOWER_BOUND" < /dev/null 2>&1 | tee "$log"
      total=$(get_val "Total time" "$log")
      comp=$(get_val "Compute time" "$log")
      comm=$(get_val "Comm time" "$log")
      gflops=$(get_val "GFlops" "$log")
      echo "${name},${np},${threads},balanced,${total},${comp},${comm},${gflops},${log}" >> "$CSV"

    done
  done
done < "$LIST"

echo "Benchmarks finished. Results -> $CSV"
ls -lh results | sed -n '1,200p'
