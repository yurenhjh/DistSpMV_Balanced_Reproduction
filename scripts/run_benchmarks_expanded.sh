#!/usr/bin/env bash
# Expanded benchmarks: np=(1 2 4 8), OMP_NUM_THREADS=(1 2 4)
# Robust: skips non-coordinate matrices, appends results, logs errors and continues on failure.

cd "$(dirname "$0")/.."
mkdir -p results
LIST=matrices/matrix_list.csv
NP_LIST=(1 2 4 8)
THREADS_LIST=(1 2 4)
NITER=5
LOWER_BOUND=0.8
CSV=results/bench_results_expanded.csv
ERRLOG=results/bench_expanded_errors.log

# initialize CSV if not present
if [ ! -f "$CSV" ]; then
  echo "matrix,np,threads,mode,total_time,compute_time,comm_time,gflops,log" > "$CSV"
fi

get_val(){
  local key="$1"; local log="$2"
  local v; v=$(grep -E "${key}" "$log" | head -n1 | sed -E 's/.*: *([0-9.+-eE]+).*/\1/') || true
  if [ -z "$v" ]; then echo NA; else echo "$v"; fi
}

while IFS=, read -r collection name; do
  mtx="matrices/${name}.mtx"
  if [ ! -f "$mtx" ]; then echo "Missing $mtx, skipping" >> "$ERRLOG"; continue; fi
  # skip non-coordinate MatrixMarket (e.g., 'array')
  header=$(sed -n '1p' "$mtx" 2>/dev/null || echo "")
  typ=$(echo "$header" | sed -E 's/.*matrix[[:space:]]+([^ ]+).*/\1/')
  if [ "$typ" != "coordinate" ]; then echo "Skipping $name (type=$typ)" >> "$ERRLOG"; continue; fi

  for np in "${NP_LIST[@]}"; do
    for threads in "${THREADS_LIST[@]}"; do
      export OMP_NUM_THREADS="$threads"

      # naive
      if ! grep -q "^${name},${np},${threads},naive," "$CSV" 2>/dev/null; then
        log="results/${name}_naive_np${np}_t${threads}.log"
        echo "Run naive: $name np=$np threads=$threads"
        if ! mpirun -np "$np" ./dist_spmv_naive "$mtx" "$NITER" < /dev/null 2>&1 | tee "$log"; then
          echo "ERROR naive $name np=$np threads=$threads" >> "$ERRLOG"
        fi
        total=$(get_val "Total time" "$log")
        comp=$(get_val "Compute time" "$log")
        comm=$(get_val "Comm time" "$log")
        gflops=$(get_val "GFlops" "$log")
        echo "${name},${np},${threads},naive,${total},${comp},${comm},${gflops},${log}" >> "$CSV"
      else
        echo "Skipping naive for $name (already done)"
      fi

      # METIS + naive
      if ! grep -q "^${name},${np},${threads},metis," "$CSV" 2>/dev/null; then
        log="results/${name}_metis_np${np}_t${threads}.log"
        echo "Run METIS+naive: $name np=$np threads=$threads"
        if ! mpirun -np "$np" ./dist_spmv_naive "$mtx" "$NITER" METIS < /dev/null 2>&1 | tee "$log"; then
          echo "ERROR metis $name np=$np threads=$threads" >> "$ERRLOG"
        fi
        total=$(get_val "Total time" "$log")
        comp=$(get_val "Compute time" "$log")
        comm=$(get_val "Comm time" "$log")
        gflops=$(get_val "GFlops" "$log")
        echo "${name},${np},${threads},metis,${total},${comp},${comm},${gflops},${log}" >> "$CSV"
      else
        echo "Skipping metis for $name (already done)"
      fi

      # balanced
      if ! grep -q "^${name},${np},${threads},balanced," "$CSV" 2>/dev/null; then
        log="results/${name}_balanced_np${np}_t${threads}.log"
        echo "Run balanced: $name np=$np threads=$threads"
        if ! mpirun -np "$np" ./dist_spmv_balanced "$mtx" "$NITER" "$LOWER_BOUND" < /dev/null 2>&1 | tee "$log"; then
          echo "ERROR balanced $name np=$np threads=$threads" >> "$ERRLOG"
        fi
        total=$(get_val "Total time" "$log")
        comp=$(get_val "Compute time" "$log")
        comm=$(get_val "Comm time" "$log")
        gflops=$(get_val "GFlops" "$log")
        echo "${name},${np},${threads},balanced,${total},${comp},${comm},${gflops},${log}" >> "$CSV"
      else
        echo "Skipping balanced for $name (already done)"
      fi

      # small pause to avoid hammering the system
      sleep 0.05

    done
  done

done < "$LIST"

echo "Expanded benchmarks finished. Results -> $CSV ; errors -> $ERRLOG"
ls -lh results | sed -n '1,200p'
