#!/usr/bin/env bash
set -euo pipefail

# Run compilation and experiments, collect metrics into CSVs, save logs to logs/
# Usage: ./run_experiments.sh

ROOT_DIR=$(cd "$(dirname "$0")" && pwd)
cd "$ROOT_DIR"

mkdir -p logs figures data

echo "Compiling binaries..."
mpicc -O3 -fopenmp -o dist_spmv_naive dist_spmv_naive.c -lm
mpicc -O3 -fopenmp -o dist_spmv_metis_naive dist_spmv_metis_naive.c -lmetis -lm
mpicc -O3 -fopenmp -o dist_spmv_balanced dist_spmv_balanced.c -lmetis -lm

export OMP_NUM_THREADS=4

MAIN_CSV="results_main.csv"
ABL_CSV="results_ablation.csv"
BAL_CSV="results_balance.csv"

# Headers
echo "Matrix,Algorithm,np,GFlops,TotalTime,ComputeTime,CommTime,Ynorm" > "$MAIN_CSV"
echo "frac,np,GFlops,TotalTime,ComputeTime,CommTime,LocalNNZ_p0,LocalNNZ_p1,RemoteNNZ_p0,RemoteNNZ_p1" > "$ABL_CSV"
echo "Matrix,Algorithm,np,NNZ_p0,NNZ_p1,LIR" > "$BAL_CSV"

matrices=("cant.mtx" "ecology1.mtx" "bcsstk30.mtx")
algs=("naive" "metis_naive" "balanced")
# Increase iteration count to reduce measurement noise on VM
niter=100

# helpers
extract_num() {
  # Portable extraction of first numeric token (int/float/scientific) from a matched line
  local file="$1"
  local pattern="$2"
  local line
  line=$(grep -E "$pattern" "$file" | tail -n1 || true)
  if [ -z "$line" ]; then echo "NA"; return; fi
  local num
  num=$(echo "$line" | awk 'match($0, /[+-]?[0-9]+(\.[0-9]+)?([eE][+-]?[0-9]+)?/) {print substr($0, RSTART, RLENGTH); exit}')
  echo "${num:-NA}"
}

extract_all_ints() {
  # Portable extraction of all integer tokens from a matched line
  local file="$1"
  local pattern="$2"
  local line
  line=$(grep -E "$pattern" "$file" | tail -n1 || true)
  if [ -z "$line" ]; then echo ""; return; fi
  echo "$line" | awk '{
    s=$0
    while (match(s, /-?[0-9]+/)) {
      printf "%s ", substr(s, RSTART, RLENGTH)
      s = substr(s, RSTART+RLENGTH)
    }
  }'
}

# Try to drop page cache if the script has permission; do nothing otherwise.
drop_caches_if_possible() {
  if [ -w /proc/sys/vm/drop_caches ]; then
    sync
    echo 3 > /proc/sys/vm/drop_caches
  else
    echo "Warning: cannot write /proc/sys/vm/drop_caches; skipping drop_caches"
  fi
}

# Main experiment
for mat in "${matrices[@]}"; do
  for alg in "${algs[@]}"; do
    for np in 1 2; do
      echo "Running $alg $mat np=$np..."
      log="logs/${mat%.*}_${alg}_np${np}.log"
      drop_caches_if_possible
      if [ "$alg" = "naive" ]; then
        mpirun -np $np ./dist_spmv_naive matrices/$mat $niter > "$log" 2>&1
      elif [ "$alg" = "metis_naive" ]; then
        mpirun -np $np ./dist_spmv_metis_naive matrices/$mat $niter > "$log" 2>&1
      else
        mpirun -np $np ./dist_spmv_balanced matrices/$mat $niter 1.0 balanced > "$log" 2>&1
      fi

      gflops=$(extract_num "$log" "GFlops")
      total=$(extract_num "$log" "Total time|总时间")
      compute=$(extract_num "$log" "Compute time|计算时间")
      comm=$(extract_num "$log" "Comm time|通信时间")
      ynorm=$(extract_num "$log" "Y-norm")

      echo "$mat,$alg,$np,$gflops,$total,$compute,$comm,$ynorm" >> "$MAIN_CSV"
      echo "Done."
    done
  done
done

# Ablation (cant, np=2)
frac_list=(0.5 0.6 0.7 0.8 0.9 1.0 1.2)
for frac in "${frac_list[@]}"; do
  log="logs/cant_ablation_frac${frac}_np2.log"
  echo "Running balanced cant np=2 frac=$frac..."
  drop_caches_if_possible
  mpirun -np 2 ./dist_spmv_balanced matrices/cant.mtx $niter $frac balanced > "$log" 2>&1
  gflops=$(extract_num "$log" "GFlops")
  total=$(extract_num "$log" "Total time|总时间")
  compute=$(extract_num "$log" "Compute time|计算时间")
  comm=$(extract_num "$log" "Comm time|通信时间")
  local_vals=( $(extract_all_ints "$log" "Per-process local_nnz") )
  remote_vals=( $(extract_all_ints "$log" "Per-process remote_nnz") )
  local_p0=${local_vals[0]:-NA}
  local_p1=${local_vals[1]:-NA}
  remote_p0=${remote_vals[0]:-NA}
  remote_p1=${remote_vals[1]:-NA}
  echo "$frac,2,$gflops,$total,$compute,$comm,$local_p0,$local_p1,$remote_p0,$remote_p1" >> "$ABL_CSV"
  echo "Done."
done

# Load-balance experiment (np=2)
for mat in "${matrices[@]}"; do
  for alg in "${algs[@]}"; do
    np=2
    log="logs/${mat%.*}_${alg}_np${np}_balance.log"
    echo "Running balance check $alg $mat np=$np..."
    drop_caches_if_possible
    if [ "$alg" = "naive" ]; then
      mpirun -np $np ./dist_spmv_naive matrices/$mat $niter > "$log" 2>&1
      nums=( $(extract_all_ints "$log" "Per-process nnz") )
      p0=${nums[0]:-NA}
      p1=${nums[1]:-NA}
    elif [ "$alg" = "metis_naive" ]; then
      mpirun -np $np ./dist_spmv_metis_naive matrices/$mat $niter > "$log" 2>&1
      nums=( $(extract_all_ints "$log" "Per-process nnz") )
      p0=${nums[0]:-NA}
      p1=${nums[1]:-NA}
    else
      mpirun -np $np ./dist_spmv_balanced matrices/$mat $niter 1.0 balanced > "$log" 2>&1
      local_vals=( $(extract_all_ints "$log" "Per-process local_nnz") )
      remote_vals=( $(extract_all_ints "$log" "Per-process remote_nnz") )
      # sum to get total nnz per process
      p0=$(( ${local_vals[0]:-0} + ${remote_vals[0]:-0} ))
      p1=$(( ${local_vals[1]:-0} + ${remote_vals[1]:-0} ))
      # also append remote nnz info to balance CSV via comment in log parsing
      echo "# RemoteNNZ: p0=${remote_vals[0]:-0} p1=${remote_vals[1]:-0}" >> "$log"
    fi
    if [ "$p0" = "NA" ] || [ "$p1" = "NA" ]; then
      lir="NA"
    else
      lir=$(awk -v a="$p0" -v b="$p1" 'BEGIN{avg=(a+b)/2; if(avg==0){print "NA"; exit} m=(a>b?a:b); printf("%.6f", m/avg - 1)}')
    fi
    echo "$mat,$alg,2,$p0,$p1,$lir" >> "$BAL_CSV"
    echo "Done."
  done
done

echo "All experiments finished."
