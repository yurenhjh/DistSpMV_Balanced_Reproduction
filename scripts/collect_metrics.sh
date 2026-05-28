#!/usr/bin/env bash
set -euo pipefail

mkdir -p data
out=data/results_perf.csv
echo "Matrix,np,mode,lower_bound_frac,total_s,compute_s,comm_s,GFlops,Y-norm" > "$out"

for f in logs/*.log; do
  [ -f "$f" ] || continue
  base=$(basename "$f")
  base=${base%.log}
  if [[ $base =~ ^(.+)_([bB]alanced|[nN]aive)_np([0-9]+)$ ]]; then
    matrix="${BASH_REMATCH[1]}"
    mode="${BASH_REMATCH[2]}"
    mode=$(echo "$mode" | tr '[:upper:]' '[:lower:]')
    np="${BASH_REMATCH[3]}"
  else
    # fallback: try to extract mode/np from file contents
    matrix="$base"
    mode=""
    np=""
  fi

  lower=$(grep -E "lower_bound" "$f" 2>/dev/null | awk '{print $NF}' | head -n1 || echo "")
  total=$(grep -E "Total time" "$f" 2>/dev/null | awk '{print $(NF-1)}' | head -n1 || echo "")
  compute=$(grep -E "Compute time" "$f" 2>/dev/null | awk '{print $(NF-1)}' | head -n1 || echo "")
  comm=$(grep -E "Comm time" "$f" 2>/dev/null | awk '{print $(NF-1)}' | head -n1 || echo "")
  gflops=$(grep -E "^GFlops" "$f" 2>/dev/null | awk '{print $NF}' | head -n1 || echo "")
  ynorm=$(grep -E "Y-norm" "$f" 2>/dev/null | awk '{print $NF}' | head -n1 || echo "")

  printf "%s,%s,%s,%s,%s,%s,%s,%s,%s\n" "$matrix" "$np" "$mode" "$lower" "$total" "$compute" "$comm" "$gflops" "$ynorm" >> "$out"
done

echo "Wrote $out"
