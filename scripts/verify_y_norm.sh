#!/usr/bin/env bash
set -euo pipefail

mkdir -p data
out=data/results_y_norm.csv
echo "Matrix,np,y_balanced,y_naive,rel_diff" > "$out"

declare -A bal naive

for f in logs/*.log; do
  [ -f "$f" ] || continue
  base=$(basename "${f%.log}")
  if [[ $base =~ ^(.+)_([bB]alanced|[nN]aive)_np([0-9]+)$ ]]; then
    matrix="${BASH_REMATCH[1]}"
    mode="${BASH_REMATCH[2]}"
    np="${BASH_REMATCH[3]}"
    key="${matrix}_np${np}"
    y=$(grep -E "Y-norm" "$f" 2>/dev/null | awk '{print $NF}' | head -n1 || echo "")
    if [[ -z "$y" ]]; then y=""; fi
    if [[ "${mode,,}" == "balanced" ]]; then
      bal["$key"]="$y"
    else
      naive["$key"]="$y"
    fi
  fi
done

for key in "${!bal[@]}"; do
  matrix_np="$key"
  yb="${bal[$key]}"
  yn="${naive[$key]:-}"
  if [[ -z "$yb" || -z "$yn" ]]; then
    rel="NA"
  else
    rel=$(awk -v a="$yb" -v b="$yn" 'BEGIN{a=a+0; b=b+0; if(a==b){print "0"; exit} da=(a-b); if(da<0) da=-da; aa=a; if(aa<0) aa=-aa; bb=b; if(bb<0) bb=-bb; maxv=(aa>bb?aa:bb); if(maxv==0){print "0"; exit} printf("%.3e", da/maxv)}')
  fi
  matrix=${matrix_np%_np*}
  np=${matrix_np##*_np}
  printf "%s,%s,%s,%s,%s\n" "$matrix" "$np" "$yb" "$yn" "$rel" >> "$out"
done

echo "Wrote $out"
