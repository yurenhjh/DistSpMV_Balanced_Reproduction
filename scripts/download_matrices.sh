#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p matrices
LIST=matrices/matrix_list.csv
while IFS=, read -r collection name; do
  echo "\n== $name ($collection) =="
  out_mtx="matrices/${name}.mtx"
  if [ -f "$out_mtx" ]; then echo "exists: $out_mtx"; continue; fi
  success=0
  # try possible URL patterns
  urls=(
    "https://suitesparse-collection-website.herokuapp.com/MM/${collection}/${name}.tar.gz"
    "https://suitesparse-collection-website.herokuapp.com/MM/${collection}/${name}.mtx.gz"
    "https://suitesparse-collection-website.herokuapp.com/MM/${collection}/${name}.mtx"
    "https://suitesparse-collection-website.herokuapp.com/MM/${collection}/${name}/${name}.mtx"
    "https://sparse.tamu.edu/MM/${collection}/${name}.mtx"
    "https://sparse.tamu.edu/MM/${collection}/${name}.mtx.gz"
    "https://sparse.tamu.edu/MM/${collection}/${name}.tar.gz"
    "https://sparse.tamu.edu/MM/${collection}/${name}/${name}.mtx"
    "https://sparse.tamu.edu/MM/${collection}/${name}/${name}.mtx.gz"
    "https://sparse.tamu.edu/mmfiles/${name}.mtx"
    "https://sparse.tamu.edu/mmfiles/${name}.mtx.gz"
  )
  for url in "${urls[@]}"; do
    echo "trying: $url"
    if curl -fSL --retry 2 -o "matrices/${name}.tmp" "$url" 2>/dev/null; then
      # ensure file exists and is non-empty
      if [ ! -s "matrices/${name}.tmp" ]; then
        echo "empty download from $url, skipping"
        rm -f "matrices/${name}.tmp" || true
        continue
      fi
      # detect HTML/text responses (some 404s return 200 with HTML)
      mimetype=$(file -b --mime-type "matrices/${name}.tmp" || true)
      if echo "$mimetype" | grep -qiE 'html|xml|text'; then
        echo "skipped HTML response from $url"
        rm -f "matrices/${name}.tmp" || true
        continue
      fi
      # infer suffix from URL or fallback to mime-type
      suffix=$(echo "$url" | sed -n 's/.*\(\.mtx\|\.mtx\.gz\|\.tar\.gz\)$/\1/p' || true)
      if [ -z "$suffix" ]; then
        case "$mimetype" in
          application/x-tar) suffix=".tar.gz" ;;
          application/gzip) suffix=".mtx.gz" ;;
          application/octet-stream|binary/octet-stream) suffix=".mtx" ;;
          *) suffix=".mtx" ;;
        esac
      fi
      dst="matrices/${name}${suffix}"
      mv "matrices/${name}.tmp" "$dst"
      fetched="$dst"
      echo "downloaded: $fetched"
      if [[ "$fetched" == *.mtx ]]; then mv "$fetched" "$out_mtx"; fi
      if [[ "$fetched" == *.mtx.gz ]]; then gzip -d -f "$fetched"; mv "${fetched%.gz}" "$out_mtx"; fi
      if [[ "$fetched" == *.tar.gz ]]; then tmpd=$(mktemp -d); tar -xzf "$fetched" -C "$tmpd"; found=0; for f in "$tmpd"/*.mtx "$tmpd"/*/*.mtx; do if [ -f "$f" ]; then mv "$f" "$out_mtx"; found=1; break; fi; done; rm -rf "$tmpd"; fi
      success=1; break
    else
      rm -f "matrices/${name}.tmp" || true
    fi
  done
  if [ "$success" -eq 0 ]; then echo "Failed to download $name"; fi
done < "$LIST"

echo "\nDownload finished. Listing matrices:"
ls -lh matrices | sed -n '1,200p'
