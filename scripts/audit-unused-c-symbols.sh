#!/usr/bin/env bash
set -euo pipefail

root=${1:-src}
max_refs=${MAX_REFS:-4}

if ! command -v rg >/dev/null 2>&1; then
  echo "error: rg is required" >&2
  exit 2
fi

tmp=$(mktemp)
trap 'rm -f "$tmp"' EXIT

rg -n --glob '*.c' \
  '^[A-Za-z_][A-Za-z0-9_ *]+[[:space:]]+[A-Za-z_][A-Za-z0-9_]*[[:space:]]*\([^;]*\)[[:space:]]*\{' \
  "$root" |
awk -F: '
  {
    line=$0
    sub(/^[^:]+:[0-9]+:/, "", line)
    sub(/\(.*/, "", line)
    gsub(/^[[:space:]]+|[[:space:]]+$/, "", line)
    n=split(line, parts, /[[:space:]*]+/)
    name=parts[n]
    if (name != "main") {
      print name "\t" $1 "\t" $2
    }
  }
' | sort -u > "$tmp"

printf '%-36s %-8s %s\n' "symbol" "refs" "definition"
printf '%-36s %-8s %s\n' "------" "----" "----------"

while IFS=$'\t' read -r name file line; do
  refs=$(rg -n --glob '*.[ch]' "\\b${name}\\b" "$root" | wc -l)
  if [ "$refs" -le "$max_refs" ]; then
    printf '%-36s %-8s %s:%s\n' "$name" "$refs" "$file" "$line"
  fi
done < "$tmp"
