#!/usr/bin/env bash
set -euo pipefail

help() {
    cat <<EOF
Usage: ${0##*/} FILE [SIZE]

Keep the first and last SIZE/2 bytes if FILE exceeds SIZE.
Insert a blank-line-separated ...[TRUNCATED]... marker between them.
The marker's bytes are additional to SIZE.

SIZE defaults to 2MB (2097152 bytes). Use a nonnegative integer with:
  no suffix or B: bytes
  K or KB: kilobytes (1024 bytes)
  M or MB: megabytes (1024*1024 bytes)
  G or GB: gigabytes (1024*1024*1024 bytes)
Suffixes are case-insensitive. A space before the suffix is allowed.
Examples: ${0##*/} file.log 100K
          ${0##*/} file.log "100 KB"
          ${0##*/} file.log 2MB
          ${0##*/} file.log 2 GB

  -h, --help  Show this help.
EOF
}

if (( $# == 1 )) && [[ "$1" == -h || "$1" == --help ]]; then
    help
    exit 0
fi
if (( $# < 1 || $# > 3 )); then
    help >&2
    exit 1
fi

file=$1
size_arg=${2-2MB}
if (( $# == 3 )); then
    size_arg="$2 $3"
fi

if [[ ! -f "$file" ]]; then
    printf 'Error: file does not exist or is not a regular file: %s\n' "$file" >&2
    exit 1
fi
size_pattern='^([0-9]+)[[:space:]]*([bB]|[kKmMgG][bB]?)?$'
if [[ ! "$size_arg" =~ $size_pattern ]]; then
    printf 'Error: invalid size: %s. Use --help for supported formats.\n' "$size_arg" >&2
    exit 1
fi
max_size=${BASH_REMATCH[1]}
suffix=${BASH_REMATCH[2]^^}
case "$suffix" in
    K|KB) multiplier=1024 ;;
    M|MB) multiplier=1048576 ;;
    G|GB) multiplier=1073741824 ;;
    *) multiplier=1 ;;
esac

# Remove leading zeros so Bash treats the value as decimal.
max_size=${max_size#"${max_size%%[!0]*}"}
max_size=${max_size:-0}
size=$(stat -c %s -- "$file")

# Compare before multiplication to avoid overflow for arbitrarily large limits.
threshold=$((size / multiplier))
if (( ${#max_size} > ${#threshold} )) ||
   { (( ${#max_size} == ${#threshold} )) && [[ "$max_size" > "$threshold" ]]; }; then
    exit 0
fi
max_size=$((max_size * multiplier))
if (( size <= max_size )); then
    exit 0
fi

half=$((max_size / 2))
tmp=$(mktemp -- "${file}.tmp.XXXXXXXX")
trap 'rm -f -- "$tmp"' EXIT

head -c "$half" -- "$file" > "$tmp"
printf '\n\n...[TRUNCATED]...\n\n' >> "$tmp"
tail -c "$half" -- "$file" >> "$tmp"
chmod --reference="$file" -- "$tmp"
mv -f -- "$tmp" "$file"
