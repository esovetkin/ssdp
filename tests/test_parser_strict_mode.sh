#!/usr/bin/env bash
set -euo pipefail

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
ssdp_bin="$(cd "$top_builddir/src" && pwd)/ssdp"

cat >"$tmp/strict.ssdp" <<'EOF'
   # leading-space comment must be ignored in strict mode
make_scalar x=a val=1 val=2 extra=3
write_array a0=a file=out.txt
EOF

(
	cd "$tmp"
	SSDP_PARSER=bison SSDP_PARSER_STRICT=1 "$ssdp_bin" -q -f strict.ssdp >/dev/null 2>strict.err
)

if ! rg -q "duplicate argument val" "$tmp/strict.err"; then
	echo "expected duplicate-argument warning"
	exit 1
fi

if ! rg -q "unknown argument extra" "$tmp/strict.err"; then
	echo "expected unknown-argument warning"
	exit 1
fi

if rg -q "Command # not defined" "$tmp/strict.err"; then
	echo "leading-space comment was not ignored in strict mode"
	exit 1
fi

if [ ! -s "$tmp/out.txt" ]; then
	echo "strict mode script did not produce output"
	exit 1
fi

echo "test_parser_strict_mode.sh: PASS"
