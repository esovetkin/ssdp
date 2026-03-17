#!/usr/bin/env bash
set -euo pipefail

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
ssdp_bin="$(cd "$top_builddir/src" && pwd)/ssdp"

cat >"$tmp/bad.ssdp" <<'EOF'
make_scalar x=a val=1
make_scalar x=b val=
EOF

(
	cd "$tmp"
	SSDP_PARSER=bison "$ssdp_bin" -q -f bad.ssdp >/dev/null 2>diag.err || true
)

if ! rg -q "Parse error at .*bad\\.ssdp:2:" "$tmp/diag.err"; then
	echo "missing line/column parse diagnostic in stderr"
	exit 1
fi

echo "test_parser_diagnostics.sh: PASS"
