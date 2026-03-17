#!/usr/bin/env bash
set -euo pipefail

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
ssdp_bin="$(cd "$top_builddir/src" && pwd)/ssdp"

cat >"$tmp/parity.ssdp" <<'EOF'
# comment line should be ignored
make_scalar val=1 x='a'
make_scalar x=b val=2
array_eval c=c b=b op=+ a=a
write_array file=out.txt a0=c
EOF

(
	cd "$tmp"
	SSDP_PARSER=legacy "$ssdp_bin" -q -f parity.ssdp >legacy.out 2>legacy.err
	mv out.txt out_legacy.txt
)

(
	cd "$tmp"
	SSDP_PARSER=bison "$ssdp_bin" -q -f parity.ssdp >bison.out 2>bison.err
	mv out.txt out_bison.txt
)

cmp "$tmp/out_legacy.txt" "$tmp/out_bison.txt"
cmp "$tmp/legacy.err" "$tmp/bison.err"

cat >"$tmp/control.ssdp" <<'EOF'
unknown_command x=1
exit
make_scalar x=late val=10
write_array a0=late file=late.txt
EOF

(
	cd "$tmp"
	SSDP_PARSER=legacy "$ssdp_bin" -q -f control.ssdp >legacy2.out 2>legacy2.err
)

(
	cd "$tmp"
	SSDP_PARSER=bison "$ssdp_bin" -q -f control.ssdp >bison2.out 2>bison2.err
)

cmp "$tmp/legacy2.err" "$tmp/bison2.err"
if [ -e "$tmp/late.txt" ]; then
	echo "late.txt should not exist after exit"
	exit 1
fi

echo "test_parser_modes.sh: PASS"
