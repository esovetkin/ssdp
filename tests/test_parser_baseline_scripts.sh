#!/usr/bin/env bash
set -euo pipefail

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
ssdp_bin="$(cd "$top_builddir/src" && pwd)/ssdp"
src_root="$(cd "$srcdir" && pwd)"

run_and_compare() {
	local case_dir="$1"
	local script_name="$2"
	shift 2
	local expected_outputs=("$@")
	local legacy_dir="${case_dir}.legacy"
	local bison_dir="${case_dir}.bison"
	local f

	cp -a "$case_dir" "$legacy_dir"
	cp -a "$case_dir" "$bison_dir"

	(
		cd "$legacy_dir"
		SSDP_PARSER=legacy "$ssdp_bin" -q -f "$script_name" >legacy.out 2>legacy.err
	)
	(
		cd "$bison_dir"
		SSDP_PARSER=bison "$ssdp_bin" -q -f "$script_name" >bison.out 2>bison.err
	)
	cmp "$legacy_dir/legacy.err" "$bison_dir/bison.err"

	for f in "${expected_outputs[@]}"; do
		if [ -e "$legacy_dir/$f" ] || [ -e "$bison_dir/$f" ]; then
			if [ ! -e "$legacy_dir/$f" ] || [ ! -e "$bison_dir/$f" ]; then
				echo "mismatch in output file presence: $f"
				exit 1
			fi
			cmp "$legacy_dir/$f" "$bison_dir/$f"
		fi
	done
}

run_sim0() {
	local src="$src_root/scripts"
	local dst="$tmp/sim0"
	mkdir -p "$dst/data"
	cp "$src/sim0.ssdp" "$dst/"
	cp "$src/data/sim0.txt" "$dst/data/"
	run_and_compare "$dst" "sim0.ssdp" sim0_static.txt sim0_static_integral.txt sim0_route.txt
}

run_sim1() {
	local src="$src_root/scripts"
	local dst="$tmp/sim1"
	mkdir -p "$dst/data/rasters"
	cp "$src/sim1.ssdp" "$dst/"
	cp "$src/data/input.h5" "$dst/data/"
	cp "$src/data/raster.list" "$dst/data/"

	while IFS= read -r rel; do
		[ -z "$rel" ] && continue
		cp "$src/$rel" "$dst/$rel"
	done < "$src/data/raster.list"

	run_and_compare "$dst" "sim1.ssdp" sim1.png
}

run_approx_test0() {
	local src="$src_root/scripts/approx"
	local dst="$tmp/approx"
	mkdir -p "$dst"
	cp "$src/test0.ssdp" "$dst/"
	cp "$src/bbox.txt" "$dst/"
	cp "$src/raster_max.list" "$dst/"

	while IFS= read -r rel; do
		[ -z "$rel" ] && continue
		mkdir -p "$dst/$(dirname "$rel")"
		cp "$src/$rel" "$dst/$rel"
	done < "$src/raster_max.list"

	run_and_compare "$dst" "test0.ssdp"
}

run_sim0
run_sim1
run_approx_test0

echo "test_parser_baseline_scripts.sh: PASS"
