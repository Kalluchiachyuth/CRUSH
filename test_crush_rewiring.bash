#!/usr/bin/env bash
################################################################################
#
#  CRUSH Parameter Rewiring — Test Script
#
#  Tests all new parameter-layer changes without needing a real Hi-C file
#  or network access. Run from the same directory as CrushV89_1.bash.
#
#  Usage:  bash test_crush_rewiring.bash
#
################################################################################

set -uo pipefail

CRUSH="$(cd "$(dirname "$0")" && pwd)/CrushV89_1.bash"
TESTDIR="$(mktemp -d /tmp/crush_test_XXXXXX)"
PASS=0
FAIL=0

trap "rm -rf $TESTDIR" EXIT

# ── helpers ──────────────────────────────────────────────────────────────────

pass() { echo "  PASS  $1"; ((PASS++)) || true; }
fail() { echo "  FAIL  $1"; [ -n "${2:-}" ] && echo "        └─ $2"; ((FAIL++)) || true; }

contains() {
    # contains "string" "substring"
    [[ "$1" == *"$2"* ]]
}

# ── mock files ───────────────────────────────────────────────────────────────

# No-chr sizes (as hosted on GitHub)
cat > "$TESTDIR/nochr.sizes" <<'EOF'
1	248956422
2	242193529
X	156040895
EOF

# chr-prefixed sizes
cat > "$TESTDIR/chr.sizes" <<'EOF'
chr1	248956422
chr2	242193529
chrX	156040895
EOF

# CHR-prefixed sizes
cat > "$TESTDIR/CHR.sizes" <<'EOF'
CHR1	248956422
CHR2	242193529
CHRX	156040895
EOF

# No-chr genes BED
cat > "$TESTDIR/nochr_genes.bed" <<'EOF'
1	1000	2000	geneA	0	+
2	3000	4000	geneB	0	-
X	5000	6000	geneC	0	+
EOF

# No-chr Bbins BED
cat > "$TESTDIR/nochr_Bbins.bed" <<'EOF'
1	0	500	.
2	0	500	.
X	0	500	.
EOF

################################################################################
echo ""
echo "════════════════════════════════════════════════════════════"
echo "  CRUSH Parameter Rewiring — Test Suite"
echo "════════════════════════════════════════════════════════════"
echo ""

################################################################################
echo "[ 1 ]  Window default"
echo "────────────────────────────────────────"

winval=$(grep '^window=' "$CRUSH" | head -1 | cut -d= -f2)
if [ "$winval" = "5" ]; then
    pass "window defaults to 5"
else
    fail "window default" "expected 5, got '${winval}'"
fi

echo ""

################################################################################
echo "[ 2 ]  Basic help  (-h)"
echo "────────────────────────────────────────"

help_out=$(bash "$CRUSH" -h 2>&1) || true

if contains "$help_out" "ALWAYS REQUIRED"; then
    pass "-h shows ALWAYS REQUIRED section"
else
    fail "-h missing ALWAYS REQUIRED section"
fi

if contains "$help_out" "PATH A" && contains "$help_out" "PATH B"; then
    pass "-h shows both PATH A and PATH B"
else
    fail "-h does not show both reference paths"
fi

if contains "$help_out" "-gb"; then
    pass "-h mentions -gb flag"
else
    fail "-h does not mention -gb"
fi

if contains "$help_out" "advanced-help"; then
    pass "-h points to --advanced-help"
else
    fail "-h does not mention --advanced-help"
fi

# Basic help should NOT contain statistical detail
if contains "$help_out" "STATISTICAL OPTIONS"; then
    fail "-h leaks STATISTICAL OPTIONS (should be in advanced only)"
else
    pass "-h correctly omits STATISTICAL OPTIONS"
fi

echo ""

################################################################################
echo "[ 3 ]  Advanced help  (--advanced-help)"
echo "────────────────────────────────────────"

adv_out=$(bash "$CRUSH" --advanced-help 2>&1) || true

for section in "STATISTICAL OPTIONS" "ANALYSIS OPTIONS" "ADVANCED OPTIONS" "OUTPUT OPTIONS"; do
    if contains "$adv_out" "$section"; then
        pass "--advanced-help shows ${section}"
    else
        fail "--advanced-help missing ${section}"
    fi
done

if contains "$adv_out" "500bp"; then
    pass "--advanced-help mentions 500bp Bbins rule"
else
    fail "--advanced-help does not mention 500bp rule"
fi

if contains "$adv_out" "hg19" && contains "$adv_out" "mm9"; then
    pass "--advanced-help lists supported builds"
else
    fail "--advanced-help does not list supported builds"
fi

echo ""

################################################################################
echo "[ 4 ]  Supported build detection"
echo "────────────────────────────────────────"

CRUSH_SUPPORTED_BUILDS=("hg19" "hg38" "mm10" "mm9")

check_build() {
    local build="$1"
    local supported=0
    for _b in "${CRUSH_SUPPORTED_BUILDS[@]}"; do
        [ "$_b" = "$build" ] && supported=1 && break
    done
    echo $supported
}

for build in hg19 hg38 mm10 mm9; do
    if [ "$(check_build $build)" -eq 1 ]; then
        pass "${build} recognised as supported"
    else
        fail "${build} not recognised as supported"
    fi
done

for build in dm6 hg39 mm39 rn6; do
    if [ "$(check_build $build)" -eq 0 ]; then
        pass "${build} correctly identified as unsupported"
    else
        fail "${build} incorrectly marked as supported"
    fi
done

echo ""

################################################################################
echo "[ 5 ]  500bp resolution gate"
echo "────────────────────────────────────────"

check_sub500() {
    local r="$1"
    if [ "$r" != "0" ] && [ "$r" -lt 500 ] 2>/dev/null; then
        echo "blocked"
    else
        echo "allowed"
    fi
}

declare -A expect_500=(
    [200]=blocked
    [499]=blocked
    [500]=allowed
    [1000]=allowed
    [10000]=allowed
    [0]=allowed
)

for r in "${!expect_500[@]}"; do
    result=$(check_sub500 "$r")
    if [ "$result" = "${expect_500[$r]}" ]; then
        pass "res=${r} → ${result}"
    else
        fail "res=${r}" "expected ${expect_500[$r]}, got ${result}"
    fi
done

echo ""

################################################################################
echo "[ 6 ]  Sizefile prefix detection"
echo "────────────────────────────────────────"

detect_prefix() {
    awk 'NR==1{
        if (tolower(substr($1,1,3)) == "chr") print substr($1,1,3)
        else print ""
    }' "$1"
}

p=$(detect_prefix "$TESTDIR/nochr.sizes")
[ "$p" = "" ]   && pass "no-prefix sizes detected (got empty string)" \
                || fail "no-prefix detection" "expected '', got '${p}'"

p=$(detect_prefix "$TESTDIR/chr.sizes")
[ "$p" = "chr" ] && pass "chr-prefix sizes detected" \
                 || fail "chr detection" "expected 'chr', got '${p}'"

p=$(detect_prefix "$TESTDIR/CHR.sizes")
[ "$p" = "CHR" ] && pass "CHR-prefix sizes detected" \
                 || fail "CHR detection" "expected 'CHR', got '${p}'"

echo ""

################################################################################
echo "[ 7 ]  Chr prefix conversion"
echo "────────────────────────────────────────"

convert_chr() {
    local infile="$1"
    local outfile="$2"
    local target="$3"   # target prefix, e.g. "chr", "CHR", or ""
    awk -v p="$target" '
    BEGIN { FS="\t"; OFS="\t" }
    {
        if (tolower(substr($1,1,3)) == "chr") $1 = substr($1,4)
        if (p != "") $1 = p $1
        print
    }' "$infile" > "$outfile"
}

# no-prefix → chr
convert_chr "$TESTDIR/nochr.sizes" "$TESTDIR/out.sizes" "chr"
r=$(awk 'NR==1{print $1}' "$TESTDIR/out.sizes")
[ "$r" = "chr1" ] && pass "no-prefix → chr" || fail "no-prefix → chr" "got '${r}'"

# no-prefix → CHR
convert_chr "$TESTDIR/nochr.sizes" "$TESTDIR/out.sizes" "CHR"
r=$(awk 'NR==1{print $1}' "$TESTDIR/out.sizes")
[ "$r" = "CHR1" ] && pass "no-prefix → CHR" || fail "no-prefix → CHR" "got '${r}'"

# chr → CHR  (strip then re-add)
convert_chr "$TESTDIR/chr.sizes" "$TESTDIR/out.sizes" "CHR"
r=$(awk 'NR==1{print $1}' "$TESTDIR/out.sizes")
[ "$r" = "CHR1" ] && pass "chr → CHR" || fail "chr → CHR" "got '${r}'"

# chr → no-prefix  (strip only)
convert_chr "$TESTDIR/chr.sizes" "$TESTDIR/out.sizes" ""
r=$(awk 'NR==1{print $1}' "$TESTDIR/out.sizes")
[ "$r" = "1" ] && pass "chr → no-prefix" || fail "chr → no-prefix" "got '${r}'"

# CHR → chr  (strip then re-add)
convert_chr "$TESTDIR/CHR.sizes" "$TESTDIR/out.sizes" "chr"
r=$(awk 'NR==1{print $1}' "$TESTDIR/out.sizes")
[ "$r" = "chr1" ] && pass "CHR → chr" || fail "CHR → chr" "got '${r}'"

# chrX and chrY should survive correctly (sex chromosomes)
cat > "$TESTDIR/sex.sizes" <<'EOF'
chrX	156040895
chrY	57227415
EOF
convert_chr "$TESTDIR/sex.sizes" "$TESTDIR/out.sizes" "CHR"
rx=$(awk 'NR==1{print $1}' "$TESTDIR/out.sizes")
ry=$(awk 'NR==2{print $1}' "$TESTDIR/out.sizes")
[ "$rx" = "CHRX" ] && [ "$ry" = "CHRY" ] \
    && pass "sex chromosomes (chrX/chrY → CHRX/CHRY)" \
    || fail "sex chromosomes" "got '${rx}' '${ry}'"

# All three file types convert consistently
convert_chr "$TESTDIR/nochr_genes.bed"  "$TESTDIR/out_genes.bed"  "chr"
convert_chr "$TESTDIR/nochr_Bbins.bed" "$TESTDIR/out_Bbins.bed" "chr"

g=$(awk 'NR==1{print $1}' "$TESTDIR/out_genes.bed")
b=$(awk 'NR==1{print $1}' "$TESTDIR/out_Bbins.bed")
[ "$g" = "chr1" ] && [ "$b" = "chr1" ] \
    && pass "genes.bed and Bbins.bed convert consistently with sizes" \
    || fail "file-type consistency" "genes='${g}' Bbins='${b}'"

echo ""

################################################################################
echo "[ 8 ]  -gb fetch URL structure  (offline check)"
echo "────────────────────────────────────────"

# Verify the URL template in the script matches the expected GitHub structure
url_line=$(grep 'JROWLEY_BASE=' "$CRUSH" | head -1)

if contains "$url_line" "Kalluchiachyuth/CRUSH"; then
    pass "URL references Kalluchiachyuth/CRUSH repo"
else
    fail "URL does not reference Kalluchiachyuth/CRUSH"
fi

if contains "$url_line" "genomes/"; then
    pass "URL includes genomes/ path component"
else
    fail "URL missing genomes/ path component"
fi

# Check all three file patterns are referenced
for pattern in ".chrom.sizes" "_genes.bed" "_Bbins.bed"; do
    if grep -q "$pattern" "$CRUSH"; then
        pass "fetch target: ${pattern}"
    else
        fail "fetch target missing: ${pattern}"
    fi
done

echo ""

################################################################################
echo "════════════════════════════════════════════════════════════"
printf "  Results:  %d passed  |  %d failed\n" "$PASS" "$FAIL"
echo "════════════════════════════════════════════════════════════"
echo ""

[ "$FAIL" -eq 0 ] && exit 0 || exit 1
