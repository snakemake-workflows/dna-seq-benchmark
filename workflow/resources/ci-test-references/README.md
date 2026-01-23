# Reference files for CI-testing

## `new_tabbed_revel_grch38.1pct.tsv.gz`
The CI tests do not allow to use large files, thus the REVEL table was downsampled for CI testing. This was done by following the commands for `grch38` described in the [VEP Plugin docs](https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#revel):
```bash
curl https://zenodo.org/records/7072866/files/revel-v1.3_all_chromosomes.zip -o revel-v1.3_all_chromosomes.zip
unzip revel-v1.3_all_chromosomes.zip
cat revel_with_transcript_ids | tr "," "\t" > tabbed_revel.tsv
sed '1s/.*/#&/' tabbed_revel.tsv > new_tabbed_revel.tsv
bgzip new_tabbed_revel.tsv
zcat new_tabbed_revel.tsv.gz | head -n1 > h
zgrep -h -v ^#chr new_tabbed_revel.tsv.gz | awk '$3 != "." ' | sort -k1,1 -k3,3n - | cat h - | bgzip -c > new_tabbed_revel_grch38.tsv.gz
tabix -f -s 1 -b 3 -e 3 new_tabbed_revel_grch38.tsv.gz
```

To downsample this file we took every 100th row from original file, zipped and indexed it:
```bash
awk 'NR==1 || (NR-1)%100==0' new_tabbed_revel_grch38.tsv > new_tabbed_revel_grch38.1pct.tsv

bgzip new_tabbed_revel_grch38.1pct.tsv
tabix -f -s 1 -b 3 -e 3 new_tabbed_revel_grch38.1pct.tsv.gz
```

## `vep_cache_113_GRCh38_chr22.tar.gz`
```bash
curl -O https://ftp.ensembl.org/pub/release-113/variation/indexed_vep_cache/homo_sapiens_vep_113_GRCh38.tar.gz
tar xzf homo_sapiens_vep_113_GRCh38.tar.gz

mkdir vep_cache_113_GRCh38_chr22
mkdir vep_cache_113_GRCh38_chr22/homo_sapiens
mkdir vep_cache_113_GRCh38_chr22/homo_sapiens/113_GRCh38/

cp -r homo_sapiens/113_GRCh38/MT vep_cache_113_GRCh38_chr22/homo_sapiens/113_GRCh38/
cp -r homo_sapiens/113_GRCh38/22 vep_cache_113_GRCh38_chr22/homo_sapiens/113_GRCh38/
cp homo_sapiens/113_GRCh38/chr_synonyms.txt vep_cache_113_GRCh38_chr22/homo_sapiens/113_GRCh38/
cp homo_sapiens/113_GRCh38/info.txt vep_cache_113_GRCh38_chr22/homo_sapiens/113_GRCh38/

bash subsample_all_vars.sh
# mamba activate bcf (next script needs bgzip)
bash index_subsample.sh

rm -rf vep_cache_113_GRCh38_chr22/homo_sapiens/113_GRCh38/22/subsampled_vars
rm -rf vep_cache_113_GRCh38_chr22/homo_sapiens/113_GRCh38/22/all_vars.gz
rm -rf vep_cache_113_GRCh38_chr22/homo_sapiens/113_GRCh38/22/all_vars.gz.csi

mv vep_cache_113_GRCh38_chr22/homo_sapiens/113_GRCh38/22/subsampled_vars.sorted.bgz vep_cache_113_GRCh38_chr22/homo_sapiens/113_GRCh38/22/all_vars.gz
mv vep_cache_113_GRCh38_chr22/homo_sapiens/113_GRCh38/22/subsampled_vars.sorted.bgz.csi vep_cache_113_GRCh38_chr22/homo_sapiens/113_GRCh38/22/all_vars.gz.csi

tar cvf - vep_cache_113_GRCh38_chr22 | gzip -v > vep_cache_113_GRCh38_chr22.tar.gz
```
with `subsample_all_vars.sh`:
```bash
#!/usr/bin/env bash
set -euo pipefail

# Source file (change if needed)
SRC="vep_cache_113_GRCh38_chr22/homo_sapiens/113_GRCh38/22/all_vars.gz"

# Fraction to sample (1%)
FRACTION=0.001

# Output file in same directory
OUT_DIR=$(dirname "$SRC")
OUT="$OUT_DIR/subsampled_vars"

# Choose decompressor for gz or plain text
if [[ "$SRC" == *.gz ]]; then
  DECOMP="gunzip -c"
else
  DECOMP="cat"
fi

# Count lines
TOTAL_LINES=$($DECOMP "$SRC" | wc -l)
if [ "$TOTAL_LINES" -eq 0 ]; then
  echo "Input file is empty: $SRC"
  exit 1
fi

# Number of lines to sample (round to nearest integer, at least 1)
NUM_TO_SAMPLE=$(awk -v n="$TOTAL_LINES" -v f="$FRACTION" 'BEGIN{ k=int(n*f + 0.5); if(k<1) k=1; print k }')

# Temporary index file (line numbers) and cleanup
IDX_FILE=$(mktemp)
trap 'rm -f "$IDX_FILE"' EXIT

# Pick NUM_TO_SAMPLE distinct random line numbers from 1..TOTAL_LINES, then sort to preserve original order
shuf -i 1-"$TOTAL_LINES" -n "$NUM_TO_SAMPLE" | sort -n > "$IDX_FILE"

# Use awk to print only those line numbers from the original (keeps original ordering)
# Note: awk reads the index file first (NR==FNR block), then reads stdin via '-' (piped decompressed content)
$DECOMP "$SRC" | awk 'NR==FNR{keep[$1]=1; next} FNR in keep{print}' "$IDX_FILE" - > "$OUT"

echo "Wrote $NUM_TO_SAMPLE lines (approx. $(awk -v n="$TOTAL_LINES" -v k="$NUM_TO_SAMPLE" 'BEGIN{printf \"%.2f\", (k/n)*100}'))% of $TOTAL_LINES lines to: $OUT"
```
and `index_subsample.sh`:
```bash
#!/usr/bin/env bash
set -euo pipefail

SRC="${1:-vep_cache_113_GRCh38_chr22/homo_sapiens/113_GRCh38/22/subsampled_vars}"
OUT="${2:-${SRC%.gz}.sorted.bgz}"   # default: same name without .gz + .sorted.bgz

# temp dir for sort (change if you want)
TMPDIR="${TMPDIR:-/tmp}"

if [ ! -r "$SRC" ]; then
  echo "Input file not found or not readable: $SRC" >&2
  exit 2
fi

# choose decompressor (works for .gz or plain files)
if [[ "$SRC" == *.gz ]]; then
  DECOMP="gunzip -c"
else
  DECOMP="cat"
fi

echo "Reading: $SRC"
echo "Writing compressed sorted output: $OUT"
echo "Index will be created next to it: ${OUT}.csi"

# Pipeline:
# 1) Decompress (if gz)
# 2) Ensure column 6 has a numeric end (set to start col5 if '.' or empty)
# 3) Sort by col1 (chr) then numeric col5 (start)
# 4) bgzip write to $OUT
# 5) tabix index with chr=1 start=5 end=6
$DECOMP "$SRC" \
  | awk -F'\t' 'BEGIN{OFS=FS} { if($6=="" || $6==".") $6=$5; print }' \
  | LC_ALL=C sort -k1,1 -k5,5n -T "$TMPDIR" \
  | bgzip -c > "$OUT"

# create (or overwrite) CSI index; chr column 1, start 5, end 6
# Use --csi to request a CSI index (output will be ${OUT}.csi)
tabix -f --csi -s 1 -b 5 -e 6 "$OUT"

echo "Done: $OUT and ${OUT}.csi"
```