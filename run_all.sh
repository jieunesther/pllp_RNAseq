#!/usr/bin/env bash
# ============================================================================
# run_all.sh  --  run the complete pllp RNA-seq pipeline from GEO inputs
# Usage:  bash run_all.sh            (all steps)
#         bash run_all.sh 05 06      (only the listed step numbers)
# Set RSCRIPT to the Rscript binary if it is not on PATH, e.g.
#   RSCRIPT=/nas/longleaf/rhel9/apps/r/4.4.0/bin/Rscript bash run_all.sh
# Each step writes results/logs/<script>.log; the run stops at the first error.
# Reproducibility assertions accumulate in results/reproducibility_check.csv.
# ============================================================================
set -euo pipefail
cd "$(dirname "$0")"
RSCRIPT="${RSCRIPT:-Rscript}"
mkdir -p results/logs

STEPS=(
  01_dataset1_DESeq2
  02_dataset1_enrichment
  03_dataset1_WGCNA
  04_dataset2_DESeq2
  05_figure2
  06_figure5AB
  07_heatmaps_curated_genes
  08_figureS2
  09_figureS3B
  10_figureS7
)

want=("$@")
for s in "${STEPS[@]}"; do
  if [ ${#want[@]} -gt 0 ] && ! printf '%s\n' "${want[@]}" | grep -qx "${s%%_*}"; then continue; fi
  echo "== $(date '+%H:%M:%S')  $s"
  "$RSCRIPT" "scripts/$s.R" > "results/logs/$s.log" 2>&1 || { echo "   FAILED - see results/logs/$s.log"; tail -20 "results/logs/$s.log"; exit 1; }
  grep -h "\[DIFF\]" "results/logs/$s.log" || true
done

echo
echo "== Reproducibility summary (results/reproducibility_check.csv)"
if [ -f results/reproducibility_check.csv ]; then
  awk -F, 'NR>1 {gsub(/"/,"",$5); n++; if ($5=="TRUE") ok++} END {printf "   %d of %d published values reproduced\n", ok, n}' results/reproducibility_check.csv
  grep -v ',TRUE,' results/reproducibility_check.csv | grep -v '^"script"' || true
fi
echo "== Done $(date '+%H:%M:%S')"
