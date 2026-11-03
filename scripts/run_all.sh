#!/usr/bin/env bash
# Run the whole primer pipeline in order. Safe to re-run: every step overwrites its own output.
#
#   ./scripts/run_all.sh            everything
#   ./scripts/run_all.sh --no-tx    curated-transcript work only (stages 1-4 + Check 1/4)
#
# --no-tx skips the three transcriptome-wide checks (Checks 2, 3, 3b) and stage 5, because
# verify_sdh3_backup.py rescans the assembly. Stages 1-4 run entirely on committed data.
#
# Requires: python3 with numpy. Full run also needs blastn/makeblastdb on PATH and data/external/
# populated (see data/external/README.md).
set -euo pipefail
cd "$(dirname "$0")/.."
export PYTHONDONTWRITEBYTECODE=1

run() { echo; echo "### $1"; python3 "$1"; }

run scripts/step1_inputs/parse_genbank.py
run scripts/step1_inputs/compare_genbank_coords.py || echo "    skipped: needs all 6 GenBank records (see reports/01)"
run scripts/step1_inputs/extract_cornus_utrs.py
# run scripts/step2_survey/analyze_variable_positions.py  # omitted: alignment CSV moved to archive
run scripts/step2_survey/unique_20mer_spans.py
run scripts/step3_design/design_final_pairs.py
run scripts/step3_design/validate_final_primers.py
run scripts/step4_specificity/run_primer_blast_check.py
run scripts/step4_specificity/curated_exhaustive_scan.py

if [ "${1:-}" != "--no-tx" ]; then
  run scripts/step5_sdh3_backup/design_sdh3_backup.py
  run scripts/step5_sdh3_backup/verify_sdh3_backup.py
  run scripts/step4_specificity/run_transcriptome_check.py
  run scripts/step4_specificity/tx_exhaustive_scan.py
  run scripts/step4_specificity/tx_sensitivity_check.py
fi

echo
echo "Done. Deliverables in data/results/, evidence in data/results/specificity/, narrative in reports/."
