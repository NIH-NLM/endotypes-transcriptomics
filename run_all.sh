#!/usr/bin/env bash
# The machine path. The same notebooks are the human path in JupyterLab.
#
#   ./run_all.sh                    steps 00-21: GSE65391 (RNA_array), federated
#   ./run_all.sh --cross-platform   also 22-31: GSE232381 (bulk_RNA_seq) and
#                                   GSE135779 (scRNA_seq), a 1.2 GB download
#   KERNEL=ir-endotypes-transcriptomics ./run_all.sh
#                                   use a named R kernel instead of "ir"
#
# Notebooks execute in place, so each committed notebook carries its own output.
# Every step reads cohorts/ (committed) and data/, and writes data/run_artifacts/
# (gitignored, regenerable) and figures/ (committed).
set -euo pipefail
cd "$(dirname "$0")/ipynb"

KERNEL="${KERNEL:-ir}"
STEPS=(00_RNA_array_download_GSE65391 01_RNA_array_metadata_and_stage
       02_RNA_array_probes_to_genes 03_RNA_array_study_batch_check
       04_RNA_array_assign_sites 05_RNA_array_plant_site_effects
       06_RNA_array_federated_combat 07_RNA_array_federated_location_scale
       08_RNA_array_compare_corrections 09_RNA_array_split_patients
       10_RNA_array_interferon_score 11_RNA_array_federated_genes_and_pcs
       12_RNA_array_federated_consensus 13_RNA_array_assign_endotypes
       14_RNA_array_heatmap_consensus 15_RNA_array_heatmap_genes_by_samples
       16_RNA_array_heatmap_visit_stability 17_RNA_array_federation_audit
       18_RNA_array_clustering_alone 19_RNA_array_shared_gene_list
       20_RNA_array_federation_benefit 21_RNA_array_heatmap_alone_vs_federated)
if [[ "${1:-}" == "--cross-platform" ]]; then
  STEPS+=(22_bulk_RNA_seq_download_GSE232381 23_bulk_RNA_seq_prepare_GSE232381
          24_scRNA_seq_download_GSE135779 25_scRNA_seq_pseudobulk_GSE135779
          26_cross_platform_shared_genes 27_cross_platform_signature_preservation
          28_cross_platform_interferon_score 29_cross_platform_transfer_endotypes
          30_cross_platform_federated_discovery 31_cross_platform_heatmaps)
fi

for nb in "${STEPS[@]}"; do
  printf '%-44s ' "$nb"
  start=$SECONDS
  if jupyter nbconvert --to notebook --execute --inplace \
       --ExecutePreprocessor.kernel_name="$KERNEL" \
       --ExecutePreprocessor.timeout=7200 "$nb.ipynb" >/dev/null 2>"/tmp/$nb.err"; then
    echo "OK  ($((SECONDS - start)) s)"
  else
    echo "FAILED"; tail -20 "/tmp/$nb.err"; exit 1
  fi
done
echo "all steps completed"
