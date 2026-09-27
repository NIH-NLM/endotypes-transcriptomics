#!/usr/bin/env bash
# Run the notebooks in order. Every analysis step is in a notebook; this script only executes them.
#
#   ./run_all.sh                      all three studies, then federation
#   ./run_all.sh GSE65391             one study (GSE65391, GSE232381 or GSE135779)
#   ./run_all.sh federation           the federation notebooks only
#   KERNEL=ir-endotypes-transcriptomics ./run_all.sh
#                                     use a named R kernel instead of "ir"
#
# R notebooks run on the R kernel ($KERNEL); the pvclust-py notebooks run on the environment's
# python3 kernel. Each kernel must belong to the active conda environment, or the script stops.
# Notebooks execute in place, so each saved notebook carries its own output; every figure is also
# written as a PNG to figures/<GSE>/.
# GSE232381 and GSE135779 read the GSE65391 module files, so run GSE65391 first.
# GSE135779 needs data/GSE135779/41590_2020_743_MOESM3_ESM.xlsx, downloaded by hand (see README);
# notebook 20 stops if it is missing.
set -euo pipefail
cd "$(dirname "$0")/ipynb"

R_KERNEL="${KERNEL:-ir}"
PY_KERNEL="python3"
case "${1:-all}" in
  GSE65391)   PATTERNS=("0[0-9]*_GSE65391_*") ;;
  GSE232381)  PATTERNS=("1[0-9]*_GSE232381_*") ;;
  GSE135779)  PATTERNS=("2[0-9]*_GSE135779_*") ;;
  federation) PATTERNS=("3[0-9]_federation_*") ;;
  all)        PATTERNS=("0[0-9]*_GSE65391_*" "1[0-9]*_GSE232381_*" "2[0-9]*_GSE135779_*" "3[0-9]_federation_*") ;;
  *) echo "unknown target: $1 (use GSE65391, GSE232381, GSE135779, federation or all)"; exit 1 ;;
esac

# Each kernel must run the interpreter of the active conda environment. A kernel called "ir"
# registered by another project runs that project's R and package versions, which change results.
check_kernel() {
  local k="$1" exe
  exe=$(jupyter kernelspec list --json | python -c "import json,sys; print(json.load(sys.stdin)['kernelspecs']['$k']['spec']['argv'][0])" 2>/dev/null || true)
  if [[ -z "$exe" ]]; then echo "kernel '$k' not found; see: jupyter kernelspec list"; exit 1; fi
  if [[ -z "${CONDA_PREFIX:-}" || "$exe" != "$CONDA_PREFIX"/* ]]; then
    echo "kernel '$k' runs $exe"
    echo "but the active environment is ${CONDA_PREFIX:-<none>}"
    echo "activate endotypes-transcriptomics and use its kernels, for example:"
    echo "  KERNEL=ir-endotypes-transcriptomics ./run_all.sh"
    exit 1
  fi
}
check_kernel "$R_KERNEL"
check_kernel "$PY_KERNEL"

for pattern in "${PATTERNS[@]}"; do
  for nb in $(ls $pattern.ipynb | sort); do
    lang=$(python -c "import json,sys; print(json.load(open(sys.argv[1]))['metadata']['kernelspec']['language'])" "$nb")
    if [[ "$lang" == "python" ]]; then k="$PY_KERNEL"; else k="$R_KERNEL"; fi
    printf '%-60s ' "$nb"
    start=$SECONDS
    if jupyter nbconvert --to notebook --execute --inplace \
         --ExecutePreprocessor.kernel_name="$k" \
         --ExecutePreprocessor.timeout=7200 "$nb" >/dev/null 2>"/tmp/${nb%.ipynb}.err"; then
      # keep the committed kernel names standard ("ir" or "python3"), whatever ran it locally
      python -c "
import json, sys
f, lang = sys.argv[1], sys.argv[2]
d = json.load(open(f))
d['metadata']['kernelspec'] = ({'display_name': 'Python 3 (ipykernel)', 'language': 'python', 'name': 'python3'}
                               if lang == 'python' else {'display_name': 'R', 'language': 'R', 'name': 'ir'})
json.dump(d, open(f, 'w'), indent=1)" "$nb" "$lang"
      echo "OK  ($((SECONDS - start)) s)"
    else
      echo "FAILED"; tail -20 "/tmp/${nb%.ipynb}.err"; exit 1
    fi
  done
done
echo "done"
