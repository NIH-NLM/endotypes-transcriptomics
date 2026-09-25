# endotypes-transcriptomics

**Patient endotypes from blood gene expression in pediatric systemic lupus erythematosus (SLE):
federated from the moment data sit at separate sites, one step per notebook, on a microarray study
and two RNA sequencing studies.**

An **endotype** is a group of patients who share a pattern of gene expression. We look for such groups
without telling the method anything about disease activity, nephritis or interferon, and then ask
what the groups mean.

This repository follows the layout of
[endotypes-proteomics](https://github.com/NIH-NLM/endotypes-proteomics).

## Three studies, two platforms

| study | data type | tissue | patients | reference |
|---|---|---|---|---|
| **GSE65391** | `RNA_array`: Illumina HumanHT-12 microarray | whole blood | 158 children with SLE, 924 visits; 46 healthy children | Banchereau R et al. *Cell* 2016;165:551–565 |
| **GSE232381** | `bulk_RNA_seq`: NovaSeq | peripheral blood cells | 16 with lupus nephritis (10 active, 6 inactive); no age in GEO | Chen YC et al. *Heliyon* 2024;10:e32303 |
| **GSE135779** | `scRNA_seq`, used as pseudobulk | PBMC | 33 children with SLE, 11 healthy children | Nehar-Belaid D et al. *Nat Immunol* 2020;21:1094–1106 |

Every notebook's name says which data type it works on.

## The two rules

**1. From step 04 on, no data leave a site.** GSE65391 is split into three simulated sites, A, B and C,
standing in for three institutions. After that, a site sends only summaries:
- per-gene counts, sums and sums of squares;
- model parameters;
- per-cluster counts and sums;
- histograms and contingency tables.

Every message goes through `send()` (`src/federation.R`), and step 17 audits the whole record. Because
the sites are simulated, some notebooks also compute the same quantity on all samples in one place,
in cells headed **"Oracle — possible only because this is a simulation"**. Those cells only prove
exactness or measure accuracy. Nothing an oracle computes feeds a later step.

**2. Platforms are never merged at the value level.** Array intensities and sequencing counts differ in
scale, noise and dynamic range. Each study is normalised by its own platform's method:
- array log2 intensities go to limma;
- counts go through TMM and voom to log2 counts per million.

The studies meet in three ways only:
- **shared gene symbols**, keyed through NCBI Entrez identifiers;
- **shared units**, with every gene standardised within its own study;
- **scale-free statistics**: the median-z interferon score, and correlation to centroids.

## Running on ADAPTS (Lifebit)

1. Log in to lifebit.ai and authenticate with your PIV (Personal Identity Verification) card.
2. Start a JupyterLab notebook (8 vCPUs, 16 GB RAM).
3. In a terminal:

```bash
git clone https://github.com/NIH-NLM/endotypes-transcriptomics.git
cd endotypes-transcriptomics
mamba env create -f endotypes-transcriptomics.yml
mamba activate endotypes-transcriptomics
Rscript -e 'IRkernel::installspec()'
python -m bash_kernel.install
```

4. Open `ipynb/00_RNA_array_download_GSE65391.ipynb` under the **R** kernel and run the notebooks in
   order, or run them all:

```bash
./run_all.sh
```

That runs steps 00–21 (GSE65391). To add the RNA sequencing studies (steps 22–31, including a 1.2 GB
download):

```bash
./run_all.sh --cross-platform
```

## Running on a laptop

The same environment file builds on macOS (Apple Silicon and Intel) and Linux. Every package has a
conda build, so nothing is installed from inside a notebook. If you already have an `ir` kernel from
another project, register this one under its own name:

```bash
Rscript -e 'IRkernel::installspec(name = "ir-endotypes-transcriptomics", displayname = "R (endotypes-transcriptomics)")'
```

```bash
KERNEL=ir-endotypes-transcriptomics ./run_all.sh --cross-platform
```

Every download happens inside the notebooks and only once (about 170 MB for GSE65391, 1.3 GB for
GSE135779).

## Repository layout

```
cohorts/     COMMITTED  frozen site assignment and discovery/validation split
genes/       COMMITTED  interferon gene sets (six-gene score; 28-gene panel to fill in)
figures/     COMMITTED  every figure, named figNN_<data type>_<what>.png
ipynb/       COMMITTED  the notebooks, and nothing else
src/         COMMITTED  paths.R, endotypes.R, federation.R, platforms.R, figures.R
reference/   COMMITTED  the earlier Nextflow scaffold, the source of the analysis logic
data/        ignored    one folder per GEO series, data/ncbi/ (gene_info), data/run_artifacts/
```

`src/paths.R` is the only place locations are written down. `data/run_artifacts/` can be deleted at
any time; `./run_all.sh` rebuilds it, and the frozen `cohorts/` files keep every draw the same.

## Notebooks

### GSE65391 (`RNA_array`): preparation, before sites exist
| notebook | does |
|---|---|
| `00_RNA_array_download_GSE65391` | series matrix, platform annotation, NCBI gene_info |
| `01_RNA_array_metadata_and_stage` | tidy metadata; stage = SLEDAI category per visit |
| `02_RNA_array_probes_to_genes` | brightest probe per Entrez gene, current NCBI symbols; expressed genes |
| `03_RNA_array_study_batch_check` | technical replicates show the study batch is already removed |

### GSE65391: sites and federated batch correction
| notebook | does |
|---|---|
| `04_RNA_array_assign_sites` | three pseudo-sites; each site's data in its own file from here on |
| `05_RNA_array_plant_site_effects` | each site distorts its own data with a known shift and scale |
| `06_RNA_array_federated_combat` | ComBat, federated: one round of per-gene sums, then local. **Equals `sva::ComBat`** |
| `07_RNA_array_federated_location_scale` | the simpler correction without empirical Bayes |
| `08_RNA_array_compare_corrections` | accuracy against the truth, from per-site sums; PCA drawn from site density grids |

### GSE65391: federated endotype discovery
| notebook | does |
|---|---|
| `09_RNA_array_split_patients` | each site splits its own patients 70/30, by patient |
| `10_RNA_array_interferon_score` | healthy reference pooled from site sums; each site scores locally |
| `11_RNA_array_federated_genes_and_pcs` | top 1,000 genes; PCs by federated subspace iteration, no Gram matrix |
| `12_RNA_array_federated_consensus` | federated consensus k-means, k by pooled PAC |
| `13_RNA_array_assign_endotypes` | federated centroids; every sample labelled at its own site |
| `14_RNA_array_heatmap_consensus` | consensus matrices, one row per site |
| `15_RNA_array_heatmap_genes_by_samples` | genes × samples per site; endotype against stage, nephritis and published groups |
| `16_RNA_array_heatmap_visit_stability` | endotype across visits per child, per site |

### GSE65391: does federation improve the endotypes?
| notebook | does |
|---|---|
| `17_RNA_array_federation_audit` | every message sent; none has one value per patient |
| `18_RNA_array_clustering_alone` | level 1: each site alone |
| `19_RNA_array_shared_gene_list` | level 2: sites share only the gene list |
| `20_RNA_array_federation_benefit` | alone → shared genes → federated model, measured four ways |
| `21_RNA_array_heatmap_alone_vs_federated` | the same patients grouped both ways |

### GSE232381 (`bulk_RNA_seq`), GSE135779 (`scRNA_seq`), and cross-platform
| notebook | does |
|---|---|
| `22_bulk_RNA_seq_download_GSE232381` | NCBI-generated raw counts (16 of 25 samples) |
| `23_bulk_RNA_seq_prepare_GSE232381` | Entrez → symbol; TMM + voom log2 counts per million |
| `24_scRNA_seq_download_GSE135779` | per-donor count matrices (1.2 GB) |
| `25_scRNA_seq_pseudobulk_GSE135779` | counts summed per child; TMM + voom |
| `26_cross_platform_shared_genes` | genes all three measure; coverage of the model |
| `27_cross_platform_signature_preservation` | does the E1–E2 axis exist on each platform? |
| `28_cross_platform_interferon_score` | each study against its own reference |
| `29_cross_platform_transfer_endotypes` | use A: send the array model to the sequencing studies |
| `30_cross_platform_federated_discovery` | use B: five sites, two platforms, one federated clustering |
| `31_cross_platform_heatmaps` | the same genes in every study, each in its own units |

## Results

**The study's own batch effect is already removed** (step 03). 23 healthy samples run in both batches
show no consistent shift.

**Federated batch correction is exact** (steps 06–08). Federated ComBat equals `sva::ComBat` on
pooled data to 1.7e-10. A planted site effect that explained 48% of every gene's variance is removed
completely, and corrected values correlate with the truth at r = 0.986.

**Two endotypes, found federated** (steps 11–16). The federated genes, components and interferon score
equal their central counterparts exactly. k = 2 has a pooled PAC of 0.155, against 0.37 or more for
larger k.
- **E1** is higher in T and B lymphocyte genes (LEF1, BCL11B, CCR7, CD79B, CD27, TCL1A). Healthy
  children look like E1.
- **E2** is higher in neutrophil granule genes (LTF, ELANE, MPO, AZU1, DEFA3/4, CAMP, CEACAM8) and red
  cell genes.
- Interferon is raised in both, so it does not separate them.
- Over all visits, E2 rises from 42% of visits with no disease activity to 62% at high activity. At
  the first visit, one sample per child, there is no association with stage (p = 0.84).
- 32% of children with two or more visits never change endotype. The endotypes look like states of
  blood composition that children move between.

**Federation is what makes the endotypes reproducible** (steps 18–21).

| | alone | + shared gene list | + federated model |
|---|---|---|---|
| k chosen at sites A, B, C | 6, 5, 2 | 6, 5, 2 | 2, 2, 2 |
| agreement between site models (ARI) | 0.27 | 0.30 | one model |
| agreement with a central analysis (ARI), A / B / C | 0.28 / 0.30 / 0.83 | 0.31 / 0.30 / 0.83 | 1.00 / 1.00 / 1.00 |

With about 37 patients each, two of three sites cannot find the structure alone, and a shared gene list
does not help. The limit is patients, not genes, which is unlike the proteomics study. Only a shared
definition lets the sites pool their validation patients. Pooled, E2 has higher neutrophil counts
(6.4 against 2.9), white cell counts (8.1 against 5.0) and C4 (q ≤ 0.02). That matches its gene
signature, and it is shown on patients who played no part in finding it.

**Across platforms** (steps 26–31):
- **The E1–E2 axis exists in all three studies.** E2 markers rise with it (median r +0.51 to +0.77)
  and E1 markers fall (−0.59 to −0.79).
- **PBMC barely contain the granulocyte genes** (5th–22nd percentile of expression), so the transfer
  gate **refuses** to apply the array model to GSE135779 (68% gene coverage, below 70%).
- **Clustering all five sites together**, on the 677 genes every study measures, keeps k = 2 and
  leaves the array endotypes nearly unchanged (ARI 0.79–0.89).
- **Both endotypes occur on every platform.** Endotype is unrelated to lupus nephritis activity in
  GSE232381 (p = 1.0).
- **The interferon score** separates SLE from healthy children in GSE135779 (p = 2e-5) but does not
  track nephritis activity in GSE232381.

## Next

1. Fill `genes/ifn-interferonopathy-28.txt` from the primary source; steps 10 and 28 then score it.
2. Differential expression between endotypes on validation patients: limma with
   `duplicateCorrelation` for the array; DESeq2 and limma-voom for counts
   (`reference/bin/run_limma.R`, `run_deseq2.R`). Federated as per-site model coefficients.
3. Prediction of disease stage from the endotype at a visit and the change since the previous visit.
