# endotypes-transcriptomics

**Patient endotype discovery from blood gene expression in pediatric systemic lupus erythematosus (SLE)**:
demonstrating federation from separate sites. We have a microarray study, one bulk RNA sequencing
study, and one single-cell RNA sequencing study analysed as pseudobulk.

An **endotype** is a group of patients who share a pattern of gene expression.

## Three studies, three different types of datasets

Each study is one site.

| study | data type | tissue | patients | reference |
|---|---|---|---|---|
| **GSE65391** | `RNA_array`: Illumina HumanHT-12 microarray | whole blood | 158 children with SLE, 924 visits; 46 healthy children | Banchereau R et al. *Cell* 2016;165:551–565 |
| **GSE232381** | `bulk_RNA_seq`: NovaSeq | peripheral blood cells | 16 with lupus nephritis (10 active, 6 inactive); no age in GEO | Chen YC et al. *Heliyon* 2024;10:e32303 |
| **GSE135779** | `scRNA_seq`, used as pseudobulk | PBMC | 33 children with SLE, 11 healthy children | Nehar-Belaid D et al. *Nat Immunol* 2020;21:1094–1106 |

Clinical data for GSE135779 come from the paper's Supplementary Table 1b; GEO holds only age and group.

## Analysis

These are three different platforms. The RNA_array is the largest; it is an open question how best to
leverage these in a federated manner. All data in this analysis are public.

Each study is prepared on its own:
- array: deposited log2 intensities; brightest probe per gene;
- RNA sequencing counts (bulk, and pseudobulk summed per child): genes with fewer than 10 counts in
  more than 90% of samples removed (WGCNA FAQ), then PFlog1pPF normalisation (Booeshaghi AS et al.,
  bioRxiv 2022, doi:10.1101/2022.05.06.490859).

Genes are named with current NCBI symbols, reached through Entrez identifiers, in every study.

Per study:
1. **WGCNA** gene modules (Langfelder and Horvath 2008), signed network.
2. **Module–trait associations** with that study's clinical variables.
3. **Heatmaps** of module hub genes × patients, clustered on both axes with `hclust` and with k-means.
4. For the two RNA sequencing studies, **preservation** of the GSE65391 modules (Langfelder et al. 2011),
   using the GSE65391 module gene lists only.

**Federation.** Each study writes its end products to `data/run_artifacts/<GSE>/`: module gene lists,
module–trait tables (r, n, Fisher z and its standard error) and a clinical dictionary. The federation
notebooks read only these files, never expression values or patient records.

## Results

All numbers are from the rendered notebooks; the notebook that produces each is named.

| | GSE65391 · RNA_array | GSE232381 · bulk_RNA_seq | GSE135779 · scRNA_seq (pseudobulk) |
|---|---|---|---|
| samples used for modules | 157 first visits | 16 | 33 children with SLE |
| WGCNA modules (05, 14, 25) | 16 | 34 (exploratory, n = 16) | 29 |
| interferon module | pink | paleturquoise | royalblue |
| 28 interferon response genes in that module (05c, 14b, 25b) | 26 of 28 | 14 of 28; all 28 kME ≥ 0.66 | 27 of 28 |
| 11 NF-κB control genes in that module | 0 of 11 | 0 of 11 | 1 of 11 (GZMB) |
| gap statistic, patients (06a, 16a, 27a) | k = 1 | k = 1 | k = 1 |
| gap statistic, hub genes | k = 15 | k = 2 (no clear k) | k = 10 |
| pvclust, interferon hub genes (08c, 18c, 29c) | one cluster, AU 1.000, BP 1.000 | with tan, AU 0.998, BP 0.102 | with midnightblue and greenyellow, AU 0.968, BP 0.033 |
| pvclust, patients with AU ≥ 0.95 (08d, 18d, 29d) | 14 clusters of 2–5 | 14 and 2 samples | 11 children, and five pairs |
| Dynamic Tree Cut, genes against modules (08b, 18b, 29b) | ARI 1.000 | ARI 0.776 | ARI 0.484 |

- **No patient groups.** The gap statistic gives k = 1 in all three studies. pvclust supports only
  small clusters of patients, on the axis where AU is anti-conservative. Dynamic Tree Cut always returns
  groups, so its patient groups are not evidence of types.
- **The interferon genes stay together** in every study and every method. The 28-gene interferon
  response gene score (de Jesus AA et al. *J Clin Invest* 2020;130:1669–1682, Supplemental Figure 3B)
  falls in one module in the array and in the scRNA_seq pseudobulk; the 11 NF-κB-only genes, the
  control, do not.
- **Projection.** The array modules scored in the array's own samples reproduce the WGCNA eigengenes
  (r = 1, largest difference 3 × 10⁻¹⁵; 05b). Projected into the other studies, they agree with a local
  refit for 12 of 16 modules (bulk_RNA_seq, r ≥ 0.976; 15b) and 13 of 16 (scRNA_seq, r ≥ 0.945; 26b).
- **Comparison with the original papers.**
  - Banchereau et al. 2016: B3 reproduced (plasma-cell module has the strongest within-child association
    with SLEDAI, t = 8.7; 06b); B5 reproduced (red-cell module higher, NK module lower in SLE; 06c);
    B2 consistent (interferon R² = 0.20 with SLEDAI; 06); B1 approximate (65.6% against 84.8%; 06c);
    B4 same direction, not significant; B6 not testable with one visit per child.
  - Chen et al. 2024: C2 partly reproduced (23 of 27 genes in the paper's direction, Spearman ρ = 0.56
    to 0.61; the five largest effects not reproduced; none significant; 16b); C3 direction only; C4 not
    the same measure (16).
  - Nehar-Belaid et al. 2020: N1 consistent (interferon module higher in SLE, q = 0.001); N2 partly;
    N3 consistent in direction (27).

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
```

4. Download the GSE135779 supplement by hand (the publisher's site does not allow scripted downloads):
   open https://www.nature.com/articles/s41590-020-0743-0, download Supplementary Tables 1–4
   (`41590_2020_743_MOESM3_ESM.xlsx`) and save it in `data/GSE135779/`.
5. Open the notebooks in `ipynb/` under the **R** kernel and run them in order, or run them all from
   the terminal:

```bash
./run_all.sh
```

   `run_all.sh` only executes the notebooks, in order, and stops at the first failure. One study at a
   time: `./run_all.sh GSE65391` (then `GSE232381`, `GSE135779`, `federation`); GSE65391 must run first.
   Every other download happens inside the notebooks, once (about 170 MB for GSE65391, 1.3 GB for
   GSE135779).

## Running on a laptop

The same environment file builds on macOS (Apple Silicon and Intel) and Linux. WGCNA has no conda
build for Apple Silicon; the WGCNA notebooks install it from CRAN in their first cell if it is missing.
If you already have an `ir` kernel from another project, it runs that project's R and package
versions, which change the results. `run_all.sh` stops if the kernel's R is not the active
environment's. Register this one under its own name:

```bash
Rscript -e 'IRkernel::installspec(name = "ir-endotypes-transcriptomics", displayname = "R (endotypes-transcriptomics)")'
```

```bash
KERNEL=ir-endotypes-transcriptomics ./run_all.sh
```

## Notebooks

One notebook per step per study. Every test states what is compared, on which samples, with which
method, and what counts as a finding. Each study ends with its comparison against the published
paper's claims. Ward's method (`ward.D2`, Minkowski p = 2) is always drawn as an uncut tree; k appears
only for k-means, taken from the gap statistic. Notebooks marked (Python) run pvclust-py.

### GSE65391 · RNA_array
| notebook | does |
|---|---|
| `00_…_download` | series matrix, platform annotation, NCBI gene_info |
| `01_…_metadata` | clinical fields; SLEDAI stage per visit; first visit per child |
| `02_…_probes_to_genes` | probes to genes; expressed genes |
| `03_…_batch_check` | technical replicates across the two array batches |
| `04_…_soft_threshold` | WGCNA power, first visit per child |
| `05_…_modules` | WGCNA modules and hub genes |
| `05b_…_projection` | every sample scored on the first-visit modules; round-trip gate |
| `05c_…_interferon_genes` | the 28 interferon response genes and 11 NF-κB control genes against the modules |
| `06_…_module_traits` | module scores vs clinical traits; claims B2, B4, B6 |
| `06a_…_gap_statistic` | k-means k for patients and genes |
| `06b_…_modules_across_visits` | mixed model over all visits; claim B3 |
| `06c_…_sle_vs_healthy` | module scores, SLE vs healthy; claims B1, B5 |
| `07_…_heatmap_hclust` | Ward trees on patients and genes, first visits |
| `07b_…_heatmap_all_visits_ordered` | every visit, children ordered by mean SLEDAI and by Ward tree |
| `07c_…_heatmap_all_visits_clustered` | every visit, Ward tree |
| `08_…_heatmap_kmeans` | k-means with gap-statistic k |
| `08b_…_dynamic_tree_cut` | Dynamic Tree Cut on the Ward trees |
| `08c_…_pvclust_genes` (Python) | AU p-values and BP, genes |
| `08d_…_pvclust_patients` (Python) | AU p-values and BP, patients; two-way figure |

### GSE232381 · bulk_RNA_seq
| notebook | does |
|---|---|
| `10_…_download` | NCBI-generated raw counts (16 of 25 samples), series matrix |
| `11_…_metadata` | active / inactive lupus nephritis; claims C0, C1 |
| `12_…_normalise` | filter, PFlog1pPF |
| `13_…_soft_threshold` | WGCNA power (exploratory, n = 16) |
| `14_…_modules` | WGCNA modules (exploratory) |
| `14b_…_interferon_genes` | the 28 interferon response genes and 11 NF-κB control genes against the modules |
| `15_…_array_module_preservation` | are the GSE65391 modules present here? |
| `15b_…_array_module_projection` | samples scored on the GSE65391 modules |
| `16_…_module_traits` | module scores vs nephritis activity; claims C3, C4 |
| `16a_…_gap_statistic` | k-means k for samples and genes |
| `16b_…_published_genes` | the paper's named genes, DESeq2 and limma-voom; claim C2 |
| `17_…_heatmap_hclust` | Ward trees on samples and genes |
| `18_…_heatmap_kmeans` | k-means with gap-statistic k |
| `18b_…_dynamic_tree_cut` | Dynamic Tree Cut on the Ward trees |
| `18c_…_pvclust_genes` (Python) | AU p-values and BP, genes |
| `18d_…_pvclust_samples` (Python) | AU p-values and BP, samples; two-way figure |

### GSE135779 · scRNA_seq as pseudobulk
| notebook | does |
|---|---|
| `20_…_download` | per-donor count matrices; checks the supplement |
| `21_…_metadata` | GEO fields joined to Supplementary Table 1b |
| `22_…_pseudobulk` | counts summed per child |
| `23_…_normalise` | filter, PFlog1pPF |
| `24_…_soft_threshold` | WGCNA power, children with SLE |
| `25_…_modules` | WGCNA modules |
| `25b_…_interferon_genes` | the 28 interferon response genes and 11 NF-κB control genes against the modules |
| `26_…_array_module_preservation` | are the GSE65391 modules present here? claim N4 |
| `26b_…_array_module_projection` | children scored on the GSE65391 modules |
| `27_…_module_traits` | module scores vs clinical traits; SLE vs healthy; claims N1–N4 |
| `27a_…_gap_statistic` | k-means k for children and genes |
| `28_…_heatmap_hclust` | Ward trees on children and genes |
| `29_…_heatmap_kmeans` | k-means with gap-statistic k |
| `29b_…_dynamic_tree_cut` | Dynamic Tree Cut on the Ward trees |
| `29c_…_pvclust_genes` (Python) | AU p-values and BP, genes |
| `29d_…_pvclust_children` (Python) | AU p-values and BP, children; two-way figure |

### Federation
| notebook | does |
|---|---|
| `30_federation_shared_traits` | which clinical traits the studies share (proposal) |

## Repository layout

```
ipynb/      notebooks
src/        paths.R: locations only
run_all.sh  executes the notebooks in order
genes/      curated gene lists: six-gene interferon score; 28-gene interferon response gene score and
            its 25- and 3-gene parts, and 11 NF-κB control genes (de Jesus et al. 2020, supplement)
data/       downloads and run artifacts; not committed, rebuilt by the notebooks
figures/    a 300 dpi PNG of every figure, by study; not committed, rebuilt by the notebooks
```
