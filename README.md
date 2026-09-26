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
method, and what counts as a finding.

### GSE65391 · RNA_array
| notebook | does |
|---|---|
| `00_GSE65391_RNA_array_download` | series matrix, platform annotation, NCBI gene_info |
| `01_GSE65391_RNA_array_metadata` | clinical fields; SLEDAI stage per visit; first visit per child |
| `02_GSE65391_RNA_array_probes_to_genes` | probes to genes; expressed genes |
| `03_GSE65391_RNA_array_batch_check` | technical replicates across the two array batches |
| `04_GSE65391_RNA_array_soft_threshold` | WGCNA power, first visit per child |
| `05_GSE65391_RNA_array_modules` | WGCNA modules and hub genes |
| `06_GSE65391_RNA_array_module_traits` | module eigengenes vs clinical traits |
| `07_GSE65391_RNA_array_heatmap_hclust` | hclust on genes and patients |
| `08_GSE65391_RNA_array_heatmap_kmeans` | k-means on genes and patients |

### GSE232381 · bulk_RNA_seq
| notebook | does |
|---|---|
| `10_GSE232381_bulk_RNA_seq_download` | NCBI-generated raw counts (16 of 25 samples), series matrix |
| `11_GSE232381_bulk_RNA_seq_metadata` | active / inactive lupus nephritis |
| `12_GSE232381_bulk_RNA_seq_normalise` | filter, PFlog1pPF |
| `13_GSE232381_bulk_RNA_seq_soft_threshold` | WGCNA power (exploratory, n = 16) |
| `14_GSE232381_bulk_RNA_seq_modules` | WGCNA modules (exploratory) |
| `15_GSE232381_bulk_RNA_seq_array_module_preservation` | are the GSE65391 modules present here? |
| `16_GSE232381_bulk_RNA_seq_module_traits` | eigengenes vs nephritis activity |
| `17_GSE232381_bulk_RNA_seq_heatmap_hclust` | hclust on genes and samples |
| `18_GSE232381_bulk_RNA_seq_heatmap_kmeans` | k-means on genes and samples |

### GSE135779 · scRNA_seq as pseudobulk
| notebook | does |
|---|---|
| `20_GSE135779_scRNA_seq_download` | per-donor count matrices; checks the supplement |
| `21_GSE135779_scRNA_seq_metadata` | GEO fields joined to Supplementary Table 1b |
| `22_GSE135779_scRNA_seq_pseudobulk` | counts summed per child |
| `23_GSE135779_scRNA_seq_normalise` | filter, PFlog1pPF |
| `24_GSE135779_scRNA_seq_soft_threshold` | WGCNA power, children with SLE |
| `25_GSE135779_scRNA_seq_modules` | WGCNA modules |
| `26_GSE135779_scRNA_seq_array_module_preservation` | are the GSE65391 modules present here? |
| `27_GSE135779_scRNA_seq_module_traits` | eigengenes vs clinical traits; SLE vs healthy |
| `28_GSE135779_scRNA_seq_heatmap_hclust` | hclust on genes and children |
| `29_GSE135779_scRNA_seq_heatmap_kmeans` | k-means on genes and children |

### Federation
| notebook | does |
|---|---|
| `30_federation_shared_traits` | which clinical traits the studies share (proposal) |

## Repository layout

```
ipynb/      notebooks
src/        paths.R: locations only
run_all.sh  executes the notebooks in order
genes/      curated gene lists (interferon score genes)
data/       downloads and run artifacts; not committed, rebuilt by the notebooks
```
