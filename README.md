# pllp RNA-seq analysis — Wang, Childers, Martinez, Bagnat & Park (PLoS Genetics)

"Genome wide transcriptional changes underlie gradual and recurrent adaptation to
protein malnutrition in zebrafish."

Scripts that reproduce every RNA-seq result and figure panel in the paper from the
data deposited in GEO **GSE320085**. Two independent bulk RNA-seq datasets are analysed:

| | Dataset 1 (main figures) | Dataset 2 (Fig S7) |
|---|---|---|
| Design | EK, adapted *pd1116*, newly generated *pd1245* × HP / LP diet, 20 dpf, 3 tank-pools per group, **18 samples** | EK vs adapted F8 *pd1245* × HP / LP, 20 dpf, **14 samples** |
| Library / sequencing | stranded mRNA-seq, NovaSeq 6000, 50 bp PE | 3'-biased UMI mRNA-seq (Plasmidsaurus) |
| Upstream processing (done, not repeated here) | FastQC 0.11.7 → Trim Galore 0.6.5 → STAR 2.7.3a (GRCz11, Ensembl 100, unique mappers) → HTSeq 0.11.1 | fastp 0.24.0 → STAR 2.7.11 → UMIcollapse 1.1.0 → featureCounts (subread 2.1.1) |
| GEO input used by these scripts | one HTSeq count file per sample | `TLHYNJ-expression-matrix.tsv` (`*_count` columns) |
| Paper panels | Fig 2A–E, 3A–B, 5A–B, 5G, S2A–D, S3A–B, S4A; S1–S2 Tables | Fig S7D–E |

## 1. Set up

```bash
git clone https://github.com/jieunesther/pllp_RNAseq.git
cd pllp_RNAseq
# dataset 1: the 18 HTSeq count files (may be .gz, may carry a GSM prefix)
#   -> data/dataset1_htseq_counts/
# dataset 2: TLHYNJ-expression-matrix.tsv
#   -> data/dataset2_TLHYNJ/
```
Each `data/*/README.md` states exactly what the folder expects. File names must contain
the `sample_id` values listed in `metadata/dataset1_samples.csv` /
`metadata/dataset2_samples.csv`.

**R packages** (R ≥ 4.4; Bioconductor 3.20 was used for the runs recorded here):
`DESeq2`, `org.Dr.eg.db`, `AnnotationDbi`, `clusterProfiler`, `enrichplot`, `WGCNA`,
`ComplexHeatmap`, `circlize`, `VennDiagram`, `ggplot2`, `dplyr`, `stringr`,
`matrixStats`, `gplots`, `RColorBrewer`, `ragg`, `systemfonts`.
```r
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "org.Dr.eg.db", "clusterProfiler", "enrichplot", "WGCNA",
                       "ComplexHeatmap", "circlize", "VennDiagram", "matrixStats", "gplots", "ragg"))
```
`02` and `03` query the KEGG REST API, so they need internet access.

## 2. Run

```bash
bash run_all.sh                 # everything, ~15 min on one CPU
bash run_all.sh 05 06           # only selected steps
RSCRIPT=/path/to/Rscript bash run_all.sh
```
or `Rscript scripts/<name>.R` one at a time, in numeric order (paths are resolved
from the script location, so the working directory does not matter).

| Step | Script | Reads | Writes | Paper |
|---|---|---|---|---|
| 01 | `01_dataset1_DESeq2.R` | `data/dataset1_htseq_counts/` | `results/dataset1/` — raw / normalised / VST matrices, `gene_annotation.csv`, 9 pairwise DESeq2 contrasts (`DESeq2_results_*`, `DESeq2_significant_*` = **S1 Table**), interaction terms, `DEG_summary_table.csv`, overlaps, `genes_unique_to_*` (**S2 Table**) | Fig 2B, 2C, S2A–C |
| 02 | `02_dataset1_enrichment.R` | `results/dataset1/` | `results/dataset1_enrichment/` — GO BP (simplified) and KEGG compareCluster tables + `.rds` | Fig 5A, 5B, S3A |
| 03 | `03_dataset1_WGCNA.R` | `results/dataset1/vst_transformed_counts.csv` | `results/dataset1_WGCNA/` — module assignments, eigengenes, module–trait correlations / p-values, per-module GO + KEGG, `Fig2E_KEGG_top5_per_module.csv` | Fig 2D, 2E, S2D |
| 04 | `04_dataset2_DESeq2.R` | `data/dataset2_TLHYNJ/` | `results/dataset2/` — VST, PCA, *pd1245* vs EK within diet | Fig S7D |
| 05 | `05_figure2.R` | 01, 03 | `figures/Figure2/` | Fig 2A–E |
| 06 | `06_figure5AB.R` | 02 | `figures/Figure5/` | Fig 5A, 5B |
| 07 | `07_heatmaps_curated_genes.R` | 01 + `curated_genes/` | `figures/Figure3/`, `Figure5/`, `FigureS4/`; matrices in `results/heatmap_matrices/` | Fig 3A, 3B, 5G, S4A |
| 08 | `08_figureS2.R` | 01, 03 | `figures/FigureS2/` | Fig S2A–D |
| 09 | `09_figureS3B.R` | 01 + `curated_genes/` | `figures/FigureS3/` | Fig S3B |
| 10 | `10_figureS7.R` | 01, 04 | `figures/FigureS7/` | Fig S7D, S7E |

Every script prints sanity checks (input dimensions, group composition, filter
counts, spot-check values), verifies that its output files exist and are non-empty,
and writes `results/sessionInfo_<script>.txt`.

## 3. Reproducibility checks

`metadata/expected_values.csv` lists the numbers printed in the paper (DEG counts per
contrast, overlaps, PCA variance, WGCNA gene number / power / module counts, panel
gene counts). Each script compares what it computed against them and appends the
outcome to `results/reproducibility_check.csv`; `run_all.sh` prints any mismatch.

What is and is not expected to reproduce bit-for-bit:

* **DESeq2 (01, 04)** — exact. Re-running with DESeq2 1.46 gives the published
  32,520 → 22,863 gene filter, all 9 × (total / up / down) DEG counts, the 2,279-gene
  union, every Venn overlap and PC1/PC2 = 48 % / 25 % (dataset 2: 37 % / 23 %).
  Log2 fold changes and adjusted p-values agree with the original run to < 1e-8.
* **WGCNA (03)** — exact given the same input: 5,716 genes, power 18, 13 modules of
  identical size, 12 non-grey modules, 7 trait-significant modules, 13 starred cells.
* **Gene symbols** — `01` names rows by `org.Dr.eg.db` symbol (Ensembl ID kept where
  no symbol exists). Symbols for ~100 of the 22,863 genes differ between
  `org.Dr.eg.db` releases. All curated gene lists in `curated_genes/` therefore carry
  Ensembl IDs, and every DE table carries an `ensembl_id` column, so the heatmap panels
  are unaffected. `results/dataset1/gene_annotation.csv` records the mapping used.
* **GO / KEGG enrichment (02, 03)** — depends on the annotation release. GO term sets
  come from `org.Dr.eg.db` + `GO.db`; KEGG pathways are downloaded live and change
  between releases (the release used is written to `KEGG_release.txt`). The selection
  of top-10 GO terms per column (Fig 5A/B) and the KEGG pathways for
  *pd1245*-unique DEGs (Fig S3A) can therefore gain or lose a few terms relative to
  the published panels, which were generated in February 2026. The tables that were
  actually plotted are kept under `published_results/` so the published panels can be
  regenerated verbatim:
  ```bash
  PLLP_USE_PUBLISHED=1 Rscript scripts/05_figure2.R     # Fig 2E from the published KEGG table
  PLLP_USE_PUBLISHED=1 Rscript scripts/06_figure5AB.R   # Fig 5A/5B from the published GO tables
  ```
  Observed on 2026-09-29 (Bioconductor 3.20, KEGG release of that day): Fig 5A keeps
  11 terms with 3 exchanged, Fig 5B draws 17 instead of 19 terms, Fig 2E exchanges 4 of
  35 pathways (Cadherin signaling and Efferocytosis enter; Necroptosis and FoxO signaling
  leave). All other panels and every number in the paper reproduce exactly.

## 4. Methods summary (as implemented)

* Filter: genes with ≥ 10 counts in ≥ 3 samples (dataset 1: 32,520 → 22,863).
* DESeq2 (Wald): full model `~ genotype + diet + genotype:diet` for VST (blind = FALSE),
  PCA (`plotPCA`, 500 most variable genes) and the genotype × diet interaction terms;
  group model `~ condition` for all pairwise contrasts. DEG: BH-adjusted p < 0.05 and
  |log2FC| > 1.
* ORA (Fig 5A/B, S3A): `clusterProfiler::compareCluster`; GO BP by symbol with the
  22,863 expressed genes as background, then `simplify(cutoff = 0.6)`; KEGG (`dre`)
  with the Entrez-mapped expressed genes as background; BH p < 0.05.
* WGCNA: top 25 % most variable VST genes (5,716); signed-hybrid network, signed TOM,
  soft power 18 (first power with scale-free R² ≥ 0.85), `blockwiseModules`
  (minModuleSize 30, mergeCutHeight 0.25, one block); module eigengenes correlated
  (Pearson) with binary genotype / diet indicators, Student p-values; per-module GO BP
  and KEGG ORA against the organism-wide background, BH p < 0.05. Fig 2E shows, for each
  trait-significant module, the five significant pathways with the largest gene overlap.
* Curated-gene heatmaps (Fig 3A, 3B, 5G, S4A, S7E): per-gene z-score of VST across all
  samples, averaged per genotype × diet group; markers from the DESeq2 contrasts
  (`*` EK vs *pd1116*, `#` *pd1116* vs *pd1245*, padj < 0.05, within diet).

## 5. Repository layout

```
R/            common.R (paths, checks), fig_style.R (journal style: Arial 8-12 pt, ≤ 7.5 × 8.75 in, 600 dpi TIFF)
scripts/      01-10 as above
metadata/     sample sheets, expected_values.csv
curated_genes/  gene lists for Fig 3A, 3B, 5G, S3B, S4A (group, symbol, Ensembl ID)
published_results/  enrichment tables exactly as plotted in the paper (Feb 2026 run)
legacy_v1.0.0/  the three scripts of the first release (kept for the record; superseded)
data/, results/, figures/   not tracked (inputs from GEO, generated outputs)
```

## 6. Data availability

Raw FASTQ files and the processed count files for both datasets: GEO **GSE320085**.
Code: this repository, archived at Zenodo (concept DOI 10.5281/zenodo.18707422).
