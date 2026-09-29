Result tables exactly as used for the published panels (original run, February
2026, R 4.x / Bioconductor 3.20 on macOS; DESeq2, clusterProfiler, org.Dr.eg.db and
KEGG releases of that date). They are kept so that panels whose content depends on
annotation releases (GO / KEGG enrichment) can be regenerated verbatim:

  Enrichment_GO_{HP,LP,HPvsLP}.csv     compareCluster(enrichGO) + simplify(0.6)   -> Fig 5A, 5B
  Enrichment_KEGG_{HP,LP,HPvsLP}.csv   compareCluster(enrichKEGG)
  Fig5AB_plotted_terms.csv             the rows drawn in Fig 5A/5B
  Fig2E_KEGG_top5_per_module.csv       the 35 dots of Fig 2E
  WGCNA_module_trait_{correlations,pvalues}.csv   Fig 2D
  DEG_summary_table.csv                Fig 2C

The KEGG table behind Fig S3A (DEGs unique to pd1245 vs EK) was not saved in the
original run; re-running scripts/02 regenerates it with the current KEGG release.

Set PLLP_USE_PUBLISHED=1 before running scripts/06_figure5AB.R to draw Fig 5A/5B
from these tables instead of from results/dataset1_enrichment/.
Note: cluster names in Enrichment_GO_*.csv use the short form "1116_HP_Up"; the
scripts map them to the long form "1116_vs_EK_HP_Up".
