Dataset 1 (main figures): the 18 per-sample HTSeq count files of GEO GSE320085,
samples GSM9533663 - GSM9533680 (titles EK_HP_1_S1 ... pd1245_LP_25_S18).
Untar GSE320085_RAW.tar here; files may be gzipped and may carry the GSM
prefix (e.g. GSM9533663_EK_HP_1_S1.count.gz). Script 01 matches each file to
metadata/dataset1_samples.csv by GSM accession or by sample title.
Alternatively a single matrix named raw_counts.csv (row = Ensembl gene ID,
one column per sample title) is accepted.
Do not put the 14 *_FigureS7RNAseq samples here; they belong to dataset 2.
