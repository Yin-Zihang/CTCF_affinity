# CTCF affinity: threshold fitting, compartments, loops, and TADs

This repository contains the scripts and processed tables used to relate measured CTCF–DNA affinity to occupancy across cell types, A/B compartments, chromatin loops, and TAD boundaries. The examples below are written for execution from the **repository root**. Each command is a single line.

## Repository layout

| Path | Contents |
| --- | --- |
| `input/` | Files supplied to the analysis scripts. See **Input data** for the expected names and formats. |
| `scripts/compartment/` | Shared CBS input parser, A/B peak analysis, and R density plotting. |
| `scripts/loop_TAD/` | ChIA-PET, Hi-C loop, and TAD-boundary analyses. |
| `scripts/thresholds/` | Affinity-rank fitting, bootstrap, sensitivity analyses, and plotting. |
| `processed_data/compartment/` | `GM12878_ab_peak_affinity.tsv`, `HUVEC_ab_peak_affinity.tsv`, and `K562_ab_peak_affinity.tsv`. |
| `processed_data/loop_TAD/` | `ChIApet_affinity_normalization.tsv`, `HiCLoop_affinity_normalization.tsv`, and `TAD_enrichment.tsv`. |
| `processed_data/threshold/` | `EMSAgroupvsFrequency.csv`. Add the fitting outputs listed below when depositing the complete threshold analysis. |

The `CBS_ab_compartments.py` module imports `CBS_compartment_common.py` by name. Keep these two files together in `scripts/compartment/`. If you renamed the original scripts yourself, check that the import statement in `CBS_ab_compartments.py` reads `importlib.import_module('CBS_compartment_common')`; renaming the file alone does not update that statement.

## 1. A/B compartment analysis

Each peak is assigned by its midpoint to an A or B compartment according to the sign of the eigenvector. The script retains peaks containing at least one matched CBS center. `EMSA` in the peak-level output is the **sum** of the matched CBS affinities; `ChIP` is the mean of the two replicate RPKM values.

GM12878:

`python scripts/compartment/CBS_ab_compartments.py --cell-type GM12878 --eigenvector input/GSE63525_GM12878_insitu_DpnII_combined_30_eigenvector_KR100kb_updated.wig --sites-bed input/final_split_to_12_groups.bed --peaks input/CTCFPeak_GSE30263_Stam_GM12878_CTCF_rep1plusrep2_htseq_counts.tsv --reads-rep1 1115584 --reads-rep2 3179032 --outdir processed_data/compartment`

HUVEC:

`python scripts/compartment/CBS_ab_compartments.py --cell-type HUVEC --eigenvector input/GSE63525_HUVEC_combined_30_eigenvector_KR100kb_updated.wig --sites-bed input/final_split_to_12_groups.bed --peaks input/CTCFPeak_GSE30263_Stam_HUVEC_CTCF_rep1plusrep2_htseq_counts.tsv --reads-rep1 3149985 --reads-rep2 2094431 --outdir processed_data/compartment`

K562:

`python scripts/compartment/CBS_ab_compartments.py --cell-type K562 --eigenvector input/GSE63525_K562_combined_30_eigenvector_KR100kb_updated.wig --sites-bed input/final_split_to_12_groups.bed --peaks input/CTCFPeak_GSE30263_Stam_K562_CTCF_rep1plusrep2_htseq_counts.tsv --reads-rep1 3393447 --reads-rep2 1459969 --outdir processed_data/compartment`

### Density plots and plotted counts

The R script reads the processed peak tables. `total` in each panel counts all retained peak rows, whereas `shown` counts only rows within the requested `ChIP` and `EMSA` axis limits. It also writes `*_count_audit.tsv` next to the requested figure.

GM12878:

`Rscript scripts/compartment/CBS_compartment_density_counts.R --table processed_data/compartment/GM12878_ab_peak_affinity.tsv --output figures/compartment/GM12878_ab_density.pdf --layout grid --max-chip 200 --max-emsa 7 --nbin 1600`

HUVEC:

`Rscript scripts/compartment/CBS_compartment_density_counts.R --table processed_data/compartment/HUVEC_ab_peak_affinity.tsv --output figures/compartment/HUVEC_ab_density.pdf --layout grid --max-chip 200 --max-emsa 7 --nbin 1600`

K562:

`Rscript scripts/compartment/CBS_compartment_density_counts.R --table processed_data/compartment/K562_ab_peak_affinity.tsv --output figures/compartment/K562_ab_density.pdf --layout grid --max-chip 200 --max-emsa 7 --nbin 1600`


## 2. ChIA-PET loops, Hi-C loops, and TAD boundaries

The three scripts classify CBSs as strong (`EMSA > 4.07`), middle (`2.25 < EMSA <= 4.07`), or weak (`EMSA <= 2.25`). 

`python ChIApet_affinity_normalization_tsv.py`

`python HiCLoop_affinity_norm_tsv.py`

`python TADS_oe_tsv.py`

Each script reads its input paths from the defaults in the code and writes its TSV to the **current working directory**. The resulting tables were placed in `processed_data/loop_TAD/` for this repository. 

## 3. Occupancy–affinity threshold fitting

The threshold script uses 61 cell or tissue types and 1,000 within-bin bootstrap resamples. The horizontal fitting coordinate is the **affinity rank**, not the EMSA value. After fitting, threshold ranks are mapped back to affinity using physical rows of the ranked TSV.

Run the analysis from the repository root:

`python scripts/thresholds/13review_bootstrap_thresholds.py --work-dir . --count-file input/12_count_dup_cell_new.txt --group-dir input/11_split_to_270_groups --outdir processed_data/threshold/bootstrap_results --prefix threshold --n-cells 61 --n-groups 280 --bootstrap 1000 --seed 20260927`

Recreate plots from existing fitting tables without another bootstrap:

`python scripts/thresholds/13review_bootstrap_thresholds_plot_all.py --work-dir . --outdir processed_data/threshold/bootstrap_results --prefix threshold --plot-existing`

