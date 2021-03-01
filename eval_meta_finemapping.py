import numpy as np
import pandas as pd
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--config", action="store", type=int, required=True)
parser.add_argument("--nbins", action="store", type=int, required=True, dest="n_bins")
args = parser.parse_args()

n_bins = args.n_bins
config = args.config
path_prefix = "/humgen/atgu1/methods/mkanai/workspace/202011_meta_finemapping"
n_phenos = 100

sum_pips = [0] * n_bins
n_causals = [0] * n_bins
bin_counts = [0] * n_bins
# total number of simulated causal SNPs
n_total_causals = 0
n_total_snps = 0

for pheno in range(1, n_phenos+1):

    true_beta = pd.read_table(
        f"{path_prefix}/true_pheno/pheno{pheno}.beta.liftover.txt", 
        usecols=["variant", "beta"]
    )

    try:
        fm_results = pd.read_table(
            f"{path_prefix}/ABF/meta.pheno{pheno}.config{config}.ABF.snp.bgz", 
            usecols=["rsid", "prob"], 
            compression="gzip"
        )
    except FileNotFoundError:
        print(f".snp.bgz file not found for pheno={pheno}, config={config}. Skipping.")
        continue
    
    pips = fm_results["prob"].to_numpy()
    
    n_total_causals += len(true_beta)
    n_total_snps += len(fm_results)

    
    # calculate calibration
    bin_width = 1/n_bins
    bins = np.arange(0, 1+bin_width, bin_width)
    bins[-1] = 1.1  # to include PIP of 1 in top bin
    for i in range(len(bins) - 1):
        in_bin = np.logical_and(pips >= bins[i], pips < bins[i+1])
        sum_pip_in_bin = fm_results.loc[in_bin, "prob"].sum()
        n_causal_in_bin = fm_results.loc[in_bin, "rsid"].isin(true_beta.variant).sum()
        bin_count_in_bin = in_bin.sum()
        
        sum_pips[i] += sum_pip_in_bin
        n_causals[i] += n_causal_in_bin
        bin_counts[i] += bin_count_in_bin
        
# calculate power 
powers = []
for b in range(n_bins):
    # can also divide by n_total_causals but I think that would be a less interpretable metric
    power = sum(n_causals[b:]) / sum(n_causals) 
    powers.append(power)

# sanity check 
# n finemapped causal variants <= total n causal variants (some missing due to not passing GWAS)
assert sum(n_causals) <= n_total_causals
# make sure binning was accurate
assert sum(bin_counts) == n_total_snps
    
# outut results
bin_counts = np.array(bin_counts)
mean_pips = np.array(sum_pips) / bin_counts
prop_causals = np.array(n_causals) / bin_counts
powers = np.array(powers)

df = pd.DataFrame({
    "bin_count": bin_counts,
    "mean_pip": mean_pips,
    "prop_causal": prop_causals,
    "power": powers
})

df.to_csv(
    f"/broad/finucanelab/relzur/meta_analysis_fm/meta.config{config}.eval.txt", 
    sep="\t", 
    index_label="bin"
)
