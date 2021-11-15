"""Computes power and calibration for a single meta-analysis configuration."""

import numpy as np
import pandas as pd
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--config", action="store", type=int, required=True)
args = parser.parse_args()

config = args.config

bins = [0, 0.01, 0.1, 0.5, 0.9, 1.1]
n_bins = len(bins) - 1
path_prefix = "/humgen/atgu1/methods/mkanai/workspace/202011_meta_finemapping"
n_phenos = 100

sum_pips = [0] * n_bins
n_causals = [0] * n_bins
bin_counts = [0] * n_bins

all_pips = np.array([])
causal_status = np.array([], dtype=bool)

# total number of simulated causal SNPs, including QCed out variants,
# excluding variants in regions not fine-mapped
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
            usecols=["rsid", "prob", "region"], 
            compression="gzip"
        )
    except FileNotFoundError:
        print(f".snp.bgz file not found for pheno={pheno}, config={config}. Skipping.")
        continue

    # only count causal variants that are in fine-mapped regions
    finemapped_regions = fm_results["region"].unique().tolist()
    for v in true_beta["variant"]:
        for region in finemapped_regions:
            chromosome = region.split(":")[0]
            start, end = tuple([int(i) for i in region.split(":")[1].split("-")])
            if (v.split(":")[0] == chromosome and 
                int(v.split(":")[1]) >= start and
                int(v.split(":")[1]) <= end):
                n_total_causals += 1
                break
    
    n_total_snps += len(fm_results)
    all_pips = np.append(all_pips, fm_results["prob"].to_numpy())
    causal_status = np.append(causal_status, fm_results["rsid"].isin(true_beta.variant))

# calculate calibration
for i in range(len(bins) - 1):
    in_bin = np.logical_and(all_pips >= bins[i], all_pips < bins[i+1])
    sum_pips[i] = all_pips[in_bin].sum()
    n_causals[i] = causal_status[in_bin].sum()
    bin_counts[i] = in_bin.sum()

# calculate power 
fractions = [0.0001, 0.0005, 0.001, 0.005, 0.01]
powers = []
sorted_causal = causal_status[np.argsort(all_pips)]
for f in fractions:
    # get number of causal SNPs in top fraction of PIPs
    top_f = sorted_causal[-round(f*n_total_snps):].sum()
    powers.append(top_f/n_total_causals)
    
# sanity checks
assert causal_status.sum() <= n_total_causals
assert sum(bin_counts) == n_total_snps
    
# outut results
bin_counts = np.array(bin_counts)
mean_pips = np.array(sum_pips) / bin_counts
prop_causals = np.array(n_causals) / bin_counts
powers = np.array(powers)

df = pd.DataFrame({
    "bin_left_edge": bins[:-1],
    "bin_count": bin_counts,
    "mean_pip": mean_pips,
    "prop_causal": prop_causals,
    "power_fraction": fractions,
    "power": powers
})

df.to_csv(
    f"/broad/finucanelab/relzur/meta_analysis_fm/eval.config{config}.txt", 
    sep="\t", 
    index=False
)
