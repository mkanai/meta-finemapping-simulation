#!/usr/bin/env python

import os
import pandas as pd
import numpy as np
from scipy.stats import chi2
import pickle
import hail as hl
from hail.linalg import BlockMatrix as bm

path_pref = "/humgen/atgu1/methods/mkanai/workspace/202011_meta_finemapping"
dentist_pref = "/broad/finucanelab/relzur/meta_analysis_fm/gnomad_filtered_chr3_LD"

def compute_p_dentist(z, r, i):
    """
    compute DENTIST test statistic and associated p value
    
    Args:
        z - numpy array of z scores for the locus
        r - numpy array of ld to lead variant
        i - index of lead variant
        
    Returns: 
        dentist T statistics
        -log10 p values for those statistics
    """
    
    assert z.shape == r.shape
    
    # remove causal variant
    remove_i = np.array([idx for idx in range(len(z)) if idx != i])
    z_i = z[i]
    z = z[remove_i]
    r = r[remove_i]
    
    # as described in the dentist manuscript
    # for computational stability, remove variants in near-perfect LD with causal variant i (r2 > 0.95) 
    high_ld = r**2 > 0.95
    z = z[~high_ld]
    r = r[~high_ld]
    
    t = ((z - (z_i * r)) ** 2) / (1 - (r ** 2))
    
    # t is approximately distributed chi squared with 1 dof
    dof = 1
    log10_p = chi2.logsf(t, dof) / -np.log(10)
    
    return t, log10_p


def main(pheno, config):

    pop_dict = {
        "AFR": "afr",
        "EUR": "nfe",
        "EAS": "eas"
    }

    results_df = pd.DataFrame(
        columns=["config", "pheno", "region", "t_dentist_fill", "log10_p_dentist_fill", "t_dentist_drop", "log10_p_dentist_drop"]
    )

    studies_per_config = pd.read_csv(f"{path_pref}/config/meta_config_hg19.tsv", sep="\t", header=None)
    pops = list(map(lambda x: pop_dict[x.split("_")[0]], studies_per_config.iloc[config-1, :]))
    pops, counts = np.unique(pops, return_counts=True)

    true_causals = pd.read_table(f"{path_pref}/true_pheno/pheno{pheno}.beta.txt", usecols=["variant", "gamma"])
    abf_results = pd.read_table(
        f"{path_pref}/ABF_hg19/meta.pheno{pheno}.config{config}.ABF.snp.bgz", 
        compression="gzip"
    )
    abf_results = pd.merge(abf_results, true_causals, how="left", on="variant")
    abf_results["gamma"] = abf_results["gamma"].fillna(False)
    abf_results["z"] = abf_results["beta"] / abf_results["se"]

    for region in abf_results.region.unique():

        abf_results_region = abf_results.loc[abf_results.region == region].copy()

        # get lead varaiant
        i = np.argmin(abf_results_region.p.to_numpy())
        var_id_i = abf_results_region.iloc[i, abf_results_region.columns.get_loc("variant")]

        # get LD
        lead_variant_missing = False # indicator that lead variant is missing from a population

        for pop, count in zip(pops, counts):

            # read ld and variant table
            ld = bm.read(
                f"{dentist_pref}/gnomad_v2.1.1_chr3_{pop}_ld_filtered.bm"
            )
            var = hl.read_table(
                f"{dentist_pref}/gnomad_v2.1.1_chr3_{pop}_var_filtered.ht"
            )

            # get lead index
            lead_var = var.filter(var.var_id == var_id_i)
            if lead_var.count() == 1:
                lead_idx = lead_var.idx.collect()[0]
            else:
                # Lead SNP not in population gnomAD LD
                lead_variant_missing = True
                print(
                    f"ERROR - none or more than one lead variant {var_id_i} found in LD matrix for " 
                    f"pheno: {pheno} config: {config} region: {region} pop: {pop}"
                )
                break

            # get indices for rest of locus
            var = hl.read_table(
                f"{dentist_pref}/gnomad_v2.1.1_chr3_{pop}_var.ht"
            )
            locus_vars = abf_results_region.variant.to_list()
            var = var.filter(hl.literal(set(locus_vars)).contains(var.var_id))
            ids = var.var_id.collect()
            idxs = var.idx.collect()

            ld = ld.filter_rows(idxs)
            ld = ld[:, lead_idx].to_numpy().squeeze()

            df = pd.DataFrame({
                "variant": ids,
                f"{pop}_ld": ld
            })

            abf_results_region = pd.merge(abf_results_region, df, how="left", left_on="variant", right_on="variant")

        if lead_variant_missing:
            results_df.loc[len(results_df)] = [config, pheno, region, None, None, None, None]

        else:
            # two approaches - fill missing LD with zeros, or restrict to only variants in all population LD matrices

            # fill with zeros 
            # (I suspect this will give you issues if the true population LD actually isn't zero, 
            # the SNP is just missing from gnomAD - you'll get significant DENTIST p values)
            abf_results_region_fill = abf_results_region.fillna(0)

            # restrict to shared SNPs between different population LD matrices
            abf_results_region_drop = abf_results_region.dropna(axis="index", subset=[f"{pop}_ld" for pop in pops])

            # compute meta-analyzed LD
            r_fill = np.zeros(len(abf_results_region_fill))
            r_drop = np.zeros(len(abf_results_region_drop))
            for pop, count in zip(pops, counts):
                r_fill = r_fill + count * abf_results_region_fill[f"{pop}_ld"]
                r_drop = r_drop + count * abf_results_region_drop[f"{pop}_ld"]

            r_fill = r_fill.to_numpy() / counts.sum()
            r_drop = r_drop.to_numpy() / counts.sum()

            # get z
            z_fill = abf_results_region_fill.z.to_numpy()
            z_drop = abf_results_region_drop.z.to_numpy()

            # get i
            i_fill = np.argmin(abf_results_region_fill.p.to_numpy())
            i_drop = np.argmin(abf_results_region_drop.p.to_numpy())

            # run dentist
            t_dentist_fill, log10_p_dentist_fill = compute_p_dentist(z_fill, r_fill, i_fill)
            t_dentist_drop, log10_p_dentist_drop = compute_p_dentist(z_drop, r_drop, i_drop)        

            # append results to dataframe
            results_df.loc[len(results_df)] = [
                config, pheno, region, t_dentist_fill, log10_p_dentist_fill, t_dentist_drop, log10_p_dentist_drop
            ]

    # write results dataframe
    results_df.to_hdf(
        f"/broad/finucanelab/relzur/meta_analysis_fm/dentist_results/pheno{pheno}_config{config}_dentist.hdf5",
        "dentist",
        mode="w"
    )


if __name__ == "__main__":

    task_id = int(os.getenv("SGE_TASK_ID")) - 1
    assert task_id is not None
    task_table = pd.read_csv("/broad/finucanelab/relzur/meta_analysis_fm/dentist_task_table.csv")
    pheno = task_table.loc[task_id, "pheno"]
    config = task_table.loc[task_id, "config"]
    hl.init(log=f"/broad/hptmp/relzur/log/apply_dentist_{task_id}.log")

    main(pheno, config)


