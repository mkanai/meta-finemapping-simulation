library(dplyr)
library(magrittr)

pops = c("AFR", "EAS", "EUR", "SAS")
h2g = 0.05
alpha = -0.38
locus_length = 3e6
frac_causal_locus = 0.5

df = data.table::fread("/humgen/atgu1/methods/mkanai/workspace/202011_meta_finemapping/1000GP_Phase3/chr3/1000GP_Phase3_chr3.lof_missense.tsv", data.table = F) %>%
  dplyr::filter(is_canonical_vep & !stringr::str_detect(gene_most_severe, "\\.")) %>%
  dplyr::select(variant, most_severe, gene_most_severe, lof)

df_frq = purrr::map_dfr(pops, function(pop) {
  data.table::fread(sprintf("/humgen/atgu1/methods/mkanai/workspace/202011_meta_finemapping/1000GP_Phase3/chr3/1000GP_Phase3_chr3.%s.afreq", pop), data.table = F) %>%
  dplyr::filter(ID %in% df$variant) %>%
  dplyr::mutate(variant = ID, pop = pop, freq = 1 - ALT_FREQS) %>%
  dplyr::select(variant, pop, freq)
})

df_frq =
  dplyr::bind_rows(
    df_frq,
    dplyr::group_by(df_frq, variant) %>%
    dplyr::summarize(pop = "popmax", freq = max(freq)) %>%
    dplyr::ungroup()
  ) %>%
  tidyr::pivot_wider(id_cols = "variant", names_prefix = "freq_", names_from = "pop", values_from = "freq")
  
df =
  dplyr::left_join(df, df_frq) %>%
  dplyr::filter(freq_popmax > 0.01 & freq_popmax < 0.99) %>%
  dplyr::mutate(
    position = as.numeric(stringr::str_split_fixed(variant, ":", 4)[,2]),
    locus = position %/% locus_length
  )


compute_true_beta = function(df) {
  x = dplyr::left_join(
    df,
    tibble::tibble(
      locus = unique(df$locus),
      # causal_locus = rbinom(length(locus), 1, frac_causal_locus) > 0)
      causal_locus = ((locus %% 2) == rbinom(1, 1, 0.5)))
    ) %>%
    dplyr::group_by(locus) %>%
    dplyr::mutate(gamma = causal_locus & variant == sample(variant, 1, prob = ifelse(most_severe == "missense_variant", 1, 10))) %>%
    dplyr::ungroup() %>%
    dplyr::filter(gamma) %>%
    dplyr::arrange(position)
  
  sigma2_g = h2g / sum((2 * x$freq_popmax * (1 - x$freq_popmax)) ** (1 + alpha))
  x =
    dplyr::mutate(x,
      var_beta = sigma2_g * (2 * freq_popmax * (1 - freq_popmax)) ** alpha,
      beta = purrr::map_dbl(var_beta, ~{rnorm(1, 0, sd=sqrt(.))}),
      ref = stringr::str_split_fixed(variant, ":", 4)[, 3],
      alt = stringr::str_split_fixed(variant, ":", 4)[, 4]
    )
  x = dplyr::select(x, variant, ref, alt, beta, setdiff(colnames(x), c("variant", "ref", "alt", "beta")))
  return(x)
}

purrr::map(1:100, function(i) {
  x = compute_true_beta(df)
  write.table(x, sprintf("/humgen/atgu1/methods/mkanai/workspace/202011_meta_finemapping/true_pheno/pheno%d.beta.txt", i), quote = F, row.names = F, sep = "\t")
})
