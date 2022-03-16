library(dplyr)
library(furrr)
library(magrittr)

plan("multisession")
basedir <- "/humgen/atgu1/methods/mkanai/workspace/202011_meta_finemapping/new_simulations"

pops <- c("AFR", "EAS", "EUR")
n_cohorts_per_pop <- 10
h2g <- 0.03
alpha <- -0.38
locus_length <- 3e6
frac_causal_locus <- 0.5
weights <- c(
  "pLoF" = 41.81,
  "Missense" = 20.55,
  "Synonymous" = 3.09,
  "UTR5" = 7.68,
  "UTR3" = 3.02,
  "Promoter" = 2.70,
  "CRE" = 1.93,
  "Non-genic" = 0.44
)

df <- data.table::fread(cmd = "gunzip -dc /humgen/atgu1/methods/mkanai/workspace/202011_meta_finemapping/1000GP_Phase3/chr3/1000GP_Phase3_chr3.snv.consequence.tsv.bgz", data.table = F)

df_frq <- purrr::map_dfr(pops, function(pop) {
  data.table::fread(sprintf("/humgen/atgu1/methods/mkanai/workspace/202011_meta_finemapping/1000GP_Phase3/chr3/1000GP_Phase3_chr3.%s.afreq", pop), data.table = F) %>%
    dplyr::filter(ID %in% df$variant) %>%
    dplyr::mutate(variant = ID, pop = pop, freq = 1 - ALT_FREQS) %>%
    dplyr::select(variant, pop, freq)
}) %>%
  tidyr::pivot_wider(id_cols = "variant", names_prefix = "freq_", names_from = "pop", values_from = "freq")

df <-
  dplyr::left_join(df, df_frq) %>%
  dplyr::mutate(
    position = as.numeric(stringr::str_split_fixed(variant, ":", 4)[, 2]),
    locus = position %/% locus_length
  )


compute_true_beta <- function(df, rg, seed) {
  set.seed(seed)
  x <- dplyr::left_join(
    df,
    tibble::tibble(
      locus = unique(df$locus),
      causal_locus = ((locus %% 2) == rbinom(1, 1, 0.5))
    )
  ) %>%
    dplyr::filter(causal_locus) %>%
    dplyr::group_by(locus) %>%
    dplyr::mutate(gamma = causal_locus & variant == sample(variant, 1, prob = weights[consequence])) %>%
    dplyr::ungroup() %>%
    dplyr::filter(gamma) %>%
    dplyr::arrange(position)

  sigma2_g <- h2g / sum((2 * x$max_maf * (1 - x$max_maf))**(1 + alpha))
  x <-
    dplyr::mutate(x,
      ref = stringr::str_split_fixed(variant, ":", 4)[, 3],
      alt = stringr::str_split_fixed(variant, ":", 4)[, 4],
      var_beta = sigma2_g * (2 * max_maf * (1 - max_maf))**alpha,
      purrr::map_dfr(var_beta, function(var_beta) {
        sigma <- diag(length(pops) * n_cohorts_per_pop)
        diag(sigma) <- var_beta
        sigma[upper.tri(sigma)] <- var_beta * rg
        sigma[lower.tri(sigma)] <- var_beta * rg
        beta <- mvtnorm::rmvnorm(1, sigma = sigma) %>%
          as.data.frame()
        colnames(beta) <- paste0("beta_", rep(pops, each = n_cohorts_per_pop), "_sim", rep(seq(n_cohorts_per_pop), length(pops)))
        return(beta)
      })
    )
  return(x)
}

furrr::future_map(1:100, function(i) {
  x <- compute_true_beta(df, rg = 1, seed = i)
  write.table(x, sprintf("%s/true_pheno/pheno%d.beta.txt", basedir, i), quote = F, row.names = F, sep = "\t")
})

furrr::future_map(101:200, function(i) {
  x <- compute_true_beta(df, rg = 0.95, seed = i)
  write.table(x, sprintf("%s/true_pheno/pheno%d.beta.txt", basedir, i), quote = F, row.names = F, sep = "\t")
})

furrr::future_map(201:300, function(i) {
  x <- compute_true_beta(df, rg = 0.9, seed = i)
  write.table(x, sprintf("%s/true_pheno/pheno%d.beta.txt", basedir, i), quote = F, row.names = F, sep = "\t")
})

furrr::future_map(301:400, function(i) {
  x <- compute_true_beta(df, rg = 0.8, seed = i)
  write.table(x, sprintf("%s/true_pheno/pheno%d.beta.txt", basedir, i), quote = F, row.names = F, sep = "\t")
})

furrr::future_map(401:500, function(i) {
  x <- compute_true_beta(df, rg = 0.5, seed = i)
  write.table(x, sprintf("%s/true_pheno/pheno%d.beta.txt", basedir, i), quote = F, row.names = F, sep = "\t")
})

liftover <- data.table::fread("/humgen/atgu1/methods/mkanai/workspace/202011_meta_finemapping/imputed/TopMed_GRCh38.variants.liftover.txt", data.table = F)
purrr::map(1:500, function(i) {
  x <- data.table::fread(sprintf("%s/true_pheno/pheno%d.beta.txt", basedir, i), data.table = F) %>%
    dplyr::mutate(
      variant_b37 = variant,
      variant = liftover$original_variant[match(variant, liftover$variant)]
    )
  write.table(x, sprintf("%s/true_pheno/pheno%d.beta.liftover.txt", basedir, i), quote = F, row.names = F, sep = "\t")
})