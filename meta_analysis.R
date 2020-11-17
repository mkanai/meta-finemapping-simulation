library(argparse)
library(dplyr)
library(furrr)

parser <- ArgumentParser()
parser$add_argument("--config", type = "integer", required = TRUE)
parser$add_argument("--pheno", type = "integer", required = TRUE)
parser$add_argument("--threads", type = "integer", default = 1)

args <- parser$parse_args()

plan(multisession, workers = args$threads)

print(Sys.time())

config = read.table("/humgen/atgu1/methods/mkanai/workspace/202011_meta_finemapping/meta_config.tsv", sep = "\t")
studies = config[args$config,]
studies = stringr::str_split_fixed(studies[studies != ""], "_", 5)
studies = stringr::str_c(studies[,1], studies[,2], studies[,3], sprintf("pheno%d", args$pheno), studies[,4], studies[,5], sep="_")

df = purrr::map_dfr(studies, function(x) {
  fname = sprintf("/humgen/atgu1/methods/mkanai/workspace/202011_meta_finemapping/assoc/pheno%d/%s.PHENO.glm.linear.gz", args$pheno, x)
  data.table::fread(fname, data.table = F) %>%
    dplyr::select("#ID", A1_CT, OBS_CT, BETA, SE, P) %>%
    dplyr::mutate(study = x) %>%
    dplyr::rename(variant = "#ID", a1count = A1_CT, n = OBS_CT, beta = BETA, se = SE, pvalue = P)
})

fixed_effect_meta = function(variant, beta, se, n, a1count) {
  n_studies = length(beta)
  if (n_studies == 1) {
    return(tibble::tibble(
      variant = variant[1],
      freq = a1count / n,
      beta = beta,
      se = se,
      pvalue = 2 * pnorm(-abs(beta / se)),
      pvalue_het = NA,
      n_studies = n_studies,
      n = n
    ))
  }
  
  inv_se2 = 1 / (se ** 2)
  unnorm_beta = beta * inv_se2
  beta_meta = sum(unnorm_beta) / sum(inv_se2)
  se_meta = sqrt(1 / sum(inv_se2))
  q_meta = sum((beta - beta_meta) ** 2 * inv_se2)
  p_meta = 2 * pnorm(-abs(beta_meta / se_meta))
  p_het = pchisq(q_meta, n_studies - 1, lower.tail = F)
  
  return(tibble::tibble(
    variant = variant[1],
    freq = sum(a1count) / sum(n),
    beta = beta_meta,
    se = se_meta,
    pvalue = p_meta,
    pvalue_het = p_het,
    n_studies = n_studies,
    n = sum(n)
  ))
}

out =
  dplyr::group_split(df, variant) %>%
  furrr::future_map_dfr(function(data) {
    fixed_effect_meta(data$variant, data$beta, data$se, data$n, data$a1count)
  })

x = stringr::str_split_fixed(out$variant, ":", 4)
out =
  dplyr::mutate(out, chromosome = x[,1], position = x[,2], ref = x[,3], alt = x[,4]) %>%
  dplyr::select(variant, chromosome, position, ref, alt, dplyr::everything())

fname = sprintf("/humgen/atgu1/methods/mkanai/workspace/202011_meta_finemapping/meta_analysis/meta.pheno%d.config%d.tsv", args$pheno, args$config)
data.table::fwrite(out, fname, quote = F, row.names = F, sep = "\t", na = "NA")
system(paste("bgzip -f", fname))

print(Sys.time())