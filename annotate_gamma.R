library(argparse)
library(dplyr)

parser <- ArgumentParser()
parser$add_argument("--snp", type = "character", required = TRUE)
parser$add_argument("--meta", type = "character", required = TRUE)
parser$add_argument("--beta", type = "character", required = TRUE)
parser$add_argument("--out", type = "character", required = TRUE)

df.abf <- data.table::fread(cmd = paste("zcat", args$snp), data.table = F) %>%
  dplyr::mutate(variant = stringr::str_c("chr", variant)) %>%
  dplyr::select(!c(rsid, chromosome, position, allele1, allele2))

df.beta <-
  read.table(args$beta, T) %>%
  dplyr::mutate(true_beta = beta, gamma = ifelse(gamma == "true", TRUE, FALSE)) %>%
  dplyr::select(variant, true_beta, most_severe, gene_most_severe, lof, dplyr::starts_with("freq_"), gamma)

df.meta <- data.table::fread(args$meta, data.table = F) %>%
  dplyr::filter(variant %in% df.beta$variant & !(variant %in% df.abf$variant)) %>%
  dplyr::mutate(
    trait = df.abf$trait[1],
    maf = 1 - abs(1 - freq),
    region = NA,
    cs = NA,
    cs_99 = NA,
    lbf = NA,
    prob = NA
  ) %>%
  dplyr::select(!c(chromosome, position, ref, alt, freq))

df <-
  dplyr::bind_rows(df.abf, df.meta) %>%
  dplyr::left_join(df.beta, by = "variant") %>%
  dplyr::arrange(desc(prob)) %>%
  dplyr::mutate(rank = row_number() / n()) %>%
  dplyr::filter(gamma | cs > 0 | cs_99 > 0 | rank < 0.01 | p_het < 0.001)

write.table(df, args$out, quote = F, row.names = F, sep = "\t")
