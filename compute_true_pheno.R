library(argparse)
library(dplyr)

parser <- ArgumentParser()
parser$add_argument("--pop", type = "character", required = TRUE)
parser$add_argument("--sim", type = "integer", required = TRUE)
parser$add_argument("--chr", type = "integer", required = TRUE)
parser$add_argument("--pheno", type = "integer", required = TRUE)

args <- parser$parse_args()

df = read.table(sprintf("%s_sim%s_chr%s_pheno%s.sscore", args$pop, args$sim, args$chr, args$pheno), T, comment.char='')

h2g_sim = var(df$SCORE1_SUM)
print(h2g_sim)
out = 
  dplyr::mutate(df,
    IID = sprintf("%s_sim%s_%d", args$pop, args$sim, as.numeric(stringr::str_split_fixed(IID, "_", 2)[,2]) + 1),
    PHENO = SCORE1_SUM + rnorm(nrow(df), sd = sqrt(1 - h2g_sim))
  ) %>%
  dplyr::select(IID, PHENO)

write.table(out, sprintf("%s_sim%s_chr%s_pheno%s.pheno", args$pop, args$sim, args$chr, args$pheno), quote = F, row.names = F, sep = "\t")

out_h2g = tibble::tibble(
  pop=args$pop,
  sim=args$sim,
  chr=args$chr,
  pheno=args$pheno,
  h2g=h2g_sim
)

write.table(out_h2g, sprintf("%s_sim%s_chr%s_pheno%s.h2g.txt", args$pop, args$sim, args$chr, args$pheno), quote = F, row.names = F, sep = "\t")
