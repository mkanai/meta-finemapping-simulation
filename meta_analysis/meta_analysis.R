library(argparse)
library(dplyr)

parser <- ArgumentParser()
parser$add_argument("--config", type = "integer", required = TRUE)
parser$add_argument("--pheno", type = "integer", required = TRUE)
parser$add_argument(
  "--config-file",
  type = "character"
)
parser$add_argument(
  "--sumstats-format",
  type = "character"
)
parser$add_argument(
  "--out-format",
  type = "character"
)
parser$add_argument(
  "--plink-path",
  type = "character",
  default = "plink"
)

args <- parser$parse_args()
print(Sys.time())

config <- read.table(args$config_file, sep = "\t")
studies <- config[args$config, ]
studies <- stringr::str_split_fixed(studies[studies != ""], "_", 5)
studies <- stringr::str_c(studies[, 1], studies[, 2], studies[, 3], sprintf("pheno%d", args$pheno), studies[, 4], studies[, 5], sep = "_")
out_fname <- sprintf(args$out_format, args$pheno, args$config)

if (stringr::str_starts(args$sumstats_format, "gs://")) {
  # download files
  cmd <- c("gsutil", "-m", "cp", sprintf(args$sumstats_format, args$pheno, studies), ".")
  print(cmd)
  system(paste(cmd, collapse = " "))
  args$sumstats_format <- basename(args$sumstats_format)
}

# run plink --meta-analysis
cmd <- c(
  args$plink_path,
  "--meta-analysis",
  sprintf(args$sumstats_format, args$pheno, studies),
  "+",
  "qt",
  "report-all",
  "study",
  "--out",
  paste0(out_fname, ".tmp")
)
print(cmd)
system(paste(cmd, collapse = " "))

df <- purrr::map_dfr(studies, function(x) {
  fname <- sprintf(args$sumstats_format, args$pheno, x)
  pop <- stringr::str_split_fixed(x, "_", 5)[1, 1]
  data.table::fread(fname, data.table = F, select = c("SNP", "A1_CT", "OBS_CT")) %>%
    dplyr::mutate(pop = pop) %>%
    dplyr::rename(variant = SNP)
}) %>%
  dplyr::group_by(variant) %>%
  dplyr::summarize(
    freq = sum(A1_CT) / (2 * sum(OBS_CT)),
    n_studies = length(A1_CT),
    n_samples = sum(OBS_CT),
    n_samples_AFR = sum(OBS_CT * (pop == "AFR")),
    n_samples_EAS = sum(OBS_CT * (pop == "EAS")),
    n_samples_EUR = sum(OBS_CT * (pop == "EUR"))
  )

for (col in paste0("n_samples_", c("AFR", "EAS", "EUR"))) {
  if (all(df[[col]] == 0)) {
    df[[col]] <- NULL
  }
}

meta <- read.table(paste0(out_fname, ".tmp.meta"), header = TRUE)
out <-
  dplyr::transmute(
    meta,
    variant = SNP,
    beta = BETA,
    se = SE,
    p = P,
    p_het = Q
  ) %>%
  dplyr::left_join(df)

x <- stringr::str_split_fixed(out$variant, ":", 4)
out <-
  dplyr::mutate(out, chromosome = x[, 1], position = x[, 2], ref = x[, 3], alt = x[, 4]) %>%
  dplyr::select(variant, chromosome, position, ref, alt, dplyr::everything()) %>%
  dplyr::mutate(
    se = format(se, digits = 4),
    freq = format(freq, digits = 4)
  )

cols <- colnames(meta)[stringr::str_detect(colnames(meta), "[BS][0-9]+")]
cols_beta <- cols[seq(1, length(cols), 2)]
cols_se <- cols[seq(2, length(cols), 2)]
out_study_beta <- meta[c("SNP", cols_beta, cols_se)]
colnames(out_study_beta) <- c("variant", paste0("beta", seq(length(cols_beta))), paste0("se", seq(length(cols_beta))))
out_fname_study_beta <- stringr::str_replace(out_fname, ".tsv$", ".study_beta.tsv")

data.table::fwrite(out, out_fname, quote = F, row.names = F, sep = "\t", na = "NA")
data.table::fwrite(out_study_beta, out_fname_study_beta, quote = F, row.names = F, sep = "\t", na = "NA")
system(paste("bgzip -f", out_fname))
system(paste("bgzip -f", out_fname_study_beta))
system(paste0("rm ", out_fname, ".tmp.{meta,log}"))

print(Sys.time())