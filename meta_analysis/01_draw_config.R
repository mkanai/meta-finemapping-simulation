library(dplyr)

df <- read.table(
  "~/src/github.com/mkanai/meta-finemapping-simulation/meta_analysis/config/assoc_list_final.txt",
  F,
  sep = "\t"
) %>%
  dplyr::filter(V2 != "1-5" & V4 %in% c("Omni25", "MEGA", "GSA")) %>%
  dplyr::mutate(
    V5 = ifelse(V5 == "TopMed", "TopMed-liftover", V5),
    fname = sprintf("%s_sim%s_chr%s_%s_%s", V1, V2, V3, V4, V5),
    batch = sprintf("%s_sim%s", V1, V2)
  )

draw_configs <- function(df, n_cohorts = 10, seed = 123456) {
  set.seed(seed)
  purrr::map_dfr(seq_len(100), function(i) {
    batches <- sample(unique(df$batch), n_cohorts)
    idx <-
      fnames <-
      dplyr::filter(df, batch %in% batches) %>%
      dplyr::group_by(batch) %>%
      dplyr::sample_n(1) %>%
      dplyr::ungroup() %>%
      dplyr::pull(fname)
    return(data.frame(t(sort(fnames))))
  }) %>%
    dplyr::distinct()
}

configs <- dplyr::filter(df, V1 == "EUR" & V5 == "1000GP3") %>%
  draw_configs()

write.table(
  configs,
  "~/src/github.com/mkanai/meta-finemapping-simulation/meta_analysis/config/EUR_random_configs10_1000GP3.tsv",
  quote = F,
  col.names = F,
  row.names = F,
  sep = "\t"
)

configs <- dplyr::filter(df, V1 == "EUR" & V4 == "Omni25") %>%
  draw_configs()

write.table(
  configs,
  "~/src/github.com/mkanai/meta-finemapping-simulation/meta_analysis/config/EUR_random_configs10_Omni25.tsv",
  quote = F,
  col.names = F,
  row.names = F,
  sep = "\t"
)

configs <- dplyr::filter(df, V1 == "EUR") %>%
  draw_configs()

write.table(
  configs,
  "~/src/github.com/mkanai/meta-finemapping-simulation/meta_analysis/config/EUR_random_configs10_mixed.tsv",
  quote = F,
  col.names = F,
  row.names = F,
  sep = "\t"
)

configs <- dplyr::filter(df, V5 == "1000GP3") %>%
  draw_configs() %>%
  dplyr::filter(!stringr::str_starts(X1, "^EUR"))

write.table(
  configs,
  "~/src/github.com/mkanai/meta-finemapping-simulation/meta_analysis/config/ALL_random_configs10_1000GP3.tsv",
  quote = F,
  col.names = F,
  row.names = F,
  sep = "\t"
)

configs <- dplyr::filter(df, V4 == "Omni25") %>%
  draw_configs() %>%
  dplyr::filter(!stringr::str_starts(X1, "^EUR"))

write.table(
  configs,
  "~/src/github.com/mkanai/meta-finemapping-simulation/meta_analysis/config/ALL_random_configs10_Omni25.tsv",
  quote = F,
  col.names = F,
  row.names = F,
  sep = "\t"
)

configs <- draw_configs(df) %>%
  dplyr::filter(!stringr::str_starts(X1, "^EUR"))

write.table(
  configs,
  "~/src/github.com/mkanai/meta-finemapping-simulation/meta_analysis/config/ALL_random_configs10_mixed.tsv",
  quote = F,
  col.names = F,
  row.names = F,
  sep = "\t"
)