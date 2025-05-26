logfile <- snakemake@log[[1]]
log_con <- file(logfile, open = "wt")
sink(log_con, type = "output", append = TRUE)
sink(log_con, type = "message", append = TRUE)

options(
  error = function() {
    traceback(2)
    sink(type = "message"); sink()
    close(log_con)
    quit(status = 1, save = "no")
  }
)

df_list <- lapply(c("Variants", "Preimputation", "Imputation"), function(name) {
  lapply(snakemake@input[[name]], \(x) 
      read.table(x, header=FALSE, sep="\t", quote="", comment.char="",
        stringsAsFactors=FALSE, col.names=c("CHROM", "POS", "REF", "ALT"))) |>
    dplyr::bind_rows() |>
    tibble::as.tibble() |>
    dplyr::mutate(SOURCE = name)
})

df <- dplyr::bind_rows(df_list) |>
  dplyr::mutate(
    SOURCE = factor(SOURCE, levels = c("Variants", "Preimputation", "Imputation")))

df |>
  readr::write_tsv(snakemake@output[["ImputationPositions"]])

df |>
  dplyr::group_by(SOURCE) |>
  dplyr::count() |>
  readr::write_tsv(snakemake@output[["ImputationSummary"]])

df |>
  dplyr::group_by(SOURCE, CHROM) |>
  dplyr::count() |>
  tidyr::pivot_wider(
    names_from = SOURCE,
    values_from = n,
    values_fill = 0
  ) |>
  readr::write_tsv(snakemake@output[["ImputationChromSummary"]])

save.image(file=snakemake@output[["RData"]])
sink(type = "message")
sink()
close(log_con)