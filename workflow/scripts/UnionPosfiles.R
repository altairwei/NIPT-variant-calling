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

df_list <- lapply(c("BaseVar", "RefPanel"), function(name) {
  lapply(snakemake@input[[name]], \(x) 
      read.table(x, header=FALSE, sep="\t", quote="", comment.char="",
        stringsAsFactors=FALSE, col.names=c("CHROM", "POS", "REF", "ALT"))) |>
    dplyr::bind_rows() |>
    tibble::as.tibble() |>
    dplyr::mutate(SOURCE = name)
})

df <- dplyr::bind_rows(df_list) |>
  dplyr::mutate(
    SOURCE = factor(SOURCE, levels = c("BaseVar", "RefPanel")))

df_filtered <- df |>
  # 先按 SOURCE 排序，这样 BaseVar 的记录会排在同一 (CHROM, POS) 组合的最前面
  dplyr::arrange(CHROM, POS, SOURCE) |>
  # 对 CHROM+POS 去重，保留每组的第一行（即 BaseVar 如存在则保留它）
  dplyr::distinct(CHROM, POS, .keep_all = TRUE) |>
  dplyr::ungroup() |>
  dplyr::select(!SOURCE)

readr::write_tsv(
  df_filtered,
  snakemake@output[["posfile"]],
  col_names = FALSE
)

sink(type = "message")
sink()
close(log_con)
