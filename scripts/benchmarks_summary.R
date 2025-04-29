#!/usr/bin/env Rscript

# Parse command line arguments; default folder is "benchmarks"
args <- commandArgs(trailingOnly = TRUE)
raw_output <- FALSE

# Check if the --raw flag is provided and remove it from the arguments
if ("--raw" %in% args) {
  raw_output <- TRUE
  args <- setdiff(args, "--raw")
}

folder <- if (length(args) >= 1) args[1] else "benchmarks"

# List all files in the specified folder (full path)
files <- list.files(folder, full.names = TRUE)

# Read each file (assuming tab-separated values with header)
data_list <- lapply(files, function(file) {
  read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
})

# Merge all data frames
merged_data <- do.call(rbind, data_list)

# Convert "h:m:s" column to seconds if it exists, store result in new column "elapsed_time_sec"
if ("h.m.s" %in% names(merged_data)) {
  merged_data$elapsed_time_sec <- sapply(merged_data$`h.m.s`, function(x) {
    parts <- unlist(strsplit(x, ":"))
    if (length(parts) == 3) {
      hours <- as.numeric(parts[1])
      minutes <- as.numeric(parts[2])
      seconds <- as.numeric(parts[3])
      return(hours * 3600 + minutes * 60 + seconds)
    } else {
      return(NA)
    }
  })
}

# Helper function to format seconds into h:m:s
format_time <- function(sec) {
  hours <- sec %/% 3600
  remainder <- sec %% 3600
  minutes <- remainder %/% 60
  seconds <- remainder %% 60
  sprintf("%d:%02d:%02d", hours, minutes, round(seconds))
}

# Define mapping for units based on known columns
unit_map <- list(
  s = " seconds",
  elapsed_time_sec = "",
  max_rss = " MB",
  max_vms = " MB",
  max_uss = " MB",
  max_pss = " MB",
  io_in = " bytes",
  io_out = " bytes",
  mean_load = "",
  cpu_time = " secs"
)

# Define mapping for column descriptions
desc_map <- list(
  s = "Wall clock time",
  h.m.s = "Wall clock time",
  elapsed_time_sec = "Wall clock time",
  max_rss = "Max RSS memory usage",
  max_vms = "Max VMS memory usage",
  max_uss = "Max USS memory usage",
  max_pss = "Max PSS memory usage",
  io_in = "I/O read",
  io_out = "I/O written",
  mean_load = "CPU load",
  cpu_time = "CPU time"
)

if (raw_output) {
  # Output merged raw data
  write.table(merged_data, file = stdout(), sep = "\t", row.names = FALSE, quote = FALSE)
} else {
  # Initialize a data frame for summary statistics
  stats <- data.frame(Column = character(),
                      Description = character(),
                      Average = numeric(),
                      Formatted = character(),
                      stringsAsFactors = FALSE)
  merged_data$s <- NULL
  for (col in names(merged_data)) {
    # Only process numeric columns
    if (is.numeric(merged_data[[col]])) {
      avg_val <- mean(merged_data[[col]], na.rm = TRUE)
      
      # Determine unit from mapping (default to empty string if not found)
      unit <- if (!is.null(unit_map[[col]])) unit_map[[col]] else ""
      
      # Format the value: if column is elapsed_time_sec, provide a h:m:s formatted string
      if (col == "elapsed_time_sec") {
        formatted_value <- paste0(format_time(avg_val), unit)
      } else {
        formatted_value <- paste0(sprintf("%.2f", avg_val), unit)
      }
      
      # Retrieve description from mapping if available
      description <- if (!is.null(desc_map[[col]])) desc_map[[col]] else ""
      
      stats <- rbind(stats, data.frame(Description = description,
                                       Average = formatted_value,
                                       stringsAsFactors = FALSE))
    }
  }
  
  # Print the summary table
  print(stats, row.names = FALSE)
}
