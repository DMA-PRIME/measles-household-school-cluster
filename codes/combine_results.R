# ==============================================================================
# COMBINE RESULTS FROM SLURM ARRAY JOB
# ==============================================================================
# Usage: Rscript combine_results.R <job_id>
# Example: Rscript combine_results.R 12345678
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  cat("Usage: Rscript combine_results.R <job_id>\n")
  cat("Example: Rscript combine_results.R 12345678\n")
  quit(status = 1)
}

results_dir <- file.path("results", args[1])

if (!dir.exists(results_dir)) {
  stop("Results directory not found: ", results_dir)
}

# Find and load all task files
rds_files <- list.files(results_dir, pattern = "task_.*\\.rds$", full.names = TRUE)
cat(sprintf("Found %d task files\n", length(rds_files)))

results_list <- lapply(rds_files, readRDS)
combined <- do.call(rbind, results_list)
combined$global_id <- 1:nrow(combined)

# Summary
cat(sprintf("\n=== SUMMARY (%d simulations) ===\n", nrow(combined)))
cat(sprintf("Mean infections:    %.1f\n", mean(combined$total_infected, na.rm=TRUE)))
cat(sprintf("Median infections:  %.0f\n", median(combined$total_infected, na.rm=TRUE)))
cat(sprintf("Max infections:     %d\n", max(combined$total_infected, na.rm=TRUE)))
cat(sprintf("Schools affected:   %.1f (mean)\n", mean(combined$schools_affected, na.rm=TRUE)))

# Save
output_file <- file.path(results_dir, "combined_results.rds")
saveRDS(combined, output_file)
write.csv(combined, sub("\\.rds$", ".csv", output_file), row.names = FALSE)

cat(sprintf("\nSaved: %s\n", output_file))
