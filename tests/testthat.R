# ==============================================================================
# TEST RUNNER
# ==============================================================================
# File: tests/testthat.R
# Purpose: Run all unit tests for the measles household-school cluster model
# Usage: Rscript tests/testthat.R  (from project root)
# ==============================================================================

library(testthat)

# Locate the project root (one level above tests/)
this_script <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
if (!is.null(this_script)) {
  project_root <- dirname(dirname(this_script))
} else {
  project_root <- normalizePath(".")
}

cat(sprintf("Project root: %s\n", project_root))

# Pass project root to test environment
Sys.setenv(PROJECT_ROOT = project_root)

# Run all tests
test_results <- test_dir(
  file.path(project_root, "tests", "testthat"),
  reporter = "progress",
  stop_on_failure = FALSE
)

# Print summary
cat("\n")
print(test_results)

# Exit with non-zero code if any tests failed
if (any(test_results$failed > 0) || any(test_results$error > 0)) {
  quit(status = 1)
}
