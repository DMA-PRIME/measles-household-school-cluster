# ==============================================================================
# LOAD ALL MODULES
# ==============================================================================
# File: load_all.R
# Purpose: Source all module files in the correct order
# Usage: source("load_all.R") from the directory containing all R files
# ==============================================================================

cat("=== Loading Multi-School Measles Network Simulation Modules ===

")

# ------------------------------------------------------------------------------
# HPC/SLURM SAFE PACKAGE SETUP
# ------------------------------------------------------------------------------
# Uses a writable user library (R_LIBS_USER). Do NOT write to the system library.
# Default behavior: do NOT install packages inside sbatch jobs.
# To explicitly allow installs: export ALLOW_R_INSTALL=true
ALLOW_INSTALL <- tolower(Sys.getenv("ALLOW_R_INSTALL", "false")) %in% c("1","true","yes","y")

userlib <- Sys.getenv("R_LIBS_USER")
if (!nzchar(userlib)) {
  userlib <- path.expand("~/R/4.5")
  Sys.setenv(R_LIBS_USER = userlib)
}
dir.create(userlib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(userlib, .libPaths()))

cat("R_LIBS_USER: ", Sys.getenv("R_LIBS_USER"), "
", sep = "")
cat("libPaths:
")
cat(paste0("  - ", .libPaths(), collapse = "
"), "

")

required_packages <- c("Rcpp", "dplyr", "ggplot2", "tidyr", "R6", "igraph")

missing <- required_packages[!vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing) > 0) {
  cat("Missing required R packages:
")
  cat(paste0("  - ", missing, collapse = "
"), "

")

  if (!ALLOW_INSTALL) {
    stop(
      "Refusing to auto-install packages during sbatch run (prevents writing to system libraries).\n",
      "Install once into your user library, then rerun. Example:\n\n",
      "  module load r/4.5.0\n",
      "  export R_LIBS_USER=$HOME/R/4.5\n",
      "  mkdir -p \"$R_LIBS_USER\"\n",
      "  R -q -e 'install.packages(c(", paste(sprintf("\"%s\"", missing), collapse = ", "),
      "), repos=\"https://cloud.r-project.org\")'\n\n",
      "If you intentionally want in-job installs (not recommended), submit with:\n",
      "  ALLOW_R_INSTALL=true sbatch submit_slurm.sh\n"
    )
  } else {
    cat("ALLOW_R_INSTALL is enabled; installing missing packages into R_LIBS_USER...
")
    install.packages(missing, lib = Sys.getenv("R_LIBS_USER"), repos = "https://cloud.r-project.org")
  }
}

for (pkg in required_packages) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}
cat("Required packages loaded.

")

# ------------------------------------------------------------------------------
# Locate module directory (same logic as your original file)
# ------------------------------------------------------------------------------
module_dir <- NULL

if (exists("odir")) {
  module_dir <- odir
}

if (is.null(module_dir) && sys.nframe() > 0) {
  ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(ofile) && is.character(ofile)) {
    module_dir <- dirname(ofile)
  }
}

if (is.null(module_dir) || !file.exists(file.path(module_dir, "rcpp_transmission.R"))) {
  if (file.exists("rcpp_transmission.R")) {
    module_dir <- "."
  } else if (file.exists("codes/rcpp_transmission.R")) {
    module_dir <- "codes"
  } else if (file.exists("modular/rcpp_transmission.R")) {
    module_dir <- "modular"
  }
}

if (is.null(module_dir) || !file.exists(file.path(module_dir, "rcpp_transmission.R"))) {
  stop("Cannot find module files. Please set working directory to the folder containing the R files,
or source this file from that directory.")
}

cat(sprintf("Module directory: %s

", normalizePath(module_dir)))

# ------------------------------------------------------------------------------
# Source modules in dependency order
# ------------------------------------------------------------------------------
modules <- c(
  "rcpp_transmission.R",
  "network_utils.R",
  "grade_utils.R",
  "population_utils.R",
  "household_utils.R",
  "simulation_utils.R",
  "simulation_parallel.R",  # optional
  "analysis_utils.R",
  "plotting_utils.R"
)

for (module in modules) {
  module_path <- file.path(module_dir, module)
  if (file.exists(module_path)) {
    cat(sprintf("Loading: %s ... ", module))

    load_result <- tryCatch({
      source(module_path, local = globalenv())
      "OK"
    }, error = function(e) {
      paste("FAILED:", conditionMessage(e))
    })

    cat(load_result, "
")

    if (startsWith(load_result, "FAILED") && module != "simulation_parallel.R") {
      stop(paste("Failed to load critical module:", module, "
", load_result))
    }
  } else {
    if (module == "simulation_parallel.R") {
      cat(sprintf("Note: Optional module not found: %s (parallel simulations disabled)
", module))
    } else {
      cat(sprintf("WARNING: Module not found: %s
", module_path))
    }
  }
}

# ------------------------------------------------------------------------------
# Verify critical functions exist
# ------------------------------------------------------------------------------
cat("
=== Verifying critical functions ===
")

critical_functions <- c(
  "cpp_school_transmission_contacts",
  "cpp_between_school_transmission",
  "draw_erlang",
  "school_transmission",
  "run_network_simulation"
)

missing_critical <- critical_functions[!vapply(critical_functions, function(f) {
  exists(f, where = globalenv(), mode = "function")
}, FUN.VALUE = logical(1))]

if (length(missing_critical) > 0) {
  cat("WARNING: Critical functions NOT loaded:
")
  for (f in missing_critical) cat(sprintf("  - %s
", f))
  cat("
Check the loading messages above for FAILED modules.
")
} else {
  cat("All critical functions loaded successfully!
")
}

cat("
=== Module loading complete ===

")
