# ==============================================================================
# BATCH SIMULATION SCRIPT FOR SLURM
# ==============================================================================
# Runs simulations in parallel using multiple cores.
#
# Notes:
# - This script does NOT install packages. Package loading/install policy is handled
#   in load_all.R using R_LIBS_USER.
# - Ensure your sbatch script exports R_LIBS_USER to a writable location.
#
# Usage:
#   Rscript run_batch_clean.R --task_id=1 --n_sims=100 --n_cores=12 --n_days=150 --output_dir=./results
# ==============================================================================

suppressWarnings(suppressMessages(library(parallel)))  # part of base R

# ------------------------------------------------------------------------------
# Ensure R_LIBS_USER (if provided) is on .libPaths() early
# ------------------------------------------------------------------------------
userlib <- Sys.getenv("R_LIBS_USER")
if (nzchar(userlib)) {
  dir.create(userlib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(userlib, .libPaths()))
}

# ==============================================================================
# PARSE ARGUMENTS (base R only)
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)

opt <- list(
  task_id = 1L,
  n_sims = 100L,
  n_cores = 12L,
  n_days = 365L,
  output_dir = "./results"
)

for (arg in args) {
  if (grepl("^--", arg)) {
    parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    if (length(parts) < 2) next
    key <- parts[1]
    value <- parts[2]
    if (key %in% names(opt)) {
      if (key %in% c("task_id", "n_sims", "n_cores", "n_days")) {
        opt[[key]] <- as.integer(value)
      } else {
        opt[[key]] <- value
      }
    }
  }
}

cat(sprintf("Task %d: Running %d simulations on %d cores
",
            opt$task_id, opt$n_sims, opt$n_cores))
cat(sprintf("Working directory: %s
", getwd()))
cat(sprintf("R_LIBS_USER: %s
", if (nzchar(userlib)) userlib else "<unset>"))
cat("libPaths:
")
cat(paste0("  - ", .libPaths(), collapse = "
"), "

")

# ==============================================================================
# LOAD MODULES / FUNCTIONS
# ==============================================================================
cat("Loading simulation modules...
")

if (!file.exists("load_all.R")) {
  stop("ERROR: load_all.R not found in ", getwd())
}
source("load_all.R", local = globalenv())

# ==============================================================================
# LOAD DATA
# ==============================================================================
if (!dir.exists("data")) {
  stop("ERROR: data/ directory not found in ", getwd())
}

possible_files <- c(
  "data/SC_vaccination_merged.csv",
  "data/school_data.csv",
  "data/schools.csv"
)

school_file <- possible_files[file.exists(possible_files)][1]
if (is.na(school_file)) {
  cat("ERROR: Could not find school data file. Files in data/:
")
  print(list.files("data"))
  stop("Please check your data file name")
}

cat(sprintf("Loading school data from: %s
", school_file))
schools_raw <- read.csv(school_file, stringsAsFactors = FALSE)

schools <- data.frame(
  school_id = seq_len(nrow(schools_raw)),
  school_name = as.character(schools_raw$School.Name),
  school_size = as.integer(schools_raw$Total.Students),
  vaccination_coverage = as.numeric(schools_raw$Percent.Immunized),
  lon = as.numeric(schools_raw$longitude),
  lat = as.numeric(schools_raw$latitude),
  stringsAsFactors = FALSE
)

schools <- schools[!is.na(schools$school_size) & schools$school_size > 0, , drop = FALSE]
schools$school_id <- seq_len(nrow(schools))

# Filter to region matching household cache (built for Spartanburg only, 93 schools).
# Without this filter, school IDs won't match the cache, corrupting all simulations.
region_counties <- c("Spartanburg")

county_col <- schools_raw$County
valid_rows <- !is.na(as.integer(schools_raw$Total.Students)) &
              as.integer(schools_raw$Total.Students) > 0
schools$county <- as.character(county_col[valid_rows])

schools <- schools[schools$county %in% region_counties, , drop = FALSE]
schools$school_id <- seq_len(nrow(schools))

cat(sprintf("Filtered to %d schools in: %s\n",
            nrow(schools), paste(region_counties, collapse = ", ")))

cat(sprintf("Loaded %d schools
", nrow(schools)))

# ==============================================================================
# PARAMETERS
# ==============================================================================
params <- list(

  # Disease progression (Erlang-distributed durations)

  latent_mean = 10,

  latent_shape = 8,

  infectious_mean = 8,

  infectious_shape = 8,

  prodromal_period = 4,

  

  # Within-school transmission

  c_within = 5,           # Mean contacts per day within class

  c_between = 2,          # Mean contacts per day between classes

  p_within = 0.10,       # Per-contact transmission probability (within class)

  p_between = 0.10,       # Per-contact transmission probability (between class)

  

  # Between-school transmission

  c_between_school = 0.15,  # Base contact rate between schools (scaled by network)

  

  # Household transmission

  hh_transmission_prob = 0.15,  # Daily per-contact probability within household

  

  # Infectiousness modifiers

  prodromal_infectiousness_multiplier = 1.0,

  rash_infectiousness_multiplier = 1.0,

  

  # Vaccine parameters

  vaccine_efficacy = 0.97,

  vaccine_infectiousness_reduction = 0.8,

  

  # Intervention parameters

  no_intervention = FALSE,

  isolation_delay_index = 2,

  isolation_delay_secondary = 2,

  isolation_period = 14,

  quarantine_contacts = TRUE,

  quarantine_efficacy = 0.90,

  quarantine_duration = 21,

  

  # Population structure

  avg_class_size = 25,

  age_range = c(5, 18)

)
# ==============================================================================
# SETUP NETWORK (Travel Time-Based)
# ==============================================================================
# Uses precomputed driving time matrix for between-school connectivity.
# Falls back to Haversine distance if travel time data is unavailable.
# ==============================================================================
cat("Generating network...\n")

travel_time_file     <- "data/Travel_time_matrix_Upstate.csv"
school_ref_file      <- "data/School_reference_Upstate.csv"
max_travel_time      <- 16     # minutes
network_weight_method <- "exponential"

if (file.exists(travel_time_file) && file.exists(school_ref_file)) {

  cat("Using precomputed travel time matrix for network...\n")

  # Load school reference to get OBJECTID -> School.Name mapping
  school_ref <- read.csv(school_ref_file, stringsAsFactors = FALSE)
  school_ref <- school_ref[!is.na(school_ref$School.Name) & school_ref$School.Name != "", ]

  # Handle OBJECTID vs ID column name
  if ("OBJECTID" %in% names(school_ref) && !"ID" %in% names(school_ref)) {
    school_ref$ID <- school_ref$OBJECTID
  }

  # Match schools to reference IDs by County + School.Name
  schools$match_key <- paste(schools$county, schools$school_name, sep = "_")
  school_ref$match_key <- paste(school_ref$County, school_ref$School.Name, sep = "_")

  ref_lookup <- school_ref[!duplicated(school_ref$match_key), c("match_key", "ID")]
  schools$ref_id <- ref_lookup$ID[match(schools$match_key, ref_lookup$match_key)]

  n_matched <- sum(!is.na(schools$ref_id))
  cat(sprintf("  Matched %d / %d schools to travel time matrix IDs\n",
              n_matched, nrow(schools)))

  if (n_matched >= 2) {
    matched_ref_ids <- schools$ref_id[!is.na(schools$ref_id)]

    tt_network <- generate_travel_time_network(
      travel_time_file = travel_time_file,
      school_reference_file = school_ref_file,
      selected_school_ids = matched_ref_ids,
      max_travel_time = max_travel_time,
      weight_method = network_weight_method
    )

    # Remap adjacency matrix to sequential school_id order
    # tt_network$adjacency is indexed by available_ids order
    # Build a full nrow(schools) x nrow(schools) matrix
    n_schools <- nrow(schools)
    adj_full <- matrix(0, nrow = n_schools, ncol = n_schools)

    # Create mapping: ref_id -> index in tt_network adjacency
    tt_ids <- as.numeric(rownames(tt_network$adjacency))
    tt_idx_lookup <- setNames(seq_along(tt_ids), tt_ids)

    for (i in 1:n_schools) {
      for (j in 1:n_schools) {
        if (i == j) next
        ri <- schools$ref_id[i]
        rj <- schools$ref_id[j]
        if (is.na(ri) || is.na(rj)) next
        ti <- tt_idx_lookup[as.character(ri)]
        tj <- tt_idx_lookup[as.character(rj)]
        if (!is.na(ti) && !is.na(tj)) {
          adj_full[i, j] <- tt_network$adjacency[ti, tj]
        }
      }
    }

    # Build igraph from remapped adjacency
    library(igraph)
    g <- graph_from_adjacency_matrix(adj_full, mode = "undirected",
                                      weighted = TRUE, diag = FALSE)
    if ("school_name" %in% names(schools)) V(g)$name <- schools$school_name
    V(g)$school_id <- schools$school_id

    network <- list(
      graph = g,
      adjacency = adj_full,
      n_schools = n_schools,
      n_edges = ecount(g),
      max_travel_time = max_travel_time,
      weight_method = network_weight_method
    )

    cat(sprintf("Travel time network: %d edges (max: %d min)\n",
                network$n_edges, max_travel_time))

    # Report unmatched schools
    unmatched <- which(is.na(schools$ref_id))
    if (length(unmatched) > 0) {
      cat(sprintf("  WARNING: %d schools not in travel time matrix (no between-school edges):\n",
                  length(unmatched)))
      for (idx in unmatched[1:min(5, length(unmatched))]) {
        cat(sprintf("    - %s\n", schools$school_name[idx]))
      }
      if (length(unmatched) > 5) cat(sprintf("    ... and %d more\n", length(unmatched) - 5))
    }
  } else {
    cat("  WARNING: Too few matches — falling back to Haversine distance network\n")
    network <- generate_distance_network(schools, max_distance_km = round(max_travel_time * 0.87))
  }

  # Clean up temporary columns
  schools$match_key <- NULL
  schools$ref_id <- NULL

} else {
  cat("Travel time files not found — using Haversine distance network\n")
  if (!file.exists(travel_time_file)) cat(sprintf("  Missing: %s\n", travel_time_file))
  if (!file.exists(school_ref_file))  cat(sprintf("  Missing: %s\n", school_ref_file))
  network <- generate_distance_network(schools, max_distance_km = round(max_travel_time * 0.87))
}

cat("Creating populations...
")
populations <- lapply(seq_len(nrow(schools)), function(i) {
  create_school_population(i, schools$school_size[i], params$avg_class_size, params$age_range)
})

# ==============================================================================
# HOUSEHOLD SETUP (use cache if available)
# ==============================================================================
hh_cache_file <- "data/household_assignment_cache.rds"
hh_pop <- NULL
household_assignment <- NULL

if (file.exists(hh_cache_file)) {
  cat("Loading cached household assignment...
")
  hh_result <- load_household_assignment(hh_cache_file, populations, schools, verbose = FALSE)

  populations <- hh_result$populations
  household_assignment <- hh_result$assignment_df

  populations <- assign_household_level_vaccination(populations, schools, household_assignment, verbose = FALSE)

  hh_pop <- create_household_population(hh_result$household_members, populations, params)
  hh_pop <- initialize_household_vaccination(populations, hh_pop, 0.8, 0.5, 0.95)
  cat(sprintf("Household transmission enabled: %d members
", nrow(hh_pop)))
} else {
  cat("No household cache found - running without household transmission
")
  for (i in seq_len(nrow(schools))) {
    populations[[i]] <- initialize_vaccination_school(populations[[i]], schools$vaccination_coverage[i])
  }
}

seed_school <- which.min(schools$vaccination_coverage)
cat(sprintf("Seed school: %d (%s) - %.1f%% coverage
",
            seed_school, schools$school_name[seed_school],
            schools$vaccination_coverage[seed_school] * 100))

# ==============================================================================
# RUN SIMULATIONS IN PARALLEL
# ==============================================================================
base_seed <- 12345L + (opt$task_id - 1L) * opt$n_sims
sim_seeds <- base_seed + seq_len(opt$n_sims) - 1L

# Flag for workers
use_households <- !is.null(hh_pop)

cat(sprintf("\nRunning %d simulations...\n", opt$n_sims))
start_time <- Sys.time()

results_list <- mclapply(seq_len(opt$n_sims), function(i) {
  tryCatch({

    # Create a fresh hh_pop copy with reset disease states for each simulation
    # (matches the logic in run_multiple_network_simulations)
    sim_hh_pop <- NULL
    if (use_households) {
      sim_hh_pop <- hh_pop
      sim_hh_pop$state <- ifelse(sim_hh_pop$is_vaccinated, "V", "S")
      sim_hh_pop$time_in_state <- 0L
      sim_hh_pop$time_since_prodromal <- NA_integer_
      sim_hh_pop$is_isolated <- FALSE
      sim_hh_pop$is_quarantined <- FALSE
      sim_hh_pop$breakthrough_infection <- FALSE
    }

    result <- run_network_simulation(
      schools = schools,
      network = network,
      params = params,
      seed_schools = seed_school,
      n_initial_infected = 5,
      n_days = opt$n_days,
      seed = sim_seeds[i],
      hh_pop = sim_hh_pop,
      household_assignment = household_assignment
    )

    # ------------------------------------------------------------------
    # FIXED: Use correct field names from run_network_simulation return.
    # The function returns result$school_summary (NOT outbreak_summary).
    # Using the wrong name silently returns NULL, and sum(NULL) = 0,
    # which is why all results previously showed zero infections.
    # ------------------------------------------------------------------

    total_inf <- result$total_infected
    schools_aff <- sum(result$school_summary$total_infected > 0)
    total_breakthrough <- result$total_breakthrough
    hh_inf <- result$total_hh_members_infected
    actual_days <- result$actual_days

    # Peak day from network-wide daily counts
    peak_day <- NA_integer_
    ndc <- result$network_daily_counts
    if (!is.null(ndc)) {
      active_cols <- intersect(names(ndc), c("P", "Ra"))
      if (length(active_cols) > 0) {
        active_sum <- rowSums(ndc[, active_cols, drop = FALSE])
        peak_day <- ndc$day[which.max(active_sum)]
      }
    }

    # Store full epidemic curve, tagged with task/sim IDs
    epidemic_curve <- ndc
    epidemic_curve$task_id <- opt$task_id
    epidemic_curve$sim_id <- i

    # Household epidemic curve (if available)
    hh_curve <- NULL
    if (!is.null(result$hh_daily_counts)) {
      hh_curve <- result$hh_daily_counts
      hh_curve$task_id <- opt$task_id
      hh_curve$sim_id <- i
    }

    # Per-school summary for this simulation
    school_summary <- result$school_summary
    school_summary$task_id <- opt$task_id
    school_summary$sim_id <- i

    list(
      summary = data.frame(
        task_id = opt$task_id,
        sim_id = i,
        seed = sim_seeds[i],
        total_infected = total_inf,
        schools_affected = schools_aff,
        total_breakthrough = total_breakthrough,
        hh_members_infected = hh_inf,
        actual_days = actual_days,
        peak_day = peak_day,
        stringsAsFactors = FALSE
      ),
      epidemic_curve = epidemic_curve,
      hh_curve = hh_curve,
      school_summary = school_summary
    )

  }, error = function(e) {
    list(
      summary = data.frame(
        task_id = opt$task_id, sim_id = i, seed = sim_seeds[i],
        total_infected = NA, schools_affected = NA,
        total_breakthrough = NA, hh_members_infected = NA,
        actual_days = NA, peak_day = NA,
        error = conditionMessage(e),
        stringsAsFactors = FALSE
      ),
      epidemic_curve = NULL,
      hh_curve = NULL,
      school_summary = NULL
    )
  })
}, mc.cores = opt$n_cores)

elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

# ==============================================================================
# ASSEMBLE AND SAVE RESULTS
# ==============================================================================
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Summary statistics (one row per simulation)
summary_df <- do.call(rbind, lapply(results_list, function(x) x$summary))

summary_file <- file.path(opt$output_dir, sprintf("task_%02d_summary.rds", opt$task_id))
saveRDS(summary_df, summary_file)
cat(sprintf("Summary saved: %s\n", summary_file))

# 2. Full epidemic curves (one row per day per simulation)
curves <- Filter(Negate(is.null), lapply(results_list, function(x) x$epidemic_curve))
if (length(curves) > 0) {
  curves_df <- do.call(rbind, curves)
  curves_file <- file.path(opt$output_dir, sprintf("task_%02d_curves.rds", opt$task_id))
  saveRDS(curves_df, curves_file)
  cat(sprintf("Epidemic curves saved: %s\n", curves_file))
}

# 3. Household epidemic curves (if applicable)
hh_curves <- Filter(Negate(is.null), lapply(results_list, function(x) x$hh_curve))
if (length(hh_curves) > 0) {
  hh_curves_df <- do.call(rbind, hh_curves)
  hh_file <- file.path(opt$output_dir, sprintf("task_%02d_hh_curves.rds", opt$task_id))
  saveRDS(hh_curves_df, hh_file)
  cat(sprintf("Household curves saved: %s\n", hh_file))
}

# 4. Per-school summaries (one row per school per simulation)
school_sums <- Filter(Negate(is.null), lapply(results_list, function(x) x$school_summary))
if (length(school_sums) > 0) {
  school_sums_df <- do.call(rbind, school_sums)
  school_file_out <- file.path(opt$output_dir, sprintf("task_%02d_schools.rds", opt$task_id))
  saveRDS(school_sums_df, school_file_out)
  cat(sprintf("School summaries saved: %s\n", school_file_out))
}

# ==============================================================================
# PRINT SUMMARY
# ==============================================================================
cat(sprintf("\nCompleted in %.1f sec (%.2f sec/sim)\n", elapsed, elapsed / opt$n_sims))
cat(sprintf("Mean student infections: %.1f, Median: %.0f, Max: %.0f\n",
            mean(summary_df$total_infected, na.rm = TRUE),
            median(summary_df$total_infected, na.rm = TRUE),
            suppressWarnings(max(summary_df$total_infected, na.rm = TRUE))))
cat(sprintf("Mean HH member infections: %.1f\n",
            mean(summary_df$hh_members_infected, na.rm = TRUE)))
cat(sprintf("Mean schools affected: %.1f out of %d\n",
            mean(summary_df$schools_affected, na.rm = TRUE), nrow(schools)))
cat(sprintf("Mean outbreak duration: %.0f days\n",
            mean(summary_df$actual_days, na.rm = TRUE)))
n_failed <- sum(is.na(summary_df$total_infected))
if (n_failed > 0) {
  cat(sprintf("WARNING: %d/%d simulations failed\n", n_failed, opt$n_sims))
}
