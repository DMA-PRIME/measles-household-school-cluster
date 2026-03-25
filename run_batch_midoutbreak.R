# ==============================================================================
# BATCH MID-OUTBREAK SIMULATION SCRIPT FOR SLURM
# ==============================================================================
# Starts simulations from an observed mid-outbreak state:
#   - Schools that have reported exposures
#   - Total confirmed cases to date
#   - Contacts currently in quarantine
# Then simulates forward to forecast remaining outbreak trajectory.
#
# Output files (same structure as run_batch_clean.R):
#   task_XX_summary.rds   - Per-simulation summary
#   task_XX_curves.rds    - Network-level epidemic curves
#   task_XX_hh_curves.rds - Household epidemic curves
#   task_XX_schools.rds   - Per-school summaries
#
# Usage:
#   Rscript run_batch_midoutbreak.R --task_id=1 --n_sims=100 --n_cores=12 \
#     --n_days=120 --output_dir=./results_forecast
# ==============================================================================

suppressWarnings(suppressMessages(library(parallel)))

userlib <- Sys.getenv("R_LIBS_USER")
if (nzchar(userlib)) {
  dir.create(userlib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(userlib, .libPaths()))
}

# ==============================================================================
# PARSE ARGUMENTS
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)

opt <- list(
  task_id = 1L,
  n_sims = 100L,
  n_cores = 12L,
  n_days = 365L,
  output_dir = "./results_forecast"
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

cat(sprintf("Task %d: Running %d mid-outbreak simulations on %d cores\n",
            opt$task_id, opt$n_sims, opt$n_cores))

# ==============================================================================
# LOAD MODULES
# ==============================================================================
cat("Loading simulation modules...\n")
source("load_all.R", local = globalenv())
source("midoutbreak_utils.R", local = globalenv())

# ==============================================================================
# LOAD DATA
# ==============================================================================
schools_raw <- read.csv("data/SC_vaccination_merged.csv", stringsAsFactors = FALSE)

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

# Filter to region matching household cache
region_counties <- c("Spartanburg","Greenville")
county_col <- schools_raw$County
valid_rows <- !is.na(as.integer(schools_raw$Total.Students)) &
              as.integer(schools_raw$Total.Students) > 0
schools$county <- as.character(county_col[valid_rows])
schools <- schools[schools$county %in% region_counties, , drop = FALSE]
schools$school_id <- seq_len(nrow(schools))

cat(sprintf("Loaded %d schools in %s\n", nrow(schools),
            paste(region_counties, collapse = ", ")))

# ==============================================================================
# DISEASE PARAMETERS
# ==============================================================================
params <- list(
 # Disease progression (Erlang-distributed durations)
  latent_mean = 10, latent_shape = 8,
  infectious_mean = 8,infectious_shape = 8,
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
# OBSERVED OUTBREAK STATE
# ==============================================================================
# *** EDIT THIS SECTION TO MATCH YOUR CURRENT SITUATION ***
#
# exposed_schools:     Names of schools that have reported measles exposure.
#                      Must match school names in the CSV (partial match OK).
#
# total_cases:         Total confirmed measles cases across ALL exposed schools
#                      so far (cumulative). These will be placed among
#                      unvaccinated students in the listed schools.
#
# quarantine_contacts: Number of contacts currently in quarantine as of today.
#                      Set to 0 if quarantine is not being used.
#
# fraction_active:     Estimated fraction of total_cases that are still
#                      infectious right now (E + P + Ra), not yet recovered.
#                      Default 0.15 = 15% still active, 85% recovered.
#                      Increase if outbreak is early/accelerating;
#                      decrease if outbreak is waning.
#
# hh_attack_rate:      Probability that a susceptible household member of an
#                      infected student has been infected. Default 0.30.
# ==============================================================================

observed_state <- list(
  exposed_schools = c("Fairforest Elementary","Boiling Springs Elementary", "Holly Springs-Motlow Elementary","Raibbow Lake Middle",
                       "Campobello-Gramling School", "Crestview Elementary",
                       "Libertas Academy - Boiling Springs", "Berry Shoals Elementary", "Oakland Elementary",
                       "T. E. Mabry Middle", "Landrum High","Starr Elementary","Global Academy of SC",
                       "Chapman High", "Boiling Springs High","Boiling Springs Middle","Abner Creek Middle",
			"Tyger River Elementary","Sugar Ridge Elementary","Cannons Elementary","Cannons Elementary","Cooley Springs-Fingerville Elementary",
			"Inman Intermediate","James H. Hendrix Elementary", "Jesse S. Bobo Elementary","Mayo Elementary",
			"Sugar Ridge Elementary") ,

  total_cases         = 993,    # Total confirmed cases to date
  quarantine_contacts = 52,    # Contacts in quarantine this week
  fraction_active     = 0.007,  # 15% of cases still infectious
  hh_attack_rate      = 0.90   # 30% household secondary attack rate
)

cat("\n--- OBSERVED STATE ---\n")
cat(sprintf("Exposed schools: %s\n", paste(observed_state$exposed_schools, collapse = ", ")))
cat(sprintf("Total cases: %d\n", observed_state$total_cases))
cat(sprintf("Quarantine contacts: %d\n", observed_state$quarantine_contacts))
cat(sprintf("Fraction active: %.0f%%\n", observed_state$fraction_active * 100))
cat(sprintf("HH attack rate: %.0f%%\n", observed_state$hh_attack_rate * 100))
cat("----------------------\n")

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
    n_schools <- nrow(schools)
    adj_full <- matrix(0, nrow = n_schools, ncol = n_schools)

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

  schools$match_key <- NULL
  schools$ref_id <- NULL

} else {
  cat("Travel time files not found — using Haversine distance network\n")
  if (!file.exists(travel_time_file)) cat(sprintf("  Missing: %s\n", travel_time_file))
  if (!file.exists(school_ref_file))  cat(sprintf("  Missing: %s\n", school_ref_file))
  network <- generate_distance_network(schools, max_distance_km = round(max_travel_time * 0.87))
}

# ==============================================================================
# HOUSEHOLD SETUP
# ==============================================================================
# ==============================================================================
# No longer depends on a pre-built cache. Loads the RTI synthetic population

# CSV and assigns households to students directly, so any region works.

#

# Required: data/synthetic_population_sc.csv (or whatever your synpop file is)

# Expected columns: hh_id, agep, person_id, hh_size, lon_4326, lat_4326

# ==============================================================================



synpop_file <- "data/synthetic_population_sc.csv"

hh_cache_file <- "data/household_assignment_cache.rds"

hh_pop <- NULL

household_assignment <- NULL



populations_template <- lapply(seq_len(nrow(schools)), function(i) {

  create_school_population(i, schools$school_size[i], params$avg_class_size, params$age_range)

})



if (file.exists(synpop_file)) {

  cat("Loading synthetic population...\n")

  synpop <- load_synthetic_population(synpop_file)



  # Check if a valid cache exists for this exact school set

  cache_valid <- FALSE

  if (file.exists(hh_cache_file)) {

    saved <- readRDS(hh_cache_file)

    if (!is.null(saved$n_schools) && saved$n_schools == nrow(schools)) {

      # Quick check: same number of schools suggests same region

      cache_valid <- TRUE

      cat("Found valid household cache matching current school count — using cache\n")

    } else {

      cat(sprintf("Cache mismatch (cache has %d schools, current has %d) — regenerating\n",

                  saved$n_schools, nrow(schools)))

    }

  }



  if (cache_valid) {

    hh_result <- load_household_assignment(hh_cache_file, populations_template, schools, verbose = FALSE)

  } else {

    cat("Assigning households to students (this may take a few minutes)...\n")

    hh_result <- assign_households_to_students(

      populations = populations_template,

      schools = schools,

      synpop = synpop,

      max_distance_km = 25,

      grade_tolerance = 2,

      verbose = TRUE

    )



    # Save cache for other SLURM tasks or future runs with same region

    save_household_assignment(hh_result, hh_cache_file, schools)

    cat(sprintf("Saved new household cache: %s\n", hh_cache_file))

  }



  populations_template <- hh_result$populations

  household_assignment <- hh_result$assignment_df



  populations_template <- assign_household_level_vaccination(

    populations_template, schools, household_assignment, verbose = FALSE

  )



  hh_pop <- create_household_population(hh_result$household_members, populations_template, params)

  hh_pop <- initialize_household_vaccination(populations_template, hh_pop, 0.8, 0.5, 0.95)



  cat(sprintf("Household transmission enabled: %d members in %d households\n",

              nrow(hh_pop), length(unique(household_assignment$hh_id))))



  # Clean up large synpop object

  rm(synpop); gc(verbose = FALSE)



} else if (file.exists(hh_cache_file)) {

  cat("Synpop file not found, falling back to household cache...\n")

  hh_result <- load_household_assignment(hh_cache_file, populations_template, schools, verbose = FALSE)

  populations_template <- hh_result$populations

  household_assignment <- hh_result$assignment_df

  populations_template <- assign_household_level_vaccination(

    populations_template, schools, household_assignment, verbose = FALSE

  )

  hh_pop <- create_household_population(hh_result$household_members, populations_template, params)

  hh_pop <- initialize_household_vaccination(populations_template, hh_pop, 0.8, 0.5, 0.95)

  cat(sprintf("Household transmission enabled: %d members\n", nrow(hh_pop)))

} else {

  cat("WARNING: Neither synpop file nor household cache found\n")

  cat("  Running WITHOUT household transmission\n")

}



use_households <- !is.null(hh_pop)


# ==============================================================================
# RUN SIMULATIONS IN PARALLEL
# ==============================================================================
base_seed <- 54321L + (opt$task_id - 1L) * opt$n_sims
sim_seeds <- base_seed + seq_len(opt$n_sims) - 1L

cat(sprintf("\nRunning %d mid-outbreak forecast simulations...\n", opt$n_sims))
start_time <- Sys.time()

results_list <- mclapply(seq_len(opt$n_sims), function(i) {
  tryCatch({

    # Fresh hh_pop per simulation
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

    result <- run_midoutbreak_simulation(
      schools = schools,
      network = network,
      params = params,
      observed_state = observed_state,
      n_days = opt$n_days,
      seed = sim_seeds[i],
      hh_pop = sim_hh_pop,
      household_assignment = household_assignment
    )

    # Extract results (same field names as run_network_simulation)
    total_inf <- result$total_infected
    schools_aff <- sum(result$school_summary$total_infected > 0)
    total_breakthrough <- result$total_breakthrough
    hh_inf <- result$total_hh_members_infected
    actual_days <- result$actual_days

    # New cases = total at end minus what we started with
    new_student_cases <- total_inf - observed_state$total_cases
    new_hh_cases <- hh_inf - round(
      length(which(!hh_pop$is_student & hh_pop$hh_id %in%
                     unique(unlist(lapply(result$populations[
                       resolve_school_ids(observed_state$exposed_schools, schools)
                     ], function(p) p$hh_id[!is.na(p$hh_id)]))))) *
        observed_state$hh_attack_rate * (1 - observed_state$fraction_active)
    )

    # Peak day
    peak_day <- NA_integer_
    ndc <- result$network_daily_counts
    if (!is.null(ndc)) {
      active_cols <- intersect(names(ndc), c("P", "Ra"))
      if (length(active_cols) > 0) {
        active_sum <- rowSums(ndc[, active_cols, drop = FALSE])
        peak_day <- ndc$day[which.max(active_sum)]
      }
    }

    # Epidemic curve
    epidemic_curve <- ndc
    epidemic_curve$task_id <- opt$task_id
    epidemic_curve$sim_id <- i

    hh_curve <- NULL
    if (!is.null(result$hh_daily_counts)) {
      hh_curve <- result$hh_daily_counts
      hh_curve$task_id <- opt$task_id
      hh_curve$sim_id <- i
    }

    school_summary <- result$school_summary
    school_summary$task_id <- opt$task_id
    school_summary$sim_id <- i

    list(
      summary = data.frame(
        task_id = opt$task_id,
        sim_id = i,
        seed = sim_seeds[i],
        total_infected = total_inf,
        new_student_cases = max(0, new_student_cases),
        schools_affected = schools_aff,
        total_breakthrough = total_breakthrough,
        hh_members_infected = hh_inf,
        actual_days = actual_days,
        peak_day = peak_day,
        observed_cases = observed_state$total_cases,
        observed_quarantine = observed_state$quarantine_contacts,
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
        total_infected = NA, new_student_cases = NA,
        schools_affected = NA, total_breakthrough = NA,
        hh_members_infected = NA, actual_days = NA, peak_day = NA,
        observed_cases = observed_state$total_cases,
        observed_quarantine = observed_state$quarantine_contacts,
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
# SAVE RESULTS
# ==============================================================================
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

summary_df <- do.call(rbind, lapply(results_list, function(x) x$summary))
summary_file <- file.path(opt$output_dir, sprintf("task_%02d_summary.rds", opt$task_id))
saveRDS(summary_df, summary_file)
cat(sprintf("Summary saved: %s\n", summary_file))

curves <- Filter(Negate(is.null), lapply(results_list, function(x) x$epidemic_curve))
if (length(curves) > 0) {
  curves_df <- do.call(rbind, curves)
  curves_file <- file.path(opt$output_dir, sprintf("task_%02d_curves.rds", opt$task_id))
  saveRDS(curves_df, curves_file)
  cat(sprintf("Epidemic curves saved: %s\n", curves_file))
}

hh_curves <- Filter(Negate(is.null), lapply(results_list, function(x) x$hh_curve))
if (length(hh_curves) > 0) {
  hh_curves_df <- do.call(rbind, hh_curves)
  hh_file <- file.path(opt$output_dir, sprintf("task_%02d_hh_curves.rds", opt$task_id))
  saveRDS(hh_curves_df, hh_file)
  cat(sprintf("Household curves saved: %s\n", hh_file))
}

school_sums <- Filter(Negate(is.null), lapply(results_list, function(x) x$school_summary))
if (length(school_sums) > 0) {
  school_sums_df <- do.call(rbind, school_sums)
  school_file_out <- file.path(opt$output_dir, sprintf("task_%02d_schools.rds", opt$task_id))
  saveRDS(school_sums_df, school_file_out)
  cat(sprintf("School summaries saved: %s\n", school_file_out))
}

# Also save the observed_state for reference
saveRDS(observed_state,
        file.path(opt$output_dir, "observed_state.rds"))

# ==============================================================================
# PRINT FORECAST SUMMARY
# ==============================================================================
cat(sprintf("\n=== FORECAST SUMMARY (Task %d) ===\n", opt$task_id))
cat(sprintf("Completed in %.1f sec (%.2f sec/sim)\n\n", elapsed, elapsed / opt$n_sims))

cat(sprintf("Starting from: %d observed cases, %d in quarantine\n",
            observed_state$total_cases, observed_state$quarantine_contacts))
cat(sprintf("Forecast horizon: %d days\n\n", opt$n_days))

cat(sprintf("%-40s %7s   [%5s - %5s]\n", "Metric", "Median", "2.5%", "97.5%"))
cat(strrep("-", 62), "\n")

print_forecast <- function(label, x) {
  cat(sprintf("%-40s %7.0f   [%5.0f - %5.0f]\n",
              label, median(x, na.rm = TRUE),
              quantile(x, 0.025, na.rm = TRUE),
              quantile(x, 0.975, na.rm = TRUE)))
}

print_forecast("Final total student infections",   summary_df$total_infected)
print_forecast("NEW student cases (from today)",    summary_df$new_student_cases)
print_forecast("Schools affected",                  summary_df$schools_affected)
print_forecast("Household members infected",        summary_df$hh_members_infected)
print_forecast("Remaining outbreak duration (days)", summary_df$actual_days)
print_forecast("Peak day (from today)",             summary_df$peak_day)

n_failed <- sum(is.na(summary_df$total_infected))
if (n_failed > 0) {
  cat(sprintf("\nWARNING: %d/%d simulations failed\n", n_failed, opt$n_sims))
}
