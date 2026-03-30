# ==============================================================================
# SCHOOL NETWORK MEASLES SIMULATION WITH HOUSEHOLD TRANSMISSION
# ==============================================================================
# Example script showing how to run simulations with household structure
# Updated: Added parallel processing and household assignment caching
# ==============================================================================

rm(list = ls())

library(readxl)
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)

# ==============================================================================
# LOAD MODULES
# ==============================================================================

# Find load_all.R - try multiple locations
load_all_candidates <- c(
  "load_all.R",                           # Current directory
  here::here("load_all.R"),               # Project root
  here::here("codes", "load_all.R"),      # codes/ subdirectory (preferred)
  here::here("modular", "load_all.R"),    # modular subdirectory
  file.path(getwd(), "load_all.R"),       # Explicit current dir
  file.path(getwd(), "codes", "load_all.R"),  # codes/ under cwd
  file.path(dirname(getwd()), "load_all.R")   # Parent directory
)

load_all_path <- NULL
for (candidate in load_all_candidates) {
  if (file.exists(candidate)) {
    load_all_path <- candidate
    break
  }
}

if (is.null(load_all_path)) {
  stop(paste(
    "Cannot find load_all.R!",
    "\nSearched:", paste(load_all_candidates, collapse = "\n  "),
    "\n\nPlease either:",
    "\n  1. Set working directory to the folder containing load_all.R:",
    "\n     setwd('path/to/modular')",
    "\n  2. Or run this script from that directory"
  ))
}

cat(sprintf("Loading modules from: %s\n", normalizePath(load_all_path)))
source(load_all_path, local = globalenv())

# ==============================================================================
# FILE PATHS (UPDATE THESE FOR YOUR DATA)
# ==============================================================================

# School data with vaccination rates
school_data_file <- "data/SC_vaccination_merged.csv"

# Travel time matrix for network connectivity
travel_time_file <- "data/Travel_time_matrix_Upstate.csv"

# School reference file with coordinates
school_reference_file <- "data/School_reference_Upstate.csv"

# RTI Synthetic population file
# Expected columns: hh_id, agep, person_id, hh_size, lon_4326, lat_4326
synpop_file <- "C:/Users/pandey7/OneDrive - Clemson University/Research/Research/Agent-Based Model/usa_synth_pop/sythetic_population_sc.csv"

# Household assignment cache (speeds up subsequent runs)
hh_cache_file <- "data/household_assignment_cache.rds"

# ==============================================================================
# CONFIGURATION - SIMULATION PARAMETERS
# ==============================================================================

n_simulations <- 10
n_days <- 365
seed_start <- 12345

# Region filter (set to NULL for all counties)
region_counties <- c("Spartanburg")#c("Spartanburg","Greenville","Pickens","Anderson","McCormick","Abbeville",
                    #  "Laurens","Cherokee","Union")

# Parallel processing settings
use_parallel <- TRUE
n_cores <- 10  # NULL = auto-detect optimal number of cores

# ==============================================================================
# CONFIGURATION - DISEASE PARAMETERS
# ==============================================================================

params <- list(
  # Disease progression (Erlang-distributed durations)
  latent_mean = 10,
  latent_shape = 4,
  infectious_mean = 8,
  infectious_shape = 4,
  prodromal_period = 4,
  
  # Within-school transmission
  c_within = 6,           # Mean contacts per day within class
  c_between = 2,          # Mean contacts per day between classes
  p_within = 0.17,       # Per-contact transmission probability (within class)
  p_between = 0.17,       # Per-contact transmission probability (between class)
  
  # Between-school transmission
  c_between_school = 0.2,  # Base contact rate between schools (scaled by network)
  
  # Household transmission
  hh_transmission_prob = 0.17,  # Daily per-contact probability within household
  
  # Infectiousness modifiers
  prodromal_infectiousness_multiplier = 1.0,
  rash_infectiousness_multiplier = 1.0,
  
  # Vaccine parameters
  vaccine_efficacy = 0.97,
  vaccine_infectiousness_reduction = 0.8,
  
  # Intervention parameters
  no_intervention = FALSE,
  isolation_delay_index = 2,
  isolation_delay_secondary = 3,
  isolation_period = 14,
  quarantine_contacts = TRUE,
  quarantine_efficacy = 0.90,
  quarantine_duration = 21,
  
  # Population structure
  avg_class_size = 25,
  age_range = c(5, 18)
)

# ==============================================================================
# CONFIGURATION - HOUSEHOLD VACCINATION CORRELATION
# ==============================================================================
# ==============================================================================
# CONFIGURATION - HOUSEHOLD VACCINATION CORRELATION
# ==============================================================================

# If one child is unvaccinated, probability siblings are also unvaccinated
sibling_unvax_prob <- 0.8

# If one child is unvaccinated, probability adults in household are unvaccinated
adult_unvax_prob <- 0.5

# Baseline adult vaccination coverage (when children are vaccinated)
baseline_adult_coverage <- 0.95

# ==============================================================================
# CONFIGURATION - NETWORK PARAMETERS
# ==============================================================================

network_type <- "travel_time"
max_travel_time <- 16  # minutes
weight_method <- "exponential"

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("=== LOADING DATA ===\n")

# Load vaccination data
school_data <- read.csv(school_data_file)
cat(sprintf("Vaccination data: %d schools\n", nrow(school_data)))

# Load school reference (with coordinates)
school_ref <- read.csv(school_reference_file)

school_ref <- school_ref %>% 
  filter(!is.na(School.Name) & School.Name != "") %>%
  rename(ref_id = OBJECTID)

cat(sprintf("Reference file: %d schools\n", nrow(school_ref)))

# Load synthetic population
synpop <- load_synthetic_population(synpop_file)

# ==============================================================================
# FILTER BY REGION
# ==============================================================================

cat("\n=== FILTERING BY REGION ===\n")

if (!is.null(region_counties)) {
  school_data_filtered <- school_data %>% 
    filter(County %in% region_counties)
  
  school_ref_filtered <- school_ref %>%
    filter(County %in% region_counties)
  
  cat(sprintf("Filtered to %d schools in region\n", nrow(school_data_filtered)))
}

# ==============================================================================
# MATCH SCHOOLS AND CREATE SCHOOL DATA FRAME
# ==============================================================================

cat("\n=== MATCHING SCHOOLS ===\n")

# Create match key (County + School Name)
school_data_filtered <- school_data_filtered %>%
  mutate(match_key = paste(County, School.Name, sep = "_"))

school_ref_filtered <- school_ref_filtered %>%
  mutate(match_key = paste(County, School.Name, sep = "_"))

# Join datasets
matched_schools <- school_data_filtered %>%
  inner_join(
    school_ref_filtered %>% select(match_key, ref_id),
    by = "match_key"
  ) %>%
  filter(!is.na(longitude) & !is.na(latitude))

cat(sprintf("Matched %d schools with coordinates\n", nrow(matched_schools)))


# Prepare schools data frame for simulation
schools <- matched_schools %>%
  mutate(
    school_id = row_number(),
    school_name = School.Name,
    school_size = as.numeric(Total.Students),
    vaccination_coverage = as.numeric(`Percent.Immunized`),
    lon = longitude,
    lat = latitude,
    county = County,
    grade_range_orig = Grade.Range  # Keep original for reference
  ) %>%
  select(school_id, school_name, school_size, vaccination_coverage, 
         lon, lat, county, grade_range_orig) %>%
  filter(!is.na(school_size) & school_size > 0)

cat(sprintf("Final school dataset: %d schools\n", nrow(schools)))

# Sanity check for vaccination coverage
if (max(schools$vaccination_coverage, na.rm = TRUE) < 0.1) {
  warning("Vaccination coverage values appear very low (all < 10%). Check if data needs adjustment.")
  cat("  Current range: ", round(range(schools$vaccination_coverage, na.rm = TRUE), 4), "\n")
} else if (max(schools$vaccination_coverage, na.rm = TRUE) > 1) {
  warning("Vaccination coverage values > 1 detected. Dividing by 100.")
  schools$vaccination_coverage <- schools$vaccination_coverage / 100
}

cat(sprintf("Vaccination coverage range: %.1f%% - %.1f%%\n",
            min(schools$vaccination_coverage, na.rm = TRUE) * 100,
            max(schools$vaccination_coverage, na.rm = TRUE) * 100))

# ==============================================================================
# STANDARDIZE GRADE RANGES
# ==============================================================================

cat("\n=== STANDARDIZING GRADE RANGES ===\n")

# Parse and standardize grade ranges
schools <- standardize_grade_range(schools, grade_col = "grade_range_orig")

# Add human-readable standardized grade string
schools <- add_standardized_grade_string(schools)

# Show some examples
cat("\nGrade range standardization examples:\n")
sample_idx <- sample(1:nrow(schools), min(10, nrow(schools)))
for (i in sample_idx) {
  cat(sprintf("  '%s' -> %s (%s)\n", 
              schools$grade_range_orig[i],
              schools$grade_range_std[i],
              schools$school_type[i]))
}

# ==============================================================================
# GENERATE NETWORK (WITH GRADE WEIGHTING)
# ==============================================================================

cat("\n=== GENERATING NETWORK ===\n")

# Use distance-based network (simpler, uses lon/lat from schools dataframe)
# max_distance_km roughly equivalent to max_travel_time in minutes (assuming ~1 km/min avg)
max_distance_km <- round(max_travel_time * 0.87)  # Approximate conversion

network <- generate_distance_network(
  schools = schools,
  max_distance_km = max_distance_km,
  weight_method = weight_method
)

cat(sprintf("Base network created: %d edges\n", network$n_edges))

# Apply grade-based weighting to network
# Schools with similar grade ranges have higher connection weights
# (e.g., high schools interact more through sports events)
network <- apply_grade_weighting_to_network(
  network = network,
  schools = schools,
  grade_weight = 0.5,      # How much to weight by grade overlap (0-1)
  same_type_bonus = 0.2    # Extra weight for same school type
)

cat("Grade-weighted network ready.\n")

# ==============================================================================
# CREATE POPULATIONS
# ==============================================================================

cat("\n=== CREATING SCHOOL POPULATIONS ===\n")

populations <- list()
for (i in 1:nrow(schools)) {
  populations[[i]] <- create_school_population(
    school_id = i,
    school_size = schools$school_size[i],
    avg_class_size = params$avg_class_size,
    age_range = params$age_range
  )
  # NOTE: Vaccination is assigned AFTER household assignment (below)
  # to ensure siblings share vaccination status
}

cat(sprintf("Created %d school populations\n", length(populations)))

# ==============================================================================
# ASSIGN HOUSEHOLDS TO STUDENTS (WITH CACHING)
# ==============================================================================

cat("\n=== ASSIGNING HOUSEHOLDS ===\n")

# Check if cached household assignment exists
if (file.exists(hh_cache_file)) {
  cat("Loading cached household assignment...\n")
  
  hh_result <- load_household_assignment(
    filepath = hh_cache_file,
    populations = populations,
    schools = schools,
    verbose = TRUE
  )
  
} else {
  cat("Computing household assignment (this will be cached for future runs)...\n")
  
  hh_result <- assign_households_to_students(
    populations = populations,
    schools = schools,
    synpop = synpop,
    max_distance_km = 30,
    grade_tolerance = 2
  )
  
  # Save for future runs
  save_household_assignment(
    hh_result = hh_result,
    filepath = hh_cache_file,
    schools = schools
  )
}

populations <- hh_result$populations
household_assignment <- hh_result$assignment_df
sibling_links <- hh_result$sibling_links
household_members <- hh_result$household_members

cat(sprintf("Assigned %d households\n", length(unique(household_assignment$hh_id))))
cat(sprintf("Households with siblings: %d\n", nrow(sibling_links)))
if ("is_cross_school" %in% names(sibling_links)) {
  cat(sprintf("  - Same-school siblings: %d\n", sum(!sibling_links$is_cross_school)))
  cat(sprintf("  - Cross-school siblings: %d\n", sum(sibling_links$is_cross_school)))
}

# ==============================================================================
# ASSIGN HOUSEHOLD-LEVEL VACCINATION (NEW ORDER - AFTER HOUSEHOLD ASSIGNMENT)
# ==============================================================================
# This ensures siblings share vaccination status while approximating
# target school-level coverage rates

populations <- assign_household_level_vaccination(
  populations = populations,
  schools = schools,
  assignment_df = household_assignment,
  verbose = TRUE
)

# ==============================================================================
# CREATE HOUSEHOLD POPULATION
# ==============================================================================

cat("\n=== CREATING HOUSEHOLD POPULATION ===\n")

hh_pop <- create_household_population(
  household_members = household_members,
  populations = populations,
  params = params
)

# Initialize household-correlated vaccination for non-students
# (copies student vaccination status and sets adult/young child status)
hh_pop <- initialize_household_vaccination(
  populations = populations,
  hh_pop = hh_pop,
  sibling_unvax_prob = sibling_unvax_prob,
  adult_unvax_prob = adult_unvax_prob,
  baseline_adult_coverage = baseline_adult_coverage
)

# ==============================================================================
# SELECT SEED SCHOOL(S)
# ==============================================================================

# Seed in school 
seed_schools <- schools %>%
  arrange(vaccination_coverage) %>%
  slice(7,8) %>%
  pull(school_id)

cat(sprintf("\nSeeding outbreak in school %d (%s, %.1f%% vaccinated)\n",
            seed_schools,
            schools$school_name[seed_schools],
            schools$vaccination_coverage[seed_schools] * 100))

# ==============================================================================
# RUN SIMULATIONS (PARALLEL OR SEQUENTIAL)
# ==============================================================================

cat("\n=== RUNNING SIMULATIONS ===\n")

if (use_parallel) {
  # Show parallel settings
  parallel_info <- detect_parallel_settings()
  
  cat(sprintf("\nRunning %d simulations in PARALLEL mode\n", n_simulations))
  
  results <- run_parallel_simulations(
    n_simulations = n_simulations,
    schools = schools,
    network = network,
    params = params,
    seed_schools = seed_schools,
    n_initial_infected = 5,
    n_days = n_days,
    seed_start = seed_start,
    n_cores = n_cores,
    verbose = TRUE,
    hh_pop = hh_pop,
    household_assignment = household_assignment,
    source_files = here("load_all.R")
  )
  
} else {
  cat(sprintf("\nRunning %d simulations in SEQUENTIAL mode\n", n_simulations))
  
  results <- run_fast_simulations(
    n_simulations = n_simulations,
    schools = schools,
    network = network,
    params = params,
    seed_schools = seed_schools,
    n_initial_infected = 5,
    n_days = n_days,
    seed_start = seed_start,
    verbose = TRUE,
    hh_pop = hh_pop,
    household_assignment = household_assignment
  )
}
# ==============================================================================
# ANALYZE RESULTS
# ==============================================================================
plot_outbreak_dynamics_report(results = results)
cat("\n=== ANALYSIS ===\n")

# Summary statistics
summary_stats <- results$summary_stats

cat("\n--- INFECTION SUMMARY ---\n")
cat(sprintf("Students infected (Median): %.0f [95%% CI: %.0f - %.0f]\n",
            median(summary_stats$total_infected),
            quantile(summary_stats$total_infected, 0.025),
            quantile(summary_stats$total_infected, 0.975)))

cat(sprintf("Household members infected (Median): %.0f [95%% CI: %.0f - %.0f]\n",
            median(summary_stats$hh_members_infected),
            quantile(summary_stats$hh_members_infected, 0.025),
            quantile(summary_stats$hh_members_infected, 0.975)))

total_infected <- summary_stats$total_infected + summary_stats$hh_members_infected
cat(sprintf("Total infected (Median): %.0f [95%% CI: %.0f - %.0f]\n",
            median(total_infected),
            quantile(total_infected, 0.025),
            quantile(total_infected, 0.975)))

cat(sprintf("\nSchools affected (Median): %.0f [95%% CI: %.0f - %.0f]\n",
            median(summary_stats$schools_affected),
            quantile(summary_stats$schools_affected, 0.025),
            quantile(summary_stats$schools_affected, 0.975)))

cat(sprintf("Outbreak duration (Median): %.0f days [95%% CI: %.0f - %.0f]\n",
            median(summary_stats$actual_days),
            quantile(summary_stats$actual_days, 0.025),
            quantile(summary_stats$actual_days, 0.975)))

cat(sprintf("\nComputation time: %.1f seconds (%.2f sec/sim)\n",
            results$computation_time,
            results$computation_time / n_simulations))

# ==============================================================================
# PLOTTING
# ==============================================================================

# Outbreak dynamics
p_dynamics <- plot_network_outbreak_dynamics(results)
p_dynamics


# Compare students vs household members
if (results$use_households) {
  infection_comparison <- data.frame(
    sim_id = 1:n_simulations,
    students = summary_stats$total_infected,
    hh_members = summary_stats$hh_members_infected
  ) %>%
    pivot_longer(cols = c(students, hh_members),
                 names_to = "group", values_to = "infected")
  
  p_compare <- ggplot(infection_comparison, aes(x = group, y = infected, fill = group)) +
    geom_boxplot(alpha = 0.7) +
    labs(title = "Infections by Group",
         x = "", y = "Number Infected") +
    theme_minimal() +
    theme(legend.position = "none")
  
  
}
p_compare

