# ==============================================================================
# ANALYSIS UTILITIES
# ==============================================================================
# File: analysis_utils.R
# Contains: Summary statistics, infection timing analysis, tables
# Dependencies: dplyr
# ==============================================================================

library(dplyr)

# ==============================================================================
# First Infection Time Calculation
# ==============================================================================

#' Calculate first infection time for each school in each simulation
#' @param results Results object from run_multiple_network_simulations
#' @return Data frame with first infection times by school and simulation
calculate_first_infection_times <- function(results) {
  
  # Use school daily data which has state counts by day
  school_daily <- results$all_school_daily_data
  
  # For each simulation and school, find the first day with any infection
  first_infection <- school_daily %>%
    mutate(
      ever_infected = E + P + Ra + Iso + R
    ) %>%
    filter(ever_infected > 0) %>%
    group_by(sim, school_id) %>%
    summarise(
      first_infection_day = min(day),
      school_name = first(school_name),
      school_size = first(school_size),
      vaccination_coverage = first(vaccination_coverage),
      .groups = "drop"
    )
  
  # Add schools that never got infected
  all_schools <- results$schools %>% select(school_id, school_name)
  all_sims <- unique(school_daily$sim)
  
  complete_grid <- expand.grid(
    sim = all_sims,
    school_id = all_schools$school_id,
    stringsAsFactors = FALSE
  ) %>%
    left_join(all_schools, by = "school_id")
  
  first_infection_complete <- complete_grid %>%
    left_join(
      first_infection %>% select(sim, school_id, first_infection_day),
      by = c("sim", "school_id")
    ) %>%
    left_join(
      results$schools %>% select(school_id, school_size, vaccination_coverage),
      by = "school_id"
    )
  
  return(first_infection_complete)
}


# ==============================================================================
# Infection Timing Summary Statistics
# ==============================================================================

#' Calculate summary statistics of infection timing by school
#' @param first_infection_times Output from calculate_first_infection_times
#' @return Data frame with median, mean, IQR of first infection times
summarize_infection_timing <- function(first_infection_times) {
  
  summary_stats <- first_infection_times %>%
    group_by(school_id, school_name, school_size, vaccination_coverage) %>%
    summarise(
      n_simulations = n(),
      n_infected = sum(!is.na(first_infection_day)),
      prob_infected = mean(!is.na(first_infection_day)),
      median_first_day = median(first_infection_day, na.rm = TRUE),
      mean_first_day = mean(first_infection_day, na.rm = TRUE),
      sd_first_day = sd(first_infection_day, na.rm = TRUE),
      q25_first_day = quantile(first_infection_day, 0.25, na.rm = TRUE),
      q75_first_day = quantile(first_infection_day, 0.75, na.rm = TRUE),
      min_first_day = min(first_infection_day, na.rm = TRUE),
      max_first_day = max(first_infection_day, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      median_first_day = ifelse(is.infinite(median_first_day), NA, median_first_day),
      mean_first_day = ifelse(is.nan(mean_first_day), NA, mean_first_day),
      min_first_day = ifelse(is.infinite(min_first_day), NA, min_first_day),
      max_first_day = ifelse(is.infinite(max_first_day), NA, max_first_day)
    ) %>%
    arrange(median_first_day)
  
  return(summary_stats)
}


# ==============================================================================
# Outbreak Summary Table
# ==============================================================================

#' Create a summary table of outbreak statistics
#' @param results Results from run_multiple_network_simulations()
#' @return Data frame with summary statistics
create_outbreak_summary_table <- function(results) {
  stats <- results$summary_stats
  
  summary_table <- data.frame(
    Metric = c(
      "Total Infected (Median)",
      "Total Infected (Mean)",
      "Total Infected (95% CI)",
      "Schools Affected (Median)",
      "Schools Affected (Mean)",
      "Duration (Median days)",
      "Outbreak Probability (>10 cases)"
    ),
    Value = c(
      round(median(stats$total_infected), 1),
      round(mean(stats$total_infected), 1),
      paste0(round(quantile(stats$total_infected, 0.025), 1), " - ", 
             round(quantile(stats$total_infected, 0.975), 1)),
      round(median(stats$schools_affected), 1),
      round(mean(stats$schools_affected), 1),
      round(median(stats$actual_days), 0),
      paste0(round(mean(stats$total_infected > 10) * 100, 1), "%")
    )
  )
  
  return(summary_table)
}


# ==============================================================================
# Infection Sequence Table
# ==============================================================================

#' Create summary table of infection sequence
#' @param results Results object from run_multiple_network_simulations
#' @param seed_schools Vector of seed school IDs
#' @return Data frame formatted for display
create_infection_sequence_table <- function(results, seed_schools = NULL) {
  
  first_times <- calculate_first_infection_times(results)
  timing_summary <- summarize_infection_timing(first_times)
  
  # Add seed indicator and rank
  table_data <- timing_summary %>%
    mutate(
      is_seed = ifelse(!is.null(seed_schools) & school_id %in% seed_schools, "Yes", "No"),
      infection_rank = rank(median_first_day, ties.method = "min", na.last = "keep")
    ) %>%
    arrange(infection_rank) %>%
    select(
      Rank = infection_rank,
      School = school_name,
      `Seed` = is_seed,
      `Size` = school_size,
      `Vax %` = vaccination_coverage,
      `Prob. Infected` = prob_infected,
      `Median Day` = median_first_day,
      `IQR` = q25_first_day,
      `IQR_upper` = q75_first_day
    ) %>%
    mutate(
      `Vax %` = round(`Vax %` * 100, 1),
      `Prob. Infected` = paste0(round(`Prob. Infected` * 100, 1), "%"),
      `Median Day` = round(`Median Day`, 1),
      `IQR` = paste0(round(`IQR`, 0), "-", round(`IQR_upper`, 0))
    ) %>%
    select(-IQR_upper)
  
  return(table_data)
}


# ==============================================================================
# Per-School Statistics
# ==============================================================================

#' Calculate per-school statistics across simulations
#' @param results Results from run_multiple_network_simulations()
#' @return Data frame with per-school statistics
calculate_school_statistics <- function(results) {
  
  school_stats <- results$all_school_data %>%
    group_by(school_id, school_name) %>%
    summarise(
      size = first(school_size),
      vax_pct = round(first(vaccination_coverage) * 100, 1),
      median_infected = median(total_infected),
      mean_infected = round(mean(total_infected), 1),
      sd_infected = round(sd(total_infected), 1),
      attack_rate_pct = round(median(attack_rate) * 100, 2),
      pct_outbreaks = round(mean(total_infected > 0) * 100, 1),
      .groups = "drop"
    ) %>%
    arrange(desc(median_infected))
  
  return(school_stats)
}


# ==============================================================================
# Outbreak Probability Table
# ==============================================================================

#' Create probability table for various outbreak thresholds
#' @param results Results from run_multiple_network_simulations()
#' @param thresholds Vector of case thresholds
#' @return Data frame with probabilities
create_probability_table <- function(results, thresholds = c(5, 10, 20, 50, 100, 200, 500)) {
  
  prob_table <- data.frame(
    Threshold = paste0(">", thresholds),
    Count = sapply(thresholds, function(t) sum(results$summary_stats$total_infected > t)),
    Probability = sapply(thresholds, function(t) {
      paste0(round(mean(results$summary_stats$total_infected > t) * 100, 1), "%")
    })
  )
  
  return(prob_table)
}

cat("Analysis utilities loaded.\n")


# ==============================================================================
# School Infection Summary with Credible Intervals
# ==============================================================================

#' Create summary table of infections per school with credible intervals
#' @param results Results object from run_multiple_network_simulations
#' @param ci_level Credible interval level (default 0.95 for 95% CI)
#' @return Data frame with mean, median, and CI for each school
create_school_infection_summary <- function(results, ci_level = 0.95) {
  
  lower_q <- (1 - ci_level) / 2
  upper_q <- 1 - lower_q
  
  school_summary <- results$all_school_data %>%
    group_by(school_id, school_name) %>%
    summarise(
      n_simulations = n(),
      mean_infected = mean(total_infected),
      sd_infected = sd(total_infected),
      median_infected = median(total_infected),
      ci_lower = quantile(total_infected, lower_q),
      ci_upper = quantile(total_infected, upper_q),
      min_infected = min(total_infected),
      max_infected = max(total_infected),
      prob_any_infection = mean(total_infected > 0),
      mean_attack_rate = mean(attack_rate),
      .groups = "drop"
    ) %>%
    left_join(
      results$schools %>% select(school_id, school_size, vaccination_coverage),
      by = "school_id"
    ) %>%
    mutate(
      is_seed = if (!is.null(results$seed_schools)) school_id %in% results$seed_schools else FALSE
    ) %>%
    arrange(desc(mean_infected))
  
  return(school_summary)
}


#' Format school infection summary for display/export
#' @param school_summary Output from create_school_infection_summary
#' @param digits Number of decimal places
#' @return Formatted data frame
format_school_infection_table <- function(school_summary, digits = 1) {
  
  formatted <- school_summary %>%
    mutate(
      `School` = school_name,
      `Size` = school_size,
      `Vax %` = round(vaccination_coverage * 100, 1),
      `Mean Infected` = round(mean_infected, digits),
      `95% CI` = paste0("(", round(ci_lower, digits), " - ", round(ci_upper, digits), ")"),
      `Median` = round(median_infected, digits),
      `Attack Rate %` = round(mean_attack_rate * 100, 2),
      `P(Outbreak)` = paste0(round(prob_any_infection * 100, 1), "%"),
      `Seed` = ifelse(is_seed, "Yes", "")
    ) %>%
    select(School, Size, `Vax %`, `Mean Infected`, `95% CI`, Median, 
           `Attack Rate %`, `P(Outbreak)`, Seed)
  
  return(formatted)
}
