# ==============================================================================
# SIMULATION CORE UTILITIES (WITH HOUSEHOLD TRANSMISSION)
# ==============================================================================
# File: simulation_utils.R
# Contains: Core simulation functions for running network simulations
# Updated: Added household transmission layer
# Dependencies: rcpp_transmission.R, population_utils.R, network_utils.R, 
#               household_utils.R
# ==============================================================================

library(dplyr)

# ==============================================================================
# Within-School Transmission Wrapper
# ==============================================================================

#' Within-school transmission (wrapper for Rcpp)
#' @param population School population data frame
#' @param params Simulation parameters
#' @param contact_history ContactHistory object
#' @return List with updated population and contact history
school_transmission <- function(population, params, contact_history) {
  result <- cpp_school_transmission_contacts(
    student_id      = population$student_id,
    class_id        = population$class_id,
    state           = population$state,
    is_vaccinated   = population$is_vaccinated,
    is_isolated     = population$is_isolated,
    is_quarantined  = population$is_quarantined,
    c_within        = params$c_within,
    c_between       = params$c_between,
    p_within        = params$p_within,
    p_between       = params$p_between,
    prodromal_mult  = params$prodromal_infectiousness_multiplier,
    rash_mult       = params$rash_infectiousness_multiplier,
    vaccine_reduction = params$vaccine_infectiousness_reduction,
    vaccine_efficacy  = params$vaccine_efficacy
  )
  
  # Update contact history
  if (length(result$contact_infector_ids) > 0) {
    contact_history$add_contacts(result$contact_infector_ids, result$contact_target_ids)
  }
  
  # Apply new exposures
  if (length(result$new_exposures) > 0) {
    for (exposed_idx in result$new_exposures) {
      population$state[exposed_idx] <- "E"
      population$time_in_state[exposed_idx] <- 0
      population$latent_duration[exposed_idx] <- draw_erlang(1, params$latent_mean, params$latent_shape)
      infectious_dur <- draw_erlang(1, params$infectious_mean, params$infectious_shape)
      population$infectious_duration[exposed_idx] <- infectious_dur
      population$prodromal_duration[exposed_idx] <- params$prodromal_period
      population$rash_duration[exposed_idx] <- max(1, infectious_dur - params$prodromal_period)
    }
  }
  
  # Mark breakthrough infections
  if (length(result$breakthrough_cases) > 0) {
    population$breakthrough_infection[result$breakthrough_cases] <- TRUE
  }
  
  return(list(
    population = population,
    contact_history = contact_history
  ))
}


# ==============================================================================
# Between-School Transmission
# ==============================================================================
# ==============================================================================

# CORRECTED: Between-School Transmission (density-dependent contact sampling)

# ==============================================================================

#

# FIX: Contacts are now sampled from the FULL school population (all states),

# not just S/V. This correctly models survey-derived contact rates which

# represent all inter-school encounters, not just encounters with susceptible

# individuals. Contacts that land on non-susceptible individuals (E, P, Ra, R)

# are "wasted" — they consume a contact event but produce no transmission.

#

# This matters significantly in mid-outbreak scenarios where a substantial

# fraction of a school may already be infected/recovered. With the old code,

# the effective between-school force of infection was artificially constant

# regardless of how depleted the susceptible pool was.

#

# Replace the between_school_transmission function in simulation_utils.R

# with this corrected version.

# ==============================================================================



between_school_transmission <- function(populations, network, params) {

  n_schools <- length(populations)

  adj <- network$adjacency



  # Find schools with infectious individuals (P or Ra, not isolated/quarantined)

  schools_with_infectors <- c()

  infectors_by_school <- list()



  for (school_idx in 1:n_schools) {

    pop <- populations[[school_idx]]

    infectious_mask <- (pop$state == "P" | pop$state == "Ra") &

                       !pop$is_isolated & !pop$is_quarantined

    if (sum(infectious_mask) > 0) {

      schools_with_infectors <- c(schools_with_infectors, school_idx)

      infectors_by_school[[as.character(school_idx)]] <-

        pop[infectious_mask, c("student_id", "school_id", "state", "is_vaccinated")]

    }

  }



  if (length(schools_with_infectors) == 0) {

    return(populations)

  }



  # -----------------------------------------------------------------------

  # FIX: Build contact pools from ALL non-isolated/non-quarantined students,

  # regardless of disease state. Survey-derived contact rates (Poisson mean)

  # represent total inter-school encounters with any individual, not just

  # encounters with susceptible people.

  # -----------------------------------------------------------------------

  all_members_by_school <- list()

  for (school_idx in 1:n_schools) {

    pop <- populations[[school_idx]]

    # Include everyone who is physically present at school

    # (exclude isolated and quarantined — they are not mixing)

    present_mask <- !pop$is_isolated & !pop$is_quarantined

    if (sum(present_mask) > 0) {

      all_members_by_school[[as.character(school_idx)]] <-

        pop[present_mask, c("student_id", "school_id", "state", "is_vaccinated")]

    }

  }



  # For each school with infectors, find connected schools

  for (source_school in schools_with_infectors) {

    school_infectors <- infectors_by_school[[as.character(source_school)]]



    # Find connected schools (edge weight > 0)

    connected_mask <- adj[source_school, ] > 0

    connected_school_ids <- which(connected_mask)

    connected_school_ids <- connected_school_ids[connected_school_ids != source_school]



    if (length(connected_school_ids) == 0) next



    # Build FULL target pools for connected schools

    target_pools <- list()

    edge_weights <- numeric()

    valid_schools <- integer()



    for (target_school in connected_school_ids) {

      target_key <- as.character(target_school)



      if (!is.null(all_members_by_school[[target_key]])) {

        target_pools <- append(target_pools, list(all_members_by_school[[target_key]]))

        edge_weights <- c(edge_weights, adj[source_school, target_school])

        valid_schools <- c(valid_schools, target_school)

      }

    }



    if (length(target_pools) == 0) next



    # Call Rcpp function — it samples contacts from the FULL pool.

    # Contacts landing on non-S/V individuals will "succeed" in the C++ code

    # but are filtered out below (only S/V contacts result in actual infection).

    result <- cpp_between_school_transmission(

      infector_ids = school_infectors$student_id,

      infector_school_ids = school_infectors$school_id,

      infector_states = school_infectors$state,

      infector_vaccinated = school_infectors$is_vaccinated,

      target_pools = target_pools,

      connected_schools = valid_schools,

      edge_weights = edge_weights,

      c_between_base = params$c_between_school,

      p_base = params$p_between,

      prodromal_mult = params$prodromal_infectiousness_multiplier,

      rash_mult = params$rash_infectiousness_multiplier,

      vaccine_reduction = params$vaccine_infectiousness_reduction,

      vaccine_efficacy = params$vaccine_efficacy

    )



    # Apply exposures — ONLY to targets still in S or V

    # Contacts that landed on E, P, Ra, R, Iso, Q* are correctly discarded

    if (length(result$exposed_student_ids) > 0) {

      for (i in seq_along(result$exposed_student_ids)) {

        target_school <- result$exposed_school_ids[i]

        target_student_id <- result$exposed_student_ids[i]



        pop <- populations[[target_school]]

        idx <- which(pop$student_id == target_student_id)



        if (length(idx) > 0 && pop$state[idx] %in% c("S", "V")) {

          pop$state[idx] <- "E"

          pop$time_in_state[idx] <- 0

          pop$latent_duration[idx] <- draw_erlang(1, params$latent_mean, params$latent_shape)

          infectious_dur <- draw_erlang(1, params$infectious_mean, params$infectious_shape)

          pop$infectious_duration[idx] <- infectious_dur

          pop$prodromal_duration[idx] <- params$prodromal_period

          pop$rash_duration[idx] <- max(1, infectious_dur - params$prodromal_period)

          pop$infection_source[idx] <- "between_school"



          if (result$is_breakthrough[i]) {

            pop$breakthrough_infection[idx] <- TRUE

          }



          populations[[target_school]] <- pop

        }

      }

    }

  }



  return(populations)

}
 

# ==============================================================================
# Run Single Network Simulation (WITH HOUSEHOLDS)
# ==============================================================================

#' Run a single network simulation with household transmission
#' @param schools Data frame with school information
#' @param network Network object
#' @param params Simulation parameters
#' @param seed_schools Vector of school IDs where outbreak starts
#' @param n_initial_infected Number of initial infections per seeded school
#' @param n_days Maximum simulation days
#' @param seed Random seed (optional)
#' @param hh_pop Household population data frame (optional, if NULL no HH transmission)
#' @param household_assignment Assignment data frame linking students to households
#' @return List with daily counts, summary stats, and final populations
run_network_simulation <- function(schools, network, params, 
                                    seed_schools, n_initial_infected = 1,
                                    n_days = 150, seed = NULL,
                                    hh_pop = NULL, household_assignment = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  n_schools <- nrow(schools)
  use_households <- !is.null(hh_pop)
  
  # Initialize populations
  populations <- list()
  contact_histories <- list()
  
  for (i in 1:n_schools) {
    populations[[i]] <- create_school_population(
      school_id = i,
      school_size = schools$school_size[i],
      avg_class_size = params$avg_class_size,
      age_range = params$age_range
    )
    
    populations[[i]] <- initialize_vaccination_school(
      populations[[i]], 
      schools$vaccination_coverage[i]
    )
    
    contact_histories[[i]] <- ContactHistory$new(window_size = params$isolation_delay_index)
  }
  
  # If households provided, assign household IDs to student populations
  if (use_households && !is.null(household_assignment)) {
    for (i in 1:n_schools) {
      school_assign <- household_assignment %>%
        filter(school_id == i)
      
      if (nrow(school_assign) > 0) {
        n_assign <- min(nrow(populations[[i]]), nrow(school_assign))
        
        if (!"hh_id" %in% names(populations[[i]])) {
          populations[[i]]$hh_id <- NA_integer_
          populations[[i]]$synpop_person_id <- NA_integer_
        }
        
        populations[[i]]$hh_id[1:n_assign] <- school_assign$hh_id[1:n_assign]
        populations[[i]]$synpop_person_id[1:n_assign] <- school_assign$synpop_person_id[1:n_assign]
      }
    }
  }
  
  # Seed infections
  populations <- seed_infections(populations, seed_schools, n_initial_infected, params)
  
  # If using households, sync student states to household population
  if (use_households) {
    hh_pop <- sync_student_states(populations, hh_pop, direction = "school_to_hh")
  }
  
  # Prepare daily counts storage
  state_cols <- c("S", "V", "E", "P", "Ra", "Iso", "R", "QS", "QV", "QE", "QP")
  
  daily_counts_list <- lapply(1:n_schools, function(s) {
    matrix(0, nrow = n_days, ncol = length(state_cols),
           dimnames = list(NULL, state_cols))
  })
  
  network_daily_counts <- matrix(0, nrow = n_days, ncol = length(state_cols),
                                  dimnames = list(NULL, state_cols))
  
  # Household daily counts if applicable
  if (use_households) {
    hh_daily_counts <- matrix(0, nrow = n_days, ncol = length(state_cols),
                               dimnames = list(NULL, state_cols))
  }
  
  outbreak_ended <- FALSE
  actual_days <- n_days
  
  # Main simulation loop
  for (day in 1:n_days) {
    
    if (!outbreak_ended) {
      # Step 1: Within-school transmission
      for (school_idx in 1:n_schools) {
        trans_result <- school_transmission(
          populations[[school_idx]], params, contact_histories[[school_idx]]
        )
        populations[[school_idx]] <- trans_result$population
        contact_histories[[school_idx]] <- trans_result$contact_history
      }
      
      # Step 2: Between-school transmission
      populations <- between_school_transmission(populations, network, params)
      
      # Step 3: Household transmission (if enabled)
      if (use_households) {
        # First sync school states to household
        hh_pop <- sync_student_states(populations, hh_pop, direction = "school_to_hh")
        
        # Perform household transmission
        hh_result <- household_transmission(populations, hh_pop, params)
        populations <- hh_result$populations
        hh_pop <- hh_result$hh_pop
        
        # Update non-student household member disease states
        hh_pop <- update_household_disease_states(hh_pop, params)
      }
      
      # Step 4: Update disease states
      for (school_idx in 1:n_schools) {
        populations[[school_idx]] <- update_disease_states(populations[[school_idx]], params)
        
        if (params$quarantine_contacts && !params$no_intervention) {
          populations[[school_idx]] <- apply_quarantine(
            populations[[school_idx]], params, contact_histories[[school_idx]]
          )
        }
      }
    }
    
    # Record daily counts
    for (school_idx in 1:n_schools) {
      daily_counts_list[[school_idx]][day, ] <- table(
        factor(populations[[school_idx]]$state, levels = state_cols)
      )
      network_daily_counts[day, ] <- network_daily_counts[day, ] + daily_counts_list[[school_idx]][day, ]
    }
    
    # Record household counts if applicable
    if (use_households) {
      # Non-student household members only (base R for speed)
      hh_non_student <- hh_pop[!hh_pop$is_student, ]
      if (nrow(hh_non_student) > 0) {
        hh_daily_counts[day, ] <- table(
          factor(hh_non_student$state, levels = state_cols)
        )
      }
    }
    
    # Check if outbreak ended
    if (!outbreak_ended) {
      active_states <- c("E", "P", "Ra", "QE", "QP", "Iso")
      total_active <- sum(network_daily_counts[day, active_states])
      
      # Also check household members if applicable
      if (use_households) {
        total_active <- total_active + sum(hh_daily_counts[day, active_states])
      }
      
      if (total_active == 0) {
        actual_days <- day
        outbreak_ended <- TRUE
      }
    }
  }
  
  # Calculate summary statistics
  school_summary <- data.frame(
    school_id = 1:n_schools,
    school_name = schools$school_name,
    school_size = schools$school_size,
    vaccination_coverage = schools$vaccination_coverage,
    total_infected = sapply(populations, function(p) {
      sum(p$state %in% c("P", "Ra", "Iso", "R", "QP"))
    }),
    breakthrough_infections = sapply(populations, function(p) {
      sum(p$breakthrough_infection)
    }),
    was_seeded = 1:n_schools %in% seed_schools
  )
  
  school_summary$attack_rate <- school_summary$total_infected / school_summary$school_size
  
  # Household summary if applicable
  if (use_households) {
    hh_summary <- summarize_household_infections(populations, hh_pop)
    
    # Calculate household-level statistics
    total_hh_members_infected <- sum(hh_pop$state %in% c("P", "Ra", "R", "E") & !hh_pop$is_student)
    total_hh_members <- sum(!hh_pop$is_student)
  } else {
    hh_summary <- NULL
    total_hh_members_infected <- 0
    total_hh_members <- 0
  }
  
  # Prepare output
  network_daily_df <- as.data.frame(network_daily_counts) %>%
    mutate(day = 0:(n_days - 1))
  
  school_daily_dfs <- lapply(1:n_schools, function(s) {
    as.data.frame(daily_counts_list[[s]]) %>%
      mutate(day = 0:(n_days - 1), school_id = s)
  })
  
  results <- list(
    network_daily_counts = network_daily_df,
    school_daily_counts = school_daily_dfs,
    school_summary = school_summary,
    total_infected = sum(school_summary$total_infected),
    total_breakthrough = sum(school_summary$breakthrough_infections),
    actual_days = actual_days,
    populations = populations,
    network = network,
    # Household results
    use_households = use_households,
    hh_summary = hh_summary,
    total_hh_members_infected = total_hh_members_infected,
    total_hh_members = total_hh_members,
    hh_pop = if (use_households) hh_pop else NULL
  )
  
  if (use_households) {
    results$hh_daily_counts <- as.data.frame(hh_daily_counts) %>%
      mutate(day = 0:(n_days - 1))
  }
  
  return(results)
}


# ==============================================================================
# Run Multiple Network Simulations (WITH HOUSEHOLDS)
# ==============================================================================

#' Run multiple network simulations with optional household transmission
#' @param n_simulations Number of simulations to run
#' @param schools Data frame with school information
#' @param network Network object
#' @param params Simulation parameters
#' @param seed_schools Vector of school IDs to seed
#' @param n_initial_infected Number of initial infections per seeded school
#' @param n_days Maximum simulation days
#' @param seed_start Starting seed for reproducibility
#' @param verbose Print progress (default TRUE)
#' @param hh_pop Household population data frame (optional)
#' @param household_assignment Assignment data frame linking students to households
#' @return List with combined results from all simulations
run_multiple_network_simulations <- function(n_simulations, 
                                              schools, 
                                              network, 
                                              params,
                                              seed_schools,
                                              n_initial_infected = 1,
                                              n_days = 150,
                                              seed_start = NULL,
                                              verbose = TRUE,
                                              hh_pop = NULL,
                                              household_assignment = NULL) {
  
  use_households <- !is.null(hh_pop)
  
  if (verbose) {
    cat("\n=== RUNNING MULTIPLE NETWORK SIMULATIONS ===\n")
    cat(sprintf("Number of simulations: %d\n", n_simulations))
    cat(sprintf("Number of schools: %d\n", nrow(schools)))
    cat(sprintf("Network edges: %d\n", network$n_edges))
    cat(sprintf("Seed school(s): %s\n", paste(seed_schools, collapse = ", ")))
    cat(sprintf("Household transmission: %s\n", ifelse(use_households, "ENABLED", "DISABLED")))
    if (use_households) {
      cat(sprintf("Total household members: %d\n", nrow(hh_pop)))
    }
    cat("\n")
  }
  
  all_network_daily <- list()
  all_school_summaries <- list()
  all_school_daily <- list()
  all_hh_daily <- list()
  all_hh_summaries <- list()
  
  summary_stats <- data.frame(
    sim_id = 1:n_simulations,
    total_infected = numeric(n_simulations),
    total_breakthrough = numeric(n_simulations),
    schools_affected = numeric(n_simulations),
    actual_days = numeric(n_simulations),
    hh_members_infected = numeric(n_simulations)
  )
  
  start_time <- Sys.time()
  
  for (i in 1:n_simulations) {
    sim_seed <- if (!is.null(seed_start)) seed_start + i else NULL
    
    # Create fresh copy of hh_pop for each simulation
    sim_hh_pop <- if (use_households) {
      hh_pop %>%
        mutate(
          state = ifelse(is_vaccinated, "V", "S"),
          time_in_state = 0,
          time_since_prodromal = NA_integer_,
          is_isolated = FALSE,
          is_quarantined = FALSE,
          breakthrough_infection = FALSE
        )
    } else {
      NULL
    }
    
    result <- run_network_simulation(
      schools = schools,
      network = network,
      params = params,
      seed_schools = seed_schools,
      n_initial_infected = n_initial_infected,
      n_days = n_days,
      seed = sim_seed,
      hh_pop = sim_hh_pop,
      household_assignment = household_assignment
    )
    
    result$network_daily_counts$sim_id <- i
    all_network_daily[[i]] <- result$network_daily_counts
    
    result$school_summary$sim_id <- i
    all_school_summaries[[i]] <- result$school_summary
    
    # Collect school daily counts
    for (s in seq_along(result$school_daily_counts)) {
      school_daily_df <- result$school_daily_counts[[s]]
      school_daily_df$sim <- i
      school_daily_df$school_name <- schools$school_name[s]
      school_daily_df$school_size <- schools$school_size[s]
      school_daily_df$vaccination_coverage <- schools$vaccination_coverage[s]
      all_school_daily[[length(all_school_daily) + 1]] <- school_daily_df
    }
    
    # Collect household results
    if (use_households && !is.null(result$hh_daily_counts)) {
      result$hh_daily_counts$sim_id <- i
      all_hh_daily[[i]] <- result$hh_daily_counts
      
      if (!is.null(result$hh_summary)) {
        result$hh_summary$sim_id <- i
        all_hh_summaries[[i]] <- result$hh_summary
      }
    }
    
    summary_stats$total_infected[i] <- result$total_infected
    summary_stats$total_breakthrough[i] <- result$total_breakthrough
    summary_stats$schools_affected[i] <- sum(result$school_summary$total_infected > 0)
    summary_stats$actual_days[i] <- result$actual_days
    summary_stats$hh_members_infected[i] <- result$total_hh_members_infected
    
    if (verbose && i %% 10 == 0) {
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      est_total <- elapsed / i * n_simulations
      cat(sprintf("Completed %d/%d (%.1f%%) - Est. total time: %.1f sec\n",
                  i, n_simulations, 100*i/n_simulations, est_total))
    }
  }
  
  all_network_data <- bind_rows(all_network_daily)
  all_school_data <- bind_rows(all_school_summaries)
  all_school_daily_data <- bind_rows(all_school_daily)
  
  # Combine household data
  all_hh_data <- if (use_households && length(all_hh_daily) > 0) {
    bind_rows(all_hh_daily)
  } else {
    NULL
  }
  
  all_hh_summary_data <- if (use_households && length(all_hh_summaries) > 0) {
    bind_rows(all_hh_summaries)
  } else {
    NULL
  }
  
  end_time <- Sys.time()
  total_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  if (verbose) {
    cat(sprintf("\n=== COMPLETED IN %.2f SECONDS ===\n", total_time))
    cat("\n=== SUMMARY STATISTICS ===\n")
    cat(sprintf("Students Infected (Median): %.1f (95%% CI: %.1f - %.1f)\n",
                median(summary_stats$total_infected),
                quantile(summary_stats$total_infected, 0.025),
                quantile(summary_stats$total_infected, 0.975)))
    cat(sprintf("Schools Affected (Median): %.1f (95%% CI: %.1f - %.1f)\n",
                median(summary_stats$schools_affected),
                quantile(summary_stats$schools_affected, 0.025),
                quantile(summary_stats$schools_affected, 0.975)))
    
    if (use_households) {
      cat(sprintf("Household Members Infected (Median): %.1f (95%% CI: %.1f - %.1f)\n",
                  median(summary_stats$hh_members_infected),
                  quantile(summary_stats$hh_members_infected, 0.025),
                  quantile(summary_stats$hh_members_infected, 0.975)))
    }
  }
  
  results <- list(
    all_network_data = all_network_data,
    all_school_data = all_school_data,
    all_school_daily_data = all_school_daily_data,
    summary_stats = summary_stats,
    params = params,
    schools = schools,
    network = network,
    n_simulations = n_simulations,
    n_days = n_days,
    computation_time = total_time,
    # Household results
    use_households = use_households,
    all_hh_data = all_hh_data,
    all_hh_summary_data = all_hh_summary_data
  )
  
  return(results)
}


cat("Simulation utilities loaded (with household transmission support).\n")
