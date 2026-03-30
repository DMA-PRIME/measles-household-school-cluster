# ==============================================================================
# POPULATION UTILITIES
# ==============================================================================
# File: population_utils.R
# Contains: Population creation, vaccination, infection seeding, disease states
# Updated: Added household fields for household transmission model
# Dependencies: R6 (for ContactHistory class)
# ==============================================================================

library(R6)

# ==============================================================================
# Erlang Distribution Helper
# ==============================================================================

#' Draw samples from an Erlang distribution
#' @param n Number of samples
#' @param mean Target mean
#' @param shape Shape parameter (integer)
#' @return Vector of integer samples (minimum 1)
draw_erlang <- function(n, mean, shape) {
  rate <- shape / mean
  samples <- rgamma(n, shape = shape, rate = rate)
  return(pmax(1, round(samples)))
}


# ==============================================================================
# Safe School Size Conversion
# ==============================================================================

#' Safely convert school size values that may contain characters like "<10"
#' @param x Vector of school sizes (numeric or character)
#' @param default_small Value to use for "<X" entries (default 5)
#' @return Numeric vector with NA for invalid entries
safe_school_size <- function(x, default_small = 5) {
  sapply(x, function(val) {
    if (is.na(val) || val == "" || val == "N/A" || val == "NA") {
      return(NA_real_)
    }
    if (is.numeric(val)) {
      return(val)
    }
    val_str <- as.character(val)
    if (grepl("^<", val_str)) {
      return(default_small)
    }
    if (grepl("^>", val_str)) {
      return(as.numeric(gsub(">", "", val_str)))
    }
    return(suppressWarnings(as.numeric(val_str)))
  })
}


# ==============================================================================
# Contact History R6 Class
# ==============================================================================

#' R6 class for tracking contact history
#' @description Stores infector-target pairs for quarantine implementation
ContactHistory <- R6::R6Class("ContactHistory",
  public = list(
    contact_list = NULL,
    window_size = NULL,
    
    initialize = function(window_size = 7) {
      self$window_size <- window_size
      self$contact_list <- list()
    },
    
    add_contacts = function(infector_ids, target_ids) {
      if (length(infector_ids) > 0) {
        self$contact_list[[length(self$contact_list) + 1]] <- list(
          infectors = infector_ids,
          targets = target_ids
        )
        
        if (length(self$contact_list) > self$window_size) {
          self$contact_list <- self$contact_list[(length(self$contact_list) - self$window_size + 1):length(self$contact_list)]
        }
      }
    },
    
    get_all_contacts = function() {
      all_infectors <- integer()
      all_targets <- integer()
      
      for (day_contacts in self$contact_list) {
        all_infectors <- c(all_infectors, day_contacts$infectors)
        all_targets <- c(all_targets, day_contacts$targets)
      }
      
      return(list(infector_ids = all_infectors, target_ids = all_targets))
    },
    
    clear = function() {
      self$contact_list <- list()
    }
  )
)


# ==============================================================================
# Create School Population
# ==============================================================================

#' Create population for a single school
#' @param school_id Unique identifier for the school
#' @param school_size Number of students
#' @param avg_class_size Average class size
#' @param age_range Vector of min and max age
#' @return Data frame with student population
create_school_population <- function(school_id, school_size, avg_class_size, age_range) {
  n_classes <- ceiling(school_size / avg_class_size)
  
  population <- data.frame(
    student_id = ((school_id - 1) * 10000) + (1:school_size),
    school_id = school_id,
    class_id = rep(1:n_classes, length.out = school_size),
    age = sample(age_range[1]:age_range[2], school_size, replace = TRUE)
  )
  
  # Disease state fields
  population$state <- "S"
  population$time_in_state <- 0
  population$time_since_prodromal <- NA
  population$is_vaccinated <- FALSE
  population$vaccine_failed <- FALSE   # TRUE = vaccine leaked (set at init, permanent)
  population$breakthrough_infection <- FALSE
  population$is_isolated <- FALSE
  population$is_quarantined <- FALSE
  population$is_index <- FALSE           # TRUE = first case detected in THIS school
  population$is_school_index <- FALSE    # TRUE = first case in this specific school
  population$newly_isolated <- FALSE
  population$infection_source <- NA_character_  # Track source: "seed", "within_school", "between_school", "household"
  
  # Disease duration fields
  population$latent_duration <- NA_real_
  population$infectious_duration <- NA_real_
  population$prodromal_duration <- NA_real_
  population$rash_duration <- NA_real_
  
  # Household fields (will be populated by assign_households_to_students)
  population$hh_id <- NA_integer_
  population$synpop_person_id <- NA_integer_
  population$hh_lon <- NA_real_
  population$hh_lat <- NA_real_
  
  return(population)
}


# ==============================================================================
# Initialize Vaccination
# ==============================================================================

#' Initialize vaccination for a school population
#' @param population School population data frame
#' @param vaccination_coverage Proportion vaccinated (0-1)
#' @param vaccine_efficacy Probability vaccine provides full protection (default 0.97)
#' @return Updated population data frame
initialize_vaccination_school <- function(population, vaccination_coverage, vaccine_efficacy = 0.97) {
  if (vaccination_coverage > 0) {
    n_vaccinated <- round(nrow(population) * vaccination_coverage)
    if (n_vaccinated > 0) {
      vaccinated_ids <- sample(1:nrow(population), n_vaccinated)
      population$is_vaccinated[vaccinated_ids] <- TRUE
      population$state[vaccinated_ids] <- "V"

      # Pre-assign vaccine failure: (1 - efficacy) fraction have leaky protection
      # These individuals CAN be infected (with reduced susceptibility via vaccine_reduction)
      # The remaining (efficacy) fraction are permanently immune and never enter susceptible pools
      n_vacc   <- length(vaccinated_ids)
      failures <- runif(n_vacc) >= vaccine_efficacy
      population$vaccine_failed[vaccinated_ids] <- failures
    }
  }
  return(population)
}


# ==============================================================================
# Seed Initial Infections
# ==============================================================================

#' Seed initial infections in specified schools
#' @param populations List of school populations
#' @param seed_schools Vector of school IDs to seed (must be indices 1:n_schools)
#' @param n_infected Number of initial infections per school
#' @param params Simulation parameters
#' @return Updated list of populations
seed_infections <- function(populations, seed_schools, n_infected, params) {
  for (school_idx in seed_schools) {
    pop <- populations[[school_idx]]
    
    susceptible_ids <- which(pop$state == "S")
    
    if (length(susceptible_ids) >= n_infected && n_infected > 0) {
      infected_ids <- sample(susceptible_ids, n_infected)
      index_id <- infected_ids[1]
      
      # Index case starts in P (prodromal)
      pop$state[index_id] <- "P"
      pop$time_in_state[index_id] <- 0
      pop$time_since_prodromal[index_id] <- 0
      pop$is_index[index_id] <- TRUE           # Original seeded case
      pop$is_school_index[index_id] <- TRUE    # Also first case in this school
      pop$infection_source[index_id] <- "seed"  # Track as seed case
      
      pop$latent_duration[index_id] <- 0
      infectious_dur <- draw_erlang(1, params$infectious_mean, params$infectious_shape)
      pop$infectious_duration[index_id] <- infectious_dur
      pop$prodromal_duration[index_id] <- params$prodromal_period
      pop$rash_duration[index_id] <- max(1, infectious_dur - params$prodromal_period)
      
      # Additional initial infections (if any)
      if (length(infected_ids) > 1) {
        others <- infected_ids[-1]
        pop$state[others] <- "E"
        pop$time_in_state[others] <- 0
        pop$infection_source[others] <- "seed"  # Additional seed cases
        
        n_others <- length(others)
        pop$latent_duration[others] <- draw_erlang(n_others, params$latent_mean, params$latent_shape)
        infectious_durs <- draw_erlang(n_others, params$infectious_mean, params$infectious_shape)
        pop$infectious_duration[others] <- infectious_durs
        pop$prodromal_duration[others] <- params$prodromal_period
        pop$rash_duration[others] <- pmax(1, infectious_durs - params$prodromal_period)
        # Note: these are NOT school index cases - they were infected along with the index
      }
    }
    
    populations[[school_idx]] <- pop
  }
  
  return(populations)
}


# ==============================================================================
# Update Disease States
# ==============================================================================

#' Update disease states for a population
#' @param population School population data frame
#' @param params Simulation parameters
#' @return Updated population data frame
update_disease_states <- function(population, params) {
  
  # E -> P transition
  e_to_p <- which(population$state == "E" & 
                  population$time_in_state >= population$latent_duration)
  if (length(e_to_p) > 0) {
    population$state[e_to_p] <- "P"
    population$time_in_state[e_to_p] <- 0
    population$time_since_prodromal[e_to_p] <- 0
  }
  
  # P -> Ra transition
  p_to_ra <- which(population$state == "P" & 
                   population$time_in_state >= population$prodromal_duration)
  if (length(p_to_ra) > 0) {
    # Check if school already has a detected case (someone with is_school_index = TRUE)
    school_has_index <- any(population$is_school_index)
    
    # If no index yet, the FIRST one transitioning to Ra becomes the school index
    if (!school_has_index && length(p_to_ra) > 0) {
      # Mark the first one as school index
      population$is_school_index[p_to_ra[1]] <- TRUE
    }
    
    population$state[p_to_ra] <- "Ra"
    population$time_in_state[p_to_ra] <- 0
  }
  
  # Ra -> Iso transition (with delay)
  # Index cases for each school: shorter delay (no surveillance in place)
  # Secondary cases: longer delay due to heightened surveillance after first case detected
  if (!params$no_intervention) {
    # School index cases (first detected case in each school)
    # Uses isolation_delay_index: days after rash onset
    ra_to_iso_index <- which(population$state == "Ra" & 
                             population$is_school_index & 
                             population$time_in_state >= params$isolation_delay_index)
    
    # Secondary cases (all cases after first in the school)
    # Uses isolation_delay_secondary: days after prodromal onset (heightened surveillance)
    ra_to_iso_secondary <- which(population$state == "Ra" & 
                                 !population$is_school_index & 
                                 !is.na(population$time_since_prodromal) &
                                 population$time_since_prodromal >= params$isolation_delay_secondary)
    
    ra_to_iso <- union(ra_to_iso_index, ra_to_iso_secondary)
    
    if (length(ra_to_iso) > 0) {
      population$state[ra_to_iso] <- "Iso"
      population$is_isolated[ra_to_iso] <- TRUE
      population$newly_isolated[ra_to_iso] <- TRUE
      population$time_in_state[ra_to_iso] <- 0
    }
  }
  
  # Iso -> R transition
  iso_to_r <- which(population$state == "Iso" & 
                    !is.na(population$time_since_prodromal) &
                    population$time_since_prodromal >= population$infectious_duration)
  if (length(iso_to_r) > 0) {
    population$state[iso_to_r] <- "R"
    population$is_isolated[iso_to_r] <- FALSE
    population$time_in_state[iso_to_r] <- 0
  }
  
  # Ra -> R transition (if no intervention)
  if (params$no_intervention) {
    ra_to_r <- which(population$state == "Ra" & 
                     !is.na(population$time_since_prodromal) &
                     population$time_since_prodromal >= population$infectious_duration)
    if (length(ra_to_r) > 0) {
      population$state[ra_to_r] <- "R"
      population$time_in_state[ra_to_r] <- 0
    }
  }
  
  # Quarantine transitions
  # QE -> QP
  qe_to_qp <- which(population$state == "QE" & 
                    population$time_in_state >= population$latent_duration)
  if (length(qe_to_qp) > 0) {
    population$state[qe_to_qp] <- "QP"
    population$time_in_state[qe_to_qp] <- 0
    population$time_since_prodromal[qe_to_qp] <- 0
  }
  
  # QP -> Iso
  qp_to_iso <- which(population$state == "QP" & 
                     population$time_in_state >= population$prodromal_duration)
  if (length(qp_to_iso) > 0) {
    population$state[qp_to_iso] <- "Iso"
    population$is_isolated[qp_to_iso] <- TRUE
    population$is_quarantined[qp_to_iso] <- FALSE
    population$time_in_state[qp_to_iso] <- 0
  }
  
  # QS/QV -> S/V (quarantine release)
  qs_release <- which(population$state == "QS" & 
                      population$time_in_state >= params$quarantine_duration)
  if (length(qs_release) > 0) {
    # Return to S or V depending on vaccination status
    for (idx in qs_release) {
      population$state[idx] <- ifelse(population$is_vaccinated[idx], "V", "S")
    }
    population$is_quarantined[qs_release] <- FALSE
    population$time_in_state[qs_release] <- 0
  }
  
  # QV -> V (vaccinated quarantine release - handle separately if needed)
  qv_ready <- which(population$state == "QV" & 
                    population$time_in_state >= params$quarantine_duration)
  if (length(qv_ready) > 0) {
    population$state[qv_ready] <- "V"
    population$is_quarantined[qv_ready] <- FALSE
    population$time_in_state[qv_ready] <- 0
  }
  
  # Increment time
  population$time_in_state <- population$time_in_state + 1
  
  infectious_states <- c("P", "Ra", "Iso", "QP")
  in_infectious <- which(population$state %in% infectious_states & 
                         !is.na(population$time_since_prodromal))
  if (length(in_infectious) > 0) {
    population$time_since_prodromal[in_infectious] <- 
      population$time_since_prodromal[in_infectious] + 1
  }
  
  population
}


# ==============================================================================
# Apply Quarantine
# ==============================================================================

#' Apply quarantine to contacts of isolated individuals
#' @param population School population data frame
#' @param params Simulation parameters
#' @param contact_history ContactHistory object
#' @return Updated population data frame
apply_quarantine <- function(population, params, contact_history) {
  if (params$quarantine_contacts == FALSE || params$quarantine_efficacy == 0) {
    return(population)
  }
  
  newly_isolated <- population$newly_isolated
  
  if (sum(newly_isolated) == 0) {
    return(population)
  }
  
  history_contacts <- contact_history$get_all_contacts()
  
  result <- cpp_apply_quarantine_with_history(
    student_id       = population$student_id,
    state            = population$state,
    is_quarantined   = population$is_quarantined,
    is_vaccinated    = population$is_vaccinated,
    newly_isolated   = newly_isolated,
    contact_history_infector = history_contacts$infector_ids,
    contact_history_target   = history_contacts$target_ids,
    quarantine_efficacy = params$quarantine_efficacy
  )
  
  if (length(result$quarantine_ids) > 0) {
    for (i in seq_along(result$quarantine_ids)) {
      idx <- which(population$student_id == result$quarantine_ids[i])
      if (length(idx) > 0) {
        population$state[idx] <- result$quarantine_states[i]
        population$is_quarantined[idx] <- TRUE
        population$time_in_state[idx] <- 0
      }
    }
  }
  
  population$newly_isolated <- FALSE
  
  return(population)
}

cat("Population utilities loaded (with household field support).\n")