# ==============================================================================
# MID-OUTBREAK SEEDING UTILITIES
# ==============================================================================
# File: midoutbreak_utils.R
# Purpose: Initialize simulations from a known mid-outbreak state rather than
#          from the very first case. Useful for real-time forecasting.
#
# Usage:
#   source("midoutbreak_utils.R")
#   result <- run_midoutbreak_simulation(schools, network, params,
#               observed_state, n_days, seed, hh_pop, household_assignment)
#
# The observed_state list describes what we know about the outbreak so far:
#   exposed_schools     - Character vector of school names that reported cases
#   total_cases         - Integer: cumulative confirmed cases to date
#   quarantine_contacts - Integer: contacts currently in quarantine (0 if none)
#   fraction_active     - Numeric 0-1: fraction of total_cases still infectious
#                         (default 0.15; most confirmed cases have recovered)
#   hh_attack_rate      - Numeric 0-1: secondary household attack rate applied
#                         to household members of infected students (default 0.3)
#   additional_vaccinations - Integer: number of susceptible students to vaccinate
#                         mid-outbreak as a reactive campaign (default 0)
#   hh_vax_fraction     - Numeric 0-1: fraction of susceptible HH members of
#                         newly vaccinated students to also vaccinate (default 0)
# ==============================================================================

# Null-coalescing operator (in case rlang is not loaded)
if (!exists("%||%")) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
}


# ==============================================================================
# Match school names to IDs
# ==============================================================================

resolve_school_ids <- function(exposed_schools, schools) {
  # Accept either names (character) or IDs (numeric)
  if (is.numeric(exposed_schools)) {
    ids <- exposed_schools[exposed_schools %in% schools$school_id]
    if (length(ids) == 0) stop("No valid school IDs found in exposed_schools")
    return(ids)
  }

  # Character matching: try exact, then fuzzy
  ids <- integer(0)
  unmatched <- character(0)

  for (name in exposed_schools) {
    # Exact match
    exact <- which(tolower(schools$school_name) == tolower(name))
    if (length(exact) == 1) {
      ids <- c(ids, schools$school_id[exact])
      next
    }

    # Partial match (school name contains the query or vice versa)
    partial <- which(grepl(tolower(name), tolower(schools$school_name), fixed = TRUE))
    if (length(partial) == 1) {
      ids <- c(ids, schools$school_id[partial])
      next
    }

    # Reverse partial (query contains school name)
    rev_partial <- which(sapply(tolower(schools$school_name),
                                function(sn) grepl(sn, tolower(name), fixed = TRUE)))
    if (length(rev_partial) == 1) {
      ids <- c(ids, schools$school_id[rev_partial])
      next
    }

    unmatched <- c(unmatched, name)
  }

  if (length(unmatched) > 0) {
    cat(sprintf("WARNING: Could not match %d school names:\n", length(unmatched)))
    for (nm in unmatched) cat(sprintf("  - '%s'\n", nm))
    cat("Available schools:\n")
    for (i in seq_len(min(20, nrow(schools)))) {
      cat(sprintf("  [%d] %s\n", schools$school_id[i], schools$school_name[i]))
    }
    if (nrow(schools) > 20) cat(sprintf("  ... and %d more\n", nrow(schools) - 20))
  }

  ids <- unique(ids)
  if (length(ids) == 0) stop("No exposed schools could be matched")
  return(ids)
}


# ==============================================================================
# Seed from observed mid-outbreak state
# ==============================================================================

#' Initialize populations to reflect an ongoing outbreak
#'
#' Distributes known cases and quarantine contacts across exposed schools.
#' Cases go to unvaccinated (susceptible) individuals. The split between
#' currently active (E/P/Ra) and recovered (R) is controlled by fraction_active.
#'
#' @param populations List of school population data frames
#' @param schools Data frame with school info
#' @param params Simulation parameters
#' @param observed_state List with: exposed_schools, total_cases,
#'   quarantine_contacts, fraction_active, hh_attack_rate
#' @param hh_pop Household population data frame (optional, can be NULL)
#' @return List with updated populations and hh_pop

seed_from_observed_state <- function(populations, schools, params,
                                     observed_state, hh_pop = NULL) {

  # --- Unpack observed state ---
  exposed_schools     <- observed_state$exposed_schools
  total_cases         <- observed_state$total_cases
  quarantine_contacts <- observed_state$quarantine_contacts %||% 0
  fraction_active     <- observed_state$fraction_active %||% 0.15
  hh_attack_rate      <- observed_state$hh_attack_rate %||% 0.30

  cat("\n=== SEEDING FROM OBSERVED MID-OUTBREAK STATE ===\n")

  # --- Resolve school names to IDs ---
  exposed_ids <- resolve_school_ids(exposed_schools, schools)
  cat(sprintf("Exposed schools matched: %d of %d provided\n",
              length(exposed_ids), length(exposed_schools)))
  for (sid in exposed_ids) {
    cat(sprintf("  [%d] %s (size: %d, vax: %.0f%%)\n",
                sid, schools$school_name[sid],
                schools$school_size[sid],
                schools$vaccination_coverage[sid] * 100))
  }

  # --- Calculate case distribution ---
  n_active    <- max(1, round(total_cases * fraction_active))
  n_recovered <- total_cases - n_active

  # Among active cases, split by disease stage duration ratios
  # latent(E) ~ 10d, prodromal(P) ~ 4d, rash(Ra) ~ 4d
  # But confirmed "cases" are typically past latent, so active = P + Ra
  # We also add a small E pool for very recent exposures
  active_E  <- max(0, round(n_active * 0.20))  # 20% recently exposed
  active_P  <- max(0, round(n_active * 0.35))  # 35% prodromal
  active_Ra <- n_active - active_E - active_P   # remainder rash

  cat(sprintf("\nTotal confirmed cases: %d\n", total_cases))
  cat(sprintf("  Currently active:  %d (E=%d, P=%d, Ra=%d)\n",
              n_active, active_E, active_P, active_Ra))
  cat(sprintf("  Already recovered: %d\n", n_recovered))
  cat(sprintf("  Quarantine contacts: %d\n", quarantine_contacts))

  # --- Determine per-school allocation (proportional to unvax pool size) ---
  # Get unvaccinated susceptible pool sizes per exposed school
  unvax_sizes <- sapply(exposed_ids, function(sid) {
    sum(populations[[sid]]$state == "S")
  })

  if (sum(unvax_sizes) == 0) {
    warning("No susceptible individuals found in exposed schools!")
    return(list(populations = populations, hh_pop = hh_pop))
  }

  # Proportional weights
  weights <- unvax_sizes / sum(unvax_sizes)

  # Allocate cases to schools (at least 1 per school if possible)
  allocate_to_schools <- function(n_total, weights) {
    n_schools <- length(weights)
    if (n_total <= 0) return(rep(0L, n_schools))
    if (n_total <= n_schools) {
      # Fewer cases than schools: assign 1 to random schools
      alloc <- rep(0L, n_schools)
      selected <- sample(n_schools, min(n_total, n_schools), prob = weights)
      alloc[selected] <- 1L
      return(alloc)
    }
    # Proportional allocation with rounding
    alloc <- round(n_total * weights)
    # Fix rounding errors
    diff <- n_total - sum(alloc)
    if (diff > 0) {
      add_to <- sample(n_schools, diff, prob = weights)
      for (idx in add_to) alloc[idx] <- alloc[idx] + 1L
    } else if (diff < 0) {
      remove_from <- sample(which(alloc > 0), abs(diff))
      for (idx in remove_from) alloc[idx] <- alloc[idx] - 1L
    }
    return(as.integer(alloc))
  }

  alloc_R   <- allocate_to_schools(n_recovered, weights)
  alloc_E   <- allocate_to_schools(active_E, weights)
  alloc_P   <- allocate_to_schools(active_P, weights)
  alloc_Ra  <- allocate_to_schools(active_Ra, weights)
  alloc_Q   <- allocate_to_schools(quarantine_contacts, weights)

  # --- Apply to each exposed school ---
  all_infected_hh_ids <- c()  # Track for household seeding

  for (i in seq_along(exposed_ids)) {
    sid <- exposed_ids[i]
    pop <- populations[[sid]]

    n_R  <- alloc_R[i]
    n_E  <- alloc_E[i]
    n_P  <- alloc_P[i]
    n_Ra <- alloc_Ra[i]
    n_Q  <- alloc_Q[i]

    # Pool of susceptible (unvaccinated) individuals
    susceptible_idx <- which(pop$state == "S")

    # Total cases for this school
    n_cases_school <- n_R + n_E + n_P + n_Ra

    if (length(susceptible_idx) < n_cases_school) {
      cat(sprintf("  WARNING: School %d has only %d susceptible but %d cases assigned; capping\n",
                  sid, length(susceptible_idx), n_cases_school))
      # Scale down proportionally
      scale <- length(susceptible_idx) / n_cases_school
      n_R  <- round(n_R * scale)
      n_E  <- round(n_E * scale)
      n_P  <- round(n_P * scale)
      n_Ra <- max(0, length(susceptible_idx) - n_R - n_E - n_P)
      n_cases_school <- n_R + n_E + n_P + n_Ra
    }

    if (n_cases_school == 0) next

    # Sample individuals for infection
    infected_idx <- sample(susceptible_idx, n_cases_school)

    cursor <- 0

    # --- Assign RECOVERED ---
    if (n_R > 0) {
      r_idx <- infected_idx[(cursor + 1):(cursor + n_R)]
      pop$state[r_idx] <- "R"
      pop$time_in_state[r_idx] <- sample(1:30, n_R, replace = TRUE)
      pop$infection_source[r_idx] <- "observed"
      # Recovered need completed durations (for bookkeeping)
      pop$latent_duration[r_idx] <- draw_erlang(n_R, params$latent_mean, params$latent_shape)
      inf_dur <- draw_erlang(n_R, params$infectious_mean, params$infectious_shape)
      pop$infectious_duration[r_idx] <- inf_dur
      pop$prodromal_duration[r_idx] <- params$prodromal_period
      pop$rash_duration[r_idx] <- pmax(1, inf_dur - params$prodromal_period)
      pop$time_since_prodromal[r_idx] <- inf_dur + 1  # past infectious period

      # Mark a recovered individual as the school index case.
      # In a mid-outbreak, the original index case has already been detected
      # and recovered. Setting is_school_index on a recovered individual
      # signals that this school has heightened surveillance, so ALL currently
      # active cases (P, Ra) will use the faster secondary isolation delay
      # rather than the slower index delay.
      pop$is_school_index[r_idx[1]] <- TRUE
      pop$is_index[r_idx[1]] <- TRUE

      cursor <- cursor + n_R
    }

    # --- Assign EXPOSED (latent) ---
    if (n_E > 0) {
      e_idx <- infected_idx[(cursor + 1):(cursor + n_E)]
      pop$state[e_idx] <- "E"
      lat_dur <- draw_erlang(n_E, params$latent_mean, params$latent_shape)
      pop$latent_duration[e_idx] <- lat_dur
      # Partially through latent period
      pop$time_in_state[e_idx] <- pmin(
        sample(1:max(1, round(params$latent_mean * 0.8)), n_E, replace = TRUE),
        lat_dur - 1
      )
      pop$infection_source[e_idx] <- "observed"
      inf_dur <- draw_erlang(n_E, params$infectious_mean, params$infectious_shape)
      pop$infectious_duration[e_idx] <- inf_dur
      pop$prodromal_duration[e_idx] <- params$prodromal_period
      pop$rash_duration[e_idx] <- pmax(1, inf_dur - params$prodromal_period)
      cursor <- cursor + n_E
    }

    # --- Assign PRODROMAL ---
    if (n_P > 0) {
      p_idx <- infected_idx[(cursor + 1):(cursor + n_P)]
      pop$state[p_idx] <- "P"
      pop$latent_duration[p_idx] <- 0  # Already past latent
      inf_dur <- draw_erlang(n_P, params$infectious_mean, params$infectious_shape)
      pop$infectious_duration[p_idx] <- inf_dur
      pop$prodromal_duration[p_idx] <- params$prodromal_period
      pop$rash_duration[p_idx] <- pmax(1, inf_dur - params$prodromal_period)
      # Partially through prodromal
      pop$time_in_state[p_idx] <- sample(0:(params$prodromal_period - 1), n_P, replace = TRUE)
      pop$time_since_prodromal[p_idx] <- pop$time_in_state[p_idx]
      # NOTE: is_school_index is NOT set on active cases in mid-outbreak.
      # The index case is already recovered. All active cases use secondary
      # isolation delay (faster detection due to heightened surveillance).
      pop$infection_source[p_idx] <- "observed"
      cursor <- cursor + n_P
    }

    # --- Assign RASH ---
    if (n_Ra > 0) {
      ra_idx <- infected_idx[(cursor + 1):(cursor + n_Ra)]
      pop$state[ra_idx] <- "Ra"
      pop$latent_duration[ra_idx] <- 0
      inf_dur <- draw_erlang(n_Ra, params$infectious_mean, params$infectious_shape)
      pop$infectious_duration[ra_idx] <- inf_dur
      pop$prodromal_duration[ra_idx] <- params$prodromal_period
      rash_dur <- pmax(1, inf_dur - params$prodromal_period)
      pop$rash_duration[ra_idx] <- rash_dur
      # Partially through rash period
      pop$time_in_state[ra_idx] <- pmin(
        sample(0:3, n_Ra, replace = TRUE),
        rash_dur - 1
      )
      pop$time_since_prodromal[ra_idx] <- params$prodromal_period + pop$time_in_state[ra_idx]
      # NOTE: is_school_index NOT set here — see recovered section above.
      pop$infection_source[ra_idx] <- "observed"
      cursor <- cursor + n_Ra
    }

    # --- Assign QUARANTINE CONTACTS ---
    # Quarantine uninfected individuals in this school (S -> QS, V -> QV)
    if (n_Q > 0) {
      remaining_S <- which(pop$state == "S")
      remaining_V <- which(pop$state == "V")
      available_Q <- c(remaining_S, remaining_V)

      if (length(available_Q) > 0) {
        n_Q_actual <- min(n_Q, length(available_Q))
        q_idx <- sample(available_Q, n_Q_actual)

        for (idx in q_idx) {
          pop$state[idx] <- ifelse(pop$is_vaccinated[idx], "QV", "QS")
          pop$is_quarantined[idx] <- TRUE
          # Partially through quarantine (assume placed recently)
          pop$time_in_state[idx] <- sample(0:7, 1)
        }
      }
    }

    # Collect household IDs of infected students (for household seeding)
    all_case_idx <- infected_idx[1:n_cases_school]
    hh_ids <- pop$hh_id[all_case_idx]
    all_infected_hh_ids <- c(all_infected_hh_ids, hh_ids[!is.na(hh_ids)])

    # Guard: if n_R was 0 (all cases still active), no recovered individual
    # carries is_school_index. In a confirmed mid-outbreak, surveillance IS
    # heightened, so we need is_school_index set somewhere. Set it on an E
    # case if available (E cases don't trigger isolation, so no delay penalty).
    # Otherwise fall back to the first active case.
    if (!any(pop$is_school_index) && n_cases_school > 0) {
      e_cases <- which(pop$state == "E" & pop$infection_source == "observed")
      if (length(e_cases) > 0) {
        pop$is_school_index[e_cases[1]] <- TRUE
      } else {
        # Last resort: flag first infected — this one case uses index delay
        pop$is_school_index[infected_idx[1]] <- TRUE
      }
    }

    populations[[sid]] <- pop

    cat(sprintf("  School %d (%s): R=%d, E=%d, P=%d, Ra=%d, Q=%d\n",
                sid, schools$school_name[sid], n_R, n_E, n_P, n_Ra,
                min(n_Q, length(c(which(pop$state == "QS"), which(pop$state == "QV"))))))
  }

  # --- Seed household members ---
  if (!is.null(hh_pop) && length(all_infected_hh_ids) > 0) {
    cat(sprintf("\nSeeding household members (attack rate: %.0f%%)...\n",
                hh_attack_rate * 100))

    # Find non-student, susceptible HH members in affected households
    hh_at_risk <- which(
      hh_pop$hh_id %in% all_infected_hh_ids &
      !hh_pop$is_student &
      hh_pop$state == "S"
    )

    if (length(hh_at_risk) > 0) {
      # Each at-risk member infected with probability hh_attack_rate
      infected_mask <- runif(length(hh_at_risk)) < hh_attack_rate
      hh_infected_idx <- hh_at_risk[infected_mask]
      n_hh_inf <- length(hh_infected_idx)

      if (n_hh_inf > 0) {
        # Most household infections are recovered by now, some active
        n_hh_R  <- round(n_hh_inf * (1 - fraction_active))
        n_hh_active <- n_hh_inf - n_hh_R

        # Recovered HH members
        if (n_hh_R > 0) {
          r_idx <- hh_infected_idx[1:n_hh_R]
          hh_pop$state[r_idx] <- "R"
          hh_pop$time_in_state[r_idx] <- sample(1:20, n_hh_R, replace = TRUE)
          inf_dur <- draw_erlang(n_hh_R, params$infectious_mean, params$infectious_shape)
          hh_pop$infectious_duration[r_idx] <- inf_dur
          hh_pop$latent_duration[r_idx] <- draw_erlang(n_hh_R, params$latent_mean, params$latent_shape)
          hh_pop$prodromal_duration[r_idx] <- params$prodromal_period
          hh_pop$rash_duration[r_idx] <- pmax(1, inf_dur - params$prodromal_period)
          hh_pop$time_since_prodromal[r_idx] <- inf_dur + 1
        }

        # Active HH members (E or P)
        if (n_hh_active > 0) {
          a_idx <- hh_infected_idx[(n_hh_R + 1):n_hh_inf]
          n_hh_E <- round(n_hh_active * 0.5)
          n_hh_P <- n_hh_active - n_hh_E

          if (n_hh_E > 0) {
            e_sel <- a_idx[1:n_hh_E]
            hh_pop$state[e_sel] <- "E"
            lat <- draw_erlang(n_hh_E, params$latent_mean, params$latent_shape)
            hh_pop$latent_duration[e_sel] <- lat
            hh_pop$time_in_state[e_sel] <- pmin(
              sample(1:max(1, round(params$latent_mean * 0.5)), n_hh_E, replace = TRUE),
              lat - 1
            )
            inf_dur <- draw_erlang(n_hh_E, params$infectious_mean, params$infectious_shape)
            hh_pop$infectious_duration[e_sel] <- inf_dur
            hh_pop$prodromal_duration[e_sel] <- params$prodromal_period
            hh_pop$rash_duration[e_sel] <- pmax(1, inf_dur - params$prodromal_period)
          }
          if (n_hh_P > 0) {
            p_sel <- a_idx[(n_hh_E + 1):n_hh_active]
            hh_pop$state[p_sel] <- "P"
            hh_pop$latent_duration[p_sel] <- 0
            inf_dur <- draw_erlang(n_hh_P, params$infectious_mean, params$infectious_shape)
            hh_pop$infectious_duration[p_sel] <- inf_dur
            hh_pop$prodromal_duration[p_sel] <- params$prodromal_period
            hh_pop$rash_duration[p_sel] <- pmax(1, inf_dur - params$prodromal_period)
            hh_pop$time_in_state[p_sel] <- sample(0:(params$prodromal_period - 1),
                                                   n_hh_P, replace = TRUE)
            hh_pop$time_since_prodromal[p_sel] <- hh_pop$time_in_state[p_sel]
          }
        }

        cat(sprintf("  Household members infected: %d (R=%d, active=%d)\n",
                    n_hh_inf, n_hh_R, n_hh_active))
      }
    }

    # Also quarantine some household members of quarantined students
    # (if quarantine_contacts > 0, some HH members may also be in quarantine)
    hh_at_risk_q <- which(
      hh_pop$hh_id %in% all_infected_hh_ids &
      !hh_pop$is_student &
      hh_pop$state == "S"
    )
    n_hh_q <- min(round(quarantine_contacts * 0.3), length(hh_at_risk_q))
    if (n_hh_q > 0) {
      q_sel <- sample(hh_at_risk_q, n_hh_q)
      hh_pop$state[q_sel] <- "QS"
      hh_pop$is_quarantined[q_sel] <- TRUE
      hh_pop$time_in_state[q_sel] <- sample(0:7, n_hh_q, replace = TRUE)
      cat(sprintf("  Household members quarantined: %d\n", n_hh_q))
    }
  }

  cat("\n=== MID-OUTBREAK SEEDING COMPLETE ===\n\n")

  return(list(populations = populations, hh_pop = hh_pop))
}


# ==============================================================================
# Apply Additional Mid-Outbreak Vaccination
# ==============================================================================

#' Vaccinate additional susceptible students mid-outbreak
#'
#' After seeding the observed state, this function vaccinates n_additional
#' currently-susceptible (state="S") students across all schools. This models
#' a reactive vaccination campaign.
#'
#' Vaccination is distributed proportionally to each school's remaining
#' susceptible pool. Vaccine failure is applied using the same efficacy
#' model as initial vaccination.
#'
#' @param populations List of school population data frames
#' @param params Simulation parameters (needs vaccine_efficacy)
#' @param n_additional Integer: number of additional students to vaccinate
#' @param hh_pop Household population (optional, for vaccinating HH members too)
#' @param hh_vax_fraction Fraction of susceptible HH members of newly vaccinated
#'   students to also vaccinate (default 0, set >0 for ring vaccination)
#' @return List with updated populations and hh_pop

apply_midoutbreak_vaccination <- function(populations, params, n_additional,
                                           hh_pop = NULL, hh_vax_fraction = 0) {

  if (n_additional <= 0) {
    return(list(populations = populations, hh_pop = hh_pop, n_vaccinated = 0L))
  }

  vaccine_efficacy <- params$vaccine_efficacy %||% 0.97

  cat(sprintf("\n=== APPLYING MID-OUTBREAK VACCINATION ===\n"))
  cat(sprintf("Target: %d additional students\n", n_additional))

  # Count susceptible students per school
  susceptible_per_school <- sapply(populations, function(p) sum(p$state == "S"))
  total_susceptible <- sum(susceptible_per_school)

  if (total_susceptible == 0) {
    cat("WARNING: No susceptible students remaining to vaccinate\n")
    return(list(populations = populations, hh_pop = hh_pop, n_vaccinated = 0L))
  }

  n_to_vax <- min(n_additional, total_susceptible)
  if (n_to_vax < n_additional) {
    cat(sprintf("NOTE: Only %d susceptible students available (requested %d)\n",
                total_susceptible, n_additional))
  }

  # Distribute proportionally to susceptible pool size
  weights <- susceptible_per_school / total_susceptible
  alloc <- round(n_to_vax * weights)

  # Fix rounding
  diff <- n_to_vax - sum(alloc)
  if (diff > 0) {
    candidates <- which(susceptible_per_school > alloc)
    if (length(candidates) > 0) {
      add_to <- sample(candidates, min(diff, length(candidates)), prob = weights[candidates])
      for (idx in add_to) alloc[idx] <- alloc[idx] + 1L
    }
  } else if (diff < 0) {
    candidates <- which(alloc > 0)
    remove_from <- sample(candidates, min(abs(diff), length(candidates)))
    for (idx in remove_from) alloc[idx] <- max(0L, alloc[idx] - 1L)
  }

  # Apply vaccination to each school
  newly_vaccinated_hh_ids <- c()
  total_vaccinated <- 0L

  for (s in seq_along(populations)) {
    n_vax_school <- alloc[s]
    if (n_vax_school <= 0) next

    pop <- populations[[s]]
    s_idx <- which(pop$state == "S")

    n_vax_actual <- min(n_vax_school, length(s_idx))
    if (n_vax_actual == 0) next

    vax_idx <- sample(s_idx, n_vax_actual)

    # Move S -> V
    pop$state[vax_idx] <- "V"
    pop$is_vaccinated[vax_idx] <- TRUE

    # Apply vaccine failure (leaky protection) — same as initialize_vaccination_school
    failures <- runif(n_vax_actual) >= vaccine_efficacy
    pop$vaccine_failed[vax_idx] <- failures

    # Collect household IDs for household ring vaccination
    hh_ids <- pop$hh_id[vax_idx]
    newly_vaccinated_hh_ids <- c(newly_vaccinated_hh_ids, hh_ids[!is.na(hh_ids)])

    populations[[s]] <- pop
    total_vaccinated <- total_vaccinated + n_vax_actual
  }

  cat(sprintf("Students vaccinated: %d across %d schools\n",
              total_vaccinated, sum(alloc > 0)))
  cat(sprintf("Remaining susceptible: %d\n", total_susceptible - total_vaccinated))

  # Optionally vaccinate household members of newly vaccinated students
  if (!is.null(hh_pop) && hh_vax_fraction > 0 && length(newly_vaccinated_hh_ids) > 0) {
    hh_at_risk <- which(
      hh_pop$hh_id %in% newly_vaccinated_hh_ids &
      !hh_pop$is_student &
      hh_pop$state == "S"
    )

    if (length(hh_at_risk) > 0) {
      vax_mask <- runif(length(hh_at_risk)) < hh_vax_fraction
      hh_vax_idx <- hh_at_risk[vax_mask]

      if (length(hh_vax_idx) > 0) {
        hh_pop$state[hh_vax_idx] <- "V"
        hh_pop$is_vaccinated[hh_vax_idx] <- TRUE
        failures <- runif(length(hh_vax_idx)) >= vaccine_efficacy
        hh_pop$breakthrough_infection[hh_vax_idx] <- FALSE
        cat(sprintf("Household members vaccinated: %d\n", length(hh_vax_idx)))
      }
    }
  }

  cat("=== MID-OUTBREAK VACCINATION COMPLETE ===\n\n")

  return(list(populations = populations, hh_pop = hh_pop, n_vaccinated = total_vaccinated))
}


# ==============================================================================
# Modified simulation runner for mid-outbreak scenarios
# ==============================================================================
#
# This is a variant of run_network_simulation that:
#   1. Creates populations and vaccination as normal
#   2. Applies mid-outbreak seeding instead of seed_infections
#   2b. Applies additional mid-outbreak vaccination (if specified)
#   3. Runs the same simulation loop
#
# It sources the simulation loop from simulation_utils.R functions.
# ==============================================================================

run_midoutbreak_simulation <- function(schools, network, params,
                                        observed_state,
                                        n_days = 150, seed = NULL,
                                        hh_pop = NULL,
                                        household_assignment = NULL) {

  if (!is.null(seed)) set.seed(seed)

  n_schools <- nrow(schools)
  use_households <- !is.null(hh_pop)

  # --- Step 1: Initialize populations (same as run_network_simulation) ---
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

  # Assign household IDs
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

  # --- Step 2: Apply mid-outbreak seeding (INSTEAD of seed_infections) ---
  seed_result <- seed_from_observed_state(populations, schools, params,
                                           observed_state, hh_pop)
  populations <- seed_result$populations
  hh_pop <- seed_result$hh_pop

  # --- Step 2b: Apply additional mid-outbreak vaccination (if specified) ---
  n_additional_vax <- observed_state$additional_vaccinations %||% 0
  hh_vax_fraction  <- observed_state$hh_vax_fraction %||% 0

  if (n_additional_vax > 0) {
    vax_result <- apply_midoutbreak_vaccination(
      populations, params, n_additional_vax,
      hh_pop = hh_pop, hh_vax_fraction = hh_vax_fraction
    )
    populations <- vax_result$populations
    hh_pop <- vax_result$hh_pop
  }

  # Sync student states to household population
  if (use_households) {
    hh_pop <- sync_student_states(populations, hh_pop, direction = "school_to_hh")
  }

  # --- Step 3: Simulation loop (identical to run_network_simulation) ---
  state_cols <- c("S", "V", "E", "P", "Ra", "Iso", "R", "QS", "QV", "QE", "QP")

  daily_counts_list <- lapply(1:n_schools, function(s) {
    matrix(0, nrow = n_days, ncol = length(state_cols),
           dimnames = list(NULL, state_cols))
  })

  network_daily_counts <- matrix(0, nrow = n_days, ncol = length(state_cols),
                                  dimnames = list(NULL, state_cols))

  if (use_households) {
    hh_daily_counts <- matrix(0, nrow = n_days, ncol = length(state_cols),
                               dimnames = list(NULL, state_cols))
  }

  outbreak_ended <- FALSE
  actual_days <- n_days

  for (day in 1:n_days) {

    if (!outbreak_ended) {
      # Within-school transmission
      for (school_idx in 1:n_schools) {
        trans_result <- school_transmission(
          populations[[school_idx]], params, contact_histories[[school_idx]]
        )
        populations[[school_idx]] <- trans_result$population
        contact_histories[[school_idx]] <- trans_result$contact_history
      }

      # Between-school transmission
      populations <- between_school_transmission(populations, network, params)

      # Household transmission
      if (use_households) {
        hh_pop <- sync_student_states(populations, hh_pop, direction = "school_to_hh")
        hh_result <- household_transmission(populations, hh_pop, params)
        populations <- hh_result$populations
        hh_pop <- hh_result$hh_pop
        hh_pop <- update_household_disease_states(hh_pop, params)
      }

      # Update disease states
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
      network_daily_counts[day, ] <- network_daily_counts[day, ] +
        daily_counts_list[[school_idx]][day, ]
    }

    if (use_households) {
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
      if (use_households) {
        total_active <- total_active + sum(hh_daily_counts[day, active_states])
      }
      if (total_active == 0) {
        actual_days <- day
        outbreak_ended <- TRUE
      }
    }
  }

  # --- Step 4: Assemble results (same structure as run_network_simulation) ---
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
    was_seeded = 1:n_schools %in% resolve_school_ids(observed_state$exposed_schools, schools)
  )
  school_summary$attack_rate <- school_summary$total_infected / school_summary$school_size

  if (use_households) {
    hh_summary <- summarize_household_infections(populations, hh_pop)
    total_hh_members_infected <- sum(hh_pop$state %in% c("P", "Ra", "R", "E") & !hh_pop$is_student)
    total_hh_members <- sum(!hh_pop$is_student)
  } else {
    hh_summary <- NULL
    total_hh_members_infected <- 0
    total_hh_members <- 0
  }

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
    use_households = use_households,
    hh_summary = hh_summary,
    total_hh_members_infected = total_hh_members_infected,
    total_hh_members = total_hh_members,
    hh_pop = if (use_households) hh_pop else NULL,
    # Extra: track what was seeded vs what was new
    observed_cases = observed_state$total_cases,
    observed_quarantine = observed_state$quarantine_contacts,
    additional_vaccinations = n_additional_vax
  )

  if (use_households) {
    results$hh_daily_counts <- as.data.frame(hh_daily_counts) %>%
      mutate(day = 0:(n_days - 1))
  }

  return(results)
}

cat("Mid-outbreak utilities loaded.\n")
cat("  - resolve_school_ids(): Match school names to IDs\n")
cat("  - seed_from_observed_state(): Initialize mid-outbreak state\n")
cat("  - apply_midoutbreak_vaccination(): Vaccinate additional students\n")
cat("  - run_midoutbreak_simulation(): Run simulation from observed state\n")
