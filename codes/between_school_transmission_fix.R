# ==============================================================================
# Between-School Transmission (density-dependent contact sampling)
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
