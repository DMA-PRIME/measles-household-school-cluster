# ==============================================================================
# PARALLEL SIMULATION UTILITIES
# ==============================================================================
# File: simulation_parallel.R
# Contains: Functions for running simulations in parallel on local machine
# Dependencies: parallel package (built into R), simulation_utils.R
# ==============================================================================

# ==============================================================================
# Setup and Configuration
# ==============================================================================

#' Detect available cores and recommend parallel settings
#' @param reserve_cores Number of cores to reserve for system (default 1)
#' @return List with system info and recommendations
detect_parallel_settings <- function(reserve_cores = 1) {
  
  n_cores <- parallel::detectCores(logical = TRUE)
  n_physical <- parallel::detectCores(logical = FALSE)
  
  os_type <- .Platform$OS.type
  is_windows <- os_type == "windows"
  
  recommended_cores <- max(1, n_physical - reserve_cores)
  
  cat("=== PARALLEL SETTINGS ===\n")
  cat(sprintf("Operating System: %s\n", Sys.info()["sysname"]))
  cat(sprintf("Logical cores: %d\n", n_cores))
  cat(sprintf("Physical cores: %d\n", n_physical))
  cat(sprintf("Recommended workers: %d (reserving %d for system)\n", 
              recommended_cores, reserve_cores))
  cat(sprintf("Cluster type: %s\n", ifelse(is_windows, "PSOCK", "FORK")))
  
  if (is_windows) {
    cat("\nNOTE: On Windows, parallel mode requires re-sourcing code on each worker.\n")
    cat("For many simulations (>50), parallel is faster despite the overhead.\n")
    cat("For fewer simulations, consider using run_fast_simulations() instead.\n")
  }
  
  return(list(
    n_cores = n_cores,
    n_physical = n_physical,
    recommended_workers = recommended_cores,
    is_windows = is_windows,
    cluster_type = ifelse(is_windows, "PSOCK", "FORK")
  ))
}


# ==============================================================================
# FAST Sequential Simulation (Optimized for single-core)
# ==============================================================================

#' Run simulations sequentially with maximum optimization
#' 
#' This is often FASTER than parallel on Windows for small/medium numbers of 
#' simulations (<50) because it avoids the overhead of creating PSOCK workers
#' and re-compiling C++ code on each worker.
#' 
#' @param n_simulations Number of simulations to run
#' @param schools Schools dataframe
#' @param network Network object
#' @param params Simulation parameters
#' @param seed_schools Vector of seed school IDs
#' @param n_initial_infected Number of initial infections per seed school
#' @param n_days Maximum simulation days
#' @param seed_start Starting seed for reproducibility
#' @param verbose Print progress messages
#' @param hh_pop Household population (optional)
#' @param household_assignment Household assignment dataframe (optional)
#' @param progress_every Print progress every N simulations (default 10)
#' @return List with all simulation results
#' @export
run_fast_simulations <- function(n_simulations,
                                  schools,
                                  network,
                                  params,
                                  seed_schools,
                                  n_initial_infected = 1,
                                  n_days = 150,
                                  seed_start = NULL,
                                  verbose = TRUE,
                                  hh_pop = NULL,
                                  household_assignment = NULL,
                                  progress_every = 10) {
  
  # Check that required function exists in global environment
  if (!exists("run_network_simulation", where = globalenv(), mode = "function")) {
    stop("run_network_simulation not found. Make sure to source('load_all.R') first.")
  }
  
  use_households <- !is.null(hh_pop)
  
  if (verbose) {
    cat("\n=== RUNNING FAST SEQUENTIAL SIMULATIONS ===\n")
    cat(sprintf("Number of simulations: %d\n", n_simulations))
    cat(sprintf("Schools: %d\n", nrow(schools)))
    cat(sprintf("Household transmission: %s\n", ifelse(use_households, "ENABLED", "DISABLED")))
    cat("\n")
  }
  
  start_time <- Sys.time()
  
  # Generate seeds
  sim_seeds <- if (!is.null(seed_start)) {
    seed_start + (1:n_simulations)
  } else {
    sample.int(.Machine$integer.max, n_simulations)
  }
  
  # Pre-allocate result list
  results_list <- vector("list", n_simulations)
  
  if (verbose) cat("Running simulations: ")
  
  for (i in 1:n_simulations) {
    # Reset household state for each simulation
    if (use_households) {
      hh_pop_sim <- hh_pop
      hh_pop_sim$state <- ifelse(hh_pop_sim$is_vaccinated, "V", "S")
      hh_pop_sim$time_in_state <- 0
      hh_pop_sim$time_since_prodromal <- NA_integer_
      hh_pop_sim$is_isolated <- FALSE
      hh_pop_sim$is_quarantined <- FALSE
      hh_pop_sim$breakthrough_infection <- FALSE
    } else {
      hh_pop_sim <- NULL
    }
    
    # Run simulation
    result <- run_network_simulation(
      schools = schools,
      network = network,
      params = params,
      seed_schools = seed_schools,
      n_initial_infected = n_initial_infected,
      n_days = n_days,
      seed = sim_seeds[i],
      hh_pop = hh_pop_sim,
      household_assignment = household_assignment
    )
    
    # Add simulation ID
    result$sim_id <- i
    result$network_daily_counts$sim_id <- i
    result$school_summary$sim_id <- i
    result$school_daily_combined$sim_id <- i
    
    if (use_households && !is.null(result$hh_daily_counts)) {
      result$hh_daily_counts$sim_id <- i
    }
    if (use_households && !is.null(result$hh_summary)) {
      result$hh_summary$sim_id <- i
    }
    
    # Create summary row
    result$summary_row <- data.frame(
      sim_id = i,
      total_infected = result$total_infections,
      schools_affected = result$schools_affected,
      peak_infected = result$peak_infections,
      actual_days = result$actual_days,
      hh_members_infected = ifelse(use_households, result$total_hh_members_infected, 0)
    )
    
    results_list[[i]] <- result
    
    if (verbose && i %% progress_every == 0) {
      cat(sprintf("%d ", i))
    }
  }
  
  if (verbose) cat("Done!\n")
  
  # Combine results
  if (verbose) cat("Combining results...\n")
  
  summary_stats <- do.call(rbind, lapply(results_list, function(r) r$summary_row))
  all_network_data <- do.call(rbind, lapply(results_list, function(r) r$network_daily_counts))
  all_school_data <- do.call(rbind, lapply(results_list, function(r) r$school_summary))
  all_school_daily_data <- do.call(rbind, lapply(results_list, function(r) r$school_daily_combined))
  
  all_hh_data <- NULL
  all_hh_summary_data <- NULL
  
  if (use_households) {
    hh_daily_list <- lapply(results_list, function(r) r$hh_daily_counts)
    hh_daily_list <- hh_daily_list[!sapply(hh_daily_list, is.null)]
    if (length(hh_daily_list) > 0) {
      all_hh_data <- do.call(rbind, hh_daily_list)
    }
    
    hh_summary_list <- lapply(results_list, function(r) r$hh_summary)
    hh_summary_list <- hh_summary_list[!sapply(hh_summary_list, is.null)]
    if (length(hh_summary_list) > 0) {
      all_hh_summary_data <- do.call(rbind, hh_summary_list)
    }
  }
  
  end_time <- Sys.time()
  total_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  if (verbose) {
    cat(sprintf("\n=== COMPLETED IN %.2f SECONDS ===\n", total_time))
    cat(sprintf("Average time per simulation: %.2f seconds\n", total_time / n_simulations))
    cat("\n=== SUMMARY STATISTICS ===\n")
    cat(sprintf("Students Infected (Median): %.1f (95%% CI: %.1f - %.1f)\n",
                median(summary_stats$total_infected),
                quantile(summary_stats$total_infected, 0.025),
                quantile(summary_stats$total_infected, 0.975)))
    if (use_households) {
      cat(sprintf("HH Members Infected (Median): %.1f (95%% CI: %.1f - %.1f)\n",
                  median(summary_stats$hh_members_infected),
                  quantile(summary_stats$hh_members_infected, 0.025),
                  quantile(summary_stats$hh_members_infected, 0.975)))
    }
  }
  
  return(list(
    summary_stats = summary_stats,
    all_network_data = all_network_data,
    all_school_data = all_school_data,
    all_school_daily_data = all_school_daily_data,
    all_hh_data = all_hh_data,
    all_hh_summary_data = all_hh_summary_data,
    n_simulations = n_simulations,
    computation_time = total_time,
    use_households = use_households,
    params = params
  ))
}


# ==============================================================================
# Single Simulation Worker Function
# ==============================================================================

#' Run a single simulation (worker function for parallel execution)
#' This is a self-contained function that can be sent to worker nodes
#' 
#' @param sim_id Simulation ID
#' @param schools Schools dataframe
#' @param network Network object
#' @param params Simulation parameters
#' @param seed_schools Vector of seed school IDs
#' @param n_initial_infected Number of initial infections per seed school
#' @param n_days Maximum simulation days
#' @param seed Random seed for this simulation
#' @param hh_pop Household population (optional)
#' @param household_assignment Household assignment dataframe (optional)
#' @return List with simulation results
run_single_simulation_worker <- function(sim_id, schools, network, params,
                                          seed_schools, n_initial_infected,
                                          n_days, seed,
                                          hh_pop = NULL, 
                                          household_assignment = NULL) {
  
  # Set seed for reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  use_households <- !is.null(hh_pop)
  
  # Reset household population state if using households
  if (use_households) {
    hh_pop$state <- ifelse(hh_pop$is_vaccinated, "V", "S")
    hh_pop$time_in_state <- 0
    hh_pop$time_since_prodromal <- NA_integer_
    hh_pop$is_isolated <- FALSE
    hh_pop$is_quarantined <- FALSE
    hh_pop$breakthrough_infection <- FALSE
  }
  
  # Run the simulation
  result <- run_network_simulation(
    schools = schools,
    network = network,
    params = params,
    seed_schools = seed_schools,
    n_initial_infected = n_initial_infected,
    n_days = n_days,
    seed = seed,
    hh_pop = hh_pop,
    household_assignment = household_assignment
  )
  
  # Add simulation ID to results
  result$sim_id <- sim_id
  result$network_daily_counts$sim_id <- sim_id
  result$school_summary$sim_id <- sim_id
  
  # Process school daily data
  school_daily_combined <- lapply(seq_along(result$school_daily_counts), function(s) {
    df <- result$school_daily_counts[[s]]
    df$sim_id <- sim_id
    df$school_id <- s
    df$school_name <- schools$school_name[s]
    df$school_size <- schools$school_size[s]
    df$vaccination_coverage <- schools$vaccination_coverage[s]
    df
  })
  result$school_daily_combined <- do.call(rbind, school_daily_combined)
  
  # Process household data
  if (use_households && !is.null(result$hh_daily_counts)) {
    result$hh_daily_counts$sim_id <- sim_id
    if (!is.null(result$hh_summary)) {
      result$hh_summary$sim_id <- sim_id
    }
  }
  
  # Create summary row
  result$summary_row <- data.frame(
    sim_id = sim_id,
    total_infected = result$total_infected,
    total_breakthrough = result$total_breakthrough,
    schools_affected = sum(result$school_summary$total_infected > 0),
    actual_days = result$actual_days,
    hh_members_infected = ifelse(use_households, result$total_hh_members_infected, 0)
  )
  
  return(result)
}


# ==============================================================================
# Main Parallel Simulation Function
# ==============================================================================

#' Run multiple simulations in parallel (OPTIMIZED)
#' 
#' This version avoids re-compiling C++ code on each worker by exporting
#' pre-compiled functions directly. Much faster on Windows PSOCK clusters.
#' 
#' @param n_simulations Number of simulations to run
#' @param schools Schools dataframe
#' @param network Network object
#' @param params Simulation parameters
#' @param seed_schools Vector of seed school IDs
#' @param n_initial_infected Number of initial infections per seed school
#' @param n_days Maximum simulation days
#' @param seed_start Starting seed for reproducibility
#' @param n_cores Number of parallel workers (NULL = auto-detect)
#' @param cluster_type "PSOCK" (Windows) or "FORK" (Linux/Mac), NULL = auto
#' @param verbose Print progress messages
#' @param hh_pop Household population (optional)
#' @param household_assignment Household assignment dataframe (optional)
#' @param source_files Vector of R files to source on workers (DEPRECATED - now handled automatically)
#' @return List with all simulation results
#' @export
run_parallel_simulations <- function(n_simulations,
                                      schools,
                                      network,
                                      params,
                                      seed_schools,
                                      n_initial_infected = 1,
                                      n_days = 150,
                                      seed_start = NULL,
                                      n_cores = NULL,
                                      cluster_type = NULL,
                                      verbose = TRUE,
                                      hh_pop = NULL,
                                      household_assignment = NULL,
                                      source_files = NULL) {
  
  # Check that required functions exist in global environment
  required_funcs <- c("run_network_simulation", "school_transmission", "draw_erlang")
  missing_funcs <- required_funcs[!sapply(required_funcs, function(f) {
    exists(f, where = globalenv(), mode = "function")
  })]
  
  if (length(missing_funcs) > 0) {
    # Debug: Show what functions ARE available
    all_funcs <- ls(envir = globalenv())
    sim_funcs <- all_funcs[grepl("run_|simulation|school_|transmission", all_funcs, ignore.case = TRUE)]
    
    stop(paste(
      "Required functions not found:", paste(missing_funcs, collapse = ", "),
      "\n\nAvailable simulation-related objects in global environment:",
      paste(sim_funcs, collapse = ", "),
      "\n\nMake sure to source('load_all.R') before calling this function.",
      "\nIf you already did, check for errors/warnings during loading."
    ))
  }
  
  use_households <- !is.null(hh_pop)
  
  # Auto-detect settings
  settings <- detect_parallel_settings()
  
  if (is.null(n_cores)) {
    n_cores <- settings$recommended_workers
  }
  
  if (is.null(cluster_type)) {
    cluster_type <- settings$cluster_type
  }
  
  # Don't use more cores than simulations
  n_cores <- min(n_cores, n_simulations)
  
  if (verbose) {
    cat("\n=== RUNNING PARALLEL SIMULATIONS (OPTIMIZED) ===\n")
    cat(sprintf("Number of simulations: %d\n", n_simulations))
    cat(sprintf("Parallel workers: %d\n", n_cores))
    cat(sprintf("Cluster type: %s\n", cluster_type))
    cat(sprintf("Schools: %d\n", nrow(schools)))
    cat(sprintf("Household transmission: %s\n", ifelse(use_households, "ENABLED", "DISABLED")))
    cat("\n")
  }
  
  start_time <- Sys.time()
  
  # Generate seeds for each simulation
  sim_seeds <- if (!is.null(seed_start)) {
    seed_start + (1:n_simulations)
  } else {
    sample.int(.Machine$integer.max, n_simulations)
  }
  
  # Create cluster
  if (verbose) cat("Creating parallel cluster...\n")
  
  if (cluster_type == "FORK") {
    # Fork-based parallelism (Linux/Mac) - inherits environment
    cl <- parallel::makeForkCluster(n_cores)
  } else {
    # PSOCK cluster (Windows) - need to source files on each worker
    # This includes C++ compilation which takes time, but is unavoidable
    cl <- parallel::makeCluster(n_cores, type = "PSOCK")
    
    if (verbose) cat("Initializing workers (this may take a moment for C++ compilation)...\n")
    
    # Find source files - try multiple methods
    load_all_path <- NULL
    
    # Method 1: User-provided path
    if (!is.null(source_files) && length(source_files) > 0 && file.exists(source_files[1])) {
      load_all_path <- normalizePath(source_files[1])
    }
    
    # Method 2: Try here package
    if (is.null(load_all_path)) {
      tryCatch({
        here_path <- here::here("load_all.R")
        if (file.exists(here_path)) {
          load_all_path <- normalizePath(here_path)
        }
      }, error = function(e) NULL)
    }
    
    # Method 3: Check working directory
    if (is.null(load_all_path)) {
      possible_paths <- c(
        "load_all.R",
        file.path(getwd(), "load_all.R"),
        file.path(getwd(), "modular", "load_all.R"),
        file.path(dirname(getwd()), "load_all.R")
      )
      for (p in possible_paths) {
        if (file.exists(p)) {
          load_all_path <- normalizePath(p)
          break
        }
      }
    }
    
    if (!is.null(load_all_path) && file.exists(load_all_path)) {
      if (verbose) cat(sprintf("Sourcing on workers: %s\n", load_all_path))
      
      # Export source path and load on workers
      parallel::clusterExport(cl, "load_all_path", envir = environment())
      
      # Source on workers with error handling
      worker_results <- parallel::clusterEvalQ(cl, {
        tryCatch({
          # Set working directory to the directory containing load_all.R
          setwd(dirname(load_all_path))
          suppressMessages({
            source(load_all_path, local = globalenv())
          })
          "OK"
        }, error = function(e) {
          paste("ERROR:", e$message)
        })
      })
      
      # Check if any workers failed
      failed_workers <- sapply(worker_results, function(r) startsWith(r, "ERROR"))
      if (any(failed_workers)) {
        parallel::stopCluster(cl)
        stop(paste("Failed to initialize workers. First error:", worker_results[failed_workers][1]))
      }
      
    } else {
      parallel::stopCluster(cl)
      stop(paste(
        "Could not find load_all.R file for parallel workers.",
        "\nSearched paths:", paste(possible_paths, collapse = ", "),
        "\nPlease provide the path via source_files parameter or run from the modular directory."
      ))
    }
    
    # Export data objects
    parallel::clusterExport(cl, c(
      "schools", "network", "params", "seed_schools", 
      "n_initial_infected", "n_days", "hh_pop", "household_assignment",
      "sim_seeds"
    ), envir = environment())
  }
  
  # Ensure cluster is stopped on exit
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  if (verbose) cat("Running simulations...\n")
  
  # Run simulations in parallel
  results_list <- parallel::parLapply(cl, 1:n_simulations, function(i) {
    run_single_simulation_worker(
      sim_id = i,
      schools = schools,
      network = network,
      params = params,
      seed_schools = seed_schools,
      n_initial_infected = n_initial_infected,
      n_days = n_days,
      seed = sim_seeds[i],
      hh_pop = hh_pop,
      household_assignment = household_assignment
    )
  })
  
  if (verbose) cat("Combining results...\n")
  
  # Combine results
  summary_stats <- do.call(rbind, lapply(results_list, function(r) r$summary_row))
  all_network_data <- do.call(rbind, lapply(results_list, function(r) r$network_daily_counts))
  all_school_data <- do.call(rbind, lapply(results_list, function(r) r$school_summary))
  all_school_daily_data <- do.call(rbind, lapply(results_list, function(r) r$school_daily_combined))
  
  # Combine household data
  all_hh_data <- NULL
  all_hh_summary_data <- NULL
  
  if (use_households) {
    hh_daily_list <- lapply(results_list, function(r) r$hh_daily_counts)
    hh_daily_list <- hh_daily_list[!sapply(hh_daily_list, is.null)]
    if (length(hh_daily_list) > 0) {
      all_hh_data <- do.call(rbind, hh_daily_list)
    }
    
    hh_summary_list <- lapply(results_list, function(r) r$hh_summary)
    hh_summary_list <- hh_summary_list[!sapply(hh_summary_list, is.null)]
    if (length(hh_summary_list) > 0) {
      all_hh_summary_data <- do.call(rbind, hh_summary_list)
    }
  }
  
  end_time <- Sys.time()
  total_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  if (verbose) {
    cat(sprintf("\n=== COMPLETED IN %.2f SECONDS ===\n", total_time))
    cat(sprintf("Average time per simulation: %.2f seconds\n", total_time / n_simulations))
    cat(sprintf("Speedup vs sequential (est.): %.1fx\n", n_cores * 0.8)) # Assume 80% efficiency
    cat("\n=== SUMMARY STATISTICS ===\n")
    cat(sprintf("Students Infected (Median): %.1f (95%% CI: %.1f - %.1f)\n",
                median(summary_stats$total_infected),
                quantile(summary_stats$total_infected, 0.025),
                quantile(summary_stats$total_infected, 0.975)))
    if (use_households) {
      cat(sprintf("HH Members Infected (Median): %.1f (95%% CI: %.1f - %.1f)\n",
                  median(summary_stats$hh_members_infected),
                  quantile(summary_stats$hh_members_infected, 0.025),
                  quantile(summary_stats$hh_members_infected, 0.975)))
    }
  }
  
  return(list(
    summary_stats = summary_stats,
    all_network_data = all_network_data,
    all_school_data = all_school_data,
    all_school_daily_data = all_school_daily_data,
    all_hh_data = all_hh_data,
    all_hh_summary_data = all_hh_summary_data,
    n_simulations = n_simulations,
    computation_time = total_time,
    use_households = use_households,
    params = params
  ))
}


# ==============================================================================
# Progress-Enabled Parallel Simulation (using pbapply)
# ==============================================================================

#' Run parallel simulations with progress bar (requires pbapply package)
#' 
#' @param ... Same arguments as run_parallel_simulations
#' @return Same as run_parallel_simulations
#' @export
run_parallel_simulations_progress <- function(n_simulations,
                                               schools,
                                               network,
                                               params,
                                               seed_schools,
                                               n_initial_infected = 1,
                                               n_days = 150,
                                               seed_start = NULL,
                                               n_cores = NULL,
                                               verbose = TRUE,
                                               hh_pop = NULL,
                                               household_assignment = NULL,
                                               source_files = NULL) {
  
  # Check if pbapply is available
  if (!requireNamespace("pbapply", quietly = TRUE)) {
    message("pbapply not installed. Using standard parallel without progress bar.")
    message("Install with: install.packages('pbapply')")
    return(run_parallel_simulations(
      n_simulations = n_simulations,
      schools = schools,
      network = network,
      params = params,
      seed_schools = seed_schools,
      n_initial_infected = n_initial_infected,
      n_days = n_days,
      seed_start = seed_start,
      n_cores = n_cores,
      verbose = verbose,
      hh_pop = hh_pop,
      household_assignment = household_assignment,
      source_files = source_files
    ))
  }
  
  use_households <- !is.null(hh_pop)
  settings <- detect_parallel_settings()
  
  if (is.null(n_cores)) {
    n_cores <- settings$recommended_workers
  }
  n_cores <- min(n_cores, n_simulations)
  
  if (verbose) {
    cat("\n=== RUNNING PARALLEL SIMULATIONS WITH PROGRESS ===\n")
    cat(sprintf("Simulations: %d | Workers: %d\n", n_simulations, n_cores))
  }
  
  start_time <- Sys.time()
  
  sim_seeds <- if (!is.null(seed_start)) {
    seed_start + (1:n_simulations)
  } else {
    sample.int(.Machine$integer.max, n_simulations)
  }
  
  # Create cluster
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  # Setup workers
  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      library(Rcpp)
      library(dplyr)
    })
  })
  
  if (!is.null(source_files) && length(source_files) > 0) {
    parallel::clusterExport(cl, "source_files", envir = environment())
    parallel::clusterEvalQ(cl, {
      for (f in source_files) {
        if (file.exists(f)) source(f, local = globalenv())
      }
    })
  }
  
  # Export data objects from local environment
  parallel::clusterExport(cl, c(
    "schools", "network", "params", "seed_schools", 
    "n_initial_infected", "n_days", "hh_pop", "household_assignment",
    "sim_seeds"
  ), envir = environment())
  
  # Run with progress bar
  results_list <- pbapply::pblapply(1:n_simulations, function(i) {
    run_single_simulation_worker(
      sim_id = i,
      schools = schools,
      network = network,
      params = params,
      seed_schools = seed_schools,
      n_initial_infected = n_initial_infected,
      n_days = n_days,
      seed = sim_seeds[i],
      hh_pop = hh_pop,
      household_assignment = household_assignment
    )
  }, cl = cl)
  
  # Combine results (same as run_parallel_simulations)
  summary_stats <- do.call(rbind, lapply(results_list, function(r) r$summary_row))
  all_network_data <- do.call(rbind, lapply(results_list, function(r) r$network_daily_counts))
  all_school_data <- do.call(rbind, lapply(results_list, function(r) r$school_summary))
  all_school_daily_data <- do.call(rbind, lapply(results_list, function(r) r$school_daily_combined))
  
  all_hh_data <- NULL
  all_hh_summary_data <- NULL
  
  if (use_households) {
    hh_daily_list <- lapply(results_list, function(r) r$hh_daily_counts)
    hh_daily_list <- hh_daily_list[!sapply(hh_daily_list, is.null)]
    if (length(hh_daily_list) > 0) {
      all_hh_data <- do.call(rbind, hh_daily_list)
    }
    
    hh_summary_list <- lapply(results_list, function(r) r$hh_summary)
    hh_summary_list <- hh_summary_list[!sapply(hh_summary_list, is.null)]
    if (length(hh_summary_list) > 0) {
      all_hh_summary_data <- do.call(rbind, hh_summary_list)
    }
  }
  
  end_time <- Sys.time()
  total_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  if (verbose) {
    cat(sprintf("\nCompleted in %.2f seconds (%.2f sec/sim)\n", 
                total_time, total_time / n_simulations))
  }
  
  list(
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
    n_cores = n_cores,
    use_households = use_households,
    all_hh_data = all_hh_data,
    all_hh_summary_data = all_hh_summary_data
  )
}


# ==============================================================================
# Batch Processing for Very Large Simulation Sets
# ==============================================================================

#' Run simulations in batches (useful for very large n_simulations)
#' Saves intermediate results to avoid memory issues
#' 
#' @param n_simulations Total number of simulations
#' @param batch_size Number of simulations per batch (default 100)
#' @param output_dir Directory to save batch results
#' @param ... Other arguments passed to run_parallel_simulations
#' @return Combined results from all batches
#' @export
run_batch_simulations <- function(n_simulations,
                                   batch_size = 100,
                                   output_dir = "batch_results",
                                   schools,
                                   network,
                                   params,
                                   seed_schools,
                                   n_initial_infected = 1,
                                   n_days = 150,
                                   seed_start = NULL,
                                   n_cores = NULL,
                                   verbose = TRUE,
                                   hh_pop = NULL,
                                   household_assignment = NULL,
                                   source_files = NULL) {
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  n_batches <- ceiling(n_simulations / batch_size)
  
  if (verbose) {
    cat("\n=== BATCH SIMULATION MODE ===\n")
    cat(sprintf("Total simulations: %d\n", n_simulations))
    cat(sprintf("Batch size: %d\n", batch_size))
    cat(sprintf("Number of batches: %d\n", n_batches))
    cat(sprintf("Output directory: %s\n", output_dir))
    cat("\n")
  }
  
  total_start <- Sys.time()
  batch_files <- character(n_batches)
  
  for (batch in 1:n_batches) {
    batch_start <- (batch - 1) * batch_size + 1
    batch_end <- min(batch * batch_size, n_simulations)
    batch_n <- batch_end - batch_start + 1
    
    if (verbose) {
      cat(sprintf("\n--- Batch %d/%d (sims %d-%d) ---\n", 
                  batch, n_batches, batch_start, batch_end))
    }
    
    # Adjust seed for this batch
    batch_seed <- if (!is.null(seed_start)) seed_start + batch_start - 1 else NULL
    
    # Run batch
    batch_results <- run_parallel_simulations(
      n_simulations = batch_n,
      schools = schools,
      network = network,
      params = params,
      seed_schools = seed_schools,
      n_initial_infected = n_initial_infected,
      n_days = n_days,
      seed_start = batch_seed,
      n_cores = n_cores,
      verbose = FALSE,
      hh_pop = hh_pop,
      household_assignment = household_assignment,
      source_files = source_files
    )
    
    # Adjust sim_ids to global numbering
    batch_results$summary_stats$sim_id <- batch_results$summary_stats$sim_id + batch_start - 1
    batch_results$all_network_data$sim_id <- batch_results$all_network_data$sim_id + batch_start - 1
    batch_results$all_school_data$sim_id <- batch_results$all_school_data$sim_id + batch_start - 1
    batch_results$all_school_daily_data$sim_id <- batch_results$all_school_daily_data$sim_id + batch_start - 1
    
    if (!is.null(batch_results$all_hh_data)) {
      batch_results$all_hh_data$sim_id <- batch_results$all_hh_data$sim_id + batch_start - 1
    }
    if (!is.null(batch_results$all_hh_summary_data)) {
      batch_results$all_hh_summary_data$sim_id <- batch_results$all_hh_summary_data$sim_id + batch_start - 1
    }
    
    # Save batch
    batch_file <- file.path(output_dir, sprintf("batch_%03d.rds", batch))
    saveRDS(batch_results, batch_file)
    batch_files[batch] <- batch_file
    
    if (verbose) {
      cat(sprintf("Saved: %s\n", batch_file))
    }
    
    # Clean up memory
    rm(batch_results)
    gc()
  }
  
  # Combine all batches
  if (verbose) cat("\nCombining batch results...\n")
  
  all_summary <- list()
  all_network <- list()
  all_school <- list()
  all_school_daily <- list()
  all_hh <- list()
  all_hh_summary <- list()
  
  for (i in 1:n_batches) {
    batch_results <- readRDS(batch_files[i])
    all_summary[[i]] <- batch_results$summary_stats
    all_network[[i]] <- batch_results$all_network_data
    all_school[[i]] <- batch_results$all_school_data
    all_school_daily[[i]] <- batch_results$all_school_daily_data
    
    if (!is.null(batch_results$all_hh_data)) {
      all_hh[[i]] <- batch_results$all_hh_data
    }
    if (!is.null(batch_results$all_hh_summary_data)) {
      all_hh_summary[[i]] <- batch_results$all_hh_summary_data
    }
    
    rm(batch_results)
  }
  
  combined_results <- list(
    summary_stats = do.call(rbind, all_summary),
    all_network_data = do.call(rbind, all_network),
    all_school_data = do.call(rbind, all_school),
    all_school_daily_data = do.call(rbind, all_school_daily),
    all_hh_data = if (length(all_hh) > 0) do.call(rbind, all_hh) else NULL,
    all_hh_summary_data = if (length(all_hh_summary) > 0) do.call(rbind, all_hh_summary) else NULL,
    params = params,
    schools = schools,
    network = network,
    n_simulations = n_simulations,
    n_batches = n_batches,
    batch_files = batch_files,
    computation_time = as.numeric(difftime(Sys.time(), total_start, units = "secs"))
  )
  
  # Save combined results
  combined_file <- file.path(output_dir, "combined_results.rds")
  saveRDS(combined_results, combined_file)
  
  if (verbose) {
    cat(sprintf("\nTotal time: %.2f seconds\n", combined_results$computation_time))
    cat(sprintf("Combined results saved: %s\n", combined_file))
  }
  
  return(combined_results)
}


cat("Parallel simulation utilities loaded.\n")
cat("Use detect_parallel_settings() to see recommended settings.\n")
cat("Use run_parallel_simulations() for parallel execution.\n")