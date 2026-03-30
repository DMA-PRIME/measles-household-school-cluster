# ==============================================================================
# UNIT TESTS: simulation_utils.R
# ==============================================================================
# Tests for: school_transmission, between_school_transmission,
#            run_network_simulation, run_multiple_network_simulations
#
# NOTE: These tests require Rcpp compilation of the C++ functions.
#       They are therefore integration tests that exercise the full stack.
# ==============================================================================

project_root <- Sys.getenv("PROJECT_ROOT", unset = normalizePath("."))

# Source all dependencies in the correct order
source(file.path(project_root, "codes", "rcpp_transmission.R"))
source(file.path(project_root, "codes", "population_utils.R"))
source(file.path(project_root, "codes", "network_utils.R"))
source(file.path(project_root, "codes", "household_utils.R"))
source(file.path(project_root, "codes", "simulation_utils.R"))

# ==============================================================================
# Local helper: build a minimal network list from an adjacency matrix
# ==============================================================================

make_network_from_matrix <- function(adj_matrix) {
  n <- nrow(adj_matrix)
  g <- if (sum(adj_matrix) > 0) {
    graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
  } else {
    make_empty_graph(n, directed = FALSE)
  }
  list(graph = g, adjacency = adj_matrix, n_schools = n, n_edges = ecount(g))
}

# ==============================================================================
# Shared test fixtures
# ==============================================================================

make_sim_params <- function() {
  list(
    avg_class_size                   = 25,
    age_range                        = c(5, 18),
    latent_mean                      = 8,
    latent_shape                     = 3,
    infectious_mean                  = 10,
    infectious_shape                 = 3,
    prodromal_period                 = 4,
    c_within                         = 10,
    c_between                        = 2,
    c_between_school                 = 0.5,
    p_within                         = 0.05,
    p_between                        = 0.05,
    prodromal_infectiousness_multiplier = 0.5,
    rash_infectiousness_multiplier      = 1.0,
    vaccine_efficacy                 = 0.97,
    vaccine_infectiousness_reduction = 0.5,
    isolation_delay_index            = 1,
    isolation_delay_secondary        = 3,
    quarantine_contacts              = FALSE,
    quarantine_efficacy              = 0,
    quarantine_duration              = 14,
    no_intervention                  = FALSE
  )
}

make_schools <- function(n = 3, size = 200, vax = 0.9) {
  data.frame(
    school_id           = 1:n,
    school_name         = paste0("School", 1:n),
    school_size         = size,
    vaccination_coverage = vax,
    stringsAsFactors    = FALSE
  )
}

make_3school_network <- function() {
  make_network_from_matrix(matrix(c(0, 0.5, 0.3,
                                    0.5, 0,   0.8,
                                    0.3, 0.8, 0  ), nrow = 3))
}

# ==============================================================================
# school_transmission
# ==============================================================================

test_that("school_transmission: returns list with population and contact_history", {
  set.seed(42)
  params <- make_sim_params()
  pop <- create_school_population(1, 100, 25, c(5, 11))
  pop <- initialize_vaccination_school(pop, 0.9)

  # Seed one prodromal case
  susc <- which(pop$state == "S")
  pop$state[susc[1]] <- "P"
  pop$time_in_state[susc[1]] <- 0
  pop$time_since_prodromal[susc[1]] <- 2
  pop$latent_duration[susc[1]] <- 0
  infectious_dur <- draw_erlang(1, params$infectious_mean, params$infectious_shape)
  pop$infectious_duration[susc[1]] <- infectious_dur
  pop$prodromal_duration[susc[1]] <- params$prodromal_period
  pop$rash_duration[susc[1]] <- max(1, infectious_dur - params$prodromal_period)
  pop$is_school_index[susc[1]] <- TRUE

  ch <- ContactHistory$new(window_size = 7)
  result <- school_transmission(pop, params, ch)

  expect_true(is.list(result))
  expect_true("population" %in% names(result))
  expect_true("contact_history" %in% names(result))
})

test_that("school_transmission: population has same number of rows", {
  set.seed(42)
  params <- make_sim_params()
  pop <- create_school_population(1, 100, 25, c(5, 11))
  pop <- initialize_vaccination_school(pop, 0.9)
  ch <- ContactHistory$new(window_size = 7)
  result <- school_transmission(pop, params, ch)
  expect_equal(nrow(result$population), nrow(pop))
})

test_that("school_transmission: no transmission without infectious individuals", {
  set.seed(99)
  params <- make_sim_params()
  pop <- create_school_population(1, 100, 25, c(5, 11))
  pop <- initialize_vaccination_school(pop, 0.0)
  # No infectious individuals
  ch <- ContactHistory$new(window_size = 7)
  result <- school_transmission(pop, params, ch)
  # No new exposures expected
  n_exposed <- sum(result$population$state == "E")
  expect_equal(n_exposed, 0)
})

# ==============================================================================
# between_school_transmission
# ==============================================================================

test_that("between_school_transmission: no change when no infectors present", {
  set.seed(42)
  params <- make_sim_params()
  schools <- make_schools(3)
  pops <- lapply(1:3, function(i) {
    p <- create_school_population(i, 100, 25, c(5, 18))
    initialize_vaccination_school(p, 0.9)
  })

  net <- make_3school_network()

  result_pops <- between_school_transmission(pops, net, params)
  # States should be unchanged (no infectors)
  for (i in 1:3) {
    expect_equal(result_pops[[i]]$state, pops[[i]]$state)
  }
})

test_that("between_school_transmission: returns list of populations", {
  set.seed(42)
  params <- make_sim_params()
  pops <- lapply(1:2, function(i) {
    p <- create_school_population(i, 50, 25, c(5, 18))
    initialize_vaccination_school(p, 0.9)
  })
  adj <- matrix(c(0, 0.8, 0.8, 0), nrow = 2)
  net <- make_network_from_matrix(adj)

  result_pops <- between_school_transmission(pops, net, params)
  expect_equal(length(result_pops), 2)
})

# ==============================================================================
# run_network_simulation
# ==============================================================================

test_that("run_network_simulation: returns list with required fields", {
  set.seed(1)
  params <- make_sim_params()
  schools <- make_schools(3, size = 100)
  net <- make_3school_network()

  result <- run_network_simulation(
    schools = schools, network = net, params = params,
    seed_schools = 1, n_initial_infected = 1, n_days = 30, seed = 42
  )

  expect_true(is.list(result))
  expect_true("school_daily_counts" %in% names(result))
  expect_true("school_summary" %in% names(result))
})

test_that("run_network_simulation: daily_counts has correct dimensions", {
  set.seed(1)
  params <- make_sim_params()
  schools <- make_schools(2, size = 80)
  net <- make_network_from_matrix(matrix(c(0, 0.5, 0.5, 0), 2))

  result <- run_network_simulation(
    schools = schools, network = net, params = params,
    seed_schools = 1, n_days = 20, seed = 10
  )

  # daily_counts should be a list with one entry per school
  expect_equal(length(result$daily_counts), 2)
})

test_that("run_network_simulation: seed school has at least one infected case", {
  set.seed(1)
  params <- make_sim_params()
  schools <- make_schools(2, size = 100, vax = 0.0)  # No vaccination
  net <- make_network_from_matrix(matrix(c(0, 0, 0, 0), 2))

  result <- run_network_simulation(
    schools = schools, network = net, params = params,
    seed_schools = 1, n_initial_infected = 1, n_days = 50, seed = 7
  )

  # School 1 should have > 0 total infected
  school1 <- result$school_summary[result$school_summary$school_id == 1, ]
  expect_gt(school1$total_infected, 0)
})

test_that("run_network_simulation: isolated school 2 is not affected when no network edges", {
  set.seed(1)
  params <- make_sim_params()
  schools <- make_schools(2, size = 100, vax = 0.0)
  # Zero adjacency = no between-school transmission
  net <- make_network_from_matrix(matrix(c(0, 0, 0, 0), 2))

  result <- run_network_simulation(
    schools = schools, network = net, params = params,
    seed_schools = 1, n_initial_infected = 1, n_days = 80, seed = 5
  )

  school2 <- result$school_summary[result$school_summary$school_id == 2, ]
  expect_equal(school2$total_infected, 0)
})

test_that("run_network_simulation: higher vaccination reduces total cases", {
  params <- make_sim_params()
  schools_low_vax  <- make_schools(1, size = 300, vax = 0.0)
  schools_high_vax <- make_schools(1, size = 300, vax = 0.95)
  net <- make_network_from_matrix(matrix(0, 1, 1))

  result_low  <- run_network_simulation(schools_low_vax,  net, params, 1, 1, 60, seed = 42)
  result_high <- run_network_simulation(schools_high_vax, net, params, 1, 1, 60, seed = 42)

  infected_low  <- result_low$school_summary$total_infected[1]
  infected_high <- result_high$school_summary$total_infected[1]
  expect_gte(infected_low, infected_high)
})

# ==============================================================================
# run_multiple_network_simulations
# ==============================================================================

test_that("run_multiple_network_simulations: returns results for all simulations", {
  params <- make_sim_params()
  schools <- make_schools(2, size = 80)
  net <- make_network_from_matrix(matrix(c(0, 0.5, 0.5, 0), 2))

  results <- run_multiple_network_simulations(
    n_simulations = 3, schools = schools, network = net,
    params = params, seed_schools = 1, n_initial_infected = 1,
    n_days = 20, seed_start = 1, verbose = FALSE
  )

  expect_true(is.list(results))
  expect_equal(nrow(results$summary_stats), 3)
})

test_that("run_multiple_network_simulations: all_school_data has n_sims * n_schools rows", {
  params <- make_sim_params()
  schools <- make_schools(2, size = 80)
  net <- make_network_from_matrix(matrix(c(0, 0.5, 0.5, 0), 2))

  results <- run_multiple_network_simulations(
    n_simulations = 4, schools = schools, network = net,
    params = params, seed_schools = 1, n_days = 20,
    seed_start = 1, verbose = FALSE
  )

  expect_equal(nrow(results$all_school_data), 4 * 2)
})
