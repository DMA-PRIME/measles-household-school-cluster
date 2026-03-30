# ==============================================================================
# UNIT TESTS: population_utils.R
# ==============================================================================
# Tests for: draw_erlang, safe_school_size, ContactHistory (R6 class),
#            create_school_population, initialize_vaccination_school,
#            seed_infections, update_disease_states
# ==============================================================================

# Source the module under test (R6 only; no Rcpp needed for these functions)
project_root <- Sys.getenv("PROJECT_ROOT", unset = normalizePath("."))
source(file.path(project_root, "codes", "population_utils.R"))

# ==============================================================================
# draw_erlang
# ==============================================================================

test_that("draw_erlang: returns n samples", {
  set.seed(42)
  samples <- draw_erlang(100, mean = 8, shape = 3)
  expect_equal(length(samples), 100)
})

test_that("draw_erlang: all samples are >= 1", {
  set.seed(42)
  samples <- draw_erlang(1000, mean = 8, shape = 3)
  expect_true(all(samples >= 1))
})

test_that("draw_erlang: mean is approximately correct", {
  set.seed(42)
  samples <- draw_erlang(10000, mean = 8, shape = 3)
  expect_equal(mean(samples), 8, tolerance = 0.5)
})

test_that("draw_erlang: returns integers (rounded values)", {
  set.seed(42)
  samples <- draw_erlang(100, mean = 8, shape = 3)
  expect_true(all(samples == round(samples)))
})

# ==============================================================================
# safe_school_size
# ==============================================================================

test_that("safe_school_size: numeric input passes through", {
  expect_equal(safe_school_size(500), 500)
  expect_equal(safe_school_size(c(100, 200, 300)), c(100, 200, 300))
})

test_that("safe_school_size: '<X' returns default_small", {
  expect_equal(safe_school_size("<10"), 5)
  expect_equal(safe_school_size("<10", default_small = 3), 3)
})

test_that("safe_school_size: '>X' returns the threshold number", {
  expect_equal(safe_school_size(">1000"), 1000)
})

test_that("safe_school_size: NA / empty / 'N/A' / 'NA' return NA", {
  expect_true(is.na(safe_school_size(NA)))
  expect_true(is.na(safe_school_size("")))
  expect_true(is.na(safe_school_size("N/A")))
  expect_true(is.na(safe_school_size("NA")))
})

test_that("safe_school_size: numeric-as-string is converted", {
  expect_equal(safe_school_size("250"), 250)
})

# ==============================================================================
# ContactHistory (R6 class)
# ==============================================================================

test_that("ContactHistory: initializes with empty contact list", {
  ch <- ContactHistory$new(window_size = 7)
  contacts <- ch$get_all_contacts()
  expect_equal(length(contacts$infector_ids), 0)
  expect_equal(length(contacts$target_ids),   0)
})

test_that("ContactHistory: add_contacts stores infector and target ids", {
  ch <- ContactHistory$new(window_size = 7)
  ch$add_contacts(c(1L, 2L), c(3L, 4L))
  contacts <- ch$get_all_contacts()
  expect_equal(sort(contacts$infector_ids), c(1L, 2L))
  expect_equal(sort(contacts$target_ids),   c(3L, 4L))
})

test_that("ContactHistory: add_contacts with no infectors is a no-op", {
  ch <- ContactHistory$new(window_size = 7)
  ch$add_contacts(integer(0), integer(0))
  contacts <- ch$get_all_contacts()
  expect_equal(length(contacts$infector_ids), 0)
})

test_that("ContactHistory: clear removes all stored contacts", {
  ch <- ContactHistory$new(window_size = 7)
  ch$add_contacts(c(1L), c(2L))
  ch$clear()
  contacts <- ch$get_all_contacts()
  expect_equal(length(contacts$infector_ids), 0)
})

test_that("ContactHistory: window limits stored history", {
  ch <- ContactHistory$new(window_size = 2)
  for (i in 1:5) {
    ch$add_contacts(as.integer(i), as.integer(i + 10))
  }
  # Should keep at most window_size entries
  expect_lte(length(ch$contact_list), 2)
})

# ==============================================================================
# create_school_population
# ==============================================================================

test_that("create_school_population: returns data frame with correct number of rows", {
  pop <- create_school_population(school_id = 1, school_size = 100,
                                  avg_class_size = 25, age_range = c(5, 11))
  expect_equal(nrow(pop), 100)
})

test_that("create_school_population: all students start in state 'S'", {
  pop <- create_school_population(1, 50, 25, c(5, 11))
  expect_true(all(pop$state == "S"))
})

test_that("create_school_population: school_id column is correct", {
  pop <- create_school_population(school_id = 3, school_size = 50,
                                  avg_class_size = 25, age_range = c(5, 11))
  expect_true(all(pop$school_id == 3))
})

test_that("create_school_population: ages are within range", {
  pop <- create_school_population(1, 100, 25, age_range = c(14, 18))
  expect_true(all(pop$age >= 14 & pop$age <= 18))
})

test_that("create_school_population: required columns are present", {
  pop <- create_school_population(1, 50, 25, c(5, 11))
  expected_cols <- c("student_id", "school_id", "class_id", "age", "state",
                     "is_vaccinated", "is_isolated", "is_quarantined",
                     "latent_duration", "infectious_duration",
                     "prodromal_duration", "rash_duration", "hh_id")
  expect_true(all(expected_cols %in% names(pop)))
})

test_that("create_school_population: class IDs cycle correctly", {
  pop <- create_school_population(1, 100, 25, c(5, 11))
  n_classes <- ceiling(100 / 25)
  expect_true(all(pop$class_id %in% 1:n_classes))
})

# ==============================================================================
# initialize_vaccination_school
# ==============================================================================

test_that("initialize_vaccination_school: zero coverage leaves all susceptible", {
  pop <- create_school_population(1, 100, 25, c(5, 11))
  pop_vacc <- initialize_vaccination_school(pop, vaccination_coverage = 0)
  expect_true(all(pop_vacc$state == "S"))
  expect_true(all(!pop_vacc$is_vaccinated))
})

test_that("initialize_vaccination_school: 100% coverage vaccinates all", {
  pop <- create_school_population(1, 100, 25, c(5, 11))
  pop_vacc <- initialize_vaccination_school(pop, vaccination_coverage = 1)
  expect_true(all(pop_vacc$is_vaccinated))
  expect_true(all(pop_vacc$state == "V"))
})

test_that("initialize_vaccination_school: coverage approximates target proportion", {
  set.seed(42)
  pop <- create_school_population(1, 1000, 25, c(5, 11))
  pop_vacc <- initialize_vaccination_school(pop, vaccination_coverage = 0.9)
  actual <- mean(pop_vacc$is_vaccinated)
  expect_equal(actual, 0.9, tolerance = 0.02)
})

test_that("initialize_vaccination_school: vaccine_failed flags assigned for vaccinated", {
  set.seed(42)
  pop <- create_school_population(1, 500, 25, c(5, 11))
  pop_vacc <- initialize_vaccination_school(pop, vaccination_coverage = 1,
                                            vaccine_efficacy = 0.97)
  vacc_rows <- pop_vacc[pop_vacc$is_vaccinated, ]
  # Only vaccinated individuals can have vaccine_failed TRUE
  expect_true(all(pop_vacc$vaccine_failed[!pop_vacc$is_vaccinated] == FALSE))
  # Approximately (1 - 0.97) = 3% should have vaccine failures
  expect_equal(mean(vacc_rows$vaccine_failed), 0.03, tolerance = 0.05)
})

# ==============================================================================
# seed_infections
# ==============================================================================

make_params <- function() {
  list(
    latent_mean        = 8,
    latent_shape       = 3,
    infectious_mean    = 10,
    infectious_shape   = 3,
    prodromal_period   = 4,
    no_intervention    = FALSE,
    isolation_delay_index     = 1,
    isolation_delay_secondary = 3,
    quarantine_contacts       = FALSE,
    quarantine_efficacy       = 0,
    quarantine_duration       = 14
  )
}

test_that("seed_infections: index case is placed in state 'P'", {
  set.seed(42)
  pop <- create_school_population(1, 50, 25, c(5, 11))
  params <- make_params()
  pops <- seed_infections(list(pop), seed_schools = 1, n_infected = 1, params = params)
  expect_equal(sum(pops[[1]]$state == "P"), 1)
})

test_that("seed_infections: index case has is_index = TRUE", {
  set.seed(42)
  pop <- create_school_population(1, 50, 25, c(5, 11))
  params <- make_params()
  pops <- seed_infections(list(pop), seed_schools = 1, n_infected = 1, params = params)
  index_rows <- pops[[1]][pops[[1]]$is_index, ]
  expect_equal(nrow(index_rows), 1)
})

test_that("seed_infections: additional seeds are placed in state 'E'", {
  set.seed(42)
  pop <- create_school_population(1, 50, 25, c(5, 11))
  params <- make_params()
  pops <- seed_infections(list(pop), seed_schools = 1, n_infected = 3, params = params)
  # 1 in P, 2 in E
  expect_equal(sum(pops[[1]]$state == "P"), 1)
  expect_equal(sum(pops[[1]]$state == "E"), 2)
})

test_that("seed_infections: infection_source is 'seed' for seeded cases", {
  set.seed(42)
  pop <- create_school_population(1, 50, 25, c(5, 11))
  params <- make_params()
  pops <- seed_infections(list(pop), seed_schools = 1, n_infected = 3, params = params)
  infected <- pops[[1]][pops[[1]]$state %in% c("P", "E"), ]
  expect_true(all(infected$infection_source == "seed"))
})

# ==============================================================================
# update_disease_states
# ==============================================================================

make_pop_with_state <- function(state_val, time_in_state = 0,
                                 time_since_prodromal = NA,
                                 latent_dur = 5, infectious_dur = 10,
                                 prodromal_dur = 4) {
  pop <- create_school_population(1, 1, 1, c(10, 10))
  pop$state           <- state_val
  pop$time_in_state   <- time_in_state
  pop$time_since_prodromal <- time_since_prodromal
  pop$latent_duration      <- latent_dur
  pop$infectious_duration  <- infectious_dur
  pop$prodromal_duration   <- prodromal_dur
  pop$rash_duration        <- infectious_dur - prodromal_dur
  pop
}

test_that("update_disease_states: E->P after latent period", {
  params <- make_params()
  pop <- make_pop_with_state("E", time_in_state = 5, latent_dur = 5)
  pop2 <- update_disease_states(pop, params)
  expect_equal(pop2$state, "P")
})

test_that("update_disease_states: P->Ra after prodromal period", {
  params <- make_params()
  pop <- make_pop_with_state("P", time_in_state = 4,
                              time_since_prodromal = 0, prodromal_dur = 4)
  pop2 <- update_disease_states(pop, params)
  expect_equal(pop2$state, "Ra")
})

test_that("update_disease_states: Ra->Iso for index case after delay", {
  params <- make_params()
  params$isolation_delay_index <- 1
  pop <- make_pop_with_state("Ra", time_in_state = 1,
                              time_since_prodromal = 5)
  pop$is_school_index <- TRUE
  pop2 <- update_disease_states(pop, params)
  expect_equal(pop2$state, "Iso")
})

test_that("update_disease_states: Iso->R after infectious duration", {
  params <- make_params()
  pop <- make_pop_with_state("Iso", time_in_state = 10,
                              time_since_prodromal = 10, infectious_dur = 10)
  pop$is_isolated <- TRUE
  pop2 <- update_disease_states(pop, params)
  expect_equal(pop2$state, "R")
})

test_that("update_disease_states: QE->QP after latent period", {
  params <- make_params()
  pop <- make_pop_with_state("QE", time_in_state = 5, latent_dur = 5)
  pop$is_quarantined <- TRUE
  pop2 <- update_disease_states(pop, params)
  expect_equal(pop2$state, "QP")
})

test_that("update_disease_states: QS released after quarantine_duration", {
  params <- make_params()
  params$quarantine_duration <- 14
  pop <- make_pop_with_state("QS", time_in_state = 14)
  pop$is_quarantined <- TRUE
  pop2 <- update_disease_states(pop, params)
  expect_equal(pop2$state, "S")
  expect_false(pop2$is_quarantined)
})

test_that("update_disease_states: QV released to V after quarantine_duration", {
  params <- make_params()
  params$quarantine_duration <- 14
  pop <- make_pop_with_state("QV", time_in_state = 14)
  pop$is_vaccinated  <- TRUE
  pop$is_quarantined <- TRUE
  pop2 <- update_disease_states(pop, params)
  expect_equal(pop2$state, "V")
  expect_false(pop2$is_quarantined)
})

test_that("update_disease_states: Ra->R in no_intervention mode", {
  params <- make_params()
  params$no_intervention <- TRUE
  pop <- make_pop_with_state("Ra", time_in_state = 6,
                              time_since_prodromal = 10, infectious_dur = 10)
  pop2 <- update_disease_states(pop, params)
  expect_equal(pop2$state, "R")
})

test_that("update_disease_states: time_in_state is incremented each call", {
  params <- make_params()
  pop <- make_pop_with_state("S", time_in_state = 0)
  pop2 <- update_disease_states(pop, params)
  expect_equal(pop2$time_in_state, 1)
})
