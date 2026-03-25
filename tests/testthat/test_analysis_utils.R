# ==============================================================================
# UNIT TESTS: analysis_utils.R
# ==============================================================================
# Tests for: create_outbreak_summary_table, calculate_school_statistics,
#            create_probability_table, summarize_infection_timing,
#            calculate_first_infection_times, create_school_infection_summary,
#            format_school_infection_table
# ==============================================================================

project_root <- Sys.getenv("PROJECT_ROOT", unset = normalizePath("."))
source(file.path(project_root, "analysis_utils.R"))

# ==============================================================================
# Shared mock results builder
# ==============================================================================

# Build a minimal mock 'results' object that mirrors what
# run_multiple_network_simulations() returns.
make_mock_results <- function(n_sims = 5, n_schools = 3) {
  set.seed(1)

  schools <- data.frame(
    school_id           = 1:n_schools,
    school_name         = paste0("School", 1:n_schools),
    school_size         = c(200, 300, 400),
    vaccination_coverage = c(0.90, 0.85, 0.95),
    stringsAsFactors    = FALSE
  )

  # per-simulation, per-school data
  all_school_data <- do.call(rbind, lapply(1:n_sims, function(s) {
    data.frame(
      sim                 = s,
      school_id           = 1:n_schools,
      school_name         = paste0("School", 1:n_schools),
      school_size         = c(200, 300, 400),
      vaccination_coverage = c(0.90, 0.85, 0.95),
      total_infected      = c(sample(0:15, 1), sample(0:20, 1), sample(0:10, 1)),
      attack_rate         = c(runif(1, 0, 0.1), runif(1, 0, 0.05), runif(1, 0, 0.03)),
      stringsAsFactors    = FALSE
    )
  }))

  # daily school counts (needed by calculate_first_infection_times)
  all_school_daily_data <- do.call(rbind, lapply(1:n_sims, function(s) {
    do.call(rbind, lapply(1:n_schools, function(sch) {
      n_days <- 30
      # School 1 gets infected starting day 5; school 2 day 10; school 3 never
      daily_infected <- switch(as.character(sch),
        "1" = c(rep(0, 4), rep(5, n_days - 4)),
        "2" = c(rep(0, 9), rep(3, n_days - 9)),
        "3" = rep(0, n_days)
      )
      data.frame(
        sim                 = s,
        school_id           = sch,
        school_name         = paste0("School", sch),
        day                 = 1:n_days,
        school_size         = c(200, 300, 400)[sch],
        vaccination_coverage = c(0.90, 0.85, 0.95)[sch],
        S   = pmax(0, 200 - daily_infected),
        V   = 0L,
        E   = 0L,
        P   = 0L,
        Ra  = daily_infected,
        Iso = 0L,
        R   = 0L,
        stringsAsFactors = FALSE
      )
    }))
  }))

  summary_stats <- data.frame(
    sim            = 1:n_sims,
    total_infected = sample(5:50, n_sims, replace = TRUE),
    schools_affected = sample(1:n_schools, n_sims, replace = TRUE),
    actual_days    = sample(20:80, n_sims, replace = TRUE),
    stringsAsFactors = FALSE
  )

  list(
    schools              = schools,
    all_school_data      = all_school_data,
    all_school_daily_data = all_school_daily_data,
    summary_stats        = summary_stats,
    seed_schools         = 1L
  )
}

# ==============================================================================
# create_outbreak_summary_table
# ==============================================================================

test_that("create_outbreak_summary_table: returns a data frame", {
  results <- make_mock_results()
  tbl <- create_outbreak_summary_table(results)
  expect_s3_class(tbl, "data.frame")
})

test_that("create_outbreak_summary_table: has 'Metric' and 'Value' columns", {
  results <- make_mock_results()
  tbl <- create_outbreak_summary_table(results)
  expect_true(all(c("Metric", "Value") %in% names(tbl)))
})

test_that("create_outbreak_summary_table: contains expected row metrics", {
  results <- make_mock_results()
  tbl <- create_outbreak_summary_table(results)
  expect_true(any(grepl("Total Infected", tbl$Metric)))
  expect_true(any(grepl("Duration", tbl$Metric)))
})

# ==============================================================================
# calculate_school_statistics
# ==============================================================================

test_that("calculate_school_statistics: returns one row per school", {
  results <- make_mock_results(n_sims = 10, n_schools = 3)
  stats <- calculate_school_statistics(results)
  expect_equal(nrow(stats), 3)
})

test_that("calculate_school_statistics: contains required columns", {
  results <- make_mock_results()
  stats <- calculate_school_statistics(results)
  expected_cols <- c("school_id", "school_name", "median_infected",
                     "mean_infected", "attack_rate_pct")
  expect_true(all(expected_cols %in% names(stats)))
})

test_that("calculate_school_statistics: attack_rate_pct is between 0 and 100", {
  results <- make_mock_results()
  stats <- calculate_school_statistics(results)
  expect_true(all(stats$attack_rate_pct >= 0 & stats$attack_rate_pct <= 100))
})

# ==============================================================================
# create_probability_table
# ==============================================================================

test_that("create_probability_table: returns data frame with correct columns", {
  results <- make_mock_results()
  pt <- create_probability_table(results, thresholds = c(5, 10, 20))
  expect_s3_class(pt, "data.frame")
  expect_true(all(c("Threshold", "Count", "Probability") %in% names(pt)))
})

test_that("create_probability_table: row count equals number of thresholds", {
  results <- make_mock_results()
  thresholds <- c(5, 10, 20, 50)
  pt <- create_probability_table(results, thresholds = thresholds)
  expect_equal(nrow(pt), length(thresholds))
})

test_that("create_probability_table: higher threshold has lower or equal count", {
  results <- make_mock_results()
  pt <- create_probability_table(results, thresholds = c(10, 50))
  expect_gte(pt$Count[1], pt$Count[2])
})

# ==============================================================================
# calculate_first_infection_times
# ==============================================================================

test_that("calculate_first_infection_times: returns a data frame", {
  results <- make_mock_results()
  fit <- calculate_first_infection_times(results)
  expect_s3_class(fit, "data.frame")
})

test_that("calculate_first_infection_times: has expected columns", {
  results <- make_mock_results()
  fit <- calculate_first_infection_times(results)
  expect_true(all(c("sim", "school_id", "first_infection_day") %in% names(fit)))
})

test_that("calculate_first_infection_times: school 3 (never infected) has NA first day", {
  results <- make_mock_results()
  fit <- calculate_first_infection_times(results)
  school3 <- fit[fit$school_id == 3, ]
  expect_true(all(is.na(school3$first_infection_day)))
})

test_that("calculate_first_infection_times: school 1 infected before school 2", {
  results <- make_mock_results(n_sims = 3)
  fit <- calculate_first_infection_times(results)
  med1 <- median(fit$first_infection_day[fit$school_id == 1], na.rm = TRUE)
  med2 <- median(fit$first_infection_day[fit$school_id == 2], na.rm = TRUE)
  expect_lt(med1, med2)
})

# ==============================================================================
# summarize_infection_timing
# ==============================================================================

test_that("summarize_infection_timing: returns summary per school", {
  results <- make_mock_results()
  fit <- calculate_first_infection_times(results)
  summ <- summarize_infection_timing(fit)
  expect_s3_class(summ, "data.frame")
  expect_true("median_first_day" %in% names(summ))
  expect_true("prob_infected" %in% names(summ))
})

test_that("summarize_infection_timing: prob_infected is between 0 and 1", {
  results <- make_mock_results()
  fit <- calculate_first_infection_times(results)
  summ <- summarize_infection_timing(fit)
  expect_true(all(summ$prob_infected >= 0 & summ$prob_infected <= 1))
})

# ==============================================================================
# create_school_infection_summary
# ==============================================================================

test_that("create_school_infection_summary: returns one row per school", {
  results <- make_mock_results(n_sims = 10, n_schools = 3)
  summ <- create_school_infection_summary(results)
  expect_equal(nrow(summ), 3)
})

test_that("create_school_infection_summary: ci_lower <= median <= ci_upper", {
  results <- make_mock_results(n_sims = 20, n_schools = 3)
  summ <- create_school_infection_summary(results)
  expect_true(all(summ$ci_lower <= summ$median_infected + 1e-9))
  expect_true(all(summ$median_infected <= summ$ci_upper + 1e-9))
})

test_that("create_school_infection_summary: contains is_seed column", {
  results <- make_mock_results()
  summ <- create_school_infection_summary(results)
  expect_true("is_seed" %in% names(summ))
})

# ==============================================================================
# format_school_infection_table
# ==============================================================================

test_that("format_school_infection_table: returns formatted data frame", {
  results <- make_mock_results()
  summ <- create_school_infection_summary(results)
  fmt <- format_school_infection_table(summ)
  expect_s3_class(fmt, "data.frame")
  expect_true("School" %in% names(fmt))
  expect_true("Mean Infected" %in% names(fmt))
  expect_true("95% CI" %in% names(fmt))
})

test_that("format_school_infection_table: Vax % is between 0 and 100", {
  results <- make_mock_results()
  summ <- create_school_infection_summary(results)
  fmt <- format_school_infection_table(summ)
  expect_true(all(fmt$`Vax %` >= 0 & fmt$`Vax %` <= 100))
})
