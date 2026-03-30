# ==============================================================================
# UNIT TESTS: household_utils.R
# ==============================================================================
# Tests for: haversine_distance, get_households_with_children
# ==============================================================================

project_root <- Sys.getenv("PROJECT_ROOT", unset = normalizePath("."))
source(file.path(project_root, "codes", "household_utils.R"))

# ==============================================================================
# haversine_distance
# ==============================================================================

test_that("haversine_distance: same point returns 0", {
  d <- haversine_distance(-80, 34, -80, 34)
  expect_equal(d, 0, tolerance = 1e-10)
})

test_that("haversine_distance: known distance is approximately correct", {
  # Columbia, SC to Charlotte, NC is roughly 170 km
  d <- haversine_distance(-81.035, 33.997, -80.843, 35.227)
  expect_gt(d, 130)
  expect_lt(d, 210)
})

test_that("haversine_distance: is symmetric (A->B == B->A)", {
  d1 <- haversine_distance(-80, 34, -81, 35)
  d2 <- haversine_distance(-81, 35, -80, 34)
  expect_equal(d1, d2, tolerance = 1e-10)
})

test_that("haversine_distance: distance increases monotonically with separation", {
  d_small <- haversine_distance(0, 0, 0.1, 0)
  d_large <- haversine_distance(0, 0, 1.0, 0)
  expect_lt(d_small, d_large)
})

test_that("haversine_distance: returns numeric value", {
  d <- haversine_distance(-80, 34, -79, 34)
  expect_true(is.numeric(d))
  expect_false(is.na(d))
})

# ==============================================================================
# get_households_with_children
# ==============================================================================

make_synpop <- function() {
  data.frame(
    hh_id      = c(1, 1, 1, 2, 2, 3),
    person_id  = 1:6,
    agep       = c(40, 10, 8, 35, 3, 42),   # hh1: adult+2 kids; hh2: adult+toddler; hh3: adult only
    hh_size    = c(3, 3, 3, 2, 2, 1),
    lon_4326   = rep(-80.0, 6),
    lat_4326   = rep( 34.0, 6),
    stringsAsFactors = FALSE
  )
}

test_that("get_households_with_children: returns only households with school-age children", {
  synpop <- make_synpop()
  result <- suppressMessages(get_households_with_children(synpop, min_age = 5, max_age = 18))
  # hh1 has two school-age kids (10, 8); hh2 has a toddler (3); hh3 is adult only
  expect_true(all(result$hh_id %in% c(1)))
  expect_equal(length(unique(result$hh_id)), 1)
})

test_that("get_households_with_children: all household members are included (not just kids)", {
  synpop <- make_synpop()
  result <- suppressMessages(get_households_with_children(synpop, min_age = 5, max_age = 18))
  # All 3 members of hh1 should be in the result
  expect_equal(nrow(result), 3)
})

test_that("get_households_with_children: adds n_school_age_children column", {
  synpop <- make_synpop()
  result <- suppressMessages(get_households_with_children(synpop))
  expect_true("n_school_age_children" %in% names(result))
})

test_that("get_households_with_children: n_school_age_children is correct", {
  synpop <- make_synpop()
  result <- suppressMessages(get_households_with_children(synpop))
  # hh1 has 2 school-age children
  hh1_rows <- result[result$hh_id == 1, ]
  expect_equal(unique(hh1_rows$n_school_age_children), 2)
})

test_that("get_households_with_children: no school-age children returns empty data frame", {
  synpop <- data.frame(
    hh_id = c(1, 1), person_id = 1:2,
    agep  = c(40, 3), hh_size = c(2, 2),
    lon_4326 = c(-80, -80), lat_4326 = c(34, 34),
    stringsAsFactors = FALSE
  )
  result <- suppressMessages(get_households_with_children(synpop, min_age = 5, max_age = 18))
  expect_equal(nrow(result), 0)
})

test_that("get_households_with_children: handles duplicate columns gracefully", {
  synpop <- make_synpop()
  # Duplicate a column name to simulate real-world messiness
  synpop$is_school_age <- FALSE   # pre-existing column - should be removed and recomputed
  result <- suppressMessages(get_households_with_children(synpop))
  expect_true("n_school_age_children" %in% names(result))
})
