# ==============================================================================
# UNIT TESTS: network_utils.R
# ==============================================================================
# Tests for: generate_distance_network
# ==============================================================================

project_root <- Sys.getenv("PROJECT_ROOT", unset = normalizePath("."))
source(file.path(project_root, "codes", "network_utils.R"))

# ==============================================================================
# Shared helper: validate network structure
# ==============================================================================

expect_valid_network <- function(net, n_schools) {
  expect_true(is.list(net))
  expect_true(all(c("graph", "adjacency", "n_schools", "n_edges") %in% names(net)))
  expect_equal(net$n_schools, n_schools)
  expect_equal(nrow(net$adjacency), n_schools)
  expect_equal(ncol(net$adjacency), n_schools)
  expect_true(all(net$adjacency >= 0))
  expect_true(all(net$adjacency <= 1.0 + 1e-9))  # weights normalised to [0,1]
  # Diagonal must be 0 (no self-loops)
  expect_equal(sum(diag(net$adjacency)), 0)
}

# ==============================================================================
# generate_distance_network
# ==============================================================================

test_that("generate_distance_network: returns valid network from school coords", {
  # Place schools ~10 km apart (same latitude, spacing ~0.09 degrees longitude)
  schools <- data.frame(
    school_id   = 1:4,
    school_name = paste0("School", 1:4),
    lon         = c(-80.00, -80.09, -80.18, -80.27),
    lat         = c( 34.00,  34.00,  34.00,  34.00)
  )
  net <- suppressMessages(generate_distance_network(schools, max_distance_km = 50))
  expect_valid_network(net, 4)
})

test_that("generate_distance_network: very small max_distance produces fewer edges", {
  schools <- data.frame(
    school_id = 1:3,
    lon = c(-80.00, -80.50, -81.00),
    lat = c( 34.00,  34.00,  34.00)
  )
  net_close  <- suppressMessages(generate_distance_network(schools, max_distance_km = 10))
  net_far    <- suppressMessages(generate_distance_network(schools, max_distance_km = 200))
  expect_lte(net_close$n_edges, net_far$n_edges)
})

test_that("generate_distance_network: weight_method argument is validated", {
  schools <- data.frame(school_id = 1:2, lon = c(-80, -80.1), lat = c(34, 34))
  expect_error(suppressMessages(
    generate_distance_network(schools, weight_method = "bad_method")
  ))
})
