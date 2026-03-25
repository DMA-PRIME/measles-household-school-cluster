# ==============================================================================
# UNIT TESTS: network_utils.R
# ==============================================================================
# Tests for: generate_random_network, generate_smallworld_network,
#            generate_scalefree_network, generate_custom_network,
#            generate_distance_network
# ==============================================================================

project_root <- Sys.getenv("PROJECT_ROOT", unset = normalizePath("."))
source(file.path(project_root, "network_utils.R"))

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
# generate_random_network
# ==============================================================================

test_that("generate_random_network: returns valid network for 10 schools", {
  set.seed(1)
  net <- generate_random_network(10, p_edge = 0.5)
  expect_valid_network(net, 10)
})

test_that("generate_random_network: p_edge=0 produces no edges", {
  set.seed(1)
  net <- generate_random_network(5, p_edge = 0)
  expect_equal(net$n_edges, 0)
  expect_equal(sum(net$adjacency), 0)
})

test_that("generate_random_network: p_edge=1 produces fully connected graph", {
  set.seed(1)
  n <- 6
  net <- generate_random_network(n, p_edge = 1)
  # max possible edges for undirected, no self-loops = n*(n-1)/2
  expect_equal(net$n_edges, n * (n - 1) / 2)
})

test_that("generate_random_network: adjacency is symmetric", {
  set.seed(42)
  net <- generate_random_network(8, p_edge = 0.5)
  expect_equal(net$adjacency, t(net$adjacency))
})

test_that("generate_random_network: edge weights are in (0, 1]", {
  set.seed(42)
  net <- generate_random_network(8, p_edge = 0.8)
  non_zero <- net$adjacency[net$adjacency > 0]
  if (length(non_zero) > 0) {
    expect_true(all(non_zero > 0))
    expect_true(all(non_zero <= 1))
  }
})

# ==============================================================================
# generate_smallworld_network
# ==============================================================================

test_that("generate_smallworld_network: returns valid network", {
  set.seed(1)
  net <- generate_smallworld_network(10, nei = 2, p_rewire = 0.1)
  expect_valid_network(net, 10)
})

test_that("generate_smallworld_network: adjacency is symmetric", {
  set.seed(1)
  net <- generate_smallworld_network(10, nei = 2, p_rewire = 0.0)
  expect_equal(net$adjacency, t(net$adjacency))
})

test_that("generate_smallworld_network: handles small n_schools gracefully", {
  set.seed(1)
  net <- generate_smallworld_network(4, nei = 2, p_rewire = 0.1)
  expect_valid_network(net, 4)
})

# ==============================================================================
# generate_scalefree_network
# ==============================================================================

test_that("generate_scalefree_network: returns valid network", {
  set.seed(1)
  net <- generate_scalefree_network(10, m = 2)
  expect_valid_network(net, 10)
})

test_that("generate_scalefree_network: adjacency is symmetric", {
  set.seed(1)
  net <- generate_scalefree_network(10, m = 2)
  expect_equal(net$adjacency, t(net$adjacency))
})

test_that("generate_scalefree_network: m capped at n_schools-1", {
  set.seed(1)
  # m=100 for a 5-node graph should be capped to 4
  net <- generate_scalefree_network(5, m = 100)
  expect_valid_network(net, 5)
})

# ==============================================================================
# generate_custom_network
# ==============================================================================

test_that("generate_custom_network: returns network matching supplied matrix", {
  adj <- matrix(c(0, 0.5, 0.3,
                  0.5, 0,  0.8,
                  0.3, 0.8, 0), nrow = 3, ncol = 3)
  net <- generate_custom_network(adj)
  expect_valid_network(net, 3)
  expect_equal(net$adjacency, adj)
})

test_that("generate_custom_network: all-zero matrix creates zero-edge network", {
  adj <- matrix(0, nrow = 4, ncol = 4)
  net <- generate_custom_network(adj)
  expect_equal(net$n_edges, 0)
})

test_that("generate_custom_network: n_schools equals matrix dimension", {
  adj <- matrix(0, nrow = 5, ncol = 5)
  adj[1, 2] <- 0.5; adj[2, 1] <- 0.5
  net <- generate_custom_network(adj)
  expect_equal(net$n_schools, 5)
})

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
