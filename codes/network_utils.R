# ==============================================================================
# NETWORK GENERATION UTILITIES
# ==============================================================================
# File: network_utils.R
# Contains: Functions for generating school networks
# Dependencies: igraph
# ==============================================================================

library(igraph)

# ==============================================================================
# Distance-Based Network (from schools dataframe with coordinates)
# ==============================================================================

#' Generate network based on geographic distance between schools
#' @param schools Data frame with school_id, lon, lat columns
#' @param max_distance_km Maximum distance for connection (kilometers)
#' @param weight_method Method for calculating edge weights
#' @return A list with adjacency matrix and igraph object
generate_distance_network <- function(schools, 
                                       max_distance_km = 50,
                                       weight_method = c("inverse_linear", "inverse_squared", "exponential")) {
  
  weight_method <- match.arg(weight_method)
  
  n_schools <- nrow(schools)
  cat(sprintf("Creating distance-based network for %d schools...\n", n_schools))
  
  # Haversine distance function
  haversine <- function(lon1, lat1, lon2, lat2) {
    R <- 6371  # Earth radius in km
    lon1_r <- lon1 * pi / 180
    lat1_r <- lat1 * pi / 180
    lon2_r <- lon2 * pi / 180
    lat2_r <- lat2 * pi / 180
    dlat <- lat2_r - lat1_r
    dlon <- lon2_r - lon1_r
    a <- sin(dlat/2)^2 + cos(lat1_r) * cos(lat2_r) * sin(dlon/2)^2
    c <- 2 * atan2(sqrt(a), sqrt(1-a))
    return(R * c)
  }
  
  # Create adjacency matrix
  adj_matrix <- matrix(0, nrow = n_schools, ncol = n_schools)
  
  for (i in 1:n_schools) {
    for (j in 1:n_schools) {
      if (i != j) {
        dist_km <- haversine(schools$lon[i], schools$lat[i], 
                             schools$lon[j], schools$lat[j])
        
        if (!is.na(dist_km) && dist_km > 0 && dist_km <= max_distance_km) {
          if (weight_method == "inverse_linear") {
            weight <- 1 - (dist_km / max_distance_km)
          } else if (weight_method == "inverse_squared") {
            weight <- 1 - (dist_km / max_distance_km)^2
          } else if (weight_method == "exponential") {
            weight <- exp(-dist_km / (max_distance_km / 3))
          } else {
            weight <- 1 - (dist_km / max_distance_km)
          }
          
          weight <- max(0.01, min(1, weight))
          adj_matrix[i, j] <- weight
        }
      }
    }
  }
  
  # Make symmetric
  adj_matrix <- pmax(adj_matrix, t(adj_matrix))
  
  # Create igraph object
  if (sum(adj_matrix) > 0) {
    g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
  } else {
    g <- make_empty_graph(n_schools, directed = FALSE)
  }
  
  # Add school names as vertex attribute
  if ("school_name" %in% names(schools)) {
    V(g)$name <- schools$school_name
  }
  V(g)$school_id <- schools$school_id
  
  cat(sprintf("Network created: %d edges (max distance: %d km)\n", ecount(g), max_distance_km))
  
  return(list(
    graph = g,
    adjacency = adj_matrix,
    n_schools = n_schools,
    n_edges = ecount(g),
    max_distance_km = max_distance_km,
    weight_method = weight_method
  ))
}


# ==============================================================================
# Travel Time-Based Network
# ==============================================================================

#' Generate network based on travel time between schools
#' @param travel_time_file Path to CSV with travel time matrix
#' @param school_reference_file Path to CSV with school ID to name mapping
#' @param selected_school_ids Vector of school IDs to include
#' @param max_travel_time Maximum travel time for connection (minutes)
#' @param weight_method Method for calculating edge weights
#' @return A list with adjacency matrix, igraph object, and school mapping
generate_travel_time_network <- function(travel_time_file, 
                                         school_reference_file,
                                         selected_school_ids = NULL,
                                         max_travel_time = 30,
                                         weight_method = c("inverse_linear", "inverse_squared", "exponential")) {
  
  weight_method <- match.arg(weight_method)
  
  # Load travel time matrix
  travel_matrix <- read.csv(travel_time_file, check.names = FALSE)
  
  # Load school reference
  school_ref <- read.csv(school_reference_file)
  
  # Handle different ID column names (ID or OBJECTID)
  if ("OBJECTID" %in% names(school_ref) && !"ID" %in% names(school_ref)) {
    school_ref$ID <- school_ref$OBJECTID
  }
  
  # Extract origin IDs from first column
  origin_ids <- travel_matrix[[1]]
  travel_values <- as.matrix(travel_matrix[, -1])
  rownames(travel_values) <- origin_ids
  
  # If no schools selected, use all
  if (is.null(selected_school_ids)) {
    selected_school_ids <- origin_ids
  }
  
  # Convert to character for matching
  selected_school_ids <- as.character(selected_school_ids)
  
  # Check which selected IDs are in the matrix
  available_ids <- intersect(selected_school_ids, as.character(origin_ids))
  
  if (length(available_ids) == 0) {
    stop("None of the selected school IDs are in the travel time matrix!")
  }
  
  if (length(available_ids) < length(selected_school_ids)) {
    missing_ids <- setdiff(selected_school_ids, available_ids)
    warning(paste("Some school IDs not found in matrix:", paste(missing_ids, collapse = ", ")))
  }
  
  # Subset the matrix
  row_idx <- match(available_ids, as.character(origin_ids))
  col_idx <- match(available_ids, colnames(travel_values))
  
  travel_subset <- travel_values[row_idx, col_idx]
  n_schools <- length(available_ids)
  
  # Create adjacency matrix based on travel time threshold
  adj_matrix <- matrix(0, nrow = n_schools, ncol = n_schools)
  rownames(adj_matrix) <- available_ids
  colnames(adj_matrix) <- available_ids
  
  for (i in 1:n_schools) {
    for (j in 1:n_schools) {
      if (i != j) {
        travel_time <- travel_subset[i, j]
        
        if (!is.na(travel_time) && travel_time > 0 && travel_time <= max_travel_time) {
          
          if (weight_method == "inverse_linear") {
            weight <- 1 - (travel_time / max_travel_time)
          } else if (weight_method == "inverse_squared") {
            weight <- 1 - (travel_time / max_travel_time)^2
          } else if (weight_method == "exponential") {
            weight <- exp(-travel_time / (max_travel_time / 2))
          } else {
            weight <- 1 - (travel_time / max_travel_time)
          }
          
          weight <- max(0.01, min(1, weight))
          adj_matrix[i, j] <- weight
        }
      }
    }
  }
  
  # Make matrix symmetric by taking maximum of (i,j) and (j,i)
  # This handles cases where travel time A->B differs from B->A
  adj_matrix <- pmax(adj_matrix, t(adj_matrix))
  
  # Create igraph object
  if (sum(adj_matrix) > 0) {
    g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
  } else {
    g <- make_empty_graph(n_schools, directed = FALSE)
  }
  
  # Get school names for the selected IDs
  school_names <- sapply(as.numeric(available_ids), function(id) {
    name <- school_ref$School.Name[school_ref$ID == id]
    if (length(name) == 0) return(paste0("School_", id))
    return(name)
  })
  
  # Add school names as vertex attribute
  V(g)$name <- school_names
  V(g)$school_id <- as.numeric(available_ids)
  
  # Create mapping dataframe
  school_mapping <- data.frame(
    matrix_id = as.numeric(available_ids),
    school_name = school_names,
    network_idx = 1:n_schools,
    stringsAsFactors = FALSE
  )
  
  return(list(
    graph = g,
    adjacency = adj_matrix,
    n_schools = n_schools,
    n_edges = ecount(g),
    school_mapping = school_mapping,
    travel_time_subset = travel_subset,
    max_travel_time = max_travel_time,
    weight_method = weight_method
  ))
}

cat("Network utilities loaded.\n")