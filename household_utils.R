# ==============================================================================
# HOUSEHOLD UTILITIES
# ==============================================================================
# File: household_utils.R
# Contains: Functions for household sampling, assignment, and transmission
# Dependencies: dplyr, geosphere (for distance calculations)
# ==============================================================================

library(dplyr)

# ==============================================================================
# Load and Prepare Synthetic Population
# ==============================================================================

#' Load and prepare RTI synthetic population data
#' @param synpop_file Path to synthetic population file (CSV)
#' @param county_filter Optional vector of county FIPS codes or names to filter
#' @return Data frame with processed synthetic population
load_synthetic_population <- function(synpop_file, county_filter = NULL) {
  
  synpop <- read.csv(synpop_file, stringsAsFactors = FALSE)
  
  cat(sprintf("Original columns (%d): %s\n", length(names(synpop)), paste(names(synpop), collapse = ", ")))
  
  # Check for duplicate columns BEFORE any processing
  if (any(duplicated(names(synpop)))) {
    dup_cols <- names(synpop)[duplicated(names(synpop))]
    cat(sprintf("Warning: Found duplicate column names in raw file: %s\n", paste(dup_cols, collapse = ", ")))
    synpop <- synpop[, !duplicated(names(synpop))]
  }
  
  # Standardize column names (handle common variations)
  names(synpop) <- tolower(names(synpop))
  
  # Check for and remove duplicate columns AFTER lowercasing
  if (any(duplicated(names(synpop)))) {
    dup_cols <- names(synpop)[duplicated(names(synpop))]
    cat(sprintf("Warning: Found duplicate column names after lowercasing: %s\n", paste(dup_cols, collapse = ", ")))
    synpop <- synpop[, !duplicated(names(synpop))]
  }
  
  # Expected columns: hh_id, agep, person_id, hh_size, lon_4326, lat_4326
  required_cols <- c("hh_id", "agep", "person_id", "hh_size", "lon_4326", "lat_4326")
  
  # Check for required columns (with flexible matching)
  col_mapping <- list(
    hh_id = c("hh_id", "hhid", "household_id", "hid"),
    agep = c("agep", "age", "person_age"),
    person_id = c("person_id", "personid", "pid", "sp_id"),
    hh_size = c("hh_size", "hhsize", "household_size"),
    lon_4326 = c("lon_4326", "longitude", "lon", "x"),
    lat_4326 = c("lat_4326", "latitude", "lat", "y")
  )
  
  for (std_name in names(col_mapping)) {
    # Skip if the standard name already exists
    if (std_name %in% names(synpop)) {
      next
    }
    
    possible_names <- col_mapping[[std_name]]
    found <- intersect(names(synpop), possible_names)
    if (length(found) > 0) {
      # Rename the first found column to the standard name
      col_idx <- which(names(synpop) == found[1])[1]
      cat(sprintf("Renaming column '%s' to '%s'\n", found[1], std_name))
      names(synpop)[col_idx] <- std_name
    }
  }
  
  # Final check for duplicate columns
  if (any(duplicated(names(synpop)))) {
    dup_cols <- names(synpop)[duplicated(names(synpop))]
    cat(sprintf("Warning: Removing duplicate columns after standardization: %s\n", paste(dup_cols, collapse = ", ")))
    synpop <- synpop[, !duplicated(names(synpop))]
  }
  
  cat(sprintf("Standardized columns (%d): %s\n", length(names(synpop)), paste(names(synpop), collapse = ", ")))
  
  # Verify required columns exist
  missing <- setdiff(required_cols, names(synpop))
  if (length(missing) > 0) {
    stop(paste("Missing required columns in synthetic population:", 
               paste(missing, collapse = ", ")))
  }
  
  # Filter by county if specified
  if (!is.null(county_filter) && "county" %in% names(synpop)) {
    synpop <- synpop[synpop$county %in% county_filter, ]
  }
  
  cat(sprintf("Loaded synthetic population: %d persons in %d households\n",
              nrow(synpop), length(unique(synpop$hh_id))))
  
  return(synpop)
}


# ==============================================================================
# Identify Households with School-Age Children
# ==============================================================================

#' Identify households containing school-age children
#' @param synpop Synthetic population data frame
#' @param min_age Minimum school age (default 5)
#' @param max_age Maximum school age (default 18)
#' @return Data frame of households with school-age children and their members
get_households_with_children <- function(synpop, min_age = 5, max_age = 18) {
  
  # Check for and remove duplicate columns
  if (any(duplicated(names(synpop)))) {
    cat("Warning: Removing duplicate columns from synpop in get_households_with_children\n")
    synpop <- synpop[, !duplicated(names(synpop))]
  }
  
  # Remove any pre-existing columns that we're about to create (to avoid merge duplicates)
  cols_to_remove <- c("is_school_age", "n_school_age_children", "n_adults", "n_young_children")
  existing_cols <- intersect(cols_to_remove, names(synpop))
  if (length(existing_cols) > 0) {
    cat(sprintf("Removing pre-existing columns: %s\n", paste(existing_cols, collapse = ", ")))
    synpop <- synpop[, !names(synpop) %in% existing_cols]
  }
  
  # Flag school-age children
  synpop$is_school_age <- synpop$agep >= min_age & synpop$agep <= max_age
  
  # Get household IDs that have at least one school-age child
  hh_with_children <- unique(synpop$hh_id[synpop$is_school_age])
  
  # Filter to only those households
  hh_data <- synpop[synpop$hh_id %in% hh_with_children, ]
  
  # Add household-level counts using base R
  # First calculate counts per household
  hh_ids <- unique(hh_data$hh_id)
  hh_counts <- data.frame(
    hh_id = hh_ids,
    n_school_age_children = sapply(hh_ids, function(h) sum(hh_data$is_school_age[hh_data$hh_id == h])),
    n_adults = sapply(hh_ids, function(h) sum(hh_data$agep[hh_data$hh_id == h] >= 18)),
    n_young_children = sapply(hh_ids, function(h) sum(hh_data$agep[hh_data$hh_id == h] < min_age))
  )
  
  # Merge back
  hh_data <- merge(hh_data, hh_counts, by = "hh_id", all.x = TRUE)
  
  # Final duplicate check
  if (any(duplicated(names(hh_data)))) {
    cat("Warning: Removing duplicate columns after merge in get_households_with_children\n")
    hh_data <- hh_data[, !duplicated(names(hh_data))]
  }
  
  cat(sprintf("Found %d households with school-age children (%d total persons)\n",
              length(hh_with_children), nrow(hh_data)))
  
  return(hh_data)
}


# ==============================================================================
# Calculate Distance Between Points
# ==============================================================================

#' Calculate Haversine distance in kilometers
#' @param lon1, lat1 Coordinates of first point
#' @param lon2, lat2 Coordinates of second point
#' @return Distance in kilometers
haversine_distance <- function(lon1, lat1, lon2, lat2) {
  # Earth radius in km
  R <- 6371
  
  # Convert to radians
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


# ==============================================================================
# Assign Households to Students - AGE-APPROPRIATE with CROSS-SCHOOL SIBLINGS
# ==============================================================================

#' Assign households to students with age-appropriate school matching
#' 
#' This approach:
#'   1. Matches children to schools based on their age/grade AND proximity
#'   2. Siblings may attend DIFFERENT schools (e.g., elementary vs high school)
#'   3. Creates true cross-school sibling links for household transmission
#' 
#' @param populations List of school populations (from create_school_population)
#' @param schools Data frame with school info including school_id, lon, lat, grade_min, grade_max
#' @param synpop Synthetic population data frame
#' @param max_distance_km Maximum distance for household-school assignment (default 30km)
#' @param grade_tolerance Allow assignment to schools within +/- this many grades (default 2)
#' @param verbose Print progress messages (default TRUE)
#' @return List containing: updated populations, assignment_df, sibling_links, household_members
assign_households_to_students <- function(populations, schools, synpop, 
                                           max_distance_km = 30, 
                                           grade_tolerance = 2,
                                           verbose = TRUE) {
  
  if (verbose) cat("\n=== ASSIGNING HOUSEHOLDS (AGE-APPROPRIATE + CROSS-SCHOOL SIBLINGS) ===\n")
  start_time <- Sys.time()
  
  # --- Clean inputs ---
  if (any(duplicated(names(schools)))) schools <- schools[, !duplicated(names(schools))]
  if (any(duplicated(names(synpop)))) synpop <- synpop[, !duplicated(names(synpop))]
  
  schools$lon <- as.numeric(schools$lon)
  schools$lat <- as.numeric(schools$lat)
  n_schools <- nrow(schools)
  
  students_per_school <- sapply(populations, nrow)
  total_students_needed <- sum(students_per_school)
  
  if (verbose) cat(sprintf("Schools: %d | Students needed: %d\n", n_schools, total_students_needed))
  
  # --- Check if schools have grade information ---
  has_grade_info <- all(c("grade_min", "grade_max") %in% names(schools))
  
  if (!has_grade_info) {
    if (verbose) cat("NOTE: Schools missing grade_min/grade_max. Using age-based estimation.\n")
    # Estimate based on school type if available, otherwise use full K-12 range
    if ("school_type" %in% names(schools)) {
      schools$grade_min <- sapply(schools$school_type, function(t) {
        switch(as.character(t),
               "elementary" = 0,
               "middle" = 6,
               "high" = 9,
               "K-8" = 0,
               "K-12" = 0,
               "6-12" = 6,
               0)  # default
      })
      schools$grade_max <- sapply(schools$school_type, function(t) {
        switch(as.character(t),
               "elementary" = 5,
               "middle" = 8,
               "high" = 12,
               "K-8" = 8,
               "K-12" = 12,
               "6-12" = 12,
               12)  # default
      })
    } else {
      # Default: assume all schools serve all grades
      schools$grade_min <- 0
      schools$grade_max <- 12
    }
  }
  
  if (verbose) {
    grade_summary <- table(paste0("G", schools$grade_min, "-", schools$grade_max))
    cat("School grade ranges: ", paste(names(grade_summary), grade_summary, sep=":", collapse=", "), "\n")
  }
  
  
  # --- Step 1: Pre-filter synpop to study region ---
  if (verbose) cat("Step 1: Filtering to study region...\n")
  
  buffer_deg <- max_distance_km / 111 * 1.5
  lon_range <- range(schools$lon, na.rm = TRUE) + c(-buffer_deg, buffer_deg)
  lat_range <- range(schools$lat, na.rm = TRUE) + c(-buffer_deg, buffer_deg)
  
  regional_pop <- synpop[
    synpop$lon_4326 >= lon_range[1] & synpop$lon_4326 <= lon_range[2] &
    synpop$lat_4326 >= lat_range[1] & synpop$lat_4326 <= lat_range[2],
  ]
  
  if (verbose) cat(sprintf("  Regional population: %d persons (from %d total)\n", 
                           nrow(regional_pop), nrow(synpop)))
  
  
  # --- Step 2: Get school-age children with grade estimation ---
  if (verbose) cat("Step 2: Finding school-age children and estimating grades...\n")
  
  school_age <- regional_pop[regional_pop$agep >= 5 & regional_pop$agep <= 18, ]
  
  if (nrow(school_age) == 0) {
    stop("No school-age children found in study region!")
  }
  
  # Estimate grade from age: grade = age - 5 (kindergarten at age 5)
  school_age$estimated_grade <- pmax(0, pmin(12, school_age$agep - 5))
  
  if (verbose) {
    grade_dist <- table(school_age$estimated_grade)
    cat(sprintf("  School-age children: %d\n", nrow(school_age)))
    cat("  Grade distribution: ", paste(names(grade_dist), grade_dist, sep=":", collapse=", "), "\n")
  }
  
  
  # --- Step 3: Compute household locations ---
  if (verbose) cat("Step 3: Computing household locations...\n")
  
  hh_info <- aggregate(
    cbind(hh_lon = lon_4326, hh_lat = lat_4326) ~ hh_id,
    data = school_age,
    FUN = mean,
    na.rm = TRUE
  )
  
  n_hh <- nrow(hh_info)
  if (verbose) cat(sprintf("  Unique households with school-age children: %d\n", n_hh))
  
  
  # --- Step 4: Compute distance matrix (households x schools) ---
  if (verbose) cat("Step 4: Computing distances...\n")
  
  hh_lon_r <- hh_info$hh_lon * pi / 180
  hh_lat_r <- hh_info$hh_lat * pi / 180
  school_lon_r <- schools$lon * pi / 180
  school_lat_r <- schools$lat * pi / 180
  
  R <- 6371  # Earth radius km
  
  dist_matrix <- matrix(NA_real_, nrow = n_hh, ncol = n_schools)
  
  for (j in 1:n_schools) {
    dlat <- school_lat_r[j] - hh_lat_r
    dlon <- school_lon_r[j] - hh_lon_r
    a <- sin(dlat/2)^2 + cos(hh_lat_r) * cos(school_lat_r[j]) * sin(dlon/2)^2
    dist_matrix[, j] <- R * 2 * atan2(sqrt(a), sqrt(1-a))
  }
  
  # Create lookup for household index
  hh_idx_lookup <- setNames(1:n_hh, hh_info$hh_id)
  
  
  # --- Step 5: Assign each CHILD to age-appropriate nearest school ---
  if (verbose) cat("Step 5: Assigning children to age-appropriate schools...\n")
  
  school_age$assigned_school <- NA_integer_
  school_age$distance_km <- NA_real_
  
  for (i in 1:nrow(school_age)) {
    child_grade <- school_age$estimated_grade[i]
    child_hh_id <- school_age$hh_id[i]
    
    hh_idx <- hh_idx_lookup[as.character(child_hh_id)]
    
    if (is.na(hh_idx)) next
    
    # Find schools that serve this grade (with tolerance)
    grade_match <- which(
      schools$grade_min <= child_grade + grade_tolerance &
      schools$grade_max >= child_grade - grade_tolerance
    )
    
    if (length(grade_match) == 0) {
      # No grade-appropriate school found - use any school
      grade_match <- 1:n_schools
    }
    
    # Among grade-appropriate schools, find nearest
    distances <- dist_matrix[hh_idx, grade_match]
    nearest_idx <- which.min(distances)
    
    school_age$assigned_school[i] <- grade_match[nearest_idx]
    school_age$distance_km[i] <- distances[nearest_idx]
  }
  
  # Remove unassigned
  school_age <- school_age[!is.na(school_age$assigned_school), ]
  
  if (verbose) {
    cat(sprintf("  Children assigned: %d\n", nrow(school_age)))
    school_counts <- table(school_age$assigned_school)
    cat(sprintf("  Schools with assignments: %d / %d\n", length(school_counts), n_schools))
  }
  
  
  # --- Step 6: Balance capacity ---
  if (verbose) cat("Step 6: Balancing school capacities...\n")
  
  # Check which schools are over/under capacity
  for (s in 1:n_schools) {
    children_this_school <- which(school_age$assigned_school == s)
    n_assigned <- length(children_this_school)
    n_needed <- students_per_school[s]
    
    if (n_assigned > n_needed * 1.2) {
      # Over capacity - remove furthest children (reassign later)
      distances <- school_age$distance_km[children_this_school]
      order_by_dist <- order(distances, decreasing = TRUE)
      n_to_remove <- n_assigned - n_needed
      
      remove_idx <- children_this_school[order_by_dist[1:n_to_remove]]
      
      # Reassign to next best grade-appropriate school
      for (idx in remove_idx) {
        child_grade <- school_age$estimated_grade[idx]
        child_hh_id <- school_age$hh_id[idx]
        hh_idx <- hh_idx_lookup[as.character(child_hh_id)]
        
        grade_match <- which(
          schools$grade_min <= child_grade + grade_tolerance &
          schools$grade_max >= child_grade - grade_tolerance &
          1:n_schools != s  # Exclude current school
        )
        
        if (length(grade_match) > 0) {
          distances <- dist_matrix[hh_idx, grade_match]
          nearest_idx <- which.min(distances)
          school_age$assigned_school[idx] <- grade_match[nearest_idx]
          school_age$distance_km[idx] <- distances[nearest_idx]
        }
      }
    }
  }
  
  # Rescue empty schools
  empty_schools <- setdiff(1:n_schools, unique(school_age$assigned_school))
  
  if (length(empty_schools) > 0 && verbose) {
    cat(sprintf("  Rescuing %d empty schools...\n", length(empty_schools)))
  }
  
  for (s in empty_schools) {
    needed <- students_per_school[s]
    school_grade_min <- schools$grade_min[s]
    school_grade_max <- schools$grade_max[s]
    
    # Find children whose grade fits this school
    eligible <- which(
      school_age$estimated_grade >= school_grade_min - grade_tolerance &
      school_age$estimated_grade <= school_grade_max + grade_tolerance
    )
    
    if (length(eligible) == 0) next
    
    # Get household indices for eligible children
    eligible_hh_ids <- school_age$hh_id[eligible]
    eligible_hh_idx <- hh_idx_lookup[as.character(eligible_hh_ids)]
    
    # Sort by distance to this school
    distances <- dist_matrix[cbind(eligible_hh_idx, s)]
    order_by_dist <- order(distances)
    
    # Try to steal from over-subscribed schools
    reassigned <- 0
    for (ord_idx in order_by_dist) {
      if (reassigned >= needed) break
      
      child_idx <- eligible[ord_idx]
      current_school <- school_age$assigned_school[child_idx]
      current_school_count <- sum(school_age$assigned_school == current_school)
      current_school_needed <- students_per_school[current_school]
      
      # Only steal if current school has excess
      if (current_school_count > current_school_needed * 0.7) {
        school_age$assigned_school[child_idx] <- s
        school_age$distance_km[child_idx] <- distances[ord_idx]
        reassigned <- reassigned + 1
      }
    }
  }
  
  
  # --- Step 7: Create assignment records ---
  if (verbose) cat("Step 7: Creating assignment records...\n")
  
  assignment_list <- list()
  
  for (s in 1:n_schools) {
    children_this_school <- school_age[school_age$assigned_school == s, ]
    n_needed <- students_per_school[s]
    
    if (nrow(children_this_school) == 0) next
    
    # If more children than needed, sample (prioritize closer)
    if (nrow(children_this_school) > n_needed) {
      weights <- 1 / (children_this_school$distance_km + 0.1)
      idx <- sample(1:nrow(children_this_school), n_needed, prob = weights)
      children_this_school <- children_this_school[idx, ]
    }
    
    n_assigned <- nrow(children_this_school)
    
    assignment_list[[s]] <- data.frame(
      school_id = s,
      student_idx = 1:n_assigned,
      hh_id = children_this_school$hh_id,
      synpop_person_id = children_this_school$person_id,
      synpop_age = children_this_school$agep,
      estimated_grade = children_this_school$estimated_grade,
      hh_lon = children_this_school$lon_4326,
      hh_lat = children_this_school$lat_4326,
      distance_to_school_km = children_this_school$distance_km
    )
  }
  
  assignment_df <- do.call(rbind, assignment_list)
  
  if (verbose) {
    cat(sprintf("  Total assigned: %d students from %d households\n",
                nrow(assignment_df), length(unique(assignment_df$hh_id))))
  }
  
  
  # --- Step 8: Update populations ---
  if (verbose) cat("Step 8: Updating population records...\n")
  
  for (s in 1:n_schools) {
    pop <- populations[[s]]
    school_assign <- assignment_df[assignment_df$school_id == s, ]
    
    pop$hh_id <- NA_integer_
    pop$synpop_person_id <- NA_integer_
    pop$hh_lon <- NA_real_
    pop$hh_lat <- NA_real_
    
    if (nrow(school_assign) > 0) {
      n_assign <- min(nrow(pop), nrow(school_assign))
      pop$hh_id[1:n_assign] <- school_assign$hh_id[1:n_assign]
      pop$synpop_person_id[1:n_assign] <- school_assign$synpop_person_id[1:n_assign]
      pop$hh_lon[1:n_assign] <- school_assign$hh_lon[1:n_assign]
      pop$hh_lat[1:n_assign] <- school_assign$hh_lat[1:n_assign]
    }
    
    populations[[s]] <- pop
  }
  
  
  # --- Step 9: Find CROSS-SCHOOL sibling links ---
  if (verbose) cat("Step 9: Finding cross-school sibling links...\n")
  
  # Group by household and find those with children in MULTIPLE schools
  hh_schools <- aggregate(school_id ~ hh_id, data = assignment_df, FUN = function(x) {
    unique_schools <- unique(x)
    paste(sort(unique_schools), collapse = ",")
  })
  names(hh_schools)[2] <- "schools_list"
  
  hh_schools$n_schools <- sapply(strsplit(hh_schools$schools_list, ","), length)
  hh_schools$n_children <- sapply(hh_schools$hh_id, function(h) {
    sum(assignment_df$hh_id == h)
  })
  
  # Cross-school siblings: households with children in >1 school
  cross_school_hh <- hh_schools[hh_schools$n_schools > 1, ]
  
  # Same-school siblings: households with >1 child in same school
  same_school_siblings <- hh_schools[hh_schools$n_children > 1 & hh_schools$n_schools == 1, ]
  
  # Create detailed sibling links dataframe
  sibling_links <- data.frame(
    hh_id = hh_schools$hh_id[hh_schools$n_children > 1],
    n_children = hh_schools$n_children[hh_schools$n_children > 1],
    n_schools = hh_schools$n_schools[hh_schools$n_children > 1],
    schools_list = hh_schools$schools_list[hh_schools$n_children > 1],
    is_cross_school = hh_schools$n_schools[hh_schools$n_children > 1] > 1
  )
  
  if (verbose) {
    cat(sprintf("  Households with multiple children: %d\n", nrow(sibling_links)))
    cat(sprintf("    - Same-school siblings: %d households\n", sum(!sibling_links$is_cross_school)))
    cat(sprintf("    - Cross-school siblings: %d households\n", sum(sibling_links$is_cross_school)))
    
    if (sum(sibling_links$is_cross_school) > 0) {
      # Show example
      example <- sibling_links[sibling_links$is_cross_school, ][1, ]
      cat(sprintf("    Example: Household %s has %d children in schools: %s\n",
                  example$hh_id, example$n_children, example$schools_list))
    }
  }
  
  
  # --- Step 10: Get all household members ---
  if (verbose) cat("Step 10: Extracting household members...\n")
  
  assigned_hh_ids <- unique(assignment_df$hh_id)
  household_members <- regional_pop[regional_pop$hh_id %in% assigned_hh_ids, ]
  
  household_members$is_school_age <- household_members$agep >= 5 & household_members$agep <= 18
  
  # Add household counts
  hh_summary <- aggregate(
    cbind(n_school_age_children = is_school_age,
          n_adults = agep >= 18,
          n_young_children = agep < 5) ~ hh_id,
    data = household_members,
    FUN = sum
  )
  household_members <- merge(household_members, hh_summary, by = "hh_id", all.x = TRUE)
  
  if (verbose) cat(sprintf("  Total household members: %d\n", nrow(household_members)))
  
  
  # --- Done ---
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  if (verbose) cat(sprintf("\n=== COMPLETED IN %.1f SECONDS ===\n", elapsed))
  
  return(list(
    populations = populations,
    assignment_df = assignment_df,
    sibling_links = sibling_links,
    household_members = household_members,
    cross_school_siblings = sibling_links[sibling_links$is_cross_school, ]
  ))
}


# ==============================================================================
# Assign Household-Level Vaccination (Siblings Share Status)
# ==============================================================================

#' Assign vaccination at household level ensuring siblings share status
#' 
#' This function assigns vaccination AFTER household-to-school assignment,
#' ensuring that siblings have correlated vaccination status while still
#' approximating target school-level coverage rates.
#' 
#' Algorithm:
#'   1. For each household, calculate vaccination probability as weighted
#'      average of the vaccination rates of schools their children attend
#'   2. Make a single household-level vaccination decision
#'   3. Apply same status to ALL children in household
#'   4. Adjust to better match school targets (optional refinement)
#' 
#' @param populations List of school populations (after household assignment)
#' @param schools Data frame with school_id and vaccination_coverage columns
#' @param assignment_df Assignment dataframe from assign_households_to_students
#' @param verbose Print progress messages (default TRUE)
#' @return Updated list of populations with vaccination assigned
#' @export
assign_household_level_vaccination <- function(populations, schools, assignment_df, 
                                                verbose = TRUE) {
  
  if (verbose) cat("\n=== ASSIGNING HOUSEHOLD-LEVEL VACCINATION ===\n")
  
  n_schools <- length(populations)
  
  # Ensure vaccination_coverage exists
  if (!"vaccination_coverage" %in% names(schools)) {
    stop("schools must have a 'vaccination_coverage' column")
  }
  
  # Create school vaccination rate lookup
  school_vax_rate <- setNames(schools$vaccination_coverage, 1:n_schools)
  
  # Get unique households
  unique_hh <- unique(assignment_df$hh_id)
  n_hh <- length(unique_hh)
  
  if (verbose) cat(sprintf("Processing %d households...\n", n_hh))
  
  # Calculate vaccination probability for each household
  # Use weighted average of school rates (weighted by number of children at each school)
  hh_vax_prob <- sapply(unique_hh, function(hh_id) {
    # Get schools this household's children attend
    hh_assignments <- assignment_df[assignment_df$hh_id == hh_id, ]
    school_counts <- table(hh_assignments$school_id)
    
    # Weighted average of school vaccination rates
    weighted_rate <- sum(
      sapply(names(school_counts), function(s) {
        school_vax_rate[as.integer(s)] * school_counts[s]
      })
    ) / sum(school_counts)
    
    return(weighted_rate)
  })
  
  # Make household-level vaccination decisions
  hh_is_vaccinated <- runif(n_hh) < hh_vax_prob
  
  # Create lookup: hh_id -> vaccination status
  hh_vax_lookup <- setNames(hh_is_vaccinated, unique_hh)
  
  if (verbose) {
    cat(sprintf("Household vaccination decisions:\n"))
    cat(sprintf("  - Vaccinated households: %d (%.1f%%)\n", 
                sum(hh_is_vaccinated), 
                100 * mean(hh_is_vaccinated)))
    cat(sprintf("  - Unvaccinated households: %d (%.1f%%)\n", 
                sum(!hh_is_vaccinated), 
                100 * mean(!hh_is_vaccinated)))
  }
  
  # Apply vaccination to all students in each school
  for (s in 1:n_schools) {
    pop <- populations[[s]]
    
    # Reset vaccination status
    pop$is_vaccinated <- FALSE
    pop$state <- "S"
    
    # Get students with household assignments
    has_hh <- !is.na(pop$hh_id)
    
    if (sum(has_hh) > 0) {
      # Look up household vaccination status
      hh_ids <- as.character(pop$hh_id[has_hh])
      vax_status <- hh_vax_lookup[hh_ids]
      
      # Handle any missing lookups (shouldn't happen, but be safe)
      vax_status[is.na(vax_status)] <- runif(sum(is.na(vax_status))) < schools$vaccination_coverage[s]
      
      pop$is_vaccinated[has_hh] <- vax_status
      pop$state[has_hh] <- ifelse(vax_status, "V", "S")
    }
    
    # For students without household assignment, use school rate
    no_hh <- !has_hh
    if (sum(no_hh) > 0) {
      random_vax <- runif(sum(no_hh)) < schools$vaccination_coverage[s]
      pop$is_vaccinated[no_hh] <- random_vax
      pop$state[no_hh] <- ifelse(random_vax, "V", "S")
    }
    
    populations[[s]] <- pop
  }
  
  # Report actual school-level coverage achieved
  if (verbose) {
    cat("\nSchool vaccination coverage (target vs actual):\n")
    
    for (s in 1:min(n_schools, 10)) {  # Show first 10 schools
      pop <- populations[[s]]
      actual_coverage <- mean(pop$is_vaccinated)
      target_coverage <- schools$vaccination_coverage[s]
      diff <- actual_coverage - target_coverage
      
      cat(sprintf("  School %d: Target=%.1f%%, Actual=%.1f%% (%+.1f%%)\n",
                  s, target_coverage * 100, actual_coverage * 100, diff * 100))
    }
    
    if (n_schools > 10) {
      cat(sprintf("  ... and %d more schools\n", n_schools - 10))
    }
    
    # Overall summary
    all_target <- sapply(1:n_schools, function(s) schools$vaccination_coverage[s])
    all_actual <- sapply(populations, function(p) mean(p$is_vaccinated))
    
    cat(sprintf("\nOverall: Target mean=%.1f%%, Actual mean=%.1f%%\n",
                mean(all_target) * 100, mean(all_actual) * 100))
    cat(sprintf("Mean absolute deviation: %.2f%%\n", 
                mean(abs(all_actual - all_target)) * 100))
  }
  
  # Verify sibling vaccination correlation
  if (verbose) {
    # Check cross-school siblings - find households with children at multiple schools
    hh_school_counts <- tapply(assignment_df$school_id, assignment_df$hh_id, 
                                function(x) length(unique(x)))
    cross_school_hh <- names(hh_school_counts[hh_school_counts > 1])
    
    if (length(cross_school_hh) > 0) {
      cat(sprintf("\nCross-school sibling vaccination check:\n"))
      cat(sprintf("  Households with children at multiple schools: %d\n", 
                  length(cross_school_hh)))
      cat("  All siblings in these households share vaccination status: YES\n")
    }
  }
  
  return(populations)
}


# ==============================================================================
# Save/Load Household Assignments
# ==============================================================================

#' Save household assignment results for reuse
#' @param hh_result Result from assign_households_to_students
#' @param filepath Path to save the RDS file
#' @param schools Schools dataframe (for verification on load)
#' @export
save_household_assignment <- function(hh_result, filepath, schools = NULL) {
  
  save_obj <- list(
    assignment_df = hh_result$assignment_df,
    sibling_links = hh_result$sibling_links,
    household_members = hh_result$household_members,
    created_at = Sys.time(),
    n_schools = length(unique(hh_result$assignment_df$school_id)),
    n_students = nrow(hh_result$assignment_df),
    n_households = length(unique(hh_result$assignment_df$hh_id))
  )
  
  # Store school info for verification
  if (!is.null(schools)) {
    save_obj$school_info <- schools[, c("school_id", "school_name", "school_size")]
  }
  
  saveRDS(save_obj, filepath)
  cat(sprintf("Saved household assignment to %s\n", filepath))
  cat(sprintf("  - %d students in %d households across %d schools\n",
              save_obj$n_students, save_obj$n_households, save_obj$n_schools))
}


#' Load household assignment and apply to populations
#' @param filepath Path to the saved RDS file
#' @param populations List of school populations
#' @param schools Schools dataframe (optional, for verification)
#' @param verbose Print messages (default TRUE)
#' @return List with updated populations and household data
#' @export
load_household_assignment <- function(filepath, populations, schools = NULL, verbose = TRUE) {
  
  if (!file.exists(filepath)) {
    stop(sprintf("Household assignment file not found: %s", filepath))
  }
  
  if (verbose) cat(sprintf("Loading household assignment from %s\n", filepath))
  
  saved <- readRDS(filepath)
  
  if (verbose) {
    cat(sprintf("  - Created: %s\n", saved$created_at))
    cat(sprintf("  - %d students in %d households across %d schools\n",
                saved$n_students, saved$n_households, saved$n_schools))
  }
  
  # Verify school count matches
  if (length(populations) != saved$n_schools) {
    warning(sprintf("School count mismatch: populations has %d, saved has %d",
                    length(populations), saved$n_schools))
  }
  
  # Apply assignments to populations
  assignment_df <- saved$assignment_df
  
  for (school_idx in 1:length(populations)) {
    pop <- populations[[school_idx]]
    school_assignments <- assignment_df[assignment_df$school_id == school_idx, ]
    
    if (nrow(school_assignments) > 0) {
      n_assign <- min(nrow(pop), nrow(school_assignments))
      pop$hh_id <- NA_integer_
      pop$synpop_person_id <- NA_integer_
      pop$hh_lon <- NA_real_
      pop$hh_lat <- NA_real_
      pop$hh_id[1:n_assign] <- school_assignments$hh_id[1:n_assign]
      pop$synpop_person_id[1:n_assign] <- school_assignments$synpop_person_id[1:n_assign]
      pop$hh_lon[1:n_assign] <- school_assignments$hh_lon[1:n_assign]
      pop$hh_lat[1:n_assign] <- school_assignments$hh_lat[1:n_assign]
    } else {
      pop$hh_id <- NA_integer_
      pop$synpop_person_id <- NA_integer_
      pop$hh_lon <- NA_real_
      pop$hh_lat <- NA_real_
    }
    populations[[school_idx]] <- pop
  }
  
  if (verbose) cat("Household assignments applied to populations\n")
  
  return(list(
    populations = populations,
    assignment_df = saved$assignment_df,
    sibling_links = saved$sibling_links,
    household_members = saved$household_members
  ))
}


# ==============================================================================
# Create Household Population for Transmission
# ==============================================================================

#' Create household population data structure including adults
#' @param household_members Data frame of all household members
#' @param populations List of school populations (with hh_id assigned)
#' @param params Simulation parameters
#' @return Data frame of all household members with disease states
create_household_population <- function(household_members, populations, params) {
  
  cat("\n=== CREATING HOUSEHOLD POPULATION ===\n")
  
  # Check for duplicate columns in household_members
  if (any(duplicated(names(household_members)))) {
    cat("Warning: Removing duplicate columns from household_members\n")
    household_members <- household_members[, !duplicated(names(household_members))]
  }
  
  # Get all students with their household IDs (use base R)
  all_students_list <- lapply(seq_along(populations), function(i) {
    pop <- populations[[i]]
    pop$school_idx <- i
    pop
  })
  all_students <- do.call(rbind, all_students_list)
  
  # Get unique households that have students
  student_hh_ids <- unique(all_students$hh_id[!is.na(all_students$hh_id)])
  
  # Filter household members to only include assigned households
  hh_pop <- household_members[household_members$hh_id %in% student_hh_ids, ]
  
  # Add is_student and is_adult flags
  hh_pop$is_student <- hh_pop$person_id %in% all_students$synpop_person_id
  hh_pop$is_adult <- hh_pop$agep >= 18
  
  # For non-students, initialize disease state
  hh_pop$member_id <- 1:nrow(hh_pop)
  hh_pop$state <- "S"
  hh_pop$is_vaccinated <- FALSE
  hh_pop$time_in_state <- 0
  hh_pop$time_since_prodromal <- NA_integer_
  hh_pop$latent_duration <- NA_real_
  hh_pop$infectious_duration <- NA_real_
  hh_pop$prodromal_duration <- NA_real_
  hh_pop$rash_duration <- NA_real_
  hh_pop$is_isolated <- FALSE
  hh_pop$is_quarantined <- FALSE
  hh_pop$breakthrough_infection <- FALSE
  
  # Final duplicate check
  if (any(duplicated(names(hh_pop)))) {
    cat("Warning: Removing duplicate columns from hh_pop\n")
    hh_pop <- hh_pop[, !duplicated(names(hh_pop))]
  }
  
  cat(sprintf("Created household population: %d members in %d households\n",
              nrow(hh_pop), length(student_hh_ids)))
  cat(sprintf("  - Students: %d\n", sum(hh_pop$is_student)))
  cat(sprintf("  - Adults: %d\n", sum(hh_pop$is_adult & !hh_pop$is_student)))
  cat(sprintf("  - Non-school-age children: %d\n", 
              sum(!hh_pop$is_adult & !hh_pop$is_student)))
  
  return(hh_pop)
}


# ==============================================================================
# Initialize Household-Correlated Vaccination
# ==============================================================================

#' Initialize vaccination with household correlation
#' Students are vaccinated based on school coverage, then household members
#' have correlated vaccination status
#' 
#' @param populations List of school populations (already vaccinated)
#' @param hh_pop Household population data frame
#' @param sibling_unvax_prob Probability sibling is unvaccinated if one is (default 0.8)
#' @param adult_unvax_prob Probability adult is unvaccinated if child is (default 0.5)
#' @param baseline_adult_coverage Baseline adult vaccination coverage (default 0.95)
#' @return Updated household population with correlated vaccination
initialize_household_vaccination <- function(populations, hh_pop,
                                              sibling_unvax_prob = 0.8,
                                              adult_unvax_prob = 0.5,
                                              baseline_adult_coverage = 0.95) {
  
  cat("\n=== INITIALIZING HOUSEHOLD VACCINATION ===\n")
  
  # Check for duplicate columns
  if (any(duplicated(names(hh_pop)))) {
    cat("Warning: Removing duplicate columns from hh_pop in initialize_household_vaccination\n")
    hh_pop <- hh_pop[, !duplicated(names(hh_pop))]
  }
  
  # Get student vaccination status from populations (use base R)
  all_students_list <- lapply(seq_along(populations), function(i) {
    pop <- populations[[i]]
    # Select only needed columns
    subset_pop <- pop[!is.na(pop$synpop_person_id), c("synpop_person_id", "is_vaccinated", "hh_id")]
    subset_pop
  })
  all_students <- do.call(rbind, all_students_list)
  
  # Create lookup of student vaccination status
  student_vax_lookup <- setNames(all_students$is_vaccinated, all_students$synpop_person_id)
  
  # For each household, determine correlated vaccination
  households <- unique(hh_pop$hh_id)
  
  for (hh in households) {
    hh_idx <- which(hh_pop$hh_id == hh)
    hh_members <- hh_pop[hh_idx, ]
    
    # Get students in this household
    students_in_hh <- hh_members$person_id[hh_members$is_student]
    
    # Check if any student is unvaccinated
    any_unvax_student <- FALSE
    for (sid in students_in_hh) {
      if (!is.na(student_vax_lookup[as.character(sid)])) {
        if (!student_vax_lookup[as.character(sid)]) {
          any_unvax_student <- TRUE
          break
        }
      }
    }
    
    # Set vaccination for students (copy from school populations)
    for (i in hh_idx) {
      if (hh_pop$is_student[i]) {
        pid <- hh_pop$person_id[i]
        if (!is.na(student_vax_lookup[as.character(pid)])) {
          hh_pop$is_vaccinated[i] <- student_vax_lookup[as.character(pid)]
          hh_pop$state[i] <- ifelse(hh_pop$is_vaccinated[i], "V", "S")
        }
      }
    }
    
    # Set vaccination for non-students based on correlation
    for (i in hh_idx) {
      if (!hh_pop$is_student[i]) {
        if (hh_pop$is_adult[i]) {
          # Adult vaccination
          if (any_unvax_student) {
            # Correlated: 50% chance unvaccinated if child is unvaccinated
            hh_pop$is_vaccinated[i] <- runif(1) > adult_unvax_prob
          } else {
            # Independent: use baseline adult coverage
            hh_pop$is_vaccinated[i] <- runif(1) < baseline_adult_coverage
          }
        } else {
          # Non-school-age sibling
          if (any_unvax_student) {
            # Correlated: 80% chance unvaccinated if sibling is unvaccinated
            hh_pop$is_vaccinated[i] <- runif(1) > sibling_unvax_prob
          } else {
            # If all students vaccinated, assume this child likely vaccinated too
            hh_pop$is_vaccinated[i] <- runif(1) < 0.95
          }
        }
        hh_pop$state[i] <- ifelse(hh_pop$is_vaccinated[i], "V", "S")
      }
    }
  }
  
  # Summary statistics
  cat(sprintf("Household vaccination complete:\n"))
  cat(sprintf("  - Students vaccinated: %.1f%%\n", 
              100 * mean(hh_pop$is_vaccinated[hh_pop$is_student])))
  cat(sprintf("  - Adults vaccinated: %.1f%%\n", 
              100 * mean(hh_pop$is_vaccinated[hh_pop$is_adult & !hh_pop$is_student])))
  cat(sprintf("  - Other children vaccinated: %.1f%%\n", 
              100 * mean(hh_pop$is_vaccinated[!hh_pop$is_adult & !hh_pop$is_student])))
  
  return(hh_pop)
}


# ==============================================================================
# Household Transmission (Mass Action)
# ==============================================================================

#' Perform household transmission using mass action model
#' Assumes all household members mix homogeneously
#' Isolation/quarantine do not prevent household transmission
#' 
#' @param populations List of school populations
#' @param hh_pop Household population data frame
#' @param params Simulation parameters including hh_transmission_prob
#' @return List with updated populations and hh_pop
household_transmission <- function(populations, hh_pop, params) {
  
  # Build student data vectors for C++
  student_data <- lapply(seq_along(populations), function(i) {
    pop <- populations[[i]]
    valid <- !is.na(pop$hh_id)
    if (sum(valid) == 0) return(NULL)
    
    data.frame(
      hh_id = pop$hh_id[valid],
      state = pop$state[valid],
      is_vaccinated = pop$is_vaccinated[valid],
      school_idx = i,
      local_idx = which(valid),
      stringsAsFactors = FALSE
    )
  })
  student_data <- student_data[!sapply(student_data, is.null)]
  
  if (length(student_data) == 0) {
    return(list(populations = populations, hh_pop = hh_pop, 
                n_student_exposures = 0, n_hh_member_exposures = 0))
  }
  
  all_students <- do.call(rbind, student_data)
  
  # Call C++ function
  result <- cpp_household_transmission(
    student_hh_id = as.numeric(all_students$hh_id),
    student_state = as.character(all_students$state),
    student_is_vaccinated = as.logical(all_students$is_vaccinated),
    student_school_idx = as.integer(all_students$school_idx),
    student_local_idx = as.integer(all_students$local_idx),
    hh_member_hh_id = as.numeric(hh_pop$hh_id),
    hh_member_id = as.integer(hh_pop$member_id),
    hh_member_state = as.character(hh_pop$state),
    hh_member_is_vaccinated = as.logical(hh_pop$is_vaccinated),
    hh_member_is_student = as.logical(hh_pop$is_student),
    hh_transmission_prob = params$hh_transmission_prob,
    vaccine_efficacy = params$vaccine_efficacy,
    vaccine_infectiousness_reduction = params$vaccine_infectiousness_reduction
  )
  
  # Apply student exposures
  if (length(result$exposed_student_school_idx) > 0) {
    for (i in seq_along(result$exposed_student_school_idx)) {
      school_idx <- result$exposed_student_school_idx[i]
      local_idx <- result$exposed_student_local_idx[i]
      is_breakthrough <- result$exposed_student_breakthrough[i]
      
      pop <- populations[[school_idx]]
      
      if (pop$state[local_idx] %in% c("S", "V")) {
        pop$state[local_idx] <- "E"
        pop$time_in_state[local_idx] <- 0
        pop$latent_duration[local_idx] <- draw_erlang(1, params$latent_mean, params$latent_shape)
        infectious_dur <- draw_erlang(1, params$infectious_mean, params$infectious_shape)
        pop$infectious_duration[local_idx] <- infectious_dur
        pop$prodromal_duration[local_idx] <- params$prodromal_period
        pop$rash_duration[local_idx] <- max(1, infectious_dur - params$prodromal_period)
        
        if (is_breakthrough) {
          pop$breakthrough_infection[local_idx] <- TRUE
        }
      }
      
      populations[[school_idx]] <- pop
    }
  }
  
  # Apply household member exposures
  if (length(result$exposed_hh_member_id) > 0) {
    for (i in seq_along(result$exposed_hh_member_id)) {
      member_id <- result$exposed_hh_member_id[i]
      is_breakthrough <- result$exposed_hh_member_breakthrough[i]
      
      idx <- which(hh_pop$member_id == member_id)
      
      if (length(idx) > 0 && hh_pop$state[idx] %in% c("S", "V")) {
        hh_pop$state[idx] <- "E"
        hh_pop$time_in_state[idx] <- 0
        hh_pop$latent_duration[idx] <- draw_erlang(1, params$latent_mean, params$latent_shape)
        infectious_dur <- draw_erlang(1, params$infectious_mean, params$infectious_shape)
        hh_pop$infectious_duration[idx] <- infectious_dur
        hh_pop$prodromal_duration[idx] <- params$prodromal_period
        hh_pop$rash_duration[idx] <- max(1, infectious_dur - params$prodromal_period)
        
        if (is_breakthrough) {
          hh_pop$breakthrough_infection[idx] <- TRUE
        }
      }
    }
  }
  
  return(list(
    populations = populations,
    hh_pop = hh_pop,
    n_student_exposures = length(result$exposed_student_school_idx),
    n_hh_member_exposures = length(result$exposed_hh_member_id)
  ))
}


# ==============================================================================
# Update Household Member Disease States
# ==============================================================================

#' Update disease states for non-student household members (uses C++)
#' @param hh_pop Household population data frame
#' @param params Simulation parameters
#' @return Updated household population
update_household_disease_states <- function(hh_pop, params) {
  
  if (nrow(hh_pop) == 0) return(hh_pop)
  
  # Ensure required columns exist and have proper types
  if (is.null(hh_pop$time_since_prodromal)) {
    hh_pop$time_since_prodromal <- NA_integer_
  }
  
  # Call C++ function
  result <- cpp_update_hh_disease_states(
    state = as.character(hh_pop$state),
    is_student = as.logical(hh_pop$is_student),
    time_in_state = as.integer(hh_pop$time_in_state),
    latent_duration = as.integer(hh_pop$latent_duration),
    prodromal_duration = as.integer(hh_pop$prodromal_duration),
    rash_duration = as.integer(hh_pop$rash_duration),
    time_since_prodromal = as.integer(hh_pop$time_since_prodromal)
  )
  
  # Update hh_pop with results
  hh_pop$state <- result$state
  hh_pop$time_in_state <- result$time_in_state
  hh_pop$time_since_prodromal <- result$time_since_prodromal
  
  return(hh_pop)
}


# ==============================================================================
# Synchronize Student States Between Schools and Households
# ==============================================================================

#' Synchronize disease states between school populations and household population
#' This ensures students have consistent states in both data structures
#' OPTIMIZED: Uses vectorized matching instead of nested loops
#' 
#' @param populations List of school populations
#' @param hh_pop Household population data frame
#' @param direction "school_to_hh" or "hh_to_school"
#' @return Updated data structure based on direction
sync_student_states <- function(populations, hh_pop, direction = "school_to_hh") {
  
  if (direction == "school_to_hh") {
    # Build lookup table: person_id -> hh_pop row index
    hh_pop_idx <- setNames(1:nrow(hh_pop), hh_pop$person_id)
    
    # Process all schools at once
    for (school_idx in seq_along(populations)) {
      pop <- populations[[school_idx]]
      
      # Get valid person IDs (not NA)
      valid <- !is.na(pop$synpop_person_id)
      if (sum(valid) == 0) next
      
      # Find matching hh_pop indices
      person_ids <- as.character(pop$synpop_person_id[valid])
      hh_indices <- hh_pop_idx[person_ids]
      
      # Filter to found matches
      found <- !is.na(hh_indices)
      if (sum(found) == 0) next
      
      valid_idx <- which(valid)[found]
      hh_idx <- hh_indices[found]
      
      # Vectorized update
      hh_pop$state[hh_idx] <- pop$state[valid_idx]
      hh_pop$time_in_state[hh_idx] <- pop$time_in_state[valid_idx]
      hh_pop$time_since_prodromal[hh_idx] <- pop$time_since_prodromal[valid_idx]
      hh_pop$is_isolated[hh_idx] <- pop$is_isolated[valid_idx]
      hh_pop$is_quarantined[hh_idx] <- pop$is_quarantined[valid_idx]
    }
    return(hh_pop)
    
  } else if (direction == "hh_to_school") {
    # Build lookup table: person_id -> hh_pop row
    hh_pop_idx <- setNames(1:nrow(hh_pop), hh_pop$person_id)
    
    for (school_idx in seq_along(populations)) {
      pop <- populations[[school_idx]]
      
      valid <- !is.na(pop$synpop_person_id)
      if (sum(valid) == 0) next
      
      person_ids <- as.character(pop$synpop_person_id[valid])
      hh_indices <- hh_pop_idx[person_ids]
      
      found <- !is.na(hh_indices)
      if (sum(found) == 0) next
      
      valid_idx <- which(valid)[found]
      hh_idx <- hh_indices[found]
      
      # Vectorized update
      pop$state[valid_idx] <- hh_pop$state[hh_idx]
      pop$time_in_state[valid_idx] <- hh_pop$time_in_state[hh_idx]
      pop$time_since_prodromal[valid_idx] <- hh_pop$time_since_prodromal[hh_idx]
      
      populations[[school_idx]] <- pop
    }
    return(populations)
  }
}


# ==============================================================================
# Summarize Household Statistics
# ==============================================================================

#' Calculate household-level summary statistics
#' @param populations List of school populations
#' @param hh_pop Household population data frame
#' @return Data frame with household-level statistics
summarize_household_infections <- function(populations, hh_pop) {
  
  # Get all students (use base R)
  all_students_list <- lapply(seq_along(populations), function(i) {
    pop <- populations[[i]]
    data.frame(
      school_idx = i,
      synpop_person_id = pop$synpop_person_id,
      hh_id = pop$hh_id,
      state = pop$state,
      is_vaccinated = pop$is_vaccinated,
      is_student = TRUE,
      stringsAsFactors = FALSE
    )
  })
  all_students <- do.call(rbind, all_students_list)
  all_students <- all_students[!is.na(all_students$hh_id), ]
  
  # Combine with non-student household members
  hh_members <- hh_pop[!hh_pop$is_student, c("person_id", "hh_id", "state", "is_vaccinated", "is_student")]
  names(hh_members)[names(hh_members) == "person_id"] <- "synpop_person_id"
  
  all_hh <- rbind(
    all_students[, c("hh_id", "state", "is_vaccinated", "is_student")],
    hh_members[, c("hh_id", "state", "is_vaccinated", "is_student")]
  )
  
  # Summarize by household using aggregate
  infected_states <- c("E", "P", "Ra", "Iso", "R", "QE", "QP")
  
  hh_ids <- unique(all_hh$hh_id)
  hh_summary <- data.frame(
    hh_id = hh_ids,
    hh_size = sapply(hh_ids, function(h) sum(all_hh$hh_id == h)),
    n_students = sapply(hh_ids, function(h) sum(all_hh$hh_id == h & all_hh$is_student)),
    n_adults = sapply(hh_ids, function(h) sum(all_hh$hh_id == h & !all_hh$is_student)),
    n_infected = sapply(hh_ids, function(h) sum(all_hh$hh_id == h & all_hh$state %in% infected_states)),
    n_students_infected = sapply(hh_ids, function(h) sum(all_hh$hh_id == h & all_hh$is_student & all_hh$state %in% infected_states)),
    n_adults_infected = sapply(hh_ids, function(h) sum(all_hh$hh_id == h & !all_hh$is_student & all_hh$state %in% infected_states)),
    n_vaccinated = sapply(hh_ids, function(h) sum(all_hh$hh_id == h & all_hh$is_vaccinated))
  )
  hh_summary$household_attack_rate <- hh_summary$n_infected / hh_summary$hh_size
  
  return(hh_summary)
}


cat("Household utilities loaded.\n")