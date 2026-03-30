# ==============================================================================
# GRADE UTILITIES
# ==============================================================================
# File: grade_utils.R
# Contains: Functions for parsing and standardizing grade ranges,
#           and grade-based transmission weighting
# ==============================================================================

library(dplyr)
library(stringr)

# ==============================================================================
# Grade Mapping Reference
# ==============================================================================

# Standard grade numeric mapping:
# -3 = Infant/Crib (0-1 years)
# -2 = Toddler/2K (2 years)  
# -1 = PreK/3K/PK (3 years)
#  0 = K4/4K/PreK (4 years)
#  1 = K5/5K/Kindergarten (5 years)
#  2-6 = 1st-5th grade (Elementary)
#  7-9 = 6th-8th grade (Middle)
# 10-13 = 9th-12th grade (High)

# ==============================================================================
# Parse Single Grade Token
# ==============================================================================

#' Convert a single grade token to numeric value
#' @param token Character string representing a grade (e.g., "K5", "9", "PreK")
#' @return Numeric grade value
parse_grade_token <- function(token) {
  if (is.na(token) || token == "" || token == "NA") {
    return(NA_real_)
  }
  
  # Clean the token
  token <- trimws(toupper(token))
  token <- gsub("[^A-Z0-9]", "", token)  # Remove special characters
  
  # Handle infant/baby/crib (very young)
  if (grepl("^(INFANT|BABY|BABIES|CRIB|NEWBORN)", token)) {
    return(-3)
  }
  
  # Handle toddler/2-year
  if (grepl("^(TODDLER|2YEAR|2YR)", token)) {
    return(-2)
  }
  
  # Handle various PreK/PK formats
  # 1K, 2K, 3K = Pre-kindergarten ages
  if (grepl("^1K$", token)) return(-2)
  if (grepl("^2K$", token)) return(-1)
  if (grepl("^3K$|^K3$", token)) return(0)
  if (grepl("^4K$|^K4$", token)) return(1)
  if (grepl("^5K$|^K5$", token)) return(2)  # Kindergarten = grade 1 equivalent
  
  # PK, PRE, PREK = Pre-kindergarten
  if (grepl("^(PK|PRE|PREK|PREKINDERGARTEN|PIP)$", token)) {
    return(0)
  }
  
  # K alone = Kindergarten (K5)
  if (grepl("^K$|^KINDERGARTEN$", token)) {
    return(2)
  }
  
  # Handle "Xth" or "Xst" or "Xnd" or "Xrd" format
  token <- gsub("(ST|ND|RD|TH|GRADE|GR)$", "", token)
  
  # Try to extract numeric grade
  if (grepl("^[0-9]+$", token)) {
    grade <- as.numeric(token)
    # Grades 1-12 map to 2-13 (since K=2)
    if (grade >= 1 && grade <= 12) {
      return(grade + 1)
    }
    # Could be age-based notation for very young
    if (grade < 1) {
      return(-2)  # Assume toddler
    }
  }
  
  # Handle "ADULT" 
  if (grepl("ADULT", token)) {
    return(14)  # Beyond high school
  }
  
  # If we can't parse, return NA
  return(NA_real_)
}


# ==============================================================================
# Parse Grade Range String
# ==============================================================================

#' Parse a grade range string into min and max grades
#' @param grade_range Character string (e.g., "K5-12", "6-8", "PreK-5th")
#' @return Named list with min_grade, max_grade, grade_span
parse_grade_range <- function(grade_range) {
  
  if (is.na(grade_range) || grade_range == "" || grade_range == "NA") {
    return(list(min_grade = NA, max_grade = NA, grade_span = NA))
  }
  
  # Clean the string
  gr <- trimws(as.character(grade_range))
  gr <- toupper(gr)
  
  # Remove leading apostrophe if present
  gr <- gsub("^'", "", gr)
  
  # Handle special full-text cases
  if (grepl("THROUGH|THRU|TO", gr)) {
    gr <- gsub("THROUGH|THRU|TO", "-", gr)
  }
  
  # Handle comma-separated grades (e.g., "6, 7, 8" or "6,7,8")
  if (grepl(",", gr)) {
    grades <- unlist(strsplit(gr, "[,\\s]+"))
    grades <- sapply(grades, parse_grade_token)
    grades <- grades[!is.na(grades)]
    if (length(grades) > 0) {
      return(list(
        min_grade = min(grades),
        max_grade = max(grades),
        grade_span = max(grades) - min(grades) + 1
      ))
    }
  }
  
  # Try to split on common delimiters
  # Handle cases like "K5-12", "6-8", "PreK - 5th"
  parts <- unlist(strsplit(gr, "[-–—\\s]+"))
  parts <- parts[parts != "" & !grepl("^(GRADE|GR)$", parts)]
  
  if (length(parts) == 1) {
    # Single grade (e.g., "K5" or "9-12" that didn't split properly)
    # Try splitting on transition from letter to number
    if (grepl("[A-Z][0-9]", parts[1]) && !grepl("-", grade_range)) {
      # Single grade like "K5"
      grade <- parse_grade_token(parts[1])
      return(list(
        min_grade = grade,
        max_grade = grade,
        grade_span = 1
      ))
    }
  }
  
  if (length(parts) >= 2) {
    # Take first and last as range bounds
    min_g <- parse_grade_token(parts[1])
    max_g <- parse_grade_token(parts[length(parts)])
    
    # Ensure min <= max
    if (!is.na(min_g) && !is.na(max_g)) {
      if (min_g > max_g) {
        temp <- min_g
        min_g <- max_g
        max_g <- temp
      }
      return(list(
        min_grade = min_g,
        max_grade = max_g,
        grade_span = max_g - min_g + 1
      ))
    }
  }
  
  # Last resort: try to find any numbers
  numbers <- as.numeric(unlist(regmatches(gr, gregexpr("[0-9]+", gr))))
  if (length(numbers) >= 2) {
    # Assume first and last numbers are the range
    return(list(
      min_grade = min(numbers) + 1,  # Convert to our grade scale
      max_grade = max(numbers) + 1,
      grade_span = max(numbers) - min(numbers) + 1
    ))
  } else if (length(numbers) == 1) {
    grade <- numbers[1] + 1
    return(list(
      min_grade = grade,
      max_grade = grade,
      grade_span = 1
    ))
  }
  
  return(list(min_grade = NA, max_grade = NA, grade_span = NA))
}


# ==============================================================================
# Classify School Type
# ==============================================================================

#' Classify school type based on grade range
#' @param min_grade Minimum grade (numeric)
#' @param max_grade Maximum grade (numeric)
#' @return Character string: "Preschool", "Elementary", "Middle", "High", 
#'         "Elementary-Middle", "Middle-High", "K-12", or "Other"
classify_school_type <- function(min_grade, max_grade) {
  
  if (is.na(min_grade) || is.na(max_grade)) {
    return("Unknown")
  }
  
  # Grade boundaries (in our numeric scale):
  # Preschool: -3 to 1 (infant through K4)
  # Elementary: 2 to 6 (K5 through 5th)
  # Middle: 7 to 9 (6th through 8th)
  # High: 10 to 13 (9th through 12th)
  
  has_preschool <- min_grade <= 1
  has_elementary <- (min_grade <= 6 && max_grade >= 2)
  has_middle <- (min_grade <= 9 && max_grade >= 7)
  has_high <- max_grade >= 10
  
  # Determine primary classification
  if (has_preschool && !has_elementary && !has_middle && !has_high) {
    return("Preschool")
  } else if (has_elementary && !has_middle && !has_high) {
    if (has_preschool) return("PreK-Elementary")
    return("Elementary")
  } else if (has_middle && !has_elementary && !has_high) {
    return("Middle")
  } else if (has_high && !has_elementary && !has_middle) {
    return("High")
  } else if (has_elementary && has_middle && !has_high) {
    return("Elementary-Middle")
  } else if (has_middle && has_high && !has_elementary) {
    return("Middle-High")
  } else if (has_elementary && has_middle && has_high) {
    return("K-12")
  } else if (has_elementary && has_high) {
    return("K-12")  # Spans all
  } else {
    return("Other")
  }
}


# ==============================================================================
# Standardize Grade Range Column
# ==============================================================================

#' Add standardized grade columns to schools data frame
#' @param schools Data frame with Grade.Range column
#' @param grade_col Name of the grade range column (default "Grade.Range")
#' @return Data frame with added columns: min_grade, max_grade, grade_span, 
#'         school_type, grade_category
standardize_grade_range <- function(schools, grade_col = "Grade.Range") {
  
  if (!grade_col %in% names(schools)) {
    warning(sprintf("Column '%s' not found in data frame", grade_col))
    return(schools)
  }
  
  cat("=== STANDARDIZING GRADE RANGES ===\n")
  
  # Parse each grade range
  parsed <- lapply(schools[[grade_col]], parse_grade_range)
  
  # Extract components
  schools$min_grade <- sapply(parsed, function(x) x$min_grade)
  schools$max_grade <- sapply(parsed, function(x) x$max_grade)
  schools$grade_span <- sapply(parsed, function(x) x$grade_span)
  
  # Classify school types
  schools$school_type <- mapply(classify_school_type, 
                                 schools$min_grade, 
                                 schools$max_grade)
  
  # Create simplified grade category for transmission matching
  # Elementary (K-5), Middle (6-8), High (9-12)
  schools$grade_category <- sapply(1:nrow(schools), function(i) {
    min_g <- schools$min_grade[i]
    max_g <- schools$max_grade[i]
    
    if (is.na(min_g) || is.na(max_g)) return("Unknown")
    
    categories <- c()
    if (min_g <= 6 && max_g >= 2) categories <- c(categories, "E")  # Elementary
    if (min_g <= 9 && max_g >= 7) categories <- c(categories, "M")  # Middle
    if (max_g >= 10) categories <- c(categories, "H")                # High
    
    if (length(categories) == 0) {
      if (max_g <= 1) return("P")  # Preschool only
      return("Unknown")
    }
    
    paste(categories, collapse = "-")
  })
  
  # Summary statistics
  cat("\nSchool Type Distribution:\n")
  print(table(schools$school_type))
  
  cat("\nGrade Category Distribution:\n
")
  print(table(schools$grade_category))
  
  # Check parsing success rate
  n_parsed <- sum(!is.na(schools$min_grade))
  cat(sprintf("\nSuccessfully parsed: %d / %d (%.1f%%)\n",
              n_parsed, nrow(schools), 100 * n_parsed / nrow(schools)))
  
  # Show failed parses if any
  failed <- schools[is.na(schools$min_grade), ]
  if (nrow(failed) > 0 && nrow(failed) <= 20) {
    cat("\nFailed to parse:\n")
    print(failed[[grade_col]])
  } else if (nrow(failed) > 20) {
    cat(sprintf("\n%d grade ranges failed to parse\n", nrow(failed)))
  }
  
  return(schools)
}


# ==============================================================================
# Calculate Grade Overlap
# ==============================================================================

#' Calculate grade overlap between two schools
#' @param min1, max1 Grade range for school 1
#' @param min2, max2 Grade range for school 2
#' @return Numeric overlap score (0 = no overlap, 1 = perfect overlap)
calculate_grade_overlap <- function(min1, max1, min2, max2) {
  
  if (is.na(min1) || is.na(max1) || is.na(min2) || is.na(max2)) {
    return(0.5)  # Unknown - use moderate default
  }
  
  # Calculate intersection
  overlap_min <- max(min1, min2)
  overlap_max <- min(max1, max2)
  
  if (overlap_min > overlap_max) {
    # No overlap
    return(0)
  }
  
  overlap_size <- overlap_max - overlap_min + 1
  
  # Normalize by the smaller school's grade span
  span1 <- max1 - min1 + 1
  span2 <- max2 - min2 + 1
  min_span <- min(span1, span2)
  
  overlap_score <- overlap_size / min_span
  
  return(min(1, overlap_score))  # Cap at 1
}


# ==============================================================================
# Create Grade Overlap Matrix
# ==============================================================================
  
#' Create a matrix of grade overlaps between all schools
#' @param schools Data frame with min_grade and max_grade columns
#' @return Matrix of overlap scores (n_schools x n_schools)
create_grade_overlap_matrix <- function(schools) {
  
  n <- nrow(schools)
  overlap_matrix <- matrix(0, nrow = n, ncol = n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      if (i == j) {
        overlap_matrix[i, j] <- 1
      } else {
        overlap_matrix[i, j] <- calculate_grade_overlap(
          schools$min_grade[i], schools$max_grade[i],
          schools$min_grade[j], schools$max_grade[j]
        )
      }
    }
  }
  
  return(overlap_matrix)
}


# ==============================================================================
# Generate Grade-Weighted Network
# ==============================================================================

#' Modify network weights based on grade overlap
#' Schools with similar grades have stronger connections (sports, events)
#' 
#' @param network Network object with adjacency matrix
#' @param schools Data frame with min_grade, max_grade columns
#' @param grade_weight How much to weight by grade overlap (0-1, default 0.5)
#'        0 = ignore grades, 1 = only connect same-grade schools
#' @param same_type_bonus Additional weight for same school type (default 0.2)
#' @return Updated network object with grade-weighted adjacency
apply_grade_weighting_to_network <- function(network, schools, 
                                              grade_weight = 0.5,
                                              same_type_bonus = 0.2) {
  
  cat("\n=== APPLYING GRADE-BASED WEIGHTING TO NETWORK ===\n")
  
  if (!"min_grade" %in% names(schools) || !"max_grade" %in% names(schools)) {
    warning("Schools data frame missing grade columns. Run standardize_grade_range() first.")
    return(network)
  }
  
  n <- nrow(schools)
  adj <- network$adjacency
  
  # Create grade overlap matrix
  grade_overlap <- create_grade_overlap_matrix(schools)
  
  # Create same-type bonus matrix
  type_bonus <- matrix(0, nrow = n, ncol = n)
  if ("school_type" %in% names(schools)) {
    for (i in 1:n) {
      for (j in 1:n) {
        if (i != j && schools$school_type[i] == schools$school_type[j]) {
          type_bonus[i, j] <- same_type_bonus
        }
      }
    }
  }
  
  # Combine: weighted average of original and grade-weighted
  # new_weight = original * (1 - grade_weight) + original * grade_overlap * grade_weight + type_bonus
  # Simplified: new_weight = original * (1 - grade_weight + grade_overlap * grade_weight) + type_bonus
  
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j && adj[i, j] > 0) {
        grade_factor <- 1 - grade_weight + grade_overlap[i, j] * grade_weight
        adj[i, j] <- adj[i, j] * grade_factor + type_bonus[i, j]
      }
    }
  }
  
  # Normalize to keep weights in reasonable range
  max_weight <- max(adj[adj > 0])
  if (max_weight > 1) {
    adj <- adj / max_weight
  }
  
  network$adjacency <- adj
  network$grade_weighted <- TRUE
  network$grade_weight_param <- grade_weight
  
  # Summary
  cat(sprintf("Grade weight parameter: %.2f\n", grade_weight))
  cat(sprintf("Same-type bonus: %.2f\n", same_type_bonus))
  
  # Show some examples
  cat("\nExample grade-weighted connections:\n")
  for (i in 1:min(3, n)) {
    connected <- which(adj[i, ] > 0)
    if (length(connected) > 0) {
      j <- connected[1]
      cat(sprintf("  %s (%s) <-> %s (%s): weight=%.3f, grade_overlap=%.2f\n",
                  schools$school_name[i], schools$school_type[i],
                  schools$school_name[j], schools$school_type[j],
                  adj[i, j], grade_overlap[i, j]))
    }
  }
  
  return(network)
}


# ==============================================================================
# Helper: Human-Readable Grade
# ==============================================================================

#' Convert numeric grade back to human-readable format
#' @param grade Numeric grade value
#' @return Character string (e.g., "K", "5th", "9th")
grade_to_string <- function(grade) {
  if (is.na(grade)) return(NA)
  
  if (grade <= -2) return("PreK")
  if (grade == -1) return("3K")
  if (grade == 0) return("4K")
  if (grade == 1) return("K4")
  if (grade == 2) return("K5")
  
  # Grades 1-12
  actual_grade <- grade - 1
  if (actual_grade >= 1 && actual_grade <= 12) {
    suffix <- switch(as.character(actual_grade %% 10),
                     "1" = if (actual_grade != 11) "st" else "th",
                     "2" = if (actual_grade != 12) "nd" else "th",
                     "3" = if (actual_grade != 13) "rd" else "th",
                     "th")
    return(paste0(actual_grade, suffix))
  }
  
  return(as.character(grade))
}


#' Create standardized grade range string
#' @param min_grade Minimum grade (numeric)
#' @param max_grade Maximum grade (numeric)
#' @return Character string (e.g., "K5-5th", "6th-8th", "9th-12th")
grade_range_to_string <- function(min_grade, max_grade) {
  if (is.na(min_grade) || is.na(max_grade)) return(NA)
  
  min_str <- grade_to_string(min_grade)
  max_str <- grade_to_string(max_grade)
  
  if (min_str == max_str) return(min_str)
  
  return(paste0(min_str, "-", max_str))
}


# ==============================================================================
# Add Standardized String Column
# ==============================================================================

#' Add a human-readable standardized grade range column
#' @param schools Data frame with min_grade and max_grade columns
#' @return Data frame with added grade_range_std column
add_standardized_grade_string <- function(schools) {
  
  schools$grade_range_std <- mapply(
    grade_range_to_string,
    schools$min_grade,
    schools$max_grade
  )
  
  return(schools)
}


cat("Grade utilities loaded.\n")
