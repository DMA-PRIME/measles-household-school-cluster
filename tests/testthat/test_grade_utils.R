# ==============================================================================
# UNIT TESTS: grade_utils.R
# ==============================================================================
# Tests for: parse_grade_token, parse_grade_range, classify_school_type,
#            standardize_grade_range, calculate_grade_overlap,
#            create_grade_overlap_matrix, grade_to_string,
#            grade_range_to_string, add_standardized_grade_string
# ==============================================================================

# Source the module under test
project_root <- Sys.getenv("PROJECT_ROOT", unset = normalizePath("."))
source(file.path(project_root, "grade_utils.R"))

# ==============================================================================
# parse_grade_token
# ==============================================================================

test_that("parse_grade_token: kindergarten variants return 2", {
  expect_equal(parse_grade_token("K5"),  2)
  expect_equal(parse_grade_token("5K"),  2)
  expect_equal(parse_grade_token("K"),   2)
  expect_equal(parse_grade_token("KINDERGARTEN"), 2)
})

test_that("parse_grade_token: pre-K variants return expected values", {
  expect_equal(parse_grade_token("PK"),   0)
  expect_equal(parse_grade_token("PreK"), 0)
  expect_equal(parse_grade_token("K4"),   1)
  expect_equal(parse_grade_token("4K"),   1)
  expect_equal(parse_grade_token("3K"),   0)
  expect_equal(parse_grade_token("K3"),   0)
  expect_equal(parse_grade_token("2K"),  -1)
  expect_equal(parse_grade_token("1K"),  -2)
})

test_that("parse_grade_token: numeric grades 1-12 map to scale+1", {
  expect_equal(parse_grade_token("1"),   2)
  expect_equal(parse_grade_token("5"),   6)
  expect_equal(parse_grade_token("6"),   7)
  expect_equal(parse_grade_token("8"),   9)
  expect_equal(parse_grade_token("9"),  10)
  expect_equal(parse_grade_token("12"), 13)
})

test_that("parse_grade_token: ordinal suffixes are stripped", {
  expect_equal(parse_grade_token("1st"),  2)
  expect_equal(parse_grade_token("2nd"),  3)
  expect_equal(parse_grade_token("3rd"),  4)
  expect_equal(parse_grade_token("9th"), 10)
  expect_equal(parse_grade_token("12th"), 13)
})

test_that("parse_grade_token: NA / empty inputs return NA", {
  expect_true(is.na(parse_grade_token(NA)))
  expect_true(is.na(parse_grade_token("")))
  expect_true(is.na(parse_grade_token("NA")))
})

test_that("parse_grade_token: special labels return correct values", {
  expect_equal(parse_grade_token("INFANT"),  -3)
  expect_equal(parse_grade_token("TODDLER"), -2)
  expect_equal(parse_grade_token("ADULT"),   14)
})

# ==============================================================================
# parse_grade_range
# ==============================================================================

test_that("parse_grade_range: simple numeric range", {
  result <- parse_grade_range("6-8")
  # 6 -> 7, 8 -> 9 in our scale
  expect_equal(result$min_grade, 7)
  expect_equal(result$max_grade, 9)
  expect_equal(result$grade_span, 3)
})

test_that("parse_grade_range: K5-12 range", {
  result <- parse_grade_range("K5-12")
  expect_equal(result$min_grade, 2)   # K5 = 2
  expect_equal(result$max_grade, 13)  # 12th = 13
  expect_equal(result$grade_span, 12)
})

test_that("parse_grade_range: single grade string", {
  result <- parse_grade_range("K5")
  expect_equal(result$min_grade, 2)
  expect_equal(result$max_grade, 2)
  expect_equal(result$grade_span, 1)
})

test_that("parse_grade_range: NA / empty inputs return NAs", {
  result <- parse_grade_range(NA)
  expect_true(is.na(result$min_grade))
  expect_true(is.na(result$max_grade))
  expect_true(is.na(result$grade_span))

  result2 <- parse_grade_range("")
  expect_true(is.na(result2$min_grade))
})

test_that("parse_grade_range: PreK-5th range", {
  result <- parse_grade_range("PreK-5th")
  expect_equal(result$min_grade, 0)  # PreK = 0
  expect_equal(result$max_grade, 6)  # 5th = 6
})

test_that("parse_grade_range: comma-separated grades", {
  result <- parse_grade_range("6,7,8")
  expect_equal(result$min_grade, 7)
  expect_equal(result$max_grade, 9)
})

# ==============================================================================
# classify_school_type
# ==============================================================================

test_that("classify_school_type: elementary (K5-5th)", {
  expect_equal(classify_school_type(2, 6), "Elementary")
})

test_that("classify_school_type: middle school (6th-8th)", {
  expect_equal(classify_school_type(7, 9), "Middle")
})

test_that("classify_school_type: high school (9th-12th)", {
  expect_equal(classify_school_type(10, 13), "High")
})

test_that("classify_school_type: K-12 spans all levels", {
  expect_equal(classify_school_type(2, 13), "K-12")
})

test_that("classify_school_type: elementary-middle combo", {
  expect_equal(classify_school_type(2, 9), "Elementary-Middle")
})

test_that("classify_school_type: middle-high combo", {
  expect_equal(classify_school_type(7, 13), "Middle-High")
})

test_that("classify_school_type: preschool only", {
  expect_equal(classify_school_type(-3, 1), "Preschool")
})

test_that("classify_school_type: NA inputs return 'Unknown'", {
  expect_equal(classify_school_type(NA, 6), "Unknown")
  expect_equal(classify_school_type(2, NA), "Unknown")
})

# ==============================================================================
# calculate_grade_overlap
# ==============================================================================

test_that("calculate_grade_overlap: identical ranges return 1", {
  expect_equal(calculate_grade_overlap(2, 6, 2, 6), 1)
})

test_that("calculate_grade_overlap: no overlap returns 0", {
  expect_equal(calculate_grade_overlap(2, 6, 7, 9), 0)
  expect_equal(calculate_grade_overlap(7, 9, 2, 6), 0)
})

test_that("calculate_grade_overlap: partial overlap returns value between 0 and 1", {
  score <- calculate_grade_overlap(2, 9, 7, 13)
  expect_gt(score, 0)
  expect_lte(score, 1)
})

test_that("calculate_grade_overlap: NA inputs return 0.5 default", {
  expect_equal(calculate_grade_overlap(NA, 6, 2, 6), 0.5)
  expect_equal(calculate_grade_overlap(2, NA, 2, 6), 0.5)
})

test_that("calculate_grade_overlap: score is capped at 1", {
  # Perfect overlap (subset case)
  score <- calculate_grade_overlap(2, 13, 7, 9)
  expect_lte(score, 1)
})

# ==============================================================================
# create_grade_overlap_matrix
# ==============================================================================

test_that("create_grade_overlap_matrix: diagonal is 1", {
  schools <- data.frame(
    min_grade = c(2,  7, 10),
    max_grade = c(6,  9, 13)
  )
  mat <- create_grade_overlap_matrix(schools)
  expect_equal(diag(mat), rep(1, 3))
})

test_that("create_grade_overlap_matrix: correct dimensions", {
  schools <- data.frame(min_grade = c(2, 7, 10), max_grade = c(6, 9, 13))
  mat <- create_grade_overlap_matrix(schools)
  expect_equal(dim(mat), c(3, 3))
})

test_that("create_grade_overlap_matrix: non-overlapping schools have 0 off-diagonal", {
  schools <- data.frame(min_grade = c(2, 10), max_grade = c(6, 13))
  mat <- create_grade_overlap_matrix(schools)
  expect_equal(mat[1, 2], 0)
  expect_equal(mat[2, 1], 0)
})

# ==============================================================================
# grade_to_string
# ==============================================================================

test_that("grade_to_string: converts numeric grades to strings", {
  expect_equal(grade_to_string(2),  "K5")
  expect_equal(grade_to_string(1),  "K4")
  expect_equal(grade_to_string(0),  "4K")
  expect_equal(grade_to_string(-1), "3K")
  expect_equal(grade_to_string(3),  "2nd")
  expect_equal(grade_to_string(4),  "3rd")
  expect_equal(grade_to_string(13), "12th")
})

test_that("grade_to_string: NA returns NA", {
  expect_true(is.na(grade_to_string(NA)))
})

# ==============================================================================
# grade_range_to_string
# ==============================================================================

test_that("grade_range_to_string: same min and max returns single grade", {
  expect_equal(grade_range_to_string(2, 2), "K5")
})

test_that("grade_range_to_string: range returns hyphenated string", {
  result <- grade_range_to_string(2, 6)
  expect_match(result, "-")
})

test_that("grade_range_to_string: NA inputs return NA", {
  expect_true(is.na(grade_range_to_string(NA, 6)))
  expect_true(is.na(grade_range_to_string(2, NA)))
})

# ==============================================================================
# add_standardized_grade_string
# ==============================================================================

test_that("add_standardized_grade_string: adds grade_range_std column", {
  schools <- data.frame(
    school_id = 1:3,
    min_grade = c(2,  7, 10),
    max_grade = c(6,  9, 13)
  )
  result <- add_standardized_grade_string(schools)
  expect_true("grade_range_std" %in% names(result))
  expect_equal(nrow(result), 3)
  expect_false(any(is.na(result$grade_range_std)))
})

# ==============================================================================
# standardize_grade_range
# ==============================================================================

test_that("standardize_grade_range: adds expected columns", {
  schools <- data.frame(
    school_id = 1:3,
    Grade.Range = c("K5-5", "6-8", "9-12")
  )
  result <- suppressMessages(suppressWarnings(standardize_grade_range(schools)))
  expect_true(all(c("min_grade", "max_grade", "grade_span", "school_type") %in% names(result)))
})

test_that("standardize_grade_range: warns if grade column missing", {
  schools <- data.frame(school_id = 1)
  expect_warning(standardize_grade_range(schools, grade_col = "NonExistent"))
})
