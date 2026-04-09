# ==============================================================================
# ZCTA-LEVEL VACCINATION SPATIAL COMPARISON
# ==============================================================================
rm(list = ls())
library(sf)
library(dplyr)
library(ggplot2)
library(patchwork)

# ==============================================================================
# CONFIGURATION
# ==============================================================================
school_data_file   <- "data/SC_vaccination_merged.csv"
synpop_file        <- "C:/Users/pandey7/OneDrive - Clemson University/Research/Research/Agent-Based Model/usa_synth_pop/sythetic_population_sc.csv"
school_ref_file    <- "data/School_reference_Upstate.csv"
hh_cache_file      <- "data/household_assignment_cache.rds"

region_counties <- c("Spartanburg", "Greenville", "Abbeville", "Anderson", "Cherokee",
                     "Laurens", "McCormick", "Oconee", "Pickens", "Union", "Greenwood")

params <- list(
  vaccine_efficacy = 0.97,
  avg_class_size = 25,
  age_range = c(5, 18)
)


# ==============================================================================
# LOAD SIMULATION MODULES
# ==============================================================================
cat("Loading simulation modules...\n")
source("codes/load_all.R", local = globalenv())

# ==============================================================================
# LOAD AND FILTER SCHOOL DATA
# ==============================================================================
cat("Loading school data...\n")
schools_raw <- read.csv(school_data_file, stringsAsFactors = FALSE)

schools <- data.frame(
  school_id   = seq_len(nrow(schools_raw)),
  school_name = as.character(schools_raw$School.Name),
  school_size = as.integer(schools_raw$Total.Students),
  vaccination_coverage = as.numeric(schools_raw$Percent.Immunized),
  lon = as.numeric(schools_raw$longitude),
  lat = as.numeric(schools_raw$latitude),
  county = as.character(schools_raw$County),
  stringsAsFactors = FALSE
)

if ("Grade.Range" %in% names(schools_raw)) {
  schools$grade_range_orig <- as.character(schools_raw$Grade.Range)
}

schools <- schools[!is.na(schools$school_size) & schools$school_size > 0, , drop = FALSE]
schools$school_id <- seq_len(nrow(schools))
schools <- schools[schools$county %in% region_counties, , drop = FALSE]
schools$school_id <- seq_len(nrow(schools))

if ("grade_range_orig" %in% names(schools)) {
  if (file.exists("modules/grade_utils.R") || exists("standardize_grade_range")) {
    tryCatch({
      schools <- standardize_grade_range(schools, grade_col = "grade_range_orig")
      schools <- add_standardized_grade_string(schools)
    }, error = function(e) {
      cat("Grade utils not available, using fallback\n")
    })
  }
  if (!"grade_min" %in% names(schools)) {
    parse_grade <- function(g) {
      g <- gsub("PK|KG|K", "0", g)
      parts <- regmatches(g, gregexpr("[0-9]+", g))[[1]]
      if (length(parts) >= 2) return(as.integer(parts[1:2]))
      if (length(parts) == 1) return(rep(as.integer(parts[1]), 2))
      return(c(0L, 12L))
    }
    parsed <- t(sapply(schools$grade_range_orig, parse_grade))
    schools$grade_min <- parsed[, 1]
    schools$grade_max <- parsed[, 2]
  }
}

cat(sprintf("Loaded %d schools in %s\n", nrow(schools),
            paste(region_counties, collapse = ", ")))

if (max(schools$vaccination_coverage, na.rm = TRUE) > 1) {
  schools$vaccination_coverage <- schools$vaccination_coverage / 100
}

cat(sprintf("Vaccination coverage range: %.1f%% - %.1f%%\n",
            min(schools$vaccination_coverage, na.rm = TRUE) * 100,
            max(schools$vaccination_coverage, na.rm = TRUE) * 100))

# ==============================================================================
# LOAD SHAPEFILES AND CLIP TO STUDY COUNTIES
# ==============================================================================
cat("\nLoading shapefiles...\n")
load("C:/Users/pandey7/Box/BoxPHI-PHMR Projects/Data/Community Data/SC Shape Files/sc_shape_files.RData")
zcta_sf <- st_transform(sc_zcta_sf, 4326)
county_sf <- st_transform(sc_county_sf, 4326)

# Identify ZCTA ID column
zcta_id_col <- intersect(
  names(zcta_sf),
  c("ZCTA5CE20", "ZCTA5CE10", "ZCTA5", "GEOID20", "GEOID10", "GEOID", "ZTCA", "ZIP", "ZIPCODE")
)[1]

if (is.na(zcta_id_col)) {
  cat("Available columns: ", paste(names(zcta_sf), collapse = ", "), "\n")
  stop("Cannot identify ZCTA ID column. Please set zcta_id_col manually.")
}

cat(sprintf("ZCTA shapefile: %d polygons, ID column: %s\n", nrow(zcta_sf), zcta_id_col))

# Clip ZCTAs to study region counties
study_counties <- county_sf %>%
  filter(NAME %in% region_counties)

study_boundary <- st_union(study_counties)

zcta_sf <- st_intersection(zcta_sf, study_boundary)
zcta_sf <- st_make_valid(zcta_sf)
cat(sprintf("ZCTAs clipped to %d study counties: %d\n",
            length(region_counties), nrow(zcta_sf)))

# ==============================================================================
# CREATE POPULATIONS AND ASSIGN HOUSEHOLDS
# ==============================================================================
cat("\nCreating school populations...\n")
populations <- lapply(seq_len(nrow(schools)), function(i) {
  create_school_population(i, schools$school_size[i], params$avg_class_size, params$age_range)
})

cat("Loading synthetic population...\n")
synpop <- load_synthetic_population(synpop_file)

cache_valid <- FALSE
if (file.exists(hh_cache_file)) {
  saved <- readRDS(hh_cache_file)
  if (!is.null(saved$n_schools) && saved$n_schools == nrow(schools)) {
    cache_valid <- TRUE
  }
}

if (cache_valid) {
  cat("Using cached household assignment...\n")
  hh_result <- load_household_assignment(hh_cache_file, populations, schools, verbose = TRUE)
} else {
  cat("Computing household assignment (this may take several minutes)...\n")
  hh_result <- assign_households_to_students(
    populations = populations,
    schools = schools,
    synpop = synpop,
    max_distance_km = 30,
    grade_tolerance = 2,
    verbose = TRUE
  )
  save_household_assignment(hh_result, hh_cache_file, schools)
}

populations <- hh_result$populations
assignment_df <- hh_result$assignment_df
household_members <- hh_result$household_members

cat(sprintf("Assigned %d students from %d households to %d schools\n",
            nrow(assignment_df),
            length(unique(assignment_df$hh_id)),
            nrow(schools)))

# ==============================================================================
# STEP 1: APPLY HOUSEHOLD-CORRELATED VACCINATION
# ==============================================================================
cat("\n=== APPLYING HOUSEHOLD-CORRELATED VACCINATION ===\n")

populations_hh <- assign_household_level_vaccination(
  populations, schools, assignment_df, verbose = TRUE
)

# ==============================================================================
# STEP 2: EXTRACT STUDENT-LEVEL DATA
# ==============================================================================
student_records <- list()
for (s in seq_len(nrow(schools))) {
  pop <- populations_hh[[s]]
  has_loc <- !is.na(pop$hh_id) & !is.na(pop$hh_lon) & !is.na(pop$hh_lat)
  if (sum(has_loc) > 0) {
    student_records[[s]] <- data.frame(
      hh_lon          = pop$hh_lon[has_loc],
      hh_lat          = pop$hh_lat[has_loc],
      is_vaccinated   = pop$is_vaccinated[has_loc],
      school_coverage = schools$vaccination_coverage[s],
      stringsAsFactors = FALSE
    )
  }
}
student_df <- do.call(rbind, student_records)
cat(sprintf("Students with household locations: %d\n", nrow(student_df)))

# ==============================================================================
# STEP 3: SPATIAL JOIN — ASSIGN STUDENTS TO ZCTAs
# ==============================================================================
cat("\n=== SPATIAL JOINING STUDENTS TO ZCTAs ===\n")

student_pts <- st_as_sf(student_df, coords = c("hh_lon", "hh_lat"), crs = 4326)
student_zcta <- st_join(student_pts, zcta_sf[, zcta_id_col], join = st_within)
student_zcta$zcta <- st_drop_geometry(student_zcta)[[zcta_id_col]]

cat(sprintf("Students matched to ZCTAs: %d / %d\n",
            sum(!is.na(student_zcta$zcta)), nrow(student_zcta)))

# ==============================================================================
# STEP 4: AGGREGATE BY ZCTA
# ==============================================================================
cat("\n=== COMPUTING ZCTA VACCINATION RATES ===\n")

# School-level: based on schools PHYSICALLY LOCATED in each ZCTA
school_sf <- st_as_sf(schools[!is.na(schools$lon) & !is.na(schools$lat), ],
                       coords = c("lon", "lat"), crs = 4326)
school_zcta_join <- st_join(school_sf, zcta_sf[, zcta_id_col], join = st_within)
school_zcta_join$zcta <- st_drop_geometry(school_zcta_join)[[zcta_id_col]]

school_level_zcta <- st_drop_geometry(school_zcta_join) %>%
  filter(!is.na(zcta)) %>%
  group_by(zcta) %>%
  summarise(
    n_schools = n(),
    school_level_vax_rate = weighted.mean(vaccination_coverage, school_size, na.rm = TRUE),
    .groups = "drop"
  )

# Household-correlated: based on children LIVING in each ZCTA
hh_zcta <- st_drop_geometry(student_zcta) %>%
  filter(!is.na(zcta)) %>%
  group_by(zcta) %>%
  summarise(
    n_children = n(),
    hh_correlated_vax_rate = mean(is_vaccinated),
    .groups = "drop"
  )

# Merge
comparison_df <- hh_zcta %>%
  full_join(school_level_zcta, by = "zcta") %>%
  mutate(difference = hh_correlated_vax_rate - school_level_vax_rate)

# Filter out ZCTAs with too few children
min_children <- 10
comparison_df <- comparison_df %>%
  mutate(
    hh_correlated_vax_rate = ifelse(
      !is.na(n_children) & n_children >= min_children,
      hh_correlated_vax_rate, NA
    ),
    difference = ifelse(!is.na(hh_correlated_vax_rate), difference, NA)
  )

cat(sprintf("ZCTAs dropped (< %d children): %d\n",
            min_children,
            sum(!is.na(hh_zcta$hh_correlated_vax_rate)) -
              sum(!is.na(comparison_df$hh_correlated_vax_rate))))
cat(sprintf("ZCTAs with school data: %d\n", sum(!is.na(comparison_df$school_level_vax_rate))))
cat(sprintf("ZCTAs with household data: %d\n", sum(!is.na(comparison_df$hh_correlated_vax_rate))))

# ==============================================================================
# STEP 5: MERGE WITH SPATIAL DATA
# ==============================================================================
zcta_plot <- zcta_sf %>%
  mutate(zcta = .data[[zcta_id_col]]) %>%
  left_join(comparison_df, by = "zcta")

cat(sprintf("ZCTAs with data: %d\n", sum(!is.na(zcta_plot$hh_correlated_vax_rate))))

cat("\n=== VACCINATION RATE COMPARISON ===\n")
cat(sprintf("School-level mean:         %.1f%%\n",
            mean(comparison_df$school_level_vax_rate, na.rm = TRUE) * 100))
cat(sprintf("Household-correlated mean: %.1f%%\n",
            mean(comparison_df$hh_correlated_vax_rate, na.rm = TRUE) * 100))
cat(sprintf("Mean difference:           %.1f pp\n",
            mean(comparison_df$difference, na.rm = TRUE) * 100))
cat(sprintf("ZCTAs where HH model < school model: %d / %d\n",
            sum(comparison_df$difference < 0, na.rm = TRUE),
            sum(!is.na(comparison_df$difference))))

# ==============================================================================
# STEP 6: PLOT
# ==============================================================================
cat("\n=== GENERATING FIGURES ===\n")

# --- Theme ---
theme_map <- function() {
  theme_void(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 5)),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "#666666", margin = margin(b = 10)),
      legend.position = "bottom",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      legend.key.width = unit(1.5, "cm"),
      legend.key.height = unit(0.4, "cm"),
      plot.margin = margin(10, 10, 10, 10),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

# --- Color scales (red-to-dark-green, tuned for measles 95% threshold) ---
vax_limits <- c(0.50, 1.00)

vax_fill <- scale_fill_gradientn(
  colors = c("#B2182B",   # deep red (< 80%)
             "#D6604D",   # red (80-85%)
             "#F4A582",   # salmon (85-90%)
             "#FDDBC7",   # pale (90-92%)
             "#FFF7C4",   # yellow (92-95%)
             "#B7E6A5",   # light green (95-97%)
             "#1A7A2E"),   # dark green (97%+)
  values = scales::rescale(c(0.50, 0.80, 0.85, 0.90, 0.95, 0.97, 1.00),
                            from = vax_limits),
  limits = vax_limits,
  labels = scales::percent_format(accuracy = 1),
  name = "Vaccination\nCoverage",
  na.value = "gray85"
)

vax_color <- scale_color_gradientn(
  colors = c("#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
             "#FFF7C4", "#B7E6A5", "#1A7A2E"),
  values = scales::rescale(c(0.50, 0.80, 0.85, 0.90, 0.95, 0.97, 1.00),
                            from = vax_limits),
  limits = vax_limits,
  guide = "none",
  na.value = "gray85"
)

# --- County borders overlay ---
county_borders <- geom_sf(
  data = study_counties, fill = NA, color = "black", linewidth = 0.6,
  inherit.aes = FALSE
)

# --- School points overlay ---
school_overlay <- geom_point(
  data = schools[!is.na(schools$lon) & !is.na(schools$lat), ],
  aes(x = lon, y = lat, color = vaccination_coverage),
  shape = 16, size = 1.8, inherit.aes = FALSE
)

# --- Map A: School-level ---
p_school <- ggplot(zcta_plot) +
  geom_sf(aes(fill = school_level_vax_rate), color = "gray50", linewidth = 0.15) +
  county_borders +
  school_overlay +
  vax_fill +
  vax_color +
  labs(
    title = "School-Level Vaccination Coverage",
    subtitle = "Enrollment-weighted average of school coverage per ZCTA"
  ) +
  theme_map()

p_school
# --- Map B: Household-correlated ---
p_household <- ggplot(zcta_plot) +
  geom_sf(aes(fill = hh_correlated_vax_rate), color = "gray50", linewidth = 0.15) +
  county_borders +
  school_overlay +
  vax_fill +
  vax_color +
  labs(
    title = "Household-Correlated Vaccination Coverage",
    subtitle = "ZCTA level vaccination coverage after household-level assignment"
  ) +
  theme_map()

# --- Combined side-by-side ---
p_combined <- p_school + p_household +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

p_combined



# --- Difference map ---
p_diff <- ggplot(zcta_plot) +
  geom_sf(aes(fill = difference), color = "gray50", linewidth = 0.15) +
  county_borders +
  school_overlay +
  scale_fill_distiller(
    palette = "RdBu", direction = -1,
    limits = c(-1, 1) * max(abs(zcta_plot$difference), na.rm = TRUE),
    labels = function(x) sprintf("%+.0f pp", x * 100),
    name = "Difference\n(HH \u2212 School)",
    na.value = "gray85"
  ) +
  vax_color +
  labs(
    title = "Difference: Household-Correlated vs. School-Level Coverage",
    subtitle = "Negative (red) = household clustering reduces effective coverage below school-level estimate"
  ) +
  theme_map()

p_diff

# --- Scatter plot ---
p_scatter <- ggplot(comparison_df %>% filter(!is.na(school_level_vax_rate),
                                              !is.na(hh_correlated_vax_rate)),
                    aes(x = school_level_vax_rate * 100,
                        y = hh_correlated_vax_rate * 100)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = 95, linetype = "dotted", color = "#B2182B", linewidth = 0.5) +
  geom_vline(xintercept = 95, linetype = "dotted", color = "#B2182B", linewidth = 0.5) +
  annotate("text", x = 52, y = 96, label = "95% herd immunity threshold",
           size = 3, color = "#B2182B", hjust = 0) +
  geom_point(aes(size = n_children), alpha = 0.6, color = "#D73027") +
  scale_size_continuous(range = c(1.5, 8), name = "Children\nin ZCTA") +
  labs(
    title = "ZCTA-Level Vaccination Coverage Comparison",
    subtitle = "Points below diagonal: household clustering reduces effective coverage",
    x = "School-Level Estimate (%)",
    y = "Household-Correlated Estimate (%)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  coord_equal()

p_scatter
