# ==============================================================================
# PLOTTING UTILITIES
# ==============================================================================
# File: plotting_utils.R
# Contains: All visualization functions for the simulation
# Dependencies: ggplot2, dplyr, tidyr, viridis, scales
# ==============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(scales)

# ==============================================================================
# Network Outbreak Dynamics Plot
# ==============================================================================

#' Plot network-level outbreak dynamics
#' @param results Results from run_multiple_network_simulations()
#' @param show_ci Show confidence intervals (default TRUE)
#' @return ggplot object
plot_network_outbreak_dynamics <- function(results, show_ci = TRUE) {
  
  daily_summary <- results$all_network_data %>%
    group_by(day) %>%
    summarise(
      mean_infected = mean(E + P + Ra + Iso + R),
      median_infected = median(E + P + Ra + Iso + R),
      lower = quantile(E + P + Ra + Iso + R, 0.025),
      upper = quantile(E + P + Ra + Iso + R, 0.975),
      .groups = "drop"
    )
  
  p <- ggplot(daily_summary, aes(x = day)) +
    geom_line(aes(y = median_infected), color = "darkred", linewidth = 1.2)
  
  if (show_ci) {
    p <- p + geom_ribbon(aes(ymin = lower, ymax = upper), 
                         fill = "darkred", alpha = 0.2)
  }
  
  p <- p +
    labs(
      title = "Network-wide Outbreak Dynamics",
      subtitle = paste0("Median with 95% CI across ", results$n_simulations, " simulations"),
      x = "Day",
      y = "Cumulative Infected"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  return(p)
}


# ==============================================================================
# Outbreak Size Histogram
# ==============================================================================

#' Plot histogram of outbreak sizes
#' @param results Results from run_multiple_network_simulations()
#' @param bins Number of bins
#' @return ggplot object
plot_outbreak_size_histogram <- function(results, bins = 30) {
  
  p <- ggplot(results$summary_stats, aes(x = total_infected)) +
    geom_histogram(bins = bins, fill = "steelblue", color = "white", alpha = 0.8) +
    geom_vline(aes(xintercept = median(total_infected)), 
               color = "darkred", linetype = "dashed", linewidth = 1) +
    labs(
      title = "Distribution of Outbreak Sizes",
      subtitle = paste0("Red line = Median (", round(median(results$summary_stats$total_infected), 1), ")"),
      x = "Total Infected Across Network",
      y = "Frequency"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  return(p)
}


# ==============================================================================
# School Outbreak Dynamics (All Schools)
# ==============================================================================

#' Plot outbreak size by school (all schools)
#' @param results Results from run_multiple_network_simulations()
#' @return ggplot object
plot_school_outbreak_dynamics <- function(results) {
  
  # Ensure school_name exists
  if ("school_name" %in% names(results$all_school_data)) {
    school_plot_data <- results$all_school_data
  } else {
    school_plot_data <- results$all_school_data %>%
      left_join(results$schools %>% select(school_id, school_name), by = "school_id")
  }
  
  if (!"school_name" %in% names(school_plot_data) || all(is.na(school_plot_data$school_name))) {
    school_plot_data$school_name <- paste0("School_", school_plot_data$school_id)
  }
  
  # Create labels - only add ID if there are duplicate school names
  name_counts <- table(unique(school_plot_data[, c("school_id", "school_name")])$school_name)
  duplicate_names <- names(name_counts[name_counts > 1])
  
  school_plot_data <- school_plot_data %>%
    mutate(school_label = ifelse(
      school_name %in% duplicate_names,
      paste0(school_name, " (", school_id, ")"),
      school_name
    ))
  
  p <- ggplot(school_plot_data, aes(x = reorder(school_label, total_infected, FUN = median), 
                                    y = total_infected, fill = was_seeded)) +
    geom_boxplot(alpha = 0.7, outlier.size = 1) +
    coord_flip() +
    scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "darkred"),
                      labels = c("FALSE" = "Not Seeded", "TRUE" = "Seeded"),
                      name = "Outbreak Origin") +
    labs(
      title = "Outbreak Size by School",
      x = "School",
      y = "Total Infected"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 7)
    )
  
  return(p)
}


# ==============================================================================
# Top N Schools by Median Infections
# ==============================================================================

#' Plot outbreak size for top N schools by median infections
#' @param results Results from run_multiple_network_simulations()
#' @param top_n Number of top schools to display (default 20)
#' @return ggplot object
plot_top_schools_outbreak <- function(results, top_n = 20) {
  
  # Ensure school_name exists
  if ("school_name" %in% names(results$all_school_data)) {
    school_plot_data <- results$all_school_data
  } else {
    school_plot_data <- results$all_school_data %>%
      left_join(results$schools %>% select(school_id, school_name), by = "school_id")
  }
  
  if (!"school_name" %in% names(school_plot_data) || all(is.na(school_plot_data$school_name))) {
    school_plot_data$school_name <- paste0("School_", school_plot_data$school_id)
  }
  
  # Calculate median infections per school
  school_medians <- school_plot_data %>%
    group_by(school_id, school_name) %>%
    summarise(
      median_infected = median(total_infected),
      was_seeded = first(was_seeded),
      .groups = "drop"
    ) %>%
    arrange(desc(median_infected)) %>%
    head(top_n)
  
  # Filter to top schools
  top_school_ids <- school_medians$school_id
  
  school_plot_filtered <- school_plot_data %>%
    filter(school_id %in% top_school_ids)
  
  # Create labels - only add ID if there are duplicate school names
  name_counts <- table(unique(school_plot_filtered[, c("school_id", "school_name")])$school_name)
  duplicate_names <- names(name_counts[name_counts > 1])
  
  school_plot_filtered <- school_plot_filtered %>%
    mutate(school_label = ifelse(
      school_name %in% duplicate_names,
      paste0(school_name, " (", school_id, ")"),
      school_name
    ))
  
  # Order by median
  school_order <- school_plot_filtered %>%
    group_by(school_label) %>%
    summarise(median_inf = median(total_infected), .groups = "drop") %>%
    arrange(median_inf) %>%
    pull(school_label)
  
  school_plot_filtered <- school_plot_filtered %>%
    mutate(school_label = factor(school_label, levels = school_order))
  
  p <- ggplot(school_plot_filtered, aes(x = school_label, y = total_infected, fill = was_seeded)) +
    geom_boxplot(alpha = 0.7, outlier.size = 1) +
    coord_flip() +
    scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "darkred"),
                      labels = c("FALSE" = "Not Seeded", "TRUE" = "Seeded"),
                      name = "Outbreak Origin") +
    labs(
      title = paste0("Top ", top_n, " Schools by Median Infections"),
      x = "School",
      y = "Total Infected"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 9)
    )
  
  return(p)
}


# ==============================================================================
# Disease State Dynamics Plot
# ==============================================================================

#' Plot disease state dynamics over time
#' @param results Results from run_multiple_network_simulations()
#' @return ggplot object
plot_disease_states <- function(results) {
  
  daily_summary <- results$all_network_data %>%
    group_by(day) %>%
    summarise(
      Exposed = mean(E), 
      Prodromal = mean(P), 
      Rash = mean(Ra),
      Isolated = mean(Iso), 
      Recovered = mean(R),
      .groups = "drop"
    ) %>%
    pivot_longer(-day, names_to = "State", values_to = "Count")
  
  p <- ggplot(daily_summary, aes(x = day, y = Count, color = State)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = c(
      "Exposed" = "orange", 
      "Prodromal" = "darkorange", 
      "Rash" = "red", 
      "Isolated" = "purple", 
      "Recovered" = "darkgreen"
    )) +
    labs(
      title = "Disease State Dynamics",
      x = "Day",
      y = "Mean Count"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  return(p)
}


# ==============================================================================
# Infection Sequence Plot (Median First Infection Time)
# ==============================================================================

#' Plot median first infection time by school (infection sequence)
#' @param results Results object from run_multiple_network_simulations
#' @param show_error_bars Show IQR error bars (default TRUE)
#' @param highlight_seed Highlight seed school(s) (default TRUE)
#' @param seed_schools Vector of seed school IDs
#' @param top_n Only show top N schools by infection probability (NULL for all)
#' @return ggplot object
plot_infection_sequence <- function(results, show_error_bars = TRUE, 
                                    highlight_seed = TRUE, seed_schools = NULL,
                                    top_n = NULL) {
  
  # Calculate timing statistics
  first_times <- calculate_first_infection_times(results)
  timing_summary <- summarize_infection_timing(first_times)
  
  # Filter to schools that got infected at least once
  timing_summary <- timing_summary %>%
    filter(prob_infected > 0)
  
  # Optionally limit to top N schools
  if (!is.null(top_n) && top_n < nrow(timing_summary)) {
    timing_summary <- timing_summary %>%
      arrange(desc(prob_infected), median_first_day) %>%
      head(top_n)
  }
  
  # Add seed indicator
  if (highlight_seed && !is.null(seed_schools)) {
    timing_summary <- timing_summary %>%
      mutate(is_seed = school_id %in% seed_schools)
  } else {
    timing_summary$is_seed <- FALSE
  }
  
  # Create labels - only add ID if there are duplicate school names
  name_counts <- table(timing_summary$school_name)
  duplicate_names <- names(name_counts[name_counts > 1])
  
  timing_summary <- timing_summary %>%
    mutate(school_label = ifelse(
      school_name %in% duplicate_names,
      paste0(school_name, " (", school_id, ")"),
      school_name
    ))
  
  # Order by median first infection day
  school_order <- timing_summary %>%
    arrange(median_first_day) %>%
    pull(school_label)
  
  timing_summary <- timing_summary %>%
    mutate(school_label = factor(school_label, levels = school_order))
  
  # Create plot
  p <- ggplot(timing_summary, aes(x = median_first_day, y = school_label)) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 8),
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
  
  # Add error bars for IQR
  if (show_error_bars) {
    p <- p + geom_errorbarh(
      aes(xmin = q25_first_day, xmax = q75_first_day),
      height = 0.3, color = "gray50", alpha = 0.7
    )
  }
  
  # Add points
  if (highlight_seed && any(timing_summary$is_seed)) {
    p <- p + geom_point(aes(color = is_seed, size = prob_infected)) +
      scale_color_manual(
        values = c("FALSE" = "steelblue", "TRUE" = "red"),
        labels = c("FALSE" = "Not Seeded", "TRUE" = "Seed School"),
        name = "Outbreak Origin"
      )
  } else {
    p <- p + geom_point(aes(size = prob_infected), color = "steelblue")
  }
  
  p <- p +
    scale_size_continuous(
      name = "Prob. Infected",
      range = c(2, 5),
      labels = percent
    ) +
    labs(
      title = "Infection Spread Sequence Across Schools",
      subtitle = "Median first infection day with IQR (ordered by timing)",
      x = "Median Day of First Infection",
      y = "School"
    )
  
  return(p)
}


# ==============================================================================
# Infection Timing Heatmap
# ==============================================================================

#' Plot heatmap of infection arrival times across simulations
#' @param results Results object from run_multiple_network_simulations
#' @param max_sims Maximum number of simulations to show (default 50)
#' @param seed_schools Vector of seed school IDs for ordering
#' @return ggplot object
plot_infection_timing_heatmap <- function(results, max_sims = 50, seed_schools = NULL) {
  
  first_times <- calculate_first_infection_times(results)
  
  # Limit number of simulations for readability
  sims_to_show <- unique(first_times$sim)
  if (length(sims_to_show) > max_sims) {
    sims_to_show <- sims_to_show[1:max_sims]
    first_times <- first_times %>% filter(sim %in% sims_to_show)
  }
  
  # Check for duplicate school names
  name_counts <- table(unique(first_times[, c("school_id", "school_name")])$school_name)
  duplicate_names <- names(name_counts[name_counts > 1])
  
  # Create labels - only add ID if there are duplicate school names
  first_times <- first_times %>%
    mutate(school_label = ifelse(
      school_name %in% duplicate_names,
      paste0(school_name, " (", school_id, ")"),
      school_name
    ))
  
  # Calculate median timing for ordering schools
  school_order <- first_times %>%
    group_by(school_id, school_label) %>%
    summarise(median_day = median(first_infection_day, na.rm = TRUE), .groups = "drop") %>%
    arrange(median_day) %>%
    pull(school_label)
  
  # Remove duplicates from school_order
  school_order <- unique(school_order)
  
  # Prepare data for heatmap
  heatmap_data <- first_times %>%
    mutate(
      school_label = factor(school_label, levels = school_order),
      first_infection_day = ifelse(is.na(first_infection_day), Inf, first_infection_day)
    )
  
  # Calculate upper limit for color scale
  finite_days <- heatmap_data$first_infection_day[is.finite(heatmap_data$first_infection_day)]
  upper_limit <- if (length(finite_days) > 0) quantile(finite_days, 0.95) else 100
  
  # Create heatmap
  p <- ggplot(heatmap_data, aes(x = factor(sim), y = school_label, fill = first_infection_day)) +
    geom_tile(color = "white", linewidth = 0.1) +
    scale_fill_viridis_c(
      name = "Day of\nFirst Infection",
      option = "plasma",
      na.value = "gray90",
      limits = c(0, upper_limit),
      oob = squish
    ) +
    labs(
      title = "Infection Arrival Time by School and Simulation",
      subtitle = "Schools ordered by median infection timing",
      x = "Simulation",
      y = "School"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 7),
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
  
  return(p)
}


# ==============================================================================
# Network Visualization
# ==============================================================================

#' Plot network visualization with outbreak data
#' @param network Network object
#' @param school_summary Summary data from a single simulation
#' @return Plot (base R graphics)
plot_network_with_outbreak <- function(network, school_summary) {
  g <- network$graph
  
  # Color nodes by infection level
  max_infected <- max(school_summary$total_infected, na.rm = TRUE)
  if (max_infected > 0) {
    colors <- colorRampPalette(c("lightgreen", "yellow", "red"))(100)
    node_colors <- colors[pmin(100, ceiling(school_summary$total_infected / max_infected * 99) + 1)]
  } else {
    node_colors <- rep("lightgreen", nrow(school_summary))
  }
  
  # Size nodes by school size
  node_sizes <- 5 + 15 * (school_summary$school_size / max(school_summary$school_size))
  
  # Plot
  plot(g, 
       vertex.color = node_colors,
       vertex.size = node_sizes,
       vertex.label = NA,
       edge.width = E(g)$weight * 2,
       edge.color = "gray70",
       main = "School Network (colored by infection level)")
  
  # Add legend
  legend("bottomright", 
         legend = c("Low", "Medium", "High"),
         fill = c("lightgreen", "yellow", "red"),
         title = "Infections",
         cex = 0.8)
}


# ==============================================================================
# HIGH-QUALITY REPORT FIGURES
# ==============================================================================

#' Plot top schools by median infections - Publication quality for reports
#' @param results Results from run_multiple_network_simulations()
#' @param top_n Number of top schools to display (default 20)
#' @param title Custom title (optional)
#' @param subtitle Custom subtitle (optional)
#' @return ggplot object suitable for reports
plot_top_schools_report <- function(results, top_n = 20, 
                                    title = NULL, 
                                    subtitle = NULL) {
  
  # Ensure school_name exists
  if ("school_name" %in% names(results$all_school_data)) {
    school_plot_data <- results$all_school_data
  } else {
    school_plot_data <- results$all_school_data %>%
      left_join(results$schools %>% select(school_id, school_name), by = "school_id")
  }
  
  if (!"school_name" %in% names(school_plot_data) || all(is.na(school_plot_data$school_name))) {
    school_plot_data$school_name <- paste0("School_", school_plot_data$school_id)
  }
  
  # Calculate median infections per school
  school_medians <- school_plot_data %>%
    group_by(school_id, school_name) %>%
    summarise(
      median_infected = median(total_infected),
      was_seeded = first(was_seeded),
      .groups = "drop"
    ) %>%
    arrange(desc(median_infected)) %>%
    head(top_n)
  
  # Filter to top schools
  top_school_ids <- school_medians$school_id
  
  school_plot_filtered <- school_plot_data %>%
    filter(school_id %in% top_school_ids)
  
  # Create labels - only add ID if there are duplicate school names
  name_counts <- table(unique(school_plot_filtered[, c("school_id", "school_name")])$school_name)
  duplicate_names <- names(name_counts[name_counts > 1])
  
  school_plot_filtered <- school_plot_filtered %>%
    mutate(school_label = ifelse(
      school_name %in% duplicate_names,
      paste0(school_name, " (", school_id, ")"),
      school_name
    ))
  
  # Create intuitive labels for public health audience
  school_plot_filtered <- school_plot_filtered %>%
    mutate(outbreak_type = ifelse(was_seeded, 
                                  "Initial Outbreak School", 
                                  "Secondary Spread"))
  
  # Order by median
  school_order <- school_plot_filtered %>%
    group_by(school_label) %>%
    summarise(median_inf = median(total_infected), .groups = "drop") %>%
    arrange(median_inf) %>%
    pull(school_label)
  
  school_plot_filtered <- school_plot_filtered %>%
    mutate(school_label = factor(school_label, levels = school_order))
  
  # Set default titles
  if (is.null(title)) {
    title <- paste0("Simulated Infection Burden: Top ", top_n, " Most Affected Schools")
  }
  if (is.null(subtitle)) {
    subtitle <- paste0("Distribution across simulations (box = IQR, line = median)")
  }
  
  # Create publication-quality plot
  p <- ggplot(school_plot_filtered, aes(x = school_label, y = total_infected, 
                                        fill = outbreak_type)) +
    geom_boxplot(alpha = 0.8, outlier.size = 1.5, outlier.alpha = 0.5) +
    coord_flip() +
    scale_fill_manual(
      values = c("Initial Outbreak School" = "#D62728", 
                 "Secondary Spread" = "#1F77B4"),
      name = "School Type"
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = "Total Infections per Simulation",
      caption = paste0("Note: Results reflect simulated outbreaks with isolation and quarantine interventions.\n",
                       "Initial outbreak school = where first case was introduced.")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      # Title formatting
      plot.title = element_text(face = "bold", size = 14, hjust = 0),
      plot.subtitle = element_text(size = 11, color = "gray30", hjust = 0),
      plot.caption = element_text(size = 9, color = "gray50", hjust = 0),
      
      # Axis formatting
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 11, margin = margin(t = 10)),
      
      # Legend formatting - larger and more visible
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 10),
      legend.key.size = unit(1.2, "lines"),
      legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.5),
      legend.margin = margin(t = 10, b = 5, l = 10, r = 10),
      
      # Panel formatting
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "gray85"),
      
      # Plot margins
      plot.margin = margin(t = 15, r = 15, b = 10, l = 10)
    )
  
  return(p)
}


#' Save plot in high-quality format for reports
#' @param plot ggplot object
#' @param filename Output filename (without extension)
#' @param width Width in inches (default 10)
#' @param height Height in inches (default 8)
#' @param dpi Resolution (default 300)
#' @param format Output format: "png", "pdf", "tiff", or "all" (default "png")
#' @return Invisible NULL, saves file(s) to disk
save_report_figure <- function(plot, filename, width = 10, height = 8, 
                               dpi = 300, format = "png") {
  
  if (format == "all" || format == "png") {
    ggsave(paste0(filename, ".png"), plot, width = width, height = height, 
           dpi = dpi, bg = "white")
    cat(sprintf("Saved: %s.png\n", filename))
  }
  
  if (format == "all" || format == "pdf") {
    ggsave(paste0(filename, ".pdf"), plot, width = width, height = height, 
           device = cairo_pdf)
    cat(sprintf("Saved: %s.pdf\n", filename))
  }
  
  if (format == "all" || format == "tiff") {
    ggsave(paste0(filename, ".tiff"), plot, width = width, height = height, 
           dpi = dpi, compression = "lzw")
    cat(sprintf("Saved: %s.tiff\n", filename))
  }
  
  invisible(NULL)
}


#' Plot infection sequence - Publication quality for reports
#' @param results Results object from run_multiple_network_simulations
#' @param seed_schools Vector of seed school IDs
#' @param top_n Only show top N schools by infection probability (NULL for all)
#' @param title Custom title (optional)
#' @return ggplot object suitable for reports
plot_infection_sequence_report <- function(results, seed_schools = NULL,
                                           top_n = 30, title = NULL) {
  
  # Calculate timing statistics
  first_times <- calculate_first_infection_times(results)
  timing_summary <- summarize_infection_timing(first_times)
  
  # Filter to schools that got infected at least once
  timing_summary <- timing_summary %>%
    filter(prob_infected > 0)
  
  # Optionally limit to top N schools
  if (!is.null(top_n) && top_n < nrow(timing_summary)) {
    timing_summary <- timing_summary %>%
      arrange(desc(prob_infected), median_first_day) %>%
      head(top_n)
  }
  
  # Add seed indicator with intuitive labels
  if (!is.null(seed_schools)) {
    timing_summary <- timing_summary %>%
      mutate(school_type = ifelse(school_id %in% seed_schools,
                                  "Initial Outbreak School",
                                  "Secondary Spread"))
  } else {
    timing_summary$school_type <- "Secondary Spread"
  }
  
  # Create labels - only add ID if there are duplicate school names
  name_counts <- table(timing_summary$school_name)
  duplicate_names <- names(name_counts[name_counts > 1])
  
  timing_summary <- timing_summary %>%
    mutate(school_label = ifelse(
      school_name %in% duplicate_names,
      paste0(school_name, " (", school_id, ")"),
      school_name
    ))
  
  # Order by median first infection day
  school_order <- timing_summary %>%
    arrange(median_first_day) %>%
    pull(school_label)
  
  timing_summary <- timing_summary %>%
    mutate(school_label = factor(school_label, levels = school_order))
  
  # Set default title
  if (is.null(title)) {
    title <- "Timing of Outbreak Spread Across Schools"
  }
  
  # Create plot
  p <- ggplot(timing_summary, aes(x = median_first_day, y = school_label)) +
    geom_errorbarh(
      aes(xmin = q25_first_day, xmax = q75_first_day),
      height = 0.3, color = "gray50", alpha = 0.7, linewidth = 0.8
    ) +
    geom_point(aes(color = school_type, size = prob_infected)) +
    scale_color_manual(
      values = c("Initial Outbreak School" = "#D62728", 
                 "Secondary Spread" = "#1F77B4"),
      name = "School Type"
    ) +
    scale_size_continuous(
      name = "Probability of\nInfection",
      range = c(3, 7),
      labels = scales::percent
    ) +
    labs(
      title = title,
      subtitle = "Median day of first infection with interquartile range",
      x = "Day of First Infection (from outbreak start)",
      y = NULL,
      caption = "Point size indicates probability that the school experienced any infections across simulations."
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0),
      plot.subtitle = element_text(size = 11, color = "gray30", hjust = 0),
      plot.caption = element_text(size = 9, color = "gray50", hjust = 0),
      
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 11, margin = margin(t = 10)),
      
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9),
      legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.5),
      legend.margin = margin(t = 10, b = 5, l = 10, r = 10),
      
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      
      plot.margin = margin(t = 15, r = 15, b = 10, l = 10)
    ) +
    guides(
      color = guide_legend(order = 1, override.aes = list(size = 4)),
      size = guide_legend(order = 2)
    )
  
  return(p)
}


#' Plot outbreak dynamics - Publication quality for reports
#' @param results Results from run_multiple_network_simulations()
#' @param title Custom title (optional)
#' @return ggplot object suitable for reports
plot_outbreak_dynamics_report <- function(results, title = NULL) {
  
  daily_summary <- results$all_network_data %>%
    group_by(day) %>%
    summarise(
      mean_infected = mean(E + P + Ra + Iso + R),
      median_infected = median(E + P + Ra + Iso + R),
      lower = quantile(E + P + Ra + Iso + R, 0.025),
      upper = quantile(E + P + Ra + Iso + R, 0.975),
      .groups = "drop"
    )
  
  if (is.null(title)) {
    title <- "Cumulative Infections Over Time"
  }
  
  p <- ggplot(daily_summary, aes(x = day)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#1F77B4", alpha = 0.2) +
    geom_line(aes(y = median_infected), color = "#1F77B4", linewidth = 1.2) +
    labs(
      title = title,
      subtitle = paste0("Median with 95% uncertainty interval (", 
                        results$n_simulations, " simulations)"),
      x = "Days Since First Case",
      y = "Cumulative Infections",
      caption = "Shaded area represents the range of outcomes in 95% of simulations."
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0),
      plot.subtitle = element_text(size = 11, color = "gray30", hjust = 0),
      plot.caption = element_text(size = 9, color = "gray50", hjust = 0),
      
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11),
      axis.title.y = element_text(margin = margin(r = 10)),
      axis.title.x = element_text(margin = margin(t = 10)),
      
      panel.grid.minor = element_blank(),
      
      plot.margin = margin(t = 15, r = 15, b = 10, l = 10)
    )
  
  return(p)
}

cat("Plotting utilities loaded.\n")
