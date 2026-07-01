# measles-household-school-cluster

An agent-based model (ABM) for simulating measles transmission across interconnected schools and households. The model captures within-school spread, between-school spread via a geographic network, and secondary household transmission — supporting both prospective outbreak projection and real-time mid-outbreak forecasting.

---

## Contents

- [Overview](#overview)
- [Prerequisites and Installation](#prerequisites-and-installation)
- [Data Requirements](#data-requirements)
- [Repository Structure](#repository-structure)
- [Model Description](#model-description)
  - [Population Structure](#population-structure)
  - [Disease States](#disease-states)
  - [Transmission Pathways](#transmission-pathways)
  - [Interventions](#interventions)
  - [School Network](#school-network)
- [Key Parameters](#key-parameters)
- [Cluster Simulation Framework](#cluster-simulation-framework)
  - [Loading Modules](#loading-modules)
  - [Preparing Data](#preparing-data)
  - [Sequential Simulations](#sequential-simulations)
  - [Parallel Simulations (Single Node)](#parallel-simulations-single-node)
  - [HPC / SLURM Batch Execution](#hpc--slurm-batch-execution)
  - [Combining Batch Results](#combining-batch-results)
- [Simulation Modes](#simulation-modes)
  - [Clean Outbreak](#clean-outbreak)
  - [Mid-Outbreak Forecasting](#mid-outbreak-forecasting)
- [Outputs](#outputs)
- [Running Tests](#running-tests)
- [Citation](#citation)

---

## Overview

This repository provides a modular, high-performance simulation framework for studying measles outbreaks in school-based populations. Designed for both desktop experimentation and large-scale HPC runs, the framework enables:

- **Prospective outbreak projection** — starting from a single index case at one or more seed schools
- **Real-time mid-outbreak forecasting** — initializing from a partially observed in-progress outbreak
- **Intervention evaluation** — comparing isolation, contact quarantine, and vaccination strategies
- **Network structure analysis** — exploring how school connectivity shapes outbreak dynamics

The core simulation engine is implemented in C++ via Rcpp for performance, with an R orchestration layer for data preparation, parallel execution, and analysis.

---

## Prerequisites and Installation

### R Version

R ≥ 4.1.0 is required. R 4.5.0 is recommended for HPC environments.

### Required R Packages

```r
install.packages(c(
  "Rcpp",      # C++ integration for the transmission engine
  "dplyr",     # Data manipulation
  "tidyr",     # Data reshaping
  "ggplot2",   # Plotting
  "R6",        # Reference classes (contact history)
  "igraph",    # School network construction
  "parallel",  # Parallel computation (base R, no install needed)
  "readxl",    # Reading vaccination data from Excel
  "here"       # Relative file paths in scripts
))
```

On an HPC cluster, install once into your user library before submitting jobs:

```bash
module load r/4.5.0
export R_LIBS_USER="$HOME/R/4.5"
mkdir -p "$R_LIBS_USER"
R -q -e 'install.packages(c("Rcpp","dplyr","tidyr","R6","igraph","ggplot2","readxl","here"), repos="https://cloud.r-project.org")'
```

A ready-made SLURM installation job is provided at `codes/install_r_deps.sbatch`:

```bash
sbatch codes/install_r_deps.sbatch
```

---

## Data Requirements

All data files are expected under the `data/` directory. The following files are used by the batch scripts and example scripts:

| File | Description |
|------|-------------|
| `data/SC_vaccination_merged.csv` | School-level vaccination coverage and enrollment. Required columns: `School.Name`, `County`, `Total.Students`, `Percent.Immunized`, `longitude`, `latitude`, `Grade.Range`. |
| `data/School_reference_Upstate.csv` | School reference table linking school names to unique IDs (`OBJECTID`) used in the travel time matrix. Required columns: `OBJECTID` (or `ID`), `School.Name`, `County`. |
| `data/Travel_time_matrix_Upstate.csv` | Precomputed pairwise driving-time matrix (minutes) between schools, indexed by `OBJECTID`. Used to build the travel-time network. |
| `data/household_assignment_cache.rds` | Cached output of `assign_households_to_students()`. Generated automatically on first run and reloaded on subsequent runs to avoid recomputing the expensive assignment step. |

> **Synthetic population:** Household assignment requires an RTI synthetic population file (CSV) with columns `hh_id`, `agep`, `person_id`, `hh_size`, `lon_4326`, `lat_4326`. Update the `synpop_file` path in `run_simulation_households.R` to point to your local copy. This file is not included in the repository due to size and licensing.

If the travel time matrix is unavailable, the scripts automatically fall back to a Haversine-distance network using school coordinates from the vaccination CSV.

---

## Repository Structure

```
codes/
├── load_all.R                        # Sources all modules in dependency order
├── rcpp_transmission.R               # C++ (Rcpp) within-school and between-school
│                                     #   transmission engine
├── network_utils.R                   # School network construction (distance,
│                                     #   travel time, small-world, scale-free)
├── grade_utils.R                     # Grade range parsing and school-type
│                                     #   classification
├── population_utils.R                # School population creation and SEIR
│                                     #   state updates
├── household_utils.R                 # Household loading, student assignment,
│                                     #   household transmission
├── simulation_utils.R                # Core simulation engine (single run and
│                                     #   sequential multi-run)
├── simulation_parallel.R             # Parallel execution across CPU cores
│                                     #   (run_parallel_simulations,
│                                     #    run_fast_simulations)
├── analysis_utils.R                  # Post-simulation summary statistics
├── plotting_utils.R                  # Epidemic curve and outbreak visualization
├── midoutbreak_utils.R               # Mid-outbreak scenario initialization
├── between_school_transmission_fix.R # Corrected between-school contact sampling
├── rcpp_within_school_fix.R          # Corrected within-school contact sampling
├── run_simulation_households.R       # End-to-end interactive example script
├── run_batch_clean.R                 # SLURM batch script — clean outbreak runs
├── run_batch_midoutbreak.R           # SLURM batch script — mid-outbreak runs
├── combine_results.R                 # Aggregates output from SLURM array tasks
├── vaccination_spatial_comparison.R  # ZCTA-level spatial vaccination analysis
└── install_r_deps.sbatch             # SLURM job for installing R dependencies
data/
├── SC_vaccination_merged.csv
├── School_reference_Upstate.csv
├── Travel_time_matrix_Upstate.csv
└── household_assignment_cache.rds    # Generated on first run
tests/
└── testthat/                         # Unit tests (testthat framework)
```

---

## Model Description

### Population Structure

Students aged 5–18 are assigned to schools based on age-appropriate grade range and geographic proximity to their household. Each school is divided into classes (default size 25 students). Students are linked to synthetic households, which also contain adults and children under age 5. Non-student household members can acquire secondary infections from infectious students but do not re-enter the school simulation.

Vaccination status is assigned at the household level: if one child in a household is unvaccinated, siblings have an elevated probability of also being unvaccinated (default 80%), and adults have a moderately elevated probability (default 50%), reflecting correlated health-seeking behavior.

### Disease States

The model uses a modified SEIR compartmental structure with additional states for interventions:

| State | Label | Description |
|-------|-------|-------------|
| Susceptible | **S** | No immunity; can be infected |
| Exposed | **E** | Latent infection; not yet infectious (Erlang-distributed, mean 10 days) |
| Prodromal | **P** | Early infectious stage with non-specific symptoms (~4 days) |
| Rash | **Ra** | Fully infectious with rash (~4 days) |
| Recovered | **R** | Immune; cannot be reinfected |
| Vaccinated | **V** | 97% fully immune; 3% leaky (can be infected with 80% reduced infectiousness) |
| Isolated | **Iso** | Detected and removed from school contacts; non-infectious |
| Quarantined | **Q\*** | Contact of an isolated case; held until symptom-onset confirmation |

### Transmission Pathways

**Within-school transmission** is modeled using C++ (Rcpp) for performance. Each day, every infectious student in state P or Ra (and not isolated) generates contacts within their class at rate `c_within` and across other classes at rate `c_between`. Contacts are drawn stochastically; transmission occurs if the contact is susceptible or has leaky vaccine status.

**Between-school transmission** occurs via a weighted school network. Each day, infectious students generate cross-school contacts scaled by the edge weight between schools and the base rate `c_between_school`. Contact pools at the target school are drawn from all non-isolated students.

**Household transmission** is applied daily: each infectious student (P or Ra) has a probability `hh_transmission_prob` of transmitting to each susceptible or leaky-vaccinated household member. Household infections are tracked separately from the school simulation.

### Interventions

| Intervention | Behavior |
|---|---|
| **Isolation** | Detected cases are isolated after `isolation_delay_index` days (first case) or `isolation_delay_secondary` days (subsequent cases). Isolation lasts `isolation_period` days (default 14). |
| **Contact quarantine** | Contacts of isolated individuals are quarantined for `quarantine_duration` days (default 21) with compliance rate `quarantine_efficacy` (default 90%). |
| **Vaccination** | Coverage is school-specific and loaded from data. Vaccine efficacy is modeled as leaky: 97% of vaccinated individuals are fully protected; 3% can be infected but transmit with 80% reduced infectiousness. |
| **No-intervention mode** | Setting `no_intervention = TRUE` disables both isolation and contact quarantine, useful for estimating the unmitigated outbreak potential. |

### School Network

Schools are nodes in an undirected weighted graph. Four network types are supported:

| Network Type | Construction | Use Case |
|---|---|---|
| **Travel-time-based** | Edge weight is an inverse-exponential function of driving time (minutes) | Realistic geographic connectivity; preferred for South Carolina data |
| **Distance-based** | Edge weight is an inverse-exponential function of Haversine distance (km) | Fallback when travel time data is unavailable |
| **Small-world (Watts-Strogatz)** | Clustered local connections with random long-range rewiring | Theoretical analysis of clustered networks |
| **Scale-free (Barabási-Albert)** | Preferential attachment; hubs drive most cross-school spread | Theoretical analysis of hub-dominated transmission |

A **grade-weighted overlay** can be applied on top of any base network: schools with overlapping grade ranges (indicating shared community ties, sibling links, or joint activities) receive higher edge weights. The overlap strength is controlled by `grade_weight` (0–1) and `same_type_bonus`.

---

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `latent_mean` | 10 days | Mean latent (exposed) period |
| `latent_shape` | 4 | Erlang shape for latent period |
| `infectious_mean` | 8 days | Total mean infectious period |
| `infectious_shape` | 4 | Erlang shape for infectious period |
| `prodromal_period` | 4 days | Duration of prodromal (P) stage |
| `c_within` | 6 | Within-class daily contacts per student |
| `c_between` | 2 | Between-class daily contacts per student |
| `p_within` | 0.17 | Per-contact transmission probability (within class) |
| `p_between` | 0.17 | Per-contact transmission probability (between classes) |
| `c_between_school` | 0.2 | Base between-school contact rate (scaled by network weight) |
| `hh_transmission_prob` | 0.17 | Daily per-contact household transmission probability |
| `vaccine_efficacy` | 0.97 | Fraction of vaccinated individuals fully protected |
| `vaccine_infectiousness_reduction` | 0.8 | Infectiousness reduction for leaky-vaccine breakthroughs |
| `isolation_delay_index` | 2 days | Detection delay for the first identified case |
| `isolation_delay_secondary` | 3 days | Detection delay for subsequent identified cases |
| `isolation_period` | 14 days | Duration of case isolation |
| `quarantine_efficacy` | 0.90 | Fraction of quarantined contacts who comply |
| `quarantine_duration` | 21 days | Length of contact quarantine |
| `avg_class_size` | 25 | Average number of students per class |
| `age_range` | c(5, 18) | Age range of the student population |

---

## Cluster Simulation Framework

The simulation framework is designed to scale from a single interactive R session to multi-node HPC array jobs. All execution paths share the same underlying `run_network_simulation()` engine, ensuring consistency of results.

### Loading Modules

All modules must be loaded before running simulations. The `load_all.R` script sources all dependencies in the correct order and compiles the C++ transmission code:

```r
# From the project root:
source("codes/load_all.R")
```

`load_all.R` automatically detects the module directory, installs missing packages into `R_LIBS_USER` if `ALLOW_R_INSTALL=true`, and verifies that all critical functions are available after loading.

### Preparing Data

Before running simulations, prepare the three core inputs:

```r
# 1. Schools data frame (school_id, school_name, school_size,
#    vaccination_coverage, lon, lat, county, grade_range_orig)
schools <- read.csv("data/SC_vaccination_merged.csv") |>
  dplyr::mutate(school_id = dplyr::row_number(), ...)

# 2. School network
network <- generate_travel_time_network(
  travel_time_file     = "data/Travel_time_matrix_Upstate.csv",
  school_reference_file = "data/School_reference_Upstate.csv",
  selected_school_ids  = schools$ref_id,
  max_travel_time      = 16,       # minutes
  weight_method        = "exponential"
)

# 3. School populations and household assignment
populations <- lapply(seq_len(nrow(schools)), function(i)
  create_school_population(i, schools$school_size[i],
                           params$avg_class_size, params$age_range))

hh_result <- assign_households_to_students(
  populations, schools, synpop,
  max_distance_km = 30, grade_tolerance = 2)

# Save for reuse (expensive to recompute)
save_household_assignment(hh_result, "data/household_assignment_cache.rds", schools)
```

Household assignment is the most time-consuming setup step. Once computed, it is cached as an `.rds` file and automatically reloaded on subsequent runs.

### Sequential Simulations

`run_fast_simulations()` runs simulations one at a time in the current R process. This is the simplest mode and requires no parallel infrastructure:

```r
results <- run_fast_simulations(
  n_simulations      = 50,
  schools            = schools,
  network            = network,
  params             = params,
  seed_schools       = seed_schools,   # integer vector of school IDs
  n_initial_infected = 5,
  n_days             = 365,
  seed_start         = 12345,
  verbose            = TRUE,
  hh_pop             = hh_pop,
  household_assignment = household_assignment
)
```

### Parallel Simulations (Single Node)

`run_parallel_simulations()` distributes simulations across CPU cores using R's `parallel` package. It auto-detects whether to use FORK (Linux/macOS) or PSOCK (Windows) clusters:

```r
results <- run_parallel_simulations(
  n_simulations      = 200,
  schools            = schools,
  network            = network,
  params             = params,
  seed_schools       = seed_schools,
  n_initial_infected = 5,
  n_days             = 365,
  seed_start         = 12345,
  n_cores            = NULL,           # NULL = auto-detect (recommended)
  cluster_type       = NULL,           # NULL = auto-detect platform
  verbose            = TRUE,
  hh_pop             = hh_pop,
  household_assignment = household_assignment
)
```

**Key arguments:**

| Argument | Default | Description |
|---|---|---|
| `n_simulations` | — | Total number of Monte Carlo replicates |
| `n_cores` | `NULL` | Number of parallel workers; `NULL` = auto-detect |
| `cluster_type` | `NULL` | `"FORK"` (Linux/macOS) or `"PSOCK"` (Windows); `NULL` = auto |
| `seed_start` | `NULL` | Base seed for reproducibility; each simulation receives `seed_start + i` |
| `hh_pop` | `NULL` | Household population data frame; omit to disable household transmission |

The function returns a named list with elements `results_list`, `summary_stats`, `computation_time`, and `use_households`.

A complete end-to-end example with data loading, network construction, household assignment, simulation, and plotting is provided in `codes/run_simulation_households.R`.

### HPC / SLURM Batch Execution

For large-scale runs (thousands of simulations), use the SLURM array batch scripts. Each task runs an independent subset of simulations, and results are combined afterwards.

**Step 1 — Install R dependencies (once):**

```bash
sbatch codes/install_r_deps.sbatch
```

**Step 2 — Submit a clean-outbreak array job:**

```bash
# Submit 10 tasks, each running 100 simulations on 12 cores
sbatch --array=1-10 --wrap="Rscript codes/run_batch_clean.R \
  --task_id=\$SLURM_ARRAY_TASK_ID \
  --n_sims=100 \
  --n_cores=12 \
  --n_days=365 \
  --output_dir=./results"
```

Or supply a task ID manually for testing:

```bash
Rscript codes/run_batch_clean.R \
  --task_id=1 \
  --n_sims=100 \
  --n_cores=12 \
  --n_days=365 \
  --output_dir=./results
```

**Step 3 — Submit a mid-outbreak forecast array job:**

```bash
sbatch --array=1-10 --wrap="Rscript codes/run_batch_midoutbreak.R \
  --task_id=\$SLURM_ARRAY_TASK_ID \
  --n_sims=100 \
  --n_cores=12 \
  --n_days=120 \
  --output_dir=./results_forecast"
```

**Batch script arguments:**

| Argument | Default | Description |
|---|---|---|
| `--task_id` | 1 | SLURM array task index; used to offset simulation seeds |
| `--n_sims` | 100 | Simulations to run in this task |
| `--n_cores` | 12 | CPU cores per task (should match `--cpus-per-task` in SBATCH) |
| `--n_days` | 365 | Maximum simulation duration in days |
| `--output_dir` | `./results` | Directory for output `.rds` files |

Each task writes four `.rds` files to `output_dir`:

| File | Contents |
|---|---|
| `task_XX_summary.rds` | One row per simulation: total infected, schools affected, duration, peak day, HH infections |
| `task_XX_curves.rds` | Daily network-wide counts by disease state (S, E, P, Ra, R, V, Iso) for every simulation |
| `task_XX_hh_curves.rds` | Daily household-level counts (when household transmission is enabled) |
| `task_XX_schools.rds` | Per-school totals: infected, attack rate, breakthrough infections, per simulation |

### Combining Batch Results

After all array tasks complete, combine the per-task output into a single file:

```bash
Rscript codes/combine_results.R <job_id>
# Example:
Rscript codes/combine_results.R 12345678
```

This reads all `task_XX_summary.rds` files from `results/<job_id>/`, concatenates them, and writes `combined_results.rds` and `combined_results.csv` to the same directory.

---

## Simulation Modes

### Clean Outbreak

The most common mode: an outbreak is seeded at one or more schools with a specified number of initial infections. Simulation proceeds forward for up to `n_days` days or until the outbreak is extinct.

Seed schools can be selected by ID or by sorting on vaccination coverage to explore worst-case scenarios:

```r
# Seed in the two schools with the lowest vaccination coverage
seed_schools <- schools |>
  dplyr::arrange(vaccination_coverage) |>
  dplyr::slice(1:2) |>
  dplyr::pull(school_id)
```

### Mid-Outbreak Forecasting

`run_midoutbreak_simulation()` initializes the simulation from a partially observed in-progress outbreak, then simulates the remaining trajectory. This is useful for real-time policy evaluation during an active event.

The observed state is described by a list:

```r
observed_state <- list(
  exposed_schools     = c("Springfield Elementary", "Westside Middle"),
  total_cases         = 32,             # cumulative confirmed cases
  quarantine_contacts = 85,             # contacts currently in quarantine
  fraction_active     = 0.15,           # fraction of cases still infectious
  hh_attack_rate      = 0.30,           # secondary HH attack rate applied to date
  additional_vaccinations = 0,          # emergency vaccinations to apply at t=0
  hh_vax_fraction     = 0.0            # fraction of HH members of new vaccinees to vaccinate
)

result <- run_midoutbreak_simulation(
  schools            = schools,
  network            = network,
  params             = params,
  observed_state     = observed_state,
  n_days             = 120,
  seed               = 42,
  hh_pop             = hh_pop,
  household_assignment = household_assignment
)
```

The mid-outbreak batch script `run_batch_midoutbreak.R` also generates a `task_XX_county.rds` file with county-level infection totals (students + household members combined).

---

## Outputs

Each call to `run_network_simulation()` returns a list with:

| Element | Description |
|---|---|
| `total_infected` | Total number of students infected across all schools |
| `total_breakthrough` | Students infected despite vaccination |
| `total_hh_members_infected` | Household members infected by students |
| `actual_days` | Number of days until the outbreak went extinct (≤ `n_days`) |
| `school_summary` | Data frame: one row per school with `total_infected`, `attack_rate`, `breakthrough_infections` |
| `network_daily_counts` | Data frame: one row per day with columns S, E, P, Ra, R, V, Iso for the full network |
| `hh_daily_counts` | Data frame: daily household-level counts (when `hh_pop` is provided) |

Across multiple runs (via `run_fast_simulations()` or `run_parallel_simulations()`), the `summary_stats` element provides per-simulation summary rows that can be used to compute medians, means, and 95% intervals for outbreak size, duration, schools affected, and household attack rates:

```r
summary_stats <- results$summary_stats

# Median and 95% interval for total student infections
cat(sprintf("Students infected — Median: %.0f  [95%% CI: %.0f – %.0f]\n",
  median(summary_stats$total_infected),
  quantile(summary_stats$total_infected, 0.025),
  quantile(summary_stats$total_infected, 0.975)))

# Median outbreak duration
cat(sprintf("Duration — Median: %.0f days  [95%% CI: %.0f – %.0f]\n",
  median(summary_stats$actual_days),
  quantile(summary_stats$actual_days, 0.025),
  quantile(summary_stats$actual_days, 0.975)))
```

Visualization functions in `plotting_utils.R` include:

- `plot_network_outbreak_dynamics(results)` — epidemic curve across all simulations
- `plot_outbreak_dynamics_report(results)` — multi-panel summary report
- School-level attack rate maps and school-type breakdowns

---

## Running Tests

```r
library(testthat)
test_dir("tests/testthat/")
```

---

## Citation

If you use this model in published research, please cite:

> Pandey A, *et al*. Agent-based simulation of measles transmission in school-household networks. *[Journal]*, *[Year]*. (in preparation)

Please also cite the R packages used for computation:

- Eddelbuettel D, François R (2011). Rcpp: Seamless R and C++ Integration. *Journal of Statistical Software*, 40(8), 1–18.
- Csárdi G, Nepusz T (2006). The igraph software package for complex network research. *InterJournal, Complex Systems*, 1695.
