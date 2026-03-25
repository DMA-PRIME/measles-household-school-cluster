# measles-household-school-cluster

An agent-based model (ABM) that simulates measles transmission across interconnected schools and households in South Carolina. The model captures within-school spread, between-school spread via a network of schools, and secondary transmission to household members.

---

## Overview

This model supports scenario analyses such as:

- Projecting outbreak size and duration from a single index case
- Real-time forecasting from a partially observed, in-progress outbreak
- Evaluating interventions (isolation, contact quarantine, vaccination coverage)
- Comparing transmission risk across school networks with different structures

---

## Model Structure

### Population

Students aged 5–18 are assigned to schools based on age and geographic proximity. Each school is divided into classes (default size 25). Students are further linked to synthetic households loaded from RTI population data, which include adults and children under 5 who can acquire secondary infections.

### Disease States

The model uses a modified SEIR framework:

| State | Meaning |
|-------|---------|
| **S** | Susceptible |
| **E** | Exposed (latent, ~10 days, Erlang-distributed) |
| **P** | Prodromal (infectious, ~4 days) |
| **Ra** | Rash (infectious, ~4 days) |
| **R** | Recovered |
| **V** | Vaccinated (97% fully immune; 3% leaky — can be infected with reduced infectiousness) |
| **Iso** | Isolated (non-infectious) |
| **Q\*** | Quarantine states (held until symptom-onset confirmation) |

### Transmission Pathways

**Within-school:** Each infectious student (P or Ra, not isolated) makes contacts at rates `c_within` (within class, default 6/day) and `c_between` (between classes, default 2/day). Contacts are sampled from all present students; transmission occurs only if the contact is susceptible or has leaky vaccine immunity.

**Between-school (network-based):** Schools are nodes in a weighted network. Edge weights are derived from geographic distance or travel time between schools, optionally adjusted for grade overlap. Infectious students generate cross-school contacts scaled by edge weight and base rate `c_between_school` (default 0.2/day). Contact pools are drawn from all non-isolated students at the target school.

**Household:** Infectious students transmit to non-student household members daily with probability `hh_transmission_prob` (default 0.17). Household members who develop infection are tracked but do not feed back into the school simulation.

### Interventions

- **Isolation:** Detected cases are isolated after a delay (2 days for the index case, 3 days for subsequent cases). Isolation lasts 14 days.
- **Contact quarantine:** Contacts of isolated individuals are quarantined for 21 days with configurable efficacy (default 90%).
- **Vaccination:** Coverage varies by school. Vaccine efficacy is modeled as leaky: 97% of vaccinated students receive complete immunity; the remaining 3% can be infected but transmit at 80% reduced infectiousness.

---

## School Network

Schools are connected by an undirected weighted graph built from one of several network types:

- **Distance-based** — edge weight is an inverse function of geographic distance
- **Travel-time-based** — edge weight is an inverse function of driving time
- **Small-world (Watts-Strogatz)** — captures clustered local connections with occasional long-range links
- **Scale-free (Barabási-Albert)** — hubs (large or central schools) drive most between-school transmission

A grade-weighted overlay can be applied so that schools with overlapping grade ranges (and therefore shared siblings or community ties) receive higher edge weights.

---

## Simulation Modes

### Clean outbreak

Starts from a single index case at one or more seed schools and simulates forward for a specified number of days (default 365).

### Mid-outbreak (real-time forecasting)

Initializes from a known in-progress state — total confirmed cases, schools already affected, quarantined contacts, and estimated fraction still infectious — and simulates the remaining outbreak trajectory. Useful for evaluating policy options during an active event.

---

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `latent_mean` | 10 days | Mean latent (exposed) period |
| `infectious_mean` | 8 days | Total mean infectious period |
| `c_within` | 6 | Within-class daily contacts per student |
| `c_between` | 2 | Between-class daily contacts per student |
| `p_within` / `p_between` | 0.17 | Per-contact transmission probability |
| `c_between_school` | 0.2 | Base between-school contact rate |
| `hh_transmission_prob` | 0.17 | Daily household transmission probability |
| `vaccine_efficacy` | 0.97 | Fraction of vaccinated students fully protected |
| `isolation_delay_index` | 2 days | Detection delay for the first identified case |
| `isolation_delay_secondary` | 3 days | Detection delay for subsequent cases |
| `quarantine_efficacy` | 0.90 | Fraction of quarantined contacts who comply |
| `quarantine_duration` | 21 days | Length of contact quarantine |

---

## Repository Structure

```
load_all.R                    # Sources all modules in dependency order
rcpp_transmission.R           # Compiled C++ (Rcpp) transmission functions
network_utils.R               # School network generation
grade_utils.R                 # Grade range parsing and school-type classification
population_utils.R            # School population creation, SEIR state updates
household_utils.R             # Household loading, student assignment, HH transmission
simulation_utils.R            # Core simulation engine (single and multiple runs)
simulation_parallel.R         # Parallel execution across CPU cores
analysis_utils.R              # Post-simulation summary statistics
plotting_utils.R              # Epidemic curve and outbreak visualization
midoutbreak_utils.R           # Mid-outbreak scenario initialization
between_school_transmission_fix.R  # Corrected between-school contact sampling
rcpp_within_school_fix.R      # Corrected within-school contact sampling
run_simulation_households.R   # End-to-end example script
run_batch_clean.R             # SLURM batch script for clean outbreak runs
run_batch_midoutbreak.R       # SLURM batch script for mid-outbreak runs
combine_results.R             # Aggregates output from multiple batch tasks
tests/                        # Unit tests (testthat framework)
```

---

## Quick Start

```r
# Load all modules
source("load_all.R")

# See run_simulation_households.R for a complete end-to-end example, including:
#   - Loading school and synthetic population data
#   - Building a travel-time school network
#   - Assigning students to households
#   - Running simulations in parallel
#   - Generating summary statistics and plots
```

### Running tests

```r
library(testthat)
test_dir("tests/testthat/")
```

### HPC / SLURM batch execution

```bash
# Clean outbreak from index case
Rscript run_batch_clean.R --task_id=1 --n_sims=100 --n_cores=12 \
  --n_days=365 --output_dir=./results

# Mid-outbreak forecast
Rscript run_batch_midoutbreak.R --task_id=1 --n_sims=100 --n_cores=12 \
  --n_days=120 --output_dir=./results_forecast
```

Results from parallel tasks are combined with `combine_results.R`.

---

## Outputs

Each simulation run returns:

- **Daily counts** by disease state for every school and the full network
- **School-level summaries:** total infected, attack rate, breakthrough infections
- **Household-level counts:** secondary infections among non-student household members
- **Epidemic curves** for visualization

Across multiple runs, the analysis utilities compute median, mean, and 95% intervals for outbreak size, duration, schools affected, and attack rates.
