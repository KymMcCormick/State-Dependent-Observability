# ============================================================
# State-dependent erosion simulation
# Abstract progressive subunit system
# Kym M. McCormick PhD
# Adelaide University
# June 2026
# ============================================================

# ---------------------------
# 0. Setup
# ---------------------------

set.seed(20260426)

required_packages <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "purrr",
  "patchwork",
  "rlang"
)

missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  stop(
    "The following packages are required but not installed: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(patchwork)
library(rlang)

# Create output folder if it does not already exist
if (!dir.exists("outputs")) {
  dir.create("outputs", recursive = TRUE)
}

# Inverse-logit helper
inv_logit <- function(x) {
  1 / (1 + exp(-x))
}

# Publication-style theme used throughout
theme_erosion <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      strip.background = element_rect(fill = "grey95", colour = "grey60"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
}

# ---------------------------
# User-controlled display choices
# ---------------------------

# Stages to display in selected figures.
# These are named so the labels are not hard-coded as Stage 3, Stage 6, etc.
stages_to_show <- c(
  early = 3,
  middle = 6,
  late = 10
)

# Latent-state bounds for the primary simulation
primary_state_lower <- 0
primary_state_cap <- 1

# Thresholds used for extent summaries.
# These are positions within the latent-state range, not raw values.
# With bounds 0-1, these become 0.25, 0.50, 0.75.
# With bounds 0-10, they would become 2.5, 5.0, 7.5.
state_threshold_positions <- c(
  low = 0.25,
  medium = 0.50,
  high = 0.75
)

# Thresholds used inside the erosion process.
# These are also positions within the latent-state range.
erosion_threshold_positions <- c(
  moderate = 0.25,
  high = 0.55
)

# ---------------------------
# Small helpers
# ---------------------------

check_named_unit_interval <- function(x, object_name) {
  if (is.null(names(x)) || any(names(x) == "")) {
    stop(
      object_name,
      " must be a named numeric vector.",
      call. = FALSE
    )
  }
  
  if (any(is.na(x)) || any(x < 0 | x > 1)) {
    stop(
      object_name,
      " must contain values between 0 and 1.",
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}

scale_positions_to_bounds <- function(positions, lower, upper) {
  lower + positions * (upper - lower)
}

make_extent_summary_exprs <- function(state_thresholds) {
  extent_summary_exprs <- list()
  
  for (threshold_name in names(state_thresholds)) {
    threshold_value <- state_thresholds[[threshold_name]]
    
    extent_summary_exprs[[paste0("extent_latent_", threshold_name)]] <-
      rlang::expr(mean(state_latent >= !!threshold_value))
    
    extent_summary_exprs[[paste0("extent_observed_", threshold_name)]] <-
      rlang::expr(mean(state_observed >= !!threshold_value, na.rm = TRUE))
  }
  
  extent_summary_exprs
}

# ============================================================
# 1. Primary simulation: constrained model for main figures
# ============================================================

simulate_progressive_system_main <- function(
    n_units_per_stage = 10000,
    n_subunits = 40,
    stages = 1:10,
    
    # Threshold locations as proportions of the latent-state range
    state_threshold_positions = c(
      low = 0.25,
      medium = 0.50,
      high = 0.75
    ),
    
    erosion_threshold_positions = c(
      moderate = 0.25,
      high = 0.55
    ),
    
    # Latent state process
    baseline_state = 0.10,
    stage_slope = 0.08,
    sd_unit = 0.15,
    sd_subunit = 0.10,
    sd_error = 0.08,
    state_lower = 0,
    state_cap = 1.00,
    
    # Subunit susceptibility
    susceptibility_low = 0.25,
    susceptibility_mid = 0.60,
    susceptibility_high = 1.00,
    susceptibility = NULL,
    
    # Observation regime
    observation_mode = c("complete", "background", "state_dependent"),
    
    # Background retention, independent of subunit state
    retention_intercept = 2.60,
    retention_stage_slope = -0.10,
    
    # State-dependent erosion
    penalty_moderate = 3.00,
    penalty_high = 8.00,
    erosion_intercept = 0.10,
    erosion_linear = 0.08,
    erosion_quadratic = 0.008,
    
    # Retention floor to avoid deterministic collapse
    retention_floor = 0.03
) {
  observation_mode <- match.arg(observation_mode)
  
  check_named_unit_interval(
    state_threshold_positions,
    "state_threshold_positions"
  )
  
  check_named_unit_interval(
    erosion_threshold_positions,
    "erosion_threshold_positions"
  )
  
  if (!all(c("moderate", "high") %in% names(erosion_threshold_positions))) {
    stop(
      "erosion_threshold_positions must include named values: moderate and high.",
      call. = FALSE
    )
  }
  
  if (state_cap <= state_lower) {
    stop(
      "state_cap must be greater than state_lower.",
      call. = FALSE
    )
  }
  
  state_thresholds <- scale_positions_to_bounds(
    positions = state_threshold_positions,
    lower = state_lower,
    upper = state_cap
  )
  
  erosion_thresholds <- scale_positions_to_bounds(
    positions = erosion_threshold_positions,
    lower = state_lower,
    upper = state_cap
  )
  
  threshold_moderate <- erosion_thresholds[["moderate"]]
  threshold_high <- erosion_thresholds[["high"]]
  
  # ---------------------------
  # Subunit susceptibility pattern
  # ---------------------------
  
  if (is.null(susceptibility)) {
    susceptibility <- rep(susceptibility_low, n_subunits)
    
    high_idx <- seq(1, n_subunits, length.out = 8) |>
      round() |>
      unique()
    
    mid_idx <- seq(
      ceiling(n_subunits * 0.60),
      ceiling(n_subunits * 0.75)
    )
    
    susceptibility[high_idx] <- susceptibility_high
    susceptibility[mid_idx] <- susceptibility_mid
  }
  
  stopifnot(length(susceptibility) == n_subunits)
  
  # ---------------------------
  # Unit-level table
  # ---------------------------
  
  units <- expand_grid(
    progression_stage = stages,
    unit_rep = 1:n_units_per_stage
  ) %>%
    mutate(
      unit_id = row_number(),
      unit_effect = rnorm(n(), mean = 0, sd = sd_unit)
    ) %>%
    select(unit_id, progression_stage, unit_effect)
  
  # ---------------------------
  # Unit x subunit table
  # ---------------------------
  
  subunit_data <- units %>%
    uncount(n_subunits, .id = "subunit_id") %>%
    mutate(
      observation_mode = observation_mode,
      susceptibility = susceptibility[subunit_id],
      
      subunit_effect = rnorm(n(), mean = 0, sd = sd_subunit),
      error = rnorm(n(), mean = 0, sd = sd_error),
      
      state_latent_raw =
        baseline_state +
        stage_slope * progression_stage * susceptibility +
        unit_effect +
        subunit_effect +
        error,
      
      state_latent = pmin(
        pmax(state_lower, state_latent_raw),
        state_cap
      ),
      
      # Background retention independent of subunit state
      retention_background =
        inv_logit(
          retention_intercept +
            retention_stage_slope * progression_stage
        ),
      
      # State-dependent erosion
      excess_moderate = pmax(0, state_latent - threshold_moderate),
      excess_high = pmax(0, state_latent - threshold_high),
      
      erosion_strength =
        erosion_intercept +
        erosion_linear * progression_stage +
        erosion_quadratic * progression_stage^2,
      
      state_penalty =
        -erosion_strength *
        (
          penalty_moderate * excess_moderate +
            penalty_high * excess_high
        ),
      
      retention_state_dependent =
        inv_logit(qlogis(retention_background) + state_penalty),
      
      retention_state_dependent =
        pmax(retention_state_dependent, retention_floor),
      
      # Observation process
      retention_probability = case_when(
        observation_mode == "complete" ~ 1,
        observation_mode == "background" ~ retention_background,
        observation_mode == "state_dependent" ~ retention_state_dependent
      ),
      
      observed_indicator = rbinom(
        n(),
        size = 1,
        prob = retention_probability
      ),
      
      state_observed = if_else(
        observed_indicator == 1,
        state_latent,
        NA_real_
      ),
      
      # Decomposition of loss
      loss_background = 1 - retention_background,
      loss_state_excess = retention_background - retention_state_dependent,
      loss_total_state_dependent = 1 - retention_state_dependent,
      
      proportion_state_dependent_loss =
        loss_state_excess / pmax(loss_total_state_dependent, 1e-12)
    )
  
  # ---------------------------
  # Unit-level summaries
  # ---------------------------
  
  extent_summary_exprs <- make_extent_summary_exprs(state_thresholds)
  
  unit_summary <- subunit_data %>%
    group_by(unit_id, progression_stage, observation_mode) %>%
    summarise(
      proportion_observed = mean(observed_indicator),
      
      mean_state_latent = mean(state_latent),
      mean_state_observed = mean(state_observed, na.rm = TRUE),
      
      !!!extent_summary_exprs,
      
      .groups = "drop"
    )
  
  list(
    subunit_data = subunit_data,
    unit_summary = unit_summary,
    state_thresholds = state_thresholds,
    erosion_thresholds = erosion_thresholds
  )
}

# ---------------------------
# Run primary observation regimes
# ---------------------------

sim_complete <- simulate_progressive_system_main(
  observation_mode = "complete",
  n_units_per_stage = 10000,
  state_lower = primary_state_lower,
  state_cap = primary_state_cap,
  state_threshold_positions = state_threshold_positions,
  erosion_threshold_positions = erosion_threshold_positions
)

sim_background <- simulate_progressive_system_main(
  observation_mode = "background",
  n_units_per_stage = 10000,
  state_lower = primary_state_lower,
  state_cap = primary_state_cap,
  state_threshold_positions = state_threshold_positions,
  erosion_threshold_positions = erosion_threshold_positions
)

sim_state_dependent <- simulate_progressive_system_main(
  observation_mode = "state_dependent",
  n_units_per_stage = 10000,
  state_lower = primary_state_lower,
  state_cap = primary_state_cap,
  state_threshold_positions = state_threshold_positions,
  erosion_threshold_positions = erosion_threshold_positions
)

primary_state_thresholds <- sim_state_dependent$state_thresholds

# ---------------------------
# Combine primary simulation outputs
# ---------------------------

unit_summary_all <- bind_rows(
  sim_complete$unit_summary,
  sim_background$unit_summary,
  sim_state_dependent$unit_summary
) %>%
  mutate(
    observation_mode = factor(
      observation_mode,
      levels = c("complete", "background", "state_dependent"),
      labels = c(
        "Complete observation",
        "Background loss",
        "State-dependent erosion"
      )
    )
  )

subunit_data_all <- bind_rows(
  sim_complete$subunit_data,
  sim_background$subunit_data,
  sim_state_dependent$subunit_data
) %>%
  mutate(
    observation_mode = factor(
      observation_mode,
      levels = c("complete", "background", "state_dependent"),
      labels = c(
        "Complete observation",
        "Background loss",
        "State-dependent erosion"
      )
    )
  )

extent_cols <- grep("^extent_", names(unit_summary_all), value = TRUE)

stage_summary <- unit_summary_all %>%
  group_by(observation_mode, progression_stage) %>%
  summarise(
    proportion_observed = mean(proportion_observed),
    
    mean_state_latent = mean(mean_state_latent),
    mean_state_observed = mean(mean_state_observed, na.rm = TRUE),
    
    across(
      all_of(extent_cols),
      ~ mean(.x, na.rm = TRUE)
    ),
    
    .groups = "drop"
  )

# Loss decomposition for diagnostics
loss_decomposition <- sim_state_dependent$subunit_data %>%
  group_by(progression_stage) %>%
  summarise(
    background_loss = mean(loss_background),
    state_dependent_excess_loss = mean(loss_state_excess),
    total_loss = mean(loss_total_state_dependent),
    proportion_loss_state_dependent =
      mean(proportion_state_dependent_loss),
    .groups = "drop"
  )

# Print quick diagnostics to console
message("\nPrimary simulation diagnostics: state-dependent erosion regime")

stage_summary %>%
  filter(observation_mode == "State-dependent erosion") %>%
  select(
    progression_stage,
    proportion_observed,
    mean_state_latent,
    mean_state_observed,
    all_of(extent_cols)
  ) %>%
  print(n = Inf)

message("\nThresholds used for primary extent summaries:")
print(primary_state_thresholds)

message("\nLoss decomposition diagnostics:")
print(loss_decomposition, n = Inf)

# ============================================================
# 2. Main manuscript figures
# ============================================================

# ---------------------------
# Figure 1: Erosion mechanism
# ---------------------------

state_bin_width <- 0.10

fig1a_df <- sim_state_dependent$subunit_data %>%
  filter(
    progression_stage %in% unname(stages_to_show)
  ) %>%
  mutate(
    progression_stage = factor(
      progression_stage,
      levels = unname(stages_to_show),
      labels = paste0(
        tools::toTitleCase(names(stages_to_show)),
        ""
      )
    ),
    state_bin = cut(
      state_latent,
      breaks = seq(
        primary_state_lower,
        primary_state_cap,
        by = state_bin_width
      ),
      include.lowest = TRUE
    )
  ) %>%
  group_by(progression_stage, state_bin) %>%
  summarise(
    bin_mid = mean(state_latent, na.rm = TRUE),
    retention = mean(observed_indicator),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(n >= 50)

fig1a <- ggplot(
  fig1a_df,
  aes(
    x = bin_mid,
    y = retention,
    group = progression_stage,
    linetype = progression_stage
  )
) +
  geom_line(linewidth = 0.75) +
  scale_x_continuous(
    breaks = seq(primary_state_lower, primary_state_cap, by = 0.2),
    limits = c(primary_state_lower, primary_state_cap)
  ) +
  scale_linetype_discrete(
    name = "Progression"
  ) +
  theme_erosion() +
  theme(
    legend.position = "none"
  ) +
  labs(
    x = "Accumulated component state",
    y = "Probability observed"
  )

fig1b_df <- sim_state_dependent$subunit_data %>%
  filter(
    progression_stage %in% unname(stages_to_show)
  ) %>%
  mutate(
    progression_stage = factor(
      progression_stage,
      levels = unname(stages_to_show),
      labels = paste0(
        tools::toTitleCase(names(stages_to_show)),
        ""
      )
    ),
    state_bin = cut(
      state_latent,
      breaks = seq(
        primary_state_lower,
        primary_state_cap,
        by = state_bin_width
      ),
      include.lowest = TRUE
    )
  ) %>%
  group_by(progression_stage, state_bin) %>%
  summarise(
    bin_mid = mean(state_latent, na.rm = TRUE),
    n_complete = n(),
    n_observed = sum(observed_indicator),
    .groups = "drop"
  ) %>%
  group_by(progression_stage) %>%
  mutate(
    remaining_proportion = (n_observed) / sum(n_complete)
  ) %>%
  ungroup()

fig1b <- ggplot(
  fig1b_df,
  aes(
    x = bin_mid,
    y = remaining_proportion,
    group = progression_stage,
    linetype = progression_stage
  )
) +
  geom_line(linewidth = 0.75) +
  scale_x_continuous(
    breaks = seq(primary_state_lower, primary_state_cap, by = 0.2),
    limits = c(primary_state_lower, primary_state_cap)
  ) +
  scale_linetype_discrete(
    name = "Progression"
  ) +
  theme_erosion() +
  theme(
    legend.position = "right"
  ) +
  labs(
    x = "Accumulated component state",
    y = "Observable proportion of components"
  )

fig1 <- fig1a + fig1b +
  plot_layout(
    guides = "collect",
    widths = c(1, 1)
  ) +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8.5),
    legend.key.width = unit(1.2, "cm"),
    legend.spacing.x = unit(0.3, "cm"),
    plot.tag = element_text(size = 11, face = "plain"),
    plot.margin = margin(4, 4, 4, 4)
  ) &
  guides(
    linetype = guide_legend(
      nrow = 1,
      byrow = TRUE,
      title.position = "left"
    )
  )

print(fig1)

ggsave(
  filename = "outputs/figure_1_distribution_erosion.png",
  plot = fig1,
  width = 7,
  height = 3.8,
  dpi = 300
)

ggsave(
  filename = "outputs/figure_1_distribution_erosion.pdf",
  plot = fig1,
  width = 7,
  height = 3.8
)

# ---------------------------
# Figure 2: Mean state trajectories
# ---------------------------

fig2_df <- stage_summary %>%
  filter(observation_mode == "State-dependent erosion") %>%
  select(
    progression_stage,
    mean_state_latent,
    mean_state_observed
  ) %>%
  pivot_longer(
    -progression_stage,
    names_to = "system",
    values_to = "mean_state"
  ) %>%
  mutate(
    system = recode(
      system,
      mean_state_latent = "Complete / latent",
      mean_state_observed = "Observed"
    )
  )

fig2 <- ggplot(
  fig2_df,
  aes(
    x = progression_stage,
    y = mean_state,
    linetype = system
  )
) +
  geom_line(linewidth = 1) +
  theme_erosion() +
  scale_x_continuous(
    breaks = sort(unique(stage_summary$progression_stage))
  ) +
  labs(
    x = "Progression stage",
    y = "Mean state",
    linetype = NULL
  )

print(fig2)

ggsave(
  filename = "outputs/figure_2_mean_state_trajectory.png",
  plot = fig2,
  width = 5.8,
  height = 4,
  dpi = 300
)

ggsave(
  filename = "outputs/figure_2_mean_state_trajectory.pdf",
  plot = fig2,
  width = 5.8,
  height = 4
)

# ---------------------------
# Figure 3: Threshold summaries
# ---------------------------

fig3_df <- stage_summary %>%
  filter(observation_mode == "State-dependent erosion") %>%
  select(
    progression_stage,
    all_of(extent_cols)
  ) %>%
  pivot_longer(
    -progression_stage,
    names_to = "measure",
    values_to = "extent"
  ) %>%
  mutate(
    system = if_else(
      grepl("^extent_latent_", measure),
      "Complete / latent",
      "Observed"
    ),
    threshold = sub(
      "^extent_(latent|observed)_",
      "",
      measure
    ),
    threshold = factor(
      threshold,
      levels = names(state_threshold_positions),
      labels = paste0(
        tools::toTitleCase(names(state_threshold_positions)),
        " threshold"
      )
    )
  )

fig3 <- ggplot(
  fig3_df,
  aes(
    x = progression_stage,
    y = extent,
    linetype = system
  )
) +
  geom_line(linewidth = 1) +
  facet_wrap(~ threshold) +
  theme_erosion() +
  scale_x_continuous(
    breaks = sort(unique(stage_summary$progression_stage))
  ) +
  labs(
    x = "Progression stage",
    y = "Extent",
    linetype = NULL
  )

print(fig3)

ggsave(
  filename = "outputs/figure_3_threshold_summaries.png",
  plot = fig3,
  width = 7,
  height = 3.8,
  dpi = 300
)

ggsave(
  filename = "outputs/figure_3_threshold_summaries.pdf",
  plot = fig3,
  width = 7,
  height = 3.8
)

# ---------------------------
# Figure 4: Background loss versus state-dependent erosion
# ---------------------------

high_threshold_name <- names(state_threshold_positions)[
  which.max(state_threshold_positions)
]

high_extent_observed_col <- paste0(
  "extent_observed_",
  high_threshold_name
)

fig4_df <- stage_summary %>%
  filter(observation_mode != "Complete observation") %>%
  select(
    observation_mode,
    progression_stage,
    mean_state_observed,
    proportion_observed,
    all_of(high_extent_observed_col)
  ) %>%
  rename(
    observed_high_extent = !!rlang::sym(high_extent_observed_col)
  ) %>%
  pivot_longer(
    cols = c(
      mean_state_observed,
      observed_high_extent,
      proportion_observed
    ),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(
      metric,
      mean_state_observed = "Observed mean state",
      observed_high_extent = "High-threshold extent",
      proportion_observed = "Proportion observed"
    ),
    metric = factor(
      metric,
      levels = c(
        "High-threshold extent",
        "Observed mean state",
        "Proportion observed"
      )
    )
  )

fig4 <- ggplot(
  fig4_df,
  aes(
    x = progression_stage,
    y = value,
    linetype = observation_mode
  )
) +
  geom_line(linewidth = 1) +
  facet_wrap(~ metric, scales = "free_y") +
  theme_erosion() +
  scale_x_continuous(
    breaks = sort(unique(stage_summary$progression_stage))
  ) +
  labs(
    x = "Progression stage",
    y = NULL,
    linetype = NULL
  )

print(fig4)

ggsave(
  filename = "outputs/figure_4_background_vs_state_dependent.png",
  plot = fig4,
  width = 7,
  height = 3.8,
  dpi = 300
)

ggsave(
  filename = "outputs/figure_4_background_vs_state_dependent.pdf",
  plot = fig4,
  width = 7,
  height = 3.8
)

# ============================================================
# 3. Relaxed simulation for supplementary sensitivity analysis
# ============================================================

simulate_progressive_system_sensitivity <- function(
    n_units_per_stage = 5000,
    n_subunits = 28,
    stages = 1:11,
    
    # Threshold locations as proportions of latent-state range
    state_threshold_positions = c(
      low = 0.25,
      medium = 0.50,
      high = 0.75
    ),
    
    erosion_threshold_positions = c(
      moderate = 0.25,
      high = 0.50
    ),
    
    # Latent state process
    baseline_state = 1.5,
    stage_slope = 0.4,
    sd_unit = 1.5,
    sd_subunit = 1.0,
    sd_error = 1.0,
    state_lower = 0,
    state_cap = 12,
    
    # Heterogeneity
    heterogeneity = c("structured", "homogeneous"),
    susceptibility_low = 0.20,
    susceptibility_mid = 0.55,
    susceptibility_high = 1.00,
    
    # Observation regime
    observation_mode = c("complete", "background", "state_dependent"),
    
    # Loss function and coupling strength
    loss_function = c("logistic", "linear", "step"),
    coupling_strength = c("none", "weak", "moderate", "strong"),
    
    # Background retention
    retention_intercept = 2.84,
    retention_stage_slope = -0.12,
    
    # Retention floor
    retention_floor = 0.06
) {
  heterogeneity <- match.arg(heterogeneity)
  observation_mode <- match.arg(observation_mode)
  loss_function <- match.arg(loss_function)
  coupling_strength <- match.arg(coupling_strength)
  
  check_named_unit_interval(
    state_threshold_positions,
    "state_threshold_positions"
  )
  
  check_named_unit_interval(
    erosion_threshold_positions,
    "erosion_threshold_positions"
  )
  
  if (!all(c("moderate", "high") %in% names(erosion_threshold_positions))) {
    stop(
      "erosion_threshold_positions must include named values: moderate and high.",
      call. = FALSE
    )
  }
  
  if (state_cap <= state_lower) {
    stop(
      "state_cap must be greater than state_lower.",
      call. = FALSE
    )
  }
  
  state_thresholds <- scale_positions_to_bounds(
    positions = state_threshold_positions,
    lower = state_lower,
    upper = state_cap
  )
  
  erosion_thresholds <- scale_positions_to_bounds(
    positions = erosion_threshold_positions,
    lower = state_lower,
    upper = state_cap
  )
  
  threshold_moderate <- erosion_thresholds[["moderate"]]
  threshold_high <- erosion_thresholds[["high"]]
  
  # Coupling parameters
  coupling_params <- list(
    none = list(
      penalty_moderate = 0.00,
      penalty_high = 0.00,
      erosion_intercept = 0.00,
      erosion_linear = 0.00,
      erosion_quadratic = 0.00
    ),
    weak = list(
      penalty_moderate = 0.25,
      penalty_high = 0.75,
      erosion_intercept = 0.02,
      erosion_linear = 0.02,
      erosion_quadratic = 0.001
    ),
    moderate = list(
      penalty_moderate = 0.50,
      penalty_high = 2.00,
      erosion_intercept = 0.05,
      erosion_linear = 0.04,
      erosion_quadratic = 0.004
    ),
    strong = list(
      penalty_moderate = 0.75,
      penalty_high = 3.00,
      erosion_intercept = 0.08,
      erosion_linear = 0.06,
      erosion_quadratic = 0.006
    )
  )
  
  pars <- coupling_params[[coupling_strength]]
  
  # Subunit susceptibility
  if (heterogeneity == "homogeneous") {
    susceptibility <- rep(1, n_subunits)
  }
  
  if (heterogeneity == "structured") {
    susceptibility <- rep(susceptibility_low, n_subunits)
    
    high_idx <- seq(1, n_subunits, length.out = 8) |>
      round() |>
      unique()
    
    mid_idx <- seq(
      ceiling(n_subunits * 0.60),
      ceiling(n_subunits * 0.75)
    )
    
    susceptibility[high_idx] <- susceptibility_high
    susceptibility[mid_idx] <- susceptibility_mid
  }
  
  units <- expand_grid(
    progression_stage = stages,
    unit_rep = 1:n_units_per_stage
  ) %>%
    mutate(
      unit_id = row_number(),
      unit_effect = rnorm(n(), 0, sd_unit)
    ) %>%
    select(unit_id, progression_stage, unit_effect)
  
  subunit_data <- units %>%
    uncount(n_subunits, .id = "subunit_id") %>%
    mutate(
      observation_mode = observation_mode,
      loss_function = loss_function,
      coupling_strength = coupling_strength,
      heterogeneity = heterogeneity,
      
      susceptibility = susceptibility[subunit_id],
      subunit_effect = rnorm(n(), 0, sd_subunit),
      error = rnorm(n(), 0, sd_error),
      
      state_latent_raw =
        baseline_state +
        stage_slope * progression_stage * susceptibility +
        unit_effect +
        subunit_effect +
        error,
      
      state_latent = pmin(
        pmax(state_lower, state_latent_raw),
        state_cap
      ),
      
      retention_background =
        inv_logit(
          retention_intercept +
            retention_stage_slope * progression_stage
        ),
      
      excess_moderate = pmax(0, state_latent - threshold_moderate),
      excess_high = pmax(0, state_latent - threshold_high),
      
      erosion_strength =
        pars$erosion_intercept +
        pars$erosion_linear * progression_stage +
        pars$erosion_quadratic * progression_stage^2
    )
  
  # Alternative state-dependent loss functions
  if (loss_function == "logistic") {
    subunit_data <- subunit_data %>%
      mutate(
        state_penalty =
          -erosion_strength *
          (
            pars$penalty_moderate * excess_moderate +
              pars$penalty_high * excess_high
          ),
        
        retention_state_dependent =
          inv_logit(qlogis(retention_background) + state_penalty)
      )
  }
  
  if (loss_function == "linear") {
    subunit_data <- subunit_data %>%
      mutate(
        linear_loss =
          erosion_strength *
          (
            0.04 * pars$penalty_moderate * excess_moderate +
              0.04 * pars$penalty_high * excess_high
          ),
        
        retention_state_dependent =
          retention_background - linear_loss
      )
  }
  
  if (loss_function == "step") {
    subunit_data <- subunit_data %>%
      mutate(
        step_penalty = case_when(
          state_latent >= threshold_high ~ 0.45 * pars$penalty_high / 3,
          state_latent >= threshold_moderate ~ 0.20 * pars$penalty_moderate,
          TRUE ~ 0
        ),
        
        retention_state_dependent =
          retention_background - erosion_strength * step_penalty
      )
  }
  
  subunit_data <- subunit_data %>%
    mutate(
      retention_state_dependent =
        pmin(pmax(retention_state_dependent, retention_floor), 1),
      
      retention_probability = case_when(
        observation_mode == "complete" ~ 1,
        observation_mode == "background" ~ retention_background,
        observation_mode == "state_dependent" ~ retention_state_dependent
      ),
      
      observed_indicator = rbinom(
        n(),
        size = 1,
        prob = retention_probability
      ),
      
      state_observed = if_else(
        observed_indicator == 1,
        state_latent,
        NA_real_
      )
    )
  
  extent_summary_exprs <- make_extent_summary_exprs(state_thresholds)
  
  unit_summary <- subunit_data %>%
    group_by(
      unit_id,
      progression_stage,
      observation_mode,
      loss_function,
      coupling_strength,
      heterogeneity
    ) %>%
    summarise(
      proportion_observed = mean(observed_indicator),
      
      mean_state_latent = mean(state_latent),
      mean_state_observed = mean(state_observed, na.rm = TRUE),
      
      !!!extent_summary_exprs,
      
      .groups = "drop"
    )
  
  list(
    subunit_data = subunit_data,
    unit_summary = unit_summary,
    state_thresholds = state_thresholds,
    erosion_thresholds = erosion_thresholds
  )
}

# ---------------------------
# Run sensitivity simulations
# ---------------------------

sensitivity_grid <- expand_grid(
  loss_function = c("logistic", "linear", "step"),
  coupling_strength = c("weak", "moderate", "strong"),
  heterogeneity = c("structured", "homogeneous")
)

sensitivity_results <- pmap(
  sensitivity_grid,
  function(loss_function, coupling_strength, heterogeneity) {
    simulate_progressive_system_sensitivity(
      observation_mode = "state_dependent",
      loss_function = loss_function,
      coupling_strength = coupling_strength,
      heterogeneity = heterogeneity,
      n_units_per_stage = 5000,
      state_threshold_positions = state_threshold_positions
    )$unit_summary
  }
)

sensitivity_summary <- bind_rows(sensitivity_results)

sensitivity_extent_cols <- grep(
  "^extent_",
  names(sensitivity_summary),
  value = TRUE
)

sensitivity_high_extent_latent_col <- paste0(
  "extent_latent_",
  high_threshold_name
)

sensitivity_high_extent_observed_col <- paste0(
  "extent_observed_",
  high_threshold_name
)

stage_summary_sensitivity <- sensitivity_summary %>%
  group_by(
    loss_function,
    coupling_strength,
    heterogeneity,
    progression_stage
  ) %>%
  summarise(
    proportion_observed = mean(proportion_observed),
    
    mean_state_latent = mean(mean_state_latent),
    mean_state_observed = mean(mean_state_observed, na.rm = TRUE),
    
    across(
      all_of(sensitivity_extent_cols),
      ~ mean(.x, na.rm = TRUE)
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    mean_bias = mean_state_latent - mean_state_observed,
    
    extent_high_bias =
      .data[[sensitivity_high_extent_latent_col]] -
      .data[[sensitivity_high_extent_observed_col]],
    
    loss_function = factor(
      loss_function,
      levels = c("linear", "logistic", "step")
    ),
    coupling_strength = factor(
      coupling_strength,
      levels = c("weak", "moderate", "strong")
    ),
    heterogeneity = factor(
      heterogeneity,
      levels = c("homogeneous", "structured")
    )
  )

# ---------------------------
# Supplementary Figure S1: relaxed assumptions
# ---------------------------

fig_s1 <- ggplot(
  stage_summary_sensitivity,
  aes(
    x = progression_stage,
    y = mean_bias,
    linetype = coupling_strength
  )
) +
  geom_line(linewidth = 0.9) +
  facet_grid(heterogeneity ~ loss_function) +
  theme_erosion() +
  scale_x_continuous(
    breaks = sort(unique(stage_summary_sensitivity$progression_stage))
  ) +
  labs(
    x = "Progression stage",
    y = "Bias in mean state",
    linetype = "Coupling",
    title = "Erosion-induced bias persists across relaxed simulation assumptions"
  )

print(fig_s1)

ggsave(
  filename = "outputs/figure_s1_relaxed_assumptions.png",
  plot = fig_s1,
  width = 8.2,
  height = 5.2,
  dpi = 300
)

ggsave(
  filename = "outputs/figure_s1_relaxed_assumptions.pdf",
  plot = fig_s1,
  width = 8.2,
  height = 5.2
)

# ---------------------------
# Save key data objects for reproducibility
# ---------------------------

saveRDS(
  list(
    sim_complete = sim_complete,
    sim_background = sim_background,
    sim_state_dependent = sim_state_dependent,
    primary_state_thresholds = primary_state_thresholds,
    stage_summary = stage_summary,
    loss_decomposition = loss_decomposition,
    stage_summary_sensitivity = stage_summary_sensitivity
  ),
  file = "outputs/simulation_outputs.rds"
)

message("\nDone. Figures and simulation outputs saved in the outputs/ folder.")
