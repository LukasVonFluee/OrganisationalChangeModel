###############################################################################
###############################################################################
###############################################################################

# Part of the project on social norms change in organizations
# by Faessler, von Flüe, Zehnder, Efferson

# Programmed by Lukas von Flüe

# simulation.R
# simulation script which requires functions.R and config.R to be in same wd

###############################################################################
###############################################################################
###############################################################################

# The interaction between agents is defined by the following coordination game:

# Coordination game payoff matrix with x_values:

#       SQ        Alt
##########################
# SQ    a+x_i     b+x_i
# Alt   a         d
##########################

# It is a strict coordination game because for each agent $i$, $a + x_{i} > a$, and $b + x_{i} < d$, where $b<d$.
# The x_i values are distributed over the open interval (0,(d-b)) according to a beta distribution with shape parameters alpha and beta.

# Risk dominance of SQ if:
# (x_i > (d − b)/2).


###############################################################################

source("functions.R")
source("config.R")

suppressWarnings(library(ggplot2)) # only required for testing purposes

###############################################################################
# Build pair table from config for the pairs of beta distributions within groups
###############################################################################


output_dir <- "simulation_results"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Build pair table
tp_cfg  <- normalise_types_pair(get0("types_pair", ifnotfound = NULL))
run_all <- isTRUE(get0("run_all_pairs", ifnotfound = FALSE))

if (run_all) {
  # There are 15 unique unordered pairs (including self-pairs)
  unique_pairs <- all_group_pairs(include_self = TRUE)   # list of length 15
  num_pairs    <- length(unique_pairs)                   # 15
  
  # Optional filtering via config$combo_ids
  cfg_ids <- get0("combo_ids", ifnotfound = NULL)
  if (is.null(cfg_ids)) {
    pair_ids <- seq_len(num_pairs)                       # use all 15 by default
  } else {
    pair_ids <- as.integer(cfg_ids)
    # keep only valid ids in [1, 15]
    valid <- !is.na(pair_ids) & pair_ids >= 1L & pair_ids <= num_pairs
    if (!all(valid)) {
      warning("Dropping invalid combo_ids: ",
              paste(pair_ids[!valid], collapse = ", "))
      pair_ids <- pair_ids[valid]
    }
    if (length(pair_ids) == 0L) {
      warning("No valid combo_ids provided; using all 15 pairs.")
      pair_ids <- seq_len(num_pairs)
    }
  }
  
  # Build the table directly from the precomputed list
  pair_tbl <- do.call(rbind, lapply(pair_ids, function(cid) {
    p <- unique_pairs[[cid]]
    data.frame(
      pair_id   = cid,
      pair_g1   = p[1],
      pair_g2   = p[2],
      pair_label = pair_label(p),
      stringsAsFactors = FALSE
    )
  }))
  
} else if (!is.null(tp_cfg)) {
  p <- resolve_group_pair(types_pair = tp_cfg, include_self = TRUE)
  pair_tbl <- data.frame(
    pair_id = NA_integer_, pair_g1 = p[1], pair_g2 = p[2],
    pair_label = pair_label(p), stringsAsFactors = FALSE
  )
  
} else {
  p <- resolve_group_pair(combo_id = 1L, include_self = TRUE)
  pair_tbl <- data.frame(
    pair_id = 1L, pair_g1 = p[1], pair_g2 = p[2],
    pair_label = pair_label(p), stringsAsFactors = FALSE
  )
}

# --- Quick sanity check: print the pairs that will be simulated ---
cat("\n=== Pairs of Beta distributions used in this simulation ===\n")
print(pair_tbl[, c("pair_id", "pair_g1", "pair_g2", "pair_label")])
cat("Total pairs:", nrow(pair_tbl), "\n\n")
if (run_all) {
  expected <- length(all_group_pairs(include_self = TRUE)) # 15
  cat("Expected 15 unique unordered pairs. Found:", nrow(pair_tbl), "\n\n")
}

###############################################################################
# Parameter grid
###############################################################################
base_grid <- expand.grid(
  t_max            = t_max,
  interventionTime = interventionTime,
  phi              = phi,
  N                = N,
  G                = G,
  groupSize        = groupSize,
  n_E_in           = n_E_in,
  n_E_out          = n_E_out,
  n_M              = n_M,
  probOut          = probOut,
  interventionTarget = interventionTarget,
  targetSuccess    = targetSuccess,
  targetScope      = targetScope,
  wE_inE           = wE_inE,
  wE_outE          = wE_outE,
  wM_M             = wM_M,
  choice_lambda    = choice_lambda,
  a                = a,
  b                = b,
  d                = d,
  g                = g,
  h                = h,
  stringsAsFactors = FALSE
)

# Cross with pairs so the pair is VISIBLE in the grid
param_combinations <- merge(
  base_grid, pair_tbl, all = TRUE
)

###############################################################################
# Optional: plot init distributions once per sim
plot_dir  <- file.path(output_dir, "plots")
###############################################################################

###############################################################################
###############################################################################
###############################################################################
# Run simulation
###############################################################################
###############################################################################
###############################################################################

results_list <- list()
sim_counter <- 0L

for (i in seq_len(nrow(param_combinations))) {
  params <- param_combinations[i, , drop = FALSE]
  pair <- c(params$pair_g1, params$pair_g2)
  
  sim_results <- vector("list", length = n_sim[1])
  names(sim_results) <- sprintf("sim_%03d", seq_len(n_sim[1]))
  
  for (sim in seq_len(n_sim[1])) {
    # initialise with explicit pair
    agents <- initializeAgents(
      N = params$N, G = params$G,
      a = params$a, b = params$b, d = params$d,
      types_pair = pair,
      shape_overrides = get0("shape_overrides", ifnotfound = NULL)
    )
    
    # optional plot of initial x_i
    if (isTRUE(plot_init)) {
      fname <- sprintf("init_density_%s_row%03d_sim%03d.png", params$pair_label, i, sim)
      plot_xi_densities(agents, b = params$b, d = params$d,
                        pair_used = pair, save = TRUE, outdir = plot_dir, fname = fname)
    }
    
    summary_results <- initialize_summary_table(params$t_max, params$G)
    
    for (t in seq_len(params$t_max)) {
      
      # Intervention
      if (t == params$interventionTime) {
        agents <- applyIntervention(
          agents, params$G, params$phi, params$d, params$b,
          params$interventionTarget,
          params$targetSuccess, params$targetScope
        )
      }
      
      ## 1. Belief updates based on last period's choices
      agents <- updateManagerBeliefs(
        agents, params$G, params$n_M, params$wM_M
      )
      agents <- updateEmployeeBeliefs(
        agents, params$G, params$n_E_in, params$n_E_out,
        params$wE_inE, params$wE_outE, params$probOut
      )
      
      ## 2. Choices this period, given beliefs + respond
      agents <- manager_coordination(
        agents, params$G, params$a, params$b, params$d, params$h, params$g
      )
      
      agents <- belief_to_choice(
        agents, params$a, params$b, params$d, params$h, params$g,
        choice_lambda = params$choice_lambda
      )
      
      ## 3. Employee coordination
      emp <- employee_coordination(
        agents, params$G,
        params$a, params$b, params$d, params$h, params$g,
        params$probOut
      )
      
      agents <- emp$agents
      
      ## 4. Record stats
      summary_results <- record_summary_statistics(
        t, agents, summary_results, params$G,
        emp$freq_coord_sq_total, emp$freq_coord_alt_total
      )
      
      
    }
    
    # tag summary with pair info
    summary_results$pair_g1 <- params$pair_g1
    summary_results$pair_g2 <- params$pair_g2
    summary_results$pair_label <- params$pair_label
    sim_results[[sim]] <- summary_results
  }
  
  sim_counter <- sim_counter + 1L
  key <- sprintf("params_%03d_%s", i, params$pair_label)
  results_list[[key]] <- list(
    parameters = params,
    simulations = sim_results
  )
}

saveRDS(results_list, file = file.path(output_dir, "results_list.rds"))
message("Saved results_list.rds with ", length(results_list), " parameter×pair cells.")
