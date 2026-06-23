###############################################################################
###############################################################################
###############################################################################

# Part of the project on social norms change in organizations
# by Faessler, von Flüe, Zehnder, Efferson

# Programmed by Lukas von Flüe

# functions.R
# functions used in the simulation

###############################################################################
###############################################################################
###############################################################################

###############################################################################

# GINI FUNCTION

###############################################################################

# If we want to measure inequality: Gini function. References: Gini (1912), Schmidt & Wichardt (2019).
# The argument of the function, x, will consist of an array with the agents' payoffs.

my_gini <- function(x) {
  n <- length(x)  
  
  sum_x <- sum(x) 
  
  sorted_x <- sort(x)  
  
  i <- 1:n
  numerator <- sum(2 * i * sorted_x)
  denominator <- n * sum_x
  
  gini <- numerator / denominator - (n + 1) / n
  
  return(gini)
}

###############################################################################

# Initialization functions

###############################################################################

# Helpers for groupwise Beta shapes
# Five named shapes:
#   EL  extreme-left skew  mode near 1
#   ER  extreme-right skew mode near 0
#   ML  modest-left        between EL and SYM
#   MR  modest-right       between ER and SYM
#   SYM symmetric          mode at 0.5
beta_shape_params <- function() {
  list(
    EL  = c(alpha = 80, beta = 2),   # mode ≈ 0.975
    ER  = c(alpha = 2,  beta = 80),  # mode ≈ 0.025
    ML  = c(alpha = 6,  beta = 3),   # mode ≈ 0.25
    MR  = c(alpha = 3,  beta = 6),   # mode ≈ 0.75
    SYM = c(alpha = 20, beta = 20)   # mode = 0.5
  )
}

# Fifteen non-duplicated group pairs for G = 2
# 1. EL–EL 2. EL–ER 3. EL–ML 4. EL–MR 5. EL–SYM 6. ER–ER 7. ER–ML 8. ER–MR 9. ER–SYM 10. ML–ML 11. ML–MR 12. ML–SYM 13. MR–MR 14. MR–SYM 15. SYM–SYM
# Accepts types_pair as c("EL","ER") or list(c("EL","ER"))
# Pair helpers (centralised)

# Accepts types_pair as c("EL","ER") or list(c("EL","ER"))
normalise_types_pair <- function(tp) {
  if (is.null(tp)) return(NULL)
  if (is.list(tp)) {
    stopifnot(length(tp) >= 1L)
    tp <- tp[[1]]
  }
  stopifnot(is.character(tp), length(tp) == 2L)
  tp
}


# Keep this helper to centralize the order of types
group_type_levels <- function() c("EL","ER","ML","MR","SYM")

# Return all unique unordered group-type pairs.
# include_self = TRUE -> include AA (EL-EL, ..., SYM-SYM)
all_group_pairs <- function(include_self = TRUE) {
  types <- group_type_levels()
  pairs <- list()
  for (i in seq_along(types)) {
    j_start <- if (include_self) i else (i + 1L)
    if (j_start > length(types)) next
    for (j in j_start:length(types)) {
      pairs[[length(pairs) + 1L]] <- c(types[i], types[j])
    }
  }
  pairs
}

# Canonicalize a pair to its unordered (unique) representation (EL-MR == MR-EL)
canonicalize_pair <- function(tp) {
  stopifnot(is.character(tp), length(tp) == 2L)
  types <- group_type_levels()
  ord <- order(match(tp, types))
  tp[ord]
}

# If types_pair is provided, we canonicalize it; otherwise pick from the unique list by id
resolve_group_pair <- function(combo_id = 1L, types_pair = NULL, include_self = TRUE) {
  if (!is.null(types_pair)) {
    if (is.list(types_pair)) {
      stopifnot(length(types_pair) >= 1L)
      types_pair <- types_pair[[1]]
    }
    return(canonicalize_pair(types_pair))
  }
  combos <- all_group_pairs(include_self = include_self)
  combo_id <- as.integer(combo_id)
  stopifnot(length(combo_id) == 1L, combo_id >= 1L, combo_id <= length(combos))
  combos[[combo_id]]
}

pairs_from_combo_ids <- function(ids, include_self = TRUE) {
  lapply(ids, function(id) resolve_group_pair(combo_id = id, include_self = include_self))
}

# Label uses the canonicalized order so EL-MR and MR-EL get the same label
pair_label <- function(pair) {
  p <- canonicalize_pair(pair)
  paste0("g1", p[1], "g2", p[2])
}



# Replacement for initializeAgents (groups first, then draw x_i per group)
# Signature simplified to reflect the new design. Call it from simulation.R.
# Note: this targets G = 2 as per our current plan.
initializeAgents <- function(
    N, G, a, b, d,
    combo_id = 1L,                 # choose one of the 10 non-duplicated pairs
    types_pair = NULL,             # or pass c("EL","ER"), c("ML","SYM"), etc.
    shape_overrides = NULL         # e.g., list(EL = c(alpha=100,beta=3))
) {
  if (G != 2L) stop("Current initializer targets G = 2. Extend mapping if needed.")
  total <- N + G  # keep your convention one manager per group
  
  # Build empty agent frame consistent with your structure
  agents <- data.frame(
    ID      = seq_len(total),
    x_i     = numeric(total),
    group   = integer(total),
    respond = integer(total),
    expSQ   = numeric(total),
    expAlt  = numeric(total),
    choice  = integer(total),
    payoff  = numeric(total),
    role    = rep("E", total),
    q       = numeric(total),
    stringsAsFactors = FALSE
  )
  
  # Shapes
  shapes <- beta_shape_params()
  if (!is.null(shape_overrides)) {
    for (nm in names(shape_overrides)) {
      stopifnot(nm %in% names(shapes))
      shapes[[nm]] <- shape_overrides[[nm]]
    }
  }
  pair <- resolve_group_pair(combo_id = combo_id, types_pair = types_pair)
  
  # Even split across groups including managers
  sizes <- rep(floor(total / G), G)
  remainder <- total - sum(sizes)
  if (remainder > 0) sizes[seq_len(remainder)] <- sizes[seq_len(remainder)] + 1
  
  # Assign groups
  agents$group <- rep(seq_len(G), times = sizes)
  
  # Draw x_i within each group from its Beta and scale to (0, d - b)
  idx_by_group <- split(seq_len(total), agents$group)
  for (g in seq_len(G)) {
    idx <- idx_by_group[[g]]
    key <- pair[g]
    pars <- shapes[[key]]
    xi <- rbeta(length(idx), shape1 = pars[["alpha"]], shape2 = pars[["beta"]])
    agents$x_i[idx] <- xi * (d - b)
  }
  
  # One random manager per group
  for (g in seq_len(G)) {
    idx <- which(agents$group == g)
    m <- sample(idx, 1)
    agents$role[m] <- "M"
  }
  
  
  ## NEW BLOCK: initialise choices correlated with x_i
  ## 0 = SQ, 1 = Alt
  ## P(Alt) = 1 - x_i_scaled, where x_i_scaled = x_i / (d - b)
  span <- d - b
  if (span <= 0) stop("initializeAgents: need d > b for scaling x_i.")
  x_scaled <- agents$x_i / span
  
  # Numerically guard probabilities to [0,1]
  prob_alt <- 1 - x_scaled
  prob_alt <- pmin(pmax(prob_alt, 0), 1)
  
  agents$choice <- rbinom(n = nrow(agents), size = 1, prob = prob_alt)
  
  agents
}



###############################################################################
# Plot helper for initial x_i densities (groups + pooled)
###############################################################################

# Requires ggplot2 (only when you actually call it)
plot_xi_densities <- function(agents, b, d, pair_used = NULL,
                              save = FALSE, outdir = "simulation_results/plots",
                              fname = NULL, width = 6.5, height = 4.5, dpi = 150) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plot_xi_densities(). Install it or set save=FALSE to skip plotting.")
  }
  ggplot2 <- asNamespace("ggplot2")
  
  rng_by_group <- split(seq_len(nrow(agents)), agents$group)
  dspan <- c(0, d - b)
  
  dens_df_list <- list()
  for (g in names(rng_by_group)) {
    xi <- agents$x_i[rng_by_group[[g]]]
    dd <- density(xi, from = dspan[1], to = dspan[2])
    dens_df_list[[paste0("G", g)]] <- data.frame(
      x = dd$x, density = dd$y, group = paste0("G", g), stringsAsFactors = FALSE
    )
  }
  dd_all <- density(agents$x_i, from = dspan[1], to = dspan[2])
  dens_df_list[["Pooled"]] <- data.frame(x = dd_all$x, density = dd_all$y, group = "Pooled")
  dens_df <- do.call(rbind, dens_df_list)
  
  # Legend labels using pair_used if given
  dens_df$group_label <- dens_df$group
  if (!is.null(pair_used) && length(pair_used) == 2L) {
    dens_df$group_label[dens_df$group == "G1"] <- paste0("G1 ", pair_used[1])
    dens_df$group_label[dens_df$group == "G2"] <- paste0("G2 ", pair_used[2])
  }
  dens_df$group_label[dens_df$group == "Pooled"] <- "Pooled"
  
  col_map <- c(
    setNames("#99c2a2", if (!is.null(pair_used)) paste0("G1 ", pair_used[1]) else "G1"),
    setNames("#f9caa7", if (!is.null(pair_used)) paste0("G2 ", pair_used[2]) else "G2"),
    setNames("#000000", "Pooled")
  )
  
  p <- ggplot2$ggplot(dens_df, ggplot2$aes(x = x, y = density, group = group_label)) +
    ggplot2$geom_line(ggplot2$aes(colour = group_label), linewidth = 0.8) +
    ggplot2$scale_colour_manual(values = col_map, guide = ggplot2$guide_legend(title = NULL)) +
    ggplot2$labs(
      title = paste0("Initial x_i densities ",
                     if (!is.null(pair_used)) paste0("[", paste0(pair_used, collapse = " vs "), "]") else ""),
      x = bquote( x[i]~"in"~~"(" * 0 * ", " * .(d - b) * ")" ),
      y = "Density"
    ) +
    ggplot2$theme_minimal(base_size = 12)
  
  if (isTRUE(save)) {
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
    if (is.null(fname)) {
      if (exists("pair_label")) {
        fname <- paste0("init_density_", pair_label(if (is.null(pair_used)) c("?", "?") else pair_used), ".png")
      } else {
        fname <- "init_density.png"
      }
    }
    ggplot2$ggsave(filename = file.path(outdir, fname), plot = p, width = width, height = height, dpi = dpi)
  } else {
    print(p)
  }
  invisible(p)
}



###############################################################################
###############################################################################
###############################################################################

# =====================================================================
# Plot theoretical x_i densities for all Beta-shape pairs (G = 2)
# Shows G1, G2, and pooled (50/50 mixture) for each of the 15
# unique unordered pairs.
#
# - Faceted overview (all pairs) if facet = TRUE
# - Optional separate plots for each pair if make_separate = TRUE
# =====================================================================
plot_theoretical_xi_pairs <- function(b, d,
                                      include_self   = TRUE,
                                      facet          = TRUE,
                                      make_separate  = FALSE,
                                      save           = FALSE,
                                      outdir         = "simulation_results/plots",
                                      fname          = "theoretical_xi_pairs.png",
                                      width          = 10,
                                      height         = 8,
                                      dpi            = 150,
                                      separate_prefix = "theoretical_xi_pair_",
                                      separate_width  = 6.5,
                                      separate_height = 4.5) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plot_theoretical_xi_pairs(). Install it or set save=FALSE to skip plotting.")
  }
  ggplot2 <- asNamespace("ggplot2")
  
  # Span of x_i support: (0, d - b), consistent with your initializer
  span <- d - b
  if (span <= 0) stop("plot_theoretical_xi_pairs: need d > b.")
  
  # Grid on x_i scale
  x_grid   <- seq(0, span, length.out = 1000L)
  x_scaled <- x_grid / span  # map to (0,1) for Beta pdf
  
  # Shapes and all group-type pairs
  shapes <- beta_shape_params()
  pairs  <- all_group_pairs(include_self = include_self)
  
  # Colour map: pooled = black
  col_map <- c(
    G1     = "#99c2a2",
    G2     = "#f9caa7",
    Pooled = "#000000"
  )
  
  dens_list <- list()
  
  for (k in seq_along(pairs)) {
    pair  <- pairs[[k]]              # e.g. c("EL","ER")
    pars1 <- shapes[[pair[1]]]
    pars2 <- shapes[[pair[2]]]
    
    # Beta pdf on (0,1), then rescale to (0, span) with Jacobian 1/span
    f1 <- dbeta(x_scaled, shape1 = pars1[["alpha"]], shape2 = pars1[["beta"]]) / span
    f2 <- dbeta(x_scaled, shape1 = pars2[["alpha"]], shape2 = pars2[["beta"]]) / span
    f_pooled <- 0.5 * (f1 + f2)
    
    lab <- pair_label(pair)  # e.g. "g1ELg2ER"
    
    dens_list[[length(dens_list) + 1L]] <- data.frame(
      x           = x_grid,
      density     = f1,
      group_label = "G1",
      pair_label  = lab,
      stringsAsFactors = FALSE
    )
    dens_list[[length(dens_list) + 1L]] <- data.frame(
      x           = x_grid,
      density     = f2,
      group_label = "G2",
      pair_label  = lab,
      stringsAsFactors = FALSE
    )
    dens_list[[length(dens_list) + 1L]] <- data.frame(
      x           = x_grid,
      density     = f_pooled,
      group_label = "Pooled",
      pair_label  = lab,
      stringsAsFactors = FALSE
    )
  }
  
  dens_df <- do.call(rbind, dens_list)
  
  if (save || make_separate) {
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  }
  
  # -------------------------
  # Faceted overview plot
  # -------------------------
  if (isTRUE(facet)) {
    p_all <- ggplot2$ggplot(
      dens_df,
      ggplot2$aes(x = x, y = density, colour = group_label)
    ) +
      ggplot2$geom_line(linewidth = 0.5) +
      ggplot2$facet_wrap(~ pair_label, ncol = 5) +
      ggplot2$scale_colour_manual(
        values = col_map,
        guide  = ggplot2$guide_legend(title = NULL)
      ) +
      ggplot2$labs(
        title = "Theoretical x_i densities for all group-type pairs",
        x     = bquote( x[i]~"in"~~"(" * 0 * ", " * .(span) * ")" ),
        y     = "Density"
      ) +
      ggplot2$theme_minimal(base_size = 11)
    
    if (isTRUE(save)) {
      ggplot2$ggsave(
        filename = file.path(outdir, fname),
        plot     = p_all,
        width    = width,
        height   = height,
        dpi      = dpi
      )
    } else {
      print(p_all)
    }
  }
  
  # -------------------------
  # Separate plots per pair
  # -------------------------
  if (isTRUE(make_separate)) {
    split_by_pair <- split(dens_df, dens_df$pair_label)
    
    for (lab in names(split_by_pair)) {
      df_pair <- split_by_pair[[lab]]
      
      p_pair <- ggplot2$ggplot(
        df_pair,
        ggplot2$aes(x = x, y = density, colour = group_label)
      ) +
        ggplot2$geom_line(linewidth = 0.7) +
        ggplot2$scale_colour_manual(
          values = col_map,
          guide  = ggplot2$guide_legend(title = NULL)
        ) +
        ggplot2$labs(
          title = paste0("Theoretical x_i densities [", lab, "]"),
          x     = bquote( x[i]~"in"~~"(" * 0 * ", " * .(span) * ")" ),
          y     = "Density"
        ) +
        ggplot2$theme_minimal(base_size = 11)
      
      if (isTRUE(save)) {
        fname_pair <- paste0(separate_prefix, lab, ".png")
        ggplot2$ggsave(
          filename = file.path(outdir, fname_pair),
          plot     = p_pair,
          width    = separate_width,
          height   = separate_height,
          dpi      = dpi
        )
      } else {
        print(p_pair)
      }
    }
  }
  
  invisible(NULL)
}




###############################################################################

# SUMMARY TABLE: Initialise the summary statistics table and create column for each of the G groups:

###############################################################################

initialize_summary_table <- function(t_max, G) {
  # Create an empty data frame
  summary_results <- data.frame(matrix(nrow = t_max, ncol = 0))
  
  # Store per-group choice frequencies and payoffs
  for (g in 1:G) {
    summary_results[[paste0("freq_choice_sq_group", g)]] <- rep(0, t_max)
    summary_results[[paste0("freq_choice_alt_group", g)]] <- rep(0, t_max)
    
    summary_results[[paste0("avg_payoff_sq_group", g)]] <- rep(0, t_max)
    summary_results[[paste0("avg_payoff_alt_group", g)]] <- rep(0, t_max)
    
    summary_results[[paste0("gini_coefficient_group", g)]] <- rep(NA, t_max)
    
    summary_results[[paste0("manager_choice_group", g)]] <- rep(NA, t_max)
  }
  
  # Store total coordination frequencies (across all agents)
  summary_results$freq_coord_sq_total <- rep(0, t_max)  # Placeholder for total coordination on SQ
  summary_results$freq_coord_alt_total <- rep(0, t_max)  # Placeholder for total coordination on Alt
  
  return(summary_results)
}


###############################################################################

# INTERVENTION FUNCTION: 

###############################################################################

# Target intervention based on intervention type (resistant/random/amenable) and apply intervention probabilistically based on agent's threshold values
# "targetScope" allows to either target phi*N agents (i.e. a fraction of phi of all employees) with targetScope="population", 
# or phi*(N/G) of all agents within groups with targetScope="group". 
# In other words, if the most amenable/resistant agents are targeted, we can do so either in whole population or within groups.

applyIntervention <- function(agents, G, phi, d, b,
                              interventionTarget = "M",
                              targetSuccess = 1,
                              targetScope = "population") {
  # Helper: map interventionTarget to role and resistance for employees
  role <- NA_character_
  targetResistance <- NA_integer_
  
  if (interventionTarget == "M") {
    role <- "M"
  } else if (interventionTarget == "E_Amenable") {
    role <- "E"; targetResistance <- 0L
  } else if (interventionTarget == "E_Resistant") {
    role <- "E"; targetResistance <- 1L
  } else if (interventionTarget == "E_Random") {
    role <- "E"; targetResistance <- 2L
  } else {
    stop("Unknown interventionTarget: ", interventionTarget)
  }
  
  # If targeting employees across the whole population
  if (role == "E" && targetScope == "population") {
    employee_indices <- which(agents$role == "E")
    num_targeted <- floor(length(employee_indices) * phi)
    
    if (num_targeted > 0) {
      if (targetResistance == 0L) {
        # Amenable: lowest x_i
        ordered_indices <- employee_indices[order(agents$x_i[employee_indices])]
        targeted_agents <- ordered_indices[1:num_targeted]
      } else if (targetResistance == 1L) {
        # Resistant: highest x_i
        ordered_indices <- employee_indices[order(agents$x_i[employee_indices], decreasing = TRUE)]
        targeted_agents <- ordered_indices[1:num_targeted]
      } else if (targetResistance == 2L) {
        # Random employees
        targeted_agents <- sample(employee_indices, num_targeted)
      } else {
        stop("Invalid targetResistance: ", targetResistance)
      }
    } else {
      return(agents)  # nothing to do
    }
    
    # If targeting employees within each group
  } else if (role == "E" && targetScope == "group") {
    targeted_agents <- integer(0)
    
    for (g in 1:G) {
      group_indices <- which(agents$group == g & agents$role == "E")
      num_targeted <- floor(length(group_indices) * phi)
      
      if (num_targeted > 0) {
        if (targetResistance == 0L) {
          ordered_indices <- group_indices[order(agents$x_i[group_indices])]
          targeted_agents <- c(targeted_agents, ordered_indices[1:num_targeted])
        } else if (targetResistance == 1L) {
          ordered_indices <- group_indices[order(agents$x_i[group_indices], decreasing = TRUE)]
          targeted_agents <- c(targeted_agents, ordered_indices[1:num_targeted])
        } else if (targetResistance == 2L) {
          targeted_agents <- c(targeted_agents, sample(group_indices, num_targeted))
        } else {
          stop("Invalid targetResistance: ", targetResistance)
        }
      }
    }
    
    # If targeting managers
  } else if (role == "M") {
    targeted_agents <- which(agents$role == "M")
    
  } else {
    stop("Unsupported combination in applyIntervention: role=", role,
         ", targetScope=", targetScope)
  }
  
  # Success logic for intervention
  if (targetSuccess == 0) {
    # Probability of success depends on x_i values
    prob_draw <- runif(length(targeted_agents))
    agents$respond[targeted_agents] <- ifelse(
      prob_draw <= (1 - (agents$x_i[targeted_agents] / (d - b))),
      1L, 0L
    )
  } else if (targetSuccess == 1) {
    # Guaranteed success for all targeted agents
    agents$respond[targeted_agents] <- 1L
  } else {
    stop("Invalid targetSuccess: ", targetSuccess)
  }
  
  return(agents)
}



###############################################################################

# Belief updating function for managers

###############################################################################

updateManagerBeliefs <- function(agents, G, n_M, wM_M) {
  
  # Precompute group-wise indices
  group_employee_indices <- lapply(1:G, function(g) which(agents$group == g & agents$role == "E"))
  group_manager_indices  <- sapply(1:G, function(g) which(agents$group == g & agents$role == "M"))
  
  # Loop over managers
  manager_indices <- which(agents$role == "M")
  
  for (i in manager_indices) {
    
    current_group <- agents$group[i]
    
    # 1️ Ingroup employees' mean (sampled with size n_M)
    in_indices <- group_employee_indices[[current_group]]
    if (length(in_indices) > 0) {
      sampled_in <- sample(agents$choice[in_indices], size = min(n_M, length(in_indices)), replace = FALSE)
      employee_mean_in <- mean(sampled_in)
    } else {
      employee_mean_in <- 0
    }
    
    # 2️ Outgroup managers' mean (full mean, no sampling)
    outgroup_manager_indices <- unlist(group_manager_indices[-current_group])
    manager_mean_out <- if (length(outgroup_manager_indices) > 0) mean(agents$choice[outgroup_manager_indices]) else 0
    
    # 3️ Combine
    q <- (1 - wM_M) * employee_mean_in + wM_M * manager_mean_out
    
    # Store the belief
    agents$q[i] <- q
  }
  
  return(agents)
}


###############################################################################

# Coordination game function for managers

###############################################################################

manager_coordination <- function(agents, G, a, b, d, h, g) {
  
  manager_indices <- which(agents$role == "M")
  
  # Ensure there are at least two managers for pairing
  if (length(manager_indices) >= 2) {
    
    # Shuffle managers for pairing
    shuffled_managers <- sample(manager_indices)
    
    # Split into pairs
    player_1 <- shuffled_managers[1:(length(shuffled_managers) / 2)]
    player_2 <- shuffled_managers[((length(shuffled_managers) / 2) + 1):length(shuffled_managers)]
    
    # Expected payoffs using correct beliefs
    agents$expSQ[player_1] <- ifelse(
      agents$respond[player_1] == 1,
      g,
      ((1 - agents$q[player_1]) * (a + agents$x_i[player_1])) +
        (agents$q[player_1] * (b + agents$x_i[player_1]))
    )
    
    agents$expAlt[player_1] <- ifelse(
      agents$respond[player_1] == 1,
      h,
      ((1 - agents$q[player_1]) * a) +
        (agents$q[player_1] * d)
    )
    
    agents$expSQ[player_2] <- ifelse(
      agents$respond[player_2] == 1,
      g,
      ((1 - agents$q[player_2]) * (a + agents$x_i[player_2])) +
        (agents$q[player_2] * (b + agents$x_i[player_2]))
    )
    
    agents$expAlt[player_2] <- ifelse(
      agents$respond[player_2] == 1,
      h,
      ((1 - agents$q[player_2]) * a) +
        (agents$q[player_2] * d)
    )
    
    # Determine choice based on expected payoffs
    agents$choice[player_1] <- ifelse(agents$expAlt[player_1] >= agents$expSQ[player_1], 1, 0)
    agents$choice[player_2] <- ifelse(agents$expAlt[player_2] >= agents$expSQ[player_2], 1, 0)
    
    # Compute payoffs
    payoff_player_1 <- numeric(length(player_1))
    payoff_player_2 <- numeric(length(player_2))
    
    both_sq <- agents$choice[player_1] == 0 & agents$choice[player_2] == 0
    both_alt <- agents$choice[player_1] == 1 & agents$choice[player_2] == 1
    
    # Payoffs for coordination on SQ
    payoff_player_1[both_sq] <- a + agents$x_i[player_1][both_sq]
    payoff_player_2[both_sq] <- a + agents$x_i[player_2][both_sq]
    
    # Payoffs for coordination on Alt
    payoff_player_1[both_alt] <- d
    payoff_player_2[both_alt] <- d
    
    # Miscoordination cases
    miscoord_1 <- agents$choice[player_1] == 0 & agents$choice[player_2] == 1
    miscoord_2 <- agents$choice[player_1] == 1 & agents$choice[player_2] == 0
    
    payoff_player_1[miscoord_1] <- b + agents$x_i[player_1][miscoord_1]
    payoff_player_2[miscoord_1] <- a
    
    payoff_player_1[miscoord_2] <- a
    payoff_player_2[miscoord_2] <- b + agents$x_i[player_2][miscoord_2]
    
    # Assign payoffs to agents
    agents$payoff[player_1] <- payoff_player_1
    agents$payoff[player_2] <- payoff_player_2
  }
  
  # Final step: Assign h to all agents with respond == 1
  agents$payoff[agents$respond == 1] <- h
  
  return(agents)
}



###############################################################################

# Belief updating function for employees

###############################################################################

updateEmployeeBeliefs <- function(agents, G, n_E_in, n_E_out, wE_inE, wE_outE, probOut) {
  
  # Precompute indices for each group
  group_employee_indices <- lapply(1:G, function(g) which(agents$group == g & agents$role == "E"))
  group_manager_indices <- as.integer(sapply(1:G, function(g) which(agents$group == g & agents$role == "M")))
  
  # Loop over employees only
  employee_indices <- which(agents$role == "E")
  
  for (i in employee_indices) {
    
    current_group <- agents$group[i]
    
    # 1️ Ingroup sampling (using n_E_in)
    in_indices <- setdiff(group_employee_indices[[current_group]], i)
    if (length(in_indices) > 0) {
      sampled_in <- sample(agents$choice[in_indices], size = min(n_E_in, length(in_indices)), replace = FALSE)
      employee_mean_in <- mean(sampled_in)
    } else {
      employee_mean_in <- 0
    }
    manager_index_in <- group_manager_indices[current_group]
    manager_choice_in <- if (length(manager_index_in) == 1) agents$choice[manager_index_in] else 0
    q_in <- wE_inE * employee_mean_in + (1 - wE_inE) * manager_choice_in
    
    # 2️ Outgroup sampling (using n_E_out)
    all_outgroup_samples <- c()
    outgroup_manager_choices <- c()
    for (g in setdiff(1:G, current_group)) {
      
      out_indices <- group_employee_indices[[g]]
      if (length(out_indices) > 0) {
        sampled_out <- sample(agents$choice[out_indices], size = min(n_E_out, length(out_indices)), replace = FALSE)
        all_outgroup_samples <- c(all_outgroup_samples, sampled_out)
      }
      
      # Save outgroup manager choice
      manager_index_out <- group_manager_indices[g]
      if (length(manager_index_out) == 1) {
        outgroup_manager_choices <- c(outgroup_manager_choices, agents$choice[manager_index_out])
      }
    }
    
    employee_mean_out <- if (length(all_outgroup_samples) > 0) mean(all_outgroup_samples) else 0
    manager_mean_out <- if (length(outgroup_manager_choices) > 0) mean(outgroup_manager_choices) else 0
    
    q_out <- wE_outE * employee_mean_out + (1 - wE_outE) * manager_mean_out
    
    # 3 Combine
    agents$q[i] <- (1 - probOut) * q_in + probOut * q_out
  }
  
  return(agents)
}




###############################################################################

# Choice based on belief

###############################################################################

###############################################################################
# Choice based on belief (probabilistic / logit)
###############################################################################
belief_to_choice <- function(agents, a, b, d, h, g, choice_lambda = 5) {
  employee_indices <- which(agents$role == "E")
  
  # Expected payoffs for SQ
  agents$expSQ[employee_indices] <- ifelse(
    agents$respond[employee_indices] == 1,
    g,
    ((1 - agents$q[employee_indices]) * (a + agents$x_i[employee_indices])) +
      (agents$q[employee_indices] * (b + agents$x_i[employee_indices]))
  )
  
  # Expected payoffs for Alt
  agents$expAlt[employee_indices] <- ifelse(
    agents$respond[employee_indices] == 1,
    h,
    ((1 - agents$q[employee_indices]) * a) +
      (agents$q[employee_indices] * d)
  )
  
  # --- Probabilistic choice ---
  # Stable logistic on the utility difference:
  # P(Alt) = 1 / (1 + exp(-lambda * (U_alt - U_sq)))
  # lambda = 0 -> random (0.5); lambda -> ∞ -> deterministic argmax
  if (isTRUE(is.infinite(choice_lambda))) {
    # back to deterministic threshold rule
    agents$choice[employee_indices] <- ifelse(
      agents$expAlt[employee_indices] >= agents$expSQ[employee_indices], 1, 0
    )
  } else {
    diffU <- agents$expAlt[employee_indices] - agents$expSQ[employee_indices]
    p_alt <- plogis(choice_lambda * diffU)  # numerically stable
    agents$choice[employee_indices] <- rbinom(length(employee_indices), size = 1, prob = p_alt)
  }
  
  agents
}


###############################################################################

# Coordination game function for employees

###############################################################################

employee_coordination <- function(agents, G, a, b, d, h, g, probOut) {
  
  freq_coord_sq_total <- 0
  freq_coord_alt_total <- 0
  
  # Step 1: Prepare group-wise pairs
  group_employee_indices <- lapply(1:G, function(g) which(agents$group == g & agents$role == "E"))
  
  # Make pairs in each group
  group_pairs <- list()
  used_pairs <- list()
  for (g in 1:G) {
    indices <- group_employee_indices[[g]]
    shuffled <- sample(indices)
    if (length(shuffled) %% 2 != 0) shuffled <- shuffled[-length(shuffled)]  # Drop last if odd
    pairs <- matrix(shuffled, ncol = 2, byrow = TRUE)
    group_pairs[[g]] <- pairs
    used_pairs[[g]] <- rep(FALSE, nrow(pairs))  # Track which pairs have been used in cross-group
  }
  
  # Step 2: Loop over each group’s pairs
  for (g in 1:G) {
    num_pairs <- nrow(group_pairs[[g]])
    if (num_pairs == 0) next  # Skip empty groups
    
    for (i in 1:num_pairs) {
      if (used_pairs[[g]][i]) next  # Already matched in cross-group
      
      pair_1 <- group_pairs[[g]][i, ]
      
      # Decide whether to attempt cross-group
      if (runif(1) < probOut) {
        # Try to find an available outgroup
        outgroups <- setdiff(1:G, g)
        outgroups <- sample(outgroups)  # Randomise order
        
        matched <- FALSE
        for (og in outgroups) {
          available_pairs <- which(!used_pairs[[og]])
          if (length(available_pairs) > 0) {
            # Found a pair to match with
            chosen_j <- sample(available_pairs, 1)
            pair_2 <- group_pairs[[og]][chosen_j, ]
            
            # Play coordination game between pair_1 and pair_2
            agents <- play_coordination_game(agents, pair_1, pair_2, a, b, d, h, g)
            freq_coord_sq_total <- freq_coord_sq_total + count_coord_sq(agents, pair_1, pair_2)
            freq_coord_alt_total <- freq_coord_alt_total + count_coord_alt(agents, pair_1, pair_2)
            
            # Mark both pairs as used
            used_pairs[[g]][i] <- TRUE
            used_pairs[[og]][chosen_j] <- TRUE
            matched <- TRUE
            break  # Done matching this pair
          }
        }
        
        if (!matched) {
          # Fallback: play within-group
          agents <- play_coordination_game(agents, pair_1, NULL, a, b, d, h, g)
          freq_coord_sq_total <- freq_coord_sq_total + count_coord_sq(agents, pair_1)
          freq_coord_alt_total <- freq_coord_alt_total + count_coord_alt(agents, pair_1)
        }
      } else {
        # Within-group interaction
        agents <- play_coordination_game(agents, pair_1, NULL, a, b, d, h, g)
        freq_coord_sq_total <- freq_coord_sq_total + count_coord_sq(agents, pair_1)
        freq_coord_alt_total <- freq_coord_alt_total + count_coord_alt(agents, pair_1)
      }
    }
  }
  
  # Final step: any remaining unused pairs in each group play within-group
  for (g in 1:G) {
    leftover_pairs <- group_pairs[[g]][which(!used_pairs[[g]]), , drop = FALSE]
    if (nrow(leftover_pairs) > 0) {
      for (i in 1:nrow(leftover_pairs)) {
        pair <- leftover_pairs[i, ]
        agents <- play_coordination_game(agents, pair, NULL, a, b, d, h, g)
        freq_coord_sq_total <- freq_coord_sq_total + count_coord_sq(agents, pair)
        freq_coord_alt_total <- freq_coord_alt_total + count_coord_alt(agents, pair)
      }
    }
  }
  
  # Always assign h to all agents with respond == 1
  agents$payoff[agents$respond == 1] <- h
  
  return(list(
    agents = agents,
    freq_coord_sq_total = freq_coord_sq_total,
    freq_coord_alt_total = freq_coord_alt_total
  ))
}



###############################################################################

# Play coordination game function used in employee_coordination function

###############################################################################

play_coordination_game <- function(agents, pair_1, pair_2 = NULL, a, b, d, h, g) {
  
  if (is.null(pair_2)) {
    # Within-group: pair_1 plays together
    player_1 <- pair_1[1]
    player_2 <- pair_1[2]
    
    both_sq <- agents$choice[player_1] == 0 & agents$choice[player_2] == 0
    both_alt <- agents$choice[player_1] == 1 & agents$choice[player_2] == 1
    
    if (both_sq) {
      agents$payoff[player_1] <- a + agents$x_i[player_1]
      agents$payoff[player_2] <- a + agents$x_i[player_2]
    } else if (both_alt) {
      agents$payoff[player_1] <- d
      agents$payoff[player_2] <- d
    } else {
      if (agents$choice[player_1] == 0) {
        agents$payoff[player_1] <- b + agents$x_i[player_1]
        agents$payoff[player_2] <- a
      } else {
        agents$payoff[player_1] <- a
        agents$payoff[player_2] <- b + agents$x_i[player_2]
      }
    }
    
  } else {
    # Cross-group: pair_1 plays with pair_2
    for (i in 1:2) {
      p1 <- pair_1[i]
      p2 <- pair_2[i]
      
      both_sq <- agents$choice[p1] == 0 & agents$choice[p2] == 0
      both_alt <- agents$choice[p1] == 1 & agents$choice[p2] == 1
      
      if (both_sq) {
        agents$payoff[p1] <- a + agents$x_i[p1]
        agents$payoff[p2] <- a + agents$x_i[p2]
      } else if (both_alt) {
        agents$payoff[p1] <- d
        agents$payoff[p2] <- d
      } else {
        if (agents$choice[p1] == 0) {
          agents$payoff[p1] <- b + agents$x_i[p1]
          agents$payoff[p2] <- a
        } else {
          agents$payoff[p1] <- a
          agents$payoff[p2] <- b + agents$x_i[p2]
        }
      }
    }
  }
  
  return(agents)
}


###############################################################################

# Counting helper function used in employee_coordination function

###############################################################################

count_coord_sq <- function(agents, pair_1, pair_2 = NULL) {
  if (is.null(pair_2)) {
    as.integer(agents$choice[pair_1[1]] == 0 & agents$choice[pair_1[2]] == 0)
  } else {
    sum(
      agents$choice[pair_1[1]] == 0 & agents$choice[pair_2[1]] == 0,
      agents$choice[pair_1[2]] == 0 & agents$choice[pair_2[2]] == 0
    )
  }
}

count_coord_alt <- function(agents, pair_1, pair_2 = NULL) {
  if (is.null(pair_2)) {
    as.integer(agents$choice[pair_1[1]] == 1 & agents$choice[pair_1[2]] == 1)
  } else {
    sum(
      agents$choice[pair_1[1]] == 1 & agents$choice[pair_2[1]] == 1,
      agents$choice[pair_1[2]] == 1 & agents$choice[pair_2[2]] == 1
    )
  }
}







###############################################################################

# Record key metrics over time

###############################################################################

record_summary_statistics <- function(t, agents, summary_results, G, freq_coord_sq_total, freq_coord_alt_total) {
  for (g in 1:G) {
    # Filter only employees for frequency calculations
    employee_indices <- which(agents$group == g & agents$role == "E")
    manager_index <- which(agents$group == g & agents$role == "M")
    
    # Frequency of choices within the group (employees only)
    summary_results[t, paste0("freq_choice_sq_group", g)] <- sum(agents$choice[employee_indices] == 0)
    summary_results[t, paste0("freq_choice_alt_group", g)] <- sum(agents$choice[employee_indices] == 1)
    
    # Average payoffs for those choosing SQ and Alt
    num_agents_sq <- sum(agents$choice[employee_indices] == 0)
    summary_results[t, paste0("avg_payoff_sq_group", g)] <- ifelse(num_agents_sq > 0, 
                                                                   mean(agents$payoff[employee_indices][agents$choice[employee_indices] == 0]), 
                                                                   NA)
    
    num_agents_alt <- sum(agents$choice[employee_indices] == 1)
    summary_results[t, paste0("avg_payoff_alt_group", g)] <- ifelse(num_agents_alt > 0, 
                                                                    mean(agents$payoff[employee_indices][agents$choice[employee_indices] == 1]), 
                                                                    NA)
    
    # Gini coefficient per group
    summary_results[t, paste0("gini_coefficient_group", g)] <- my_gini(agents$payoff[employee_indices])
    
    # Track manager choice
    summary_results[t, paste0("manager_choice_group", g)] <- agents$choice[manager_index]
  }
  
  # Store total coordination on SQ and Alt
  summary_results[t, "freq_coord_sq_total"] <- freq_coord_sq_total
  summary_results[t, "freq_coord_alt_total"] <- freq_coord_alt_total
  
  return(summary_results)
}


