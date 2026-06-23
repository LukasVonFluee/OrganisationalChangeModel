###############################################################################
###############################################################################
###############################################################################

# Part of the project on social norms change in organizations
# by Faessler, von Flüe, Efferson, and Zehnder

# Programmed by Lukas von Flüe

# This file: config.R
# This script defines the parameters and is used by the script simulation.R

###############################################################################
###############################################################################
###############################################################################

################################################################################
# Parameters
################################################################################

# Note, parameters are specified as arrays so that we can define multiple values
# Those values are then used in the simulation.R script for the parameter grid

################################################################################
# Number of simulations
n_sim <- c(50)
################################################################################
# Timing parameters
t_max <- c(100)  # Number of periods
interventionTime <- c(10)  # Time at which the intervention is applied
################################################################################

################################################################################
# Population size parameters:
N <- c(100)  # Total number of agents
G <- c(2) # number of groups
groupSize <- c(N/G) # group size
n_E_in <- c(5) # sample size that employees use for their ingroup
n_E_out <- c(5) # sample size that employees use for outgroup employees
n_M <- c(N/G) # sample size for managers
################################################################################

################################################################################
# Coordination game parameters:
a <- c(0.75)
b <- c(0)
d <- c(1)
g <- c(0)  # Payoff for choosing SQ for targeted agents who respond
h <- c(2)  # Payoff for choosing Alt for targeted agents who respond
################################################################################

################################################################################
# Belief weights for employees:

# Ingroup beliefs:
wE_inE <- c(0.2, 0.8) # Weight on ingroup employees
wE_inM <- 1 - wE_inE # Weight on ingroup manager

# Outgroup beliefs:
wE_outE <- c(1) # Weight on outgroup employees
wE_outM <- 1 - wE_outE # We agreed to set wE_outM to zero.

# Belief weights for managers:

wM_M <- c(0.2, 0.8) # Weight put on other manager
wM_inE <- 1 - wM_M # Weight put on ingroup employees. We agreed to set weight put on outgroup employees to zero.
################################################################################

################################################################################
# For beliefs to choice function:
# Logit "rationality"/inverse temperature (0 = random; larger = closer to deterministic)
# choice_lambda = Inf # returns original deterministic rule
choice_lambda <- c(Inf)   # we can try 1, 3, 5, 10, Inf for sensitivity

################################################################################

# Intervention parameters:
# phi <- c(0.2, 0.4, 0.6, 0.8) # Size of intervention
phi <- c(0.25, 0.5, 0.75) # use this for testing

# Collapsed intervention target:
#   "E_Amenable" -> target employees with low x_i
#   "E_Resistant" -> target employees with high x_i
#   "E_Random" -> target a random set of employees
#   "M" -> target managers
interventionTarget <- c("E_Amenable", "E_Resistant", "E_Random", "M")

targetSuccess <- c(0,1) 
# If = 0, the likelihood of success is proportional to (1 - x_i/(d - b))
# If = 1, the intervention is successful with certainty for targeted agents.

################################################################################

################################################################################
# Cross-group interactions:
probOut <- c(0, 0.25, 0.5, 0.75) # Probability that agents interact across groups for the coordination game in each period
################################################################################

# Beta distribution for x_i values
# Choose ONE of the three blocks below

# Five named shapes:
#   EL  extreme-left skew  mode near 1
#   ER  extreme-right skew mode near 0
#   ML  modest-left        between EL and SYM
#   MR  modest-right       between ER and SYM
#   SYM symmetric          mode at 0.5

# 1) RUN A SPECIFIC PAIR (enable this, disable the others)
# types_pair <- c("EL","ER")
# run_all_pairs <- FALSE
# combo_ids <- 1:10        # ignored when run_all_pairs == FALSE and types_pair is set

# 2) OR: RUN ALL PAIRS (disable 1), enable this)
types_pair <- NULL
run_all_pairs <- TRUE

# 3) OR: FALL BACK TO A SINGLE COMBO BY ID (disable 1) and 2); handled in simulation.R default)
# types_pair <- NULL
# run_all_pairs <- FALSE
# (simulation.R will default to combo_id 1 => EL–ER)

plot_init <- FALSE # will save one PNG per simulation initialisation (per param row × sim). Flip to FALSE before big runs

################################################################################
# Controlling whether intervention size applies to whole population or within groups:
targetScope <- c("population")
# targetScope <- c("population","group")
# "targetScope" allows to either target phi*N agents from the whole population (i.e. a fraction of phi of all employees) with targetScope="population", 
# or phi*(N/G) of all agents within groups with targetScope="group". 
# In other words, if the most amenable/resistant agents are targeted, we can do so either in whole population or within groups. We probably leave this at targetScope="population"