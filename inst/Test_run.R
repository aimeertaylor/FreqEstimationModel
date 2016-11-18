rm(list=ls(all=TRUE)) # Remove all variables from workspace
library(FreqEstimationModel) # Load the haplotype frequency estimation package

# TO-DO: Set working directory to source file location

# Select data to run; returns data_summary
no_markers <- 3 
data_seed <- 2  
load(sprintf('./Data_no_markers%s_seed%s.RData', no_markers, data_seed)) # Load data

# Summary of simulated data
str(data_summary)
# Note that, in the simulated data, the columns with prevalance data
# are labeled locas1... The MOI is the 'true', simulated MOI
# The subsequent columns are the 'true', simulated haplotype counts
head(data_summary$Data)

thinning_interval <- 1 # Change this variable to increase the number of iterations per chain that are not saved in memory
no_traces_preburnin <- 10000 # For more traces but managable pdfs, don't exceed 10k and increase thinning interval instead.

# MCMC variables
no_mcmc_chains <- 3 # Number of mcmcm chains to run
parallel <- FALSE # Set to true if running code in parallel (this option isn't active yet)
NGS <- FALSE # Set to true if data are in NGS format (this option isn't active yet)
log_like_zero <- FALSE # QC check; set log(p(yi|ai)) to zero (only impacts NGS) and wipes data (see below)
mcmc_variable_list <- list(no_mcmc_chains = no_mcmc_chains,
                           no_traces_preburnin = no_traces_preburnin,
                           thinning_interval = thinning_interval,
                           NGS = NGS,
                           log_like_zero = log_like_zero)


# MOI variables
if(log_like_zero){ # When log_like_zero == TRUE, code should return the prior, hence use a prior that is easy to eye-ball
  moi_prior <- 'Uniform'
} else {
  moi_prior <- 'Poisson' # Choose between 'Uniform', 'Poisson', 'Geometric' or 'nBinomial'
}
moi_max <- 8 # Maximum MOI regarded as possible by the model (haven't tested beyond 20)
moi_hyperparameter <- 3 # Specify mean MOI (Hyperparameter for the prior on the MOI)
moi_size_hyperparameter <- 0.5 # Only applies if moi_prior == 'nBinomial' (Hyperparameter for the prior on the MOI in the prior is the negative Binomial)
moi_prior_min2 <- NULL # Specify the lower bound for the MOI per individual
moi_initial <- NULL # If null, the initial vector of multiplicties of infection (MOI) is set internally, otherwise set to inputted moi_initial
moi_list <- list(moi_hyperparameter = moi_hyperparameter,
                 moi_size_hyperparameter = moi_size_hyperparameter,
                 moi_prior = moi_prior,
                 moi_max = moi_max,
                 moi_prior_min2 = moi_prior_min2,
                 moi_initial = moi_initial)


# Preprocess data
if(!NGS | log_like_zero){ # Set to true to augment partial observations rather than discard them
  augment_missing_data <- TRUE
} else {
  augment_missing_data <- FALSE
}
processed_data_list <- preprocess_data(data_summary,
                                       log_like_zero,
                                       NGS,
                                       augment_missing_data,
                                       moi_prior_min2)


# frequency variables
frequency_hyperparameter <-rep(1,processed_data_list$no_haplotypes) # The Parameter vector for the dirichlet prior on the frequency estimate
frequency_initial <- NULL # If null, external_frequency_initial set internally, otherwise set to inputted external_frequency_initial
frequency_list <- list(frequency_hyperparameter = frequency_hyperparameter,
                       frequency_initial = frequency_initial)


# Run MCMC
set.seed(1) # Set seed for reproducibility
results <- mcmc_sampling_parallel(processed_data_list,
                         moi_list,
                         frequency_list,
                         mcmc_variable_list,
                         cores_max = 3)


# Save results
arguments <- names(as.list(args(mcmc_sampling)))
arguments <- arguments[-length(arguments)]
if(!log_like_zero){
  dataset<-'Data'
} else {
  dataset<-'QC'
} # Class character: set to dataset name in order to save pdf with data set name if REAL=TRUE and to apend if REAL=FALSR
filename <- sprintf('%s_%sIterations_NoMarkers%s_NGS%s_LogLikeZero%s_seed%s',dataset, no_traces_preburnin*thinning_interval, no_markers, NGS, log_like_zero, data_seed, moi_prior)
save(list = c('results', arguments), file = sprintf('./%s.RData', filename))



# Generate pdf of results
source('./visualise_results.R')
visualise_results(results,
                  data_summary,
                  PDF = sprintf('./%s.pdf', filename),
                  child = FALSE,
                  augment_missing_data,
                  Simulated = TRUE)


# Convergence diagnostics
mcmc_frequency_chains <- mcmc.list(lapply(alply(results$genotype_freq_store_chains, no_mcmc_chains), mcmc))
gelman.diag(mcmc_frequency_chains, transform = TRUE, multivariate = FALSE) # Transformed (using log or logit) to improve the normality of the distribution







