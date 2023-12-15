rm(list=ls(all=TRUE)) # Remove all variables from workspace
library(FreqEstimationModel) # Load the haplotype frequency estimation package

# TO-DO: Set working directory to source file location

# Select data to run; returns data_summary
no_markers <- 3
data_seed <- 2
load(sprintf('./Data_no_markers%s_seed%s.RData', no_markers, data_seed)) # Load data

# Summary of simulated data
str(data_summary)

# Note that, in the simulated data, the columns with prevalence data
# are labelled locus1... The MOI is the 'true' simulated MOI.
# The subsequent columns are the 'true' simulated haplotype counts.
head(data_summary$Data)

thinning_interval <- 1 # Change this variable to increase the number of iterations per chain that are not saved in memory
no_traces_preburnin <- 10000 # For more traces but manageable pdf plots, don't exceed 10k and increase thinning interval instead.

# MCMC variables
no_mcmc_chains <- 3 # Number of MCMC chains to run
parallel <- FALSE # Set to true if running code in parallel (this option isn't active yet)
NGS <- FALSE # Set to true if data are in NGS format (this option isn't active yet)
log_like_zero <- FALSE # QC check; set log(p(yi|ai)) to zero (only impacts NGS) and wipes data (see below)
mcmc_variable_list <- list(no_mcmc_chains = no_mcmc_chains,
                           no_traces_preburnin = no_traces_preburnin,
                           thinning_interval = thinning_interval,
                           NGS = NGS,
                           log_like_zero = log_like_zero)


# MOI variables
if (log_like_zero) { # When log_like_zero == TRUE, code should return the prior, hence use a prior that is easy to eye-ball
  moi_prior <- 'Uniform'
} else {
  moi_prior <- 'Poisson' # Choose between 'Uniform', 'Poisson', 'Geometric' or 'nBinomial'
}

moi_max <- 8 # Maximum MOI regarded as possible by the model (I haven't tested beyond 20)
moi_hyperparameter <- 3 # Specify population-average MOI (parameter of the prior on the MOI)
moi_size_hyperparameter <- 0.5 # Only applies if moi_prior == 'nBinomial' (hyperparameter for the prior on the MOI if the prior is the negative Binomial)
moi_prior_min2 <- NULL # Specify the lower bound for the MOI per individual
moi_initial <- NULL # If null, the initial vector of multiplicities of infection (MOI) is set internally, otherwise set to input moi_initial
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
frequency_hyperparameter <-rep(1,processed_data_list$no_haplotypes) # The Parameter vector for the Dirichlet prior on the frequency estimate
frequency_initial <- NULL # If null, external_frequency_initial set internally, otherwise set to input external_frequency_initial
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
if(log_like_zero){
  dataset<-'QC'
} else {
  dataset <- 'Data'
}
filename <- sprintf('%s_%sIterations_NoMarkers%s_NGS%s_LogLikeZero%s_seed%s',dataset, no_traces_preburnin*thinning_interval, no_markers, NGS, log_like_zero, data_seed)
save(list = c('results', arguments), file = sprintf('./%s.RData', filename))


# ==============================================================================
# Quick example of estimates that can be extracted from results.
# N.B. Before reporting estimates, MCMC convergence should be checked; see below.
# ==============================================================================

# Install, load, and attach packages that are helpful for unpacking results
if(!require("plyr")) install.packages("plyr")
if(!require("coda")) install.packages("coda")
if(!require("abind")) install.packages("abind")

# Generate numerical approximations of posteriors by removing burn-in
burnin <- 1:(0.5*mcmc_variable_list$no_traces_preburnin) # remove first half
if(mcmc_variable_list$no_mcmc_chains > 1){
  alply_genotype_freq_store_chains_burnin <- plyr::alply(results$genotype_freq_store_chains[-burnin,,],3)
}else{
  alply_genotype_freq_store_chains_burnin <- results$genotype_freq_store_chains[-burnin,,]
}


# Use mcmc.list to represent frequency chains across parallel runs (see ?coda::mcmc.list)
mcmc_frequency_chains <- coda::mcmc.list(lapply(alply_genotype_freq_store_chains_burnin, # 3 for splitting by third dimension (nothing to do with no. of chains)
                                                coda::mcmc,
                                                start = (max(burnin)+1)*mcmc_variable_list$thinning_interval,
                                                end = mcmc_variable_list$no_traces_preburnin*mcmc_variable_list$thinning_interval,
                                                thin = mcmc_variable_list$thinning_interval))

mcmc_As <- abind::abind(plyr::alply(results$genotype_count_store_chains[-burnin,,,],4), along = 1) # Haplotype counts for all chains excluding burn in
mcmc_mois <- apply(mcmc_As, c(1,2), sum)

# Extract population-level frequency estimates:
# Haplotype names are the rownames. For example, 001 represents a haplotype across three biallelic loci where zero
# represent whatever the user encodes as zero in the input data (e.g. reference
# allele, minor allele, sensitive allele) and one represents whatever the user
# encodes as one (e.g. non-reference allele, major allele, resistant allele)
pop_freq <- cbind(mean = summary(mcmc_frequency_chains)$statistics[,"Mean"],
                      median = summary(mcmc_frequency_chains)$quantiles[,3],
                      "CI2.5%" = summary(mcmc_frequency_chains)$quantiles[,1],
                      "CI97.5%" = CI_upper<-summary(mcmc_frequency_chains)$quantiles[,5])

# Extract per-sample MOI estimates:
compute_mode <- function(x) {z <- table(x); names(z)[which.max(z)]}
qprobs <- c(0.025, 0.5, 0.975)
ind_MOI <- data.frame(colMeans(mcmc_mois),
                          as.numeric(apply(mcmc_mois, 2, compute_mode)),
                          t(apply(mcmc_mois, 2, quantile, probs = qprobs)))
colnames(ind_MOI) <- c("mean","mode","median","CI2.5%","CI97.5%")

# Extract per-sample most likely combination of phased haplotypes
mode_hap_count_chr <- apply(mcmc_As, 2, function(x) {compute_mode(apply(x, 1, paste, collapse = ""))})
ind_As<- t(sapply(mode_hap_count_chr, function(x) as.numeric(strsplit(x, split = "")[[1]])))
colnames(mode_hap_count_num) <- dimnames(mcmc_As)[[3]]

# Prevalence estimates
# The prevalence of a haplotype is the proportion of samples that contain that haplotype
# Since we can't observe haplotypes directly, we cannot compute it directly
# However, there are various ways we can estimate it using the results above

# Using the population-level frequencies at the population-level:
# We can estimate population-level prevalence using the population-frequencies and
# some population measure of the moi:
# The probability that one or more clones in the infection have the ith haplotype is one minus the probability that no clones
# in the infection have the ith haplotype:
pop_prev <- 1-(1-pop_freq)^median(mcmc_mois)

# At the sample-level using per-sample haplotype counts:
ind_prev <- t(apply(mcmc_As, 2, function(x) colMeans(x > 0)))
ind_freq <- t(apply(mcmc_As, 2, function(x) colMeans(x / rowSums(x))))

# Print posterior estimates:
ind_MOI # sample-level MOI estimates
ind_As # sample-level estimates of phased haplotype counts
ind_prev # individual prevalence estimates
ind_freq # individual frequency estimates
pop_freq # population-level frequency estimates
pop_prev # population prevalence estimates

# At the population-level using per-sample haplotype count estimates:
plot(x = pop_prev[,"median"], y = colMeans(ind_prev), pch = 20, bty = "n",
     xlim = c(0,1), ylim = c(0,1),
     main = "Probability that a haplotype is found in any given infection",
     xlab = "based on posterior population-level median estimates",
     ylab = "based on posterior per-sample haplotype count distribution")
text(x = pop_prev[,"median"], y = colMeans(ind_prev),
     labels = rownames(pop_prev), pos = 4, xpd = NA)
abline(a = 0, b = 1, lty = "dashed")

# At the population-level using per-sample haplotype frequency estimates
plot(x = pop_freq[,"mean"], y = colMeans(ind_freq), pch = 20, bty = "n",
     xlim = c(0,1), ylim = c(0,1),
     main = "Population-level haplotype frequency",
     xlab = "Parameter estimate",
     ylab = "Based on posterior per-sample haplotype frequency estimates")
text(x = pop_freq[,"mean"], y = colMeans(ind_freq),
     labels = rownames(pop_prev), pos = 4, xpd = NA)
abline(a = 0, b = 1, lty = "dashed")




# ==============================================================================
# Full summary of frequency results
# To check for convergence etc. the mcmc chains should be interrogated before
# reporting results. For convenience during my PhD, I atomised this using
# visualise_results() - see below - which I should have called interrogate
# results or something similar. To see how the interrogation is done manually,
# set child = FALSE, Simulated = TRUE, open visualise_results.R, and run the
# code chunk-by-chunk.
# ==============================================================================
# Generate pdf of results
source('./visualise_results.R') # Should have called this process results or similar
visualise_results(results,
                  data_summary,
                  PDF = sprintf('./%s.pdf', filename),
                  child = FALSE, # This relates to chapter five of my PhD thesis
                  augment_missing_data,
                  Simulated = TRUE)

# Convergence diagnostics
mcmc_frequency_chains <- coda::mcmc.list(lapply(plyr::alply(results$genotype_freq_store_chains, no_mcmc_chains), coda::mcmc))
coda::gelman.diag(mcmc_frequency_chains, transform = TRUE, multivariate = FALSE) # Transformed (using log or logit) to improve the normality of the distribution


