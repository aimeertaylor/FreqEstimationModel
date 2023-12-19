################################################################################
# Script to explore the difference between estimating allele frequences by
# summing over haplotype frequency estimates versus fitting the model to data on
# each locus separately.
#
# Beware: this code is suboptimal (lots of copy and paste)
#
# Take away: slightly greater uncertainty when summing over haplotypes
################################################################################


rm(list=ls(all=TRUE)) # Remove all variables from workspace

library(FreqEstimationModel) # Load the haplotype frequency estimation package
library("plyr")
library("coda")
library("abind")

# ==============================================================================
# Select data to run; returns data_summary
# ==============================================================================
no_markers <- 3
data_seed <- 1
load(sprintf('../data/Data_no_markers%s_seed%s.RData', no_markers, data_seed)) 


# ==============================================================================
# MCMC and model parameters
# ==============================================================================
# MCMC variables
thinning_interval <- 1 # Change this variable to increase the number of iterations per chain that are not saved in memory
no_traces_preburnin <- 10000 # For more traces but manageable pdf plots, don't exceed 10k and increase thinning interval instead.
no_mcmc_chains <- 3 # Number of MCMC chains to run
burnin <- 1:(0.5*no_traces_preburnin)  
mcmc_variable_list <- list(no_mcmc_chains = no_mcmc_chains,
                           no_traces_preburnin = no_traces_preburnin,
                           thinning_interval = thinning_interval,
                           NGS = FALSE,
                           log_like_zero = FALSE)

# MOI variables
moi_prior <- 'Poisson' # Choose between 'Uniform', 'Poisson', 'Geometric' or 'nBinomial'
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


# ==============================================================================
# Process data: three SNPs together versus separately: 
# ==============================================================================
just_data <- data_summary
just_data$Data <- data_summary$Data[,1:no_markers]
processed_data_list_all <- preprocess_data(data_summary = just_data,
                                       log_like_zero = FALSE,
                                       NGS = FALSE,
                                       augment_missing_data = TRUE,
                                       moi_prior_min2)


just_data <- data_summary
just_data$Data <- data_summary$Data[, 1, drop = FALSE]
processed_data_list_1 <- preprocess_data(data_summary = just_data,
                                           log_like_zero = FALSE,
                                           NGS = FALSE,
                                           augment_missing_data = TRUE,
                                           moi_prior_min2)

just_data <- data_summary
just_data$Data <- data_summary$Data[, 2, drop = FALSE]
processed_data_list_2 <- preprocess_data(data_summary = just_data,
                                         log_like_zero = FALSE,
                                         NGS = FALSE,
                                         augment_missing_data = TRUE,
                                         moi_prior_min2)

just_data <- data_summary
just_data$Data <- data_summary$Data[, 3, drop = FALSE]
processed_data_list_3 <- preprocess_data(data_summary = just_data,
                                         log_like_zero = FALSE,
                                         NGS = FALSE,
                                         augment_missing_data = TRUE,
                                         moi_prior_min2)



# ==============================================================================
# Run MCMC
# ==============================================================================
set.seed(1) # Set seed for reproducibility
frequency_hyperparameter <-rep(1,processed_data_list_all$no_haplotypes) # The Parameter vector for the Dirichlet prior on the frequency estimate
frequency_initial <- NULL # If null, external_frequency_initial set internally, otherwise set to input external_frequency_initial
frequency_list <- list(frequency_hyperparameter = frequency_hyperparameter, frequency_initial = frequency_initial)
results_all <- mcmc_sampling_parallel(processed_data_list_all,
                                  moi_list,
                                  frequency_list,
                                  mcmc_variable_list,
                                  cores_max = 3)



set.seed(1) # Set seed for reproducibility
frequency_hyperparameter <-rep(1,processed_data_list_1$no_haplotypes) # The Parameter vector for the Dirichlet prior on the frequency estimate
frequency_initial <- NULL # If null, external_frequency_initial set internally, otherwise set to input external_frequency_initial
frequency_list <- list(frequency_hyperparameter = frequency_hyperparameter, frequency_initial = frequency_initial)
results_1 <- mcmc_sampling_parallel(processed_data_list_1,
                                      moi_list,
                                      frequency_list,
                                      mcmc_variable_list,
                                      cores_max = 3)

set.seed(1) # Set seed for reproducibility
frequency_hyperparameter <-rep(1,processed_data_list_2$no_haplotypes) # The Parameter vector for the Dirichlet prior on the frequency estimate
frequency_initial <- NULL # If null, external_frequency_initial set internally, otherwise set to input external_frequency_initial
frequency_list <- list(frequency_hyperparameter = frequency_hyperparameter, frequency_initial = frequency_initial)
results_2 <- mcmc_sampling_parallel(processed_data_list_2,
                                      moi_list,
                                      frequency_list,
                                      mcmc_variable_list,
                                      cores_max = 3)

set.seed(1) # Set seed for reproducibility
frequency_hyperparameter <-rep(1,processed_data_list_3$no_haplotypes) # The Parameter vector for the Dirichlet prior on the frequency estimate
frequency_initial <- NULL # If null, external_frequency_initial set internally, otherwise set to input external_frequency_initial
frequency_list <- list(frequency_hyperparameter = frequency_hyperparameter, frequency_initial = frequency_initial)
results_3 <- mcmc_sampling_parallel(processed_data_list_3,
                                      moi_list,
                                      frequency_list,
                                      mcmc_variable_list,
                                      cores_max = 3)



# ==============================================================================
# Extract allele frequencies 
# ==============================================================================

# All loci together: 
alply_genotype_freq_store_chains_burnin <- plyr::alply(results_all$genotype_freq_store_chains[-burnin,,],3)
mcmc_frequencies <- do.call(rbind, alply_genotype_freq_store_chains_burnin)
hap_freq <- t(apply(mcmc_frequencies, 2, quantile, prob = c(0.025, 0.5, 0.975)))
haplotypes_str <- colnames(mcmc_frequencies) 
haplotypes_chr <- do.call(rbind, strsplit(haplotypes_str, split = ""))
ind <- apply(haplotypes_chr, 1, function(x) x == "1") # Compute the frequency of the default allele
mcmc_allele_freq <- sapply(1:nrow(ind), function(i) rowSums(mcmc_frequencies[,ind[i,]]))
allele_freq_all <- t(apply(mcmc_allele_freq, 2, quantile, prob = c(0.025, 0.5, 0.975)))


# First locus:
alply_genotype_freq_store_chains_burnin <- plyr::alply(results_1$genotype_freq_store_chains[-burnin,,],3)
mcmc_frequency_chains <- coda::mcmc.list(lapply(alply_genotype_freq_store_chains_burnin, 
                                                coda::mcmc,
                                                start = (max(burnin)+1)*mcmc_variable_list$thinning_interval,
                                                end = mcmc_variable_list$no_traces_preburnin*mcmc_variable_list$thinning_interval,
                                                thin = mcmc_variable_list$thinning_interval))
allele_freq_1 = summary(mcmc_frequency_chains)$quantiles["1",c("2.5%", "50%", "97.5%")]


# Second locus: 
alply_genotype_freq_store_chains_burnin <- plyr::alply(results_2$genotype_freq_store_chains[-burnin,,],3)
mcmc_frequency_chains <- coda::mcmc.list(lapply(alply_genotype_freq_store_chains_burnin, 
                                                coda::mcmc,
                                                start = (max(burnin)+1)*mcmc_variable_list$thinning_interval,
                                                end = mcmc_variable_list$no_traces_preburnin*mcmc_variable_list$thinning_interval,
                                                thin = mcmc_variable_list$thinning_interval))
allele_freq_2 = summary(mcmc_frequency_chains)$quantiles["1",c("2.5%", "50%", "97.5%")]


# Third locus: 
alply_genotype_freq_store_chains_burnin <- plyr::alply(results_3$genotype_freq_store_chains[-burnin,,],3)
mcmc_frequency_chains <- coda::mcmc.list(lapply(alply_genotype_freq_store_chains_burnin, 
                                                coda::mcmc,
                                                start = (max(burnin)+1)*mcmc_variable_list$thinning_interval,
                                                end = mcmc_variable_list$no_traces_preburnin*mcmc_variable_list$thinning_interval,
                                                thin = mcmc_variable_list$thinning_interval))
allele_freq_3 = summary(mcmc_frequency_chains)$quantiles["1",c("2.5%", "50%", "97.5%")]


# ==============================================================================
# Results
# ==============================================================================
allele_freq_sep <- rbind(allele_freq_1, allele_freq_2, allele_freq_3)
allele_freq_all <- cbind("50%" = allele_freq_all, 
                         "2.5%" = allele_freq_all_LCI, 
                         "97.5%" = allele_freq_all_UCI)

plot(x = allele_freq_all[,"50%"], y = allele_freq_sep[,"50%"], bty = "n", 
     ylim = c(0,1), xlim = c(0,1), 
     xlab = "Computed from haplotype frequencies", 
     ylab = "Computed from SNP-wise data")
abline(a = 0, b = 1, lty = "dashed")
segments(x0 = allele_freq_all[,"50%"], x1 = allele_freq_all[,"50%"], 
         y0 = allele_freq_sep[,"2.5%"], y1 = allele_freq_sep[,"97.5%"])
segments(x0 = allele_freq_all[,"2.5%"], x1 = allele_freq_all[,"97.5%"], 
         y0 = allele_freq_sep[,"50%"], y1 = allele_freq_sep[,"50%"])
