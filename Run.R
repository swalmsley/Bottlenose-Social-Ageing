

library(targets)

###############
# Run analysis
tar_make() 
###############


# Visualize analysis
tar_visnetwork() # visualizes analysis pipeline
tar_visnetwork(targets_only = TRUE, label = 'time') # visualizes analysis pipeline

# # Examine object from analysis
# tar_read(yearPlots, branches=1)

# Diagnostics 
View(tar_meta()) # useful tool for diagnostics
View(tar_meta(targets_only = TRUE)) # simplified



# Using the maximum response value as the number of trials.
# Make sure this isn't happening more broadly?
# Try for regular brms_fit?





# Paper insights

# With a sparse network and prior that expects high edges, dyads with fewer data will get dragged down less (will be closer to innappropriately high prior)
# Regular messing around with priors can be useful for navigating sparser vs. denser networks, zero-inflated prior is an extension of this
# Yeah I'm not convinced I want to encode the zero-inflation thing... seems too strict
# Concern with using a low prior: would need a LOT of data for a high edge weight to emerge
# Bashing SRI for being certain about dyads never seen together being 0s -- similarly something desirable for our model is that 
# dyads never seen together but with few opportunities may interact
# For a single-prior model, boosting low probabilities helps fit to sparser network

# Problem with low-centered prior is that you may not get adequate movement towards higher values when you need them -- this is something we can check of course
# Compare relative to expectations from O'Brien work, for example; are we happy maxing out at 30-40% association indices, or do we think higher is justified?
# Certainly one of the zero-inflated models we fit seems to be doing something kinda whacky; really high values, even higher than HWI...

# Note that their fixed scenario has some moderate curve AND the 0s are almost always expected to be 0s... not sure what we think of this
# Note that their x axis is much bigger though, so possible that you get a little uplift in the "never togethers" right around 0, seems like you do

# Next section about diagnostic stuff should be useful...

# Last sentence: No option is perfect: aim for incorporating previous knowledge, incorporating uncertainty, and testing for downstream pathologies





##########################
# Calculating predictions
##########################

# pred_data <- expand.grid(deltaMinAge = c(-15,0,15), 
#                          sex = 'Male',
#                          Title = NA,
#                          numSamplingPeriodsByYear=2.4,
#                          meanMinAge=15,
#                          year = NA,
#                          eigen_SD=0.1)
# pred <- data.table(epred_draws(tar_read(m5), pred_data, re_formula = NA)) # N.B. 'NULL' means it INCLUDES conditional effects!
# pred[,mean(.epred),by=deltaMinAge]
# 
# pred_data <- expand.grid(deltaMinAge = c(-15,0,15), 
#                          sex = 'Female-Juvenile',
#                          Title = NA,
#                          numSamplingPeriodsByYear=2.4,
#                          meanMinAge=15,
#                          year = NA,
#                          eigen_SD=0.1)
# pred <- data.table(epred_draws(tar_read(m10), pred_data, re_formula = NA)) # N.B. 'NULL' means it INCLUDES conditional effects!
# pred[,mean(.epred),by=deltaMinAge]

# Examples of individuals with opposite within-individual effects
# Female group size: 5024, 1757
