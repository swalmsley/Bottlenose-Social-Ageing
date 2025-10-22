

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


fit <- tar_read(ma_m2)
c <- coef(fit, prob=c(0.05, 0.95))$Title

######
# years <- years[!years %in% c(4, 18)]
###### Specify and review this!

# take random draws instead of first 100? ######


###### NO Male-Male associations in 1993? Pretty weird...
###### This needs to be investigated


######   # results from just 3 ff in 2006, potentially all associated -- sample size does matter for global network properties
# check that not also true of nodal properties... I suppose only thing it doesn't matter for is edge weights themselves

###### confirm that we're applying inv_logit to samples prior to taking mean, and so on? also more feedback on this would be nice (check bisonR)


####### check inv_logit everywhere, also post on brms forum

###### make sure SD not stressing out models

###### review poisson priors

###### mention 1000 draws in text

###### took off SD for several models - put it back on

###### check edge_sex_associations - seems to be broken in Relationships!! 


library(ggplot2)
library(dplyr)




plot_within_slopes(tar_read(ma_m1))
plot_within_slopes(tar_read(ma_m2))

plot_within_slopes(tar_read(ma_m3))
plot_within_slopes(tar_read(ma_m4))

plot_within_slopes(tar_read(ma_m5))
plot_within_slopes(tar_read(ma_m6))

plot_within_slopes(tar_read(ma_m7))
plot_within_slopes(tar_read(ma_m8))

plot_within_slopes(tar_read(ma_m9))
plot_within_slopes(tar_read(ma_m10))







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
