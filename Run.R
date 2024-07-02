

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



# Extras ------------------------------------------------------------------

# # Example contrasts
# c <- as_draws_df(tar_read(m1))$b_meanMinAge - as_draws_df(tar_read(m2))$b_meanMinAge
# c <- as_draws_df(tar_read(m3))$b_meanMinAge - as_draws_df(tar_read(m4))$b_meanMinAge
# c <- as_draws_df(tar_read(m5))$b_meanMinAge - as_draws_df(tar_read(m6))$b_meanMinAge
# c <- as_draws_df(tar_read(m7))$b_meanMinAge - as_draws_df(tar_read(m8))$b_meanMinAge
# c <- as_draws_df(tar_read(m9))$b_meanMinAge - as_draws_df(tar_read(m10))$b_meanMinAge
# PI(c)
# 
# # Example R2 calculations
# library(performance)
# r2_bayes(tar_read(m1))
# r2_bayes(tar_read(m2))
# r2_bayes(tar_read(m3))
# r2_bayes(tar_read(m4))
# r2_bayes(tar_read(m5))
# r2_bayes(tar_read(m6))
# r2_bayes(tar_read(m7))
# r2_bayes(tar_read(m8))
# r2_bayes(tar_read(m9))
# r2_bayes(tar_read(m10))
# # marginal is fixed only, conditional is both


# Examine network model fits

# mods <- tar_read(brms_fit) 
# mods[[1]]
# mods[[2]]
# mods[[3]]
# mods[[4]] 
# mods[[5]]
# mods[[6]] 
# mods[[7]]
# mods[[8]]
# mods[[9]]
# mods[[10]]
# mods[[11]]
# mods[[12]]
# mods[[13]]
# mods[[14]]
# mods[[15]]
# mods[[16]]
# mods[[17]]
# mods[[18]]
# mods[[19]]
# mods[[20]]
# mods[[21]]
# mods[[22]]
# mods[[23]]
# mods[[24]] 
# mods[[25]]
# mods[[26]]


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
