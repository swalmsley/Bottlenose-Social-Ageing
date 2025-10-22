

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

