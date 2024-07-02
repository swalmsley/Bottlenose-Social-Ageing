
library(targets)

# Source
tar_source('R') # will do all in 'R' folder

# Seed
tar_option_set(seed = 1234)

# Variables
suppressMessages(set_cmdstan_path(path='C:/Users/sjfwa/AppData/Local/R/cmdstan-2.33.1')) # cmdstan path

# Targets
list(

  
  #######################################################################
  ## Read in data
  tar_target(raw_photoData, read_data_excel('./input/LV_SS_Master_1988-2023_excel.xlsx')), ###### make a file watching target
  tar_target(raw_envData, read_data('./input/Scotian_Shelf_1988-2019_edited_HW_environment.csv')),
  
  # calculate Kappa statistic for classification of young animals
  tar_target(kappa, young_NBW_kappa()),
  
  
  #######################################################################
  ## Calculate associations and run edge weight models
  
  # process photo data
  tar_target(side_ungrouped, process_photoID(raw_photoData, 'Left', 'Gully')), # only left-sided photos from the Gully

  # carve up some groups based on timing
  tar_target(groups, carve_groups(side_ungrouped, nbw_meta)),
  tar_target(group_df_all, merge(groups, nbw_meta, by=c('Title','year'),all.x=TRUE)),
  tar_target(group_df, group_df_all[!is.na(sex),,]),
  
  
  
  ### Association edge weights
  
  # group by year for annual networks
  tar_group_by(side, side_ungrouped, year),
  
  # prepare binary association data for annual edge weight models
  tar_target(associations, build_binary_association_data(side), pattern=map(side), iteration="list"), ######

  # run edge weight models for each year
  tar_target(brms_fit, brm(
    formula = together ~ (1|dyad),
    data = associations,
    family = bernoulli(),
    prior = c(
     set_prior("normal(-1.5, 1.5)", class = "Intercept"),
     set_prior("normal(0, 1)", class = "sd", group = "dyad", coef = "Intercept")),
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.95)),
    pattern=map(associations), iteration='list'),
  
  
  #######################################################################
  ## Prepare dataframes for modelling
  
  tar_target(nodal_traits, extract_trait(brms_fit),
             pattern=map(brms_fit), iteration='list'),
  
  
  # combine trait dataframes from each year
  tar_target(combined_traits, combine_trait_dfs(nodal_traits, side)),  # assume order is preserved for now and check based on side, then make it more robust secondarily
  
  
  # build nbw metadata
  tar_target(nbw_meta_raw, build_nbw_metadata(raw_photoData, 'left')), ###### side-specific meta?
  tar_target(nbw_meta, add_age_to_meta(nbw_meta_raw)),
  
  
  # merge traits and metadata for analysis dataframe
  tar_target(model_df, merge(combined_traits, nbw_meta, by=c('year','Title'))),
  
  
  # extract biopsy details
  tar_target(biopsy_titles, biopsySexTitles(raw_photoData)),
  
  # final changes to dataframe
  tar_target(model_df_temp_1, model_df[,eigen:=eigen-0.001,]),
  tar_target(model_df_temp_2, model_df_temp_1[!(year_group %in% c(1,4,6,7,12,14,16,17,18,24)),,]), # Cut out years with very few observations 
  tar_target(model_df_temp_3, model_df_temp_2[year>=2000,,]),
  tar_target(model_df_temp, model_df_temp_3[,numYears2:=.N,by=Title]), # calculated based on span of analysis
  tar_target(degree_df_all, add_degree(side, associations, model_df_temp)),
  tar_target(degree_df, degree_df_all[!is.na(sex),,]),
  
  
  #######################################################################
  ## Trait models
  
  ## GROUP SIZE ##
  
  tar_target(m1, brm(
    formula = gs ~ meanMinAge + deltaMinAge + (deltaMinAge|Title) + (1|year),
    data = group_df[ sex=='Male' & numYears>1 & ageClass=='Adult' & year>2000,,], 
    family = 'Poisson',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.95),
    prior = c(set_prior('normal(1, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('exponential(1)', class='sd'),
              set_prior('lkj_corr_cholesky(1)', class='L')))),
  tar_target(m2, brm(
    formula = gs ~  meanMinAge + deltaMinAge + (deltaMinAge|Title) + (1|year),
    data = group_df[ sex=='Female-Juvenile' & numYears>1 & ageClass=='Adult' & year>2000,,], 
    family = 'Poisson',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.95),
    prior = c(set_prior('normal(1, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('exponential(1)', class='sd'),
              set_prior('lkj_corr_cholesky(1)', class='L')))),
  tar_target(m2b, brm(
    formula = gs ~  meanMinAge + deltaMinAge + (deltaMinAge|Title) + (1|year),
    data = group_df[ sex=='Female-Juvenile' & numYears>1 & ageClass=='Adult' & year>2000 & yearSpan>10,,], 
    family = 'Poisson',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.95),
    prior = c(set_prior('normal(1, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('exponential(1)', class='sd'),
              set_prior('lkj_corr_cholesky(1)', class='L')))),
  tar_target(m1biopsy, brm(
    formula = gs ~ meanMinAge + deltaMinAge + (deltaMinAge|Title) + (1|year),
    data = group_df[Title %in% biopsy_titles & sex=='Male' & numYears>1 & ageClass=='Adult' & year>2000,,], 
    family = 'Poisson',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.95),
    prior = c(set_prior('normal(1, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('exponential(1)', class='sd'),
              set_prior('lkj_corr_cholesky(1)', class='L')))),
  tar_target(m2biopsy, brm(
    formula = gs ~  meanMinAge + deltaMinAge + (deltaMinAge|Title) + (1|year),
    data = group_df[Title %in% biopsy_titles & sex=='Female-Juvenile' & numYears>1 & ageClass=='Adult' & year>2000,,], 
    family = 'Poisson',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(0.5, 1)', class='Intercept'),
              set_prior('normal(0, 0.5)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 0.5)', class='b', coef='meanMinAge'),
              set_prior('exponential(1)', class='sd'),
              set_prior('lkj_corr_cholesky(1)', class='L')))),
  


  ## NUMBER OF SOCIAL PARTNERS ##

  tar_target(m3, brm(
    formula = degree ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) + (1|year),
    data = degree_df[ ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Male',,],
    family = 'Poisson',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.95),
    prior = c(set_prior('normal(1.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd'),
              set_prior('lkj_corr_cholesky(1)', class='L')))),
  tar_target(m4, brm(
    formula = degree ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) + (1|year),
    data = degree_df[ ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Female-Juvenile',,],
    family = 'Poisson',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(1.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd'),
              set_prior('lkj_corr_cholesky(1)', class='L')))),
  tar_target(m4b, brm(
    formula = degree ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) + (1|year),
    data = degree_df[ ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Female-Juvenile' & yearSpan>10,,],
    family = 'Poisson',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(1.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd'),
              set_prior('lkj_corr_cholesky(1)', class='L')))),
  tar_target(m3biopsy, brm(
    formula = degree ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) + (1|year),
    data = degree_df[Title %in% biopsy_titles & ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Male',,],
    family = 'Poisson',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.95),
    prior = c(set_prior('normal(1.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd'),
              set_prior('lkj_corr_cholesky(1)', class='L')))),
  tar_target(m4biopsy, brm(
    formula = degree ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) + (1|year),
    data = degree_df[Title %in% biopsy_titles & ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Female-Juvenile',,],
    family = 'Poisson',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(1.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd'),
              set_prior('lkj_corr_cholesky(1)', class='L')))),
  
  
  ## CENTRALITY ##

  tar_target(m5, brm(
    formula = eigen | mi(eigen_SD) ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) +  (1|year),
    data = degree_df[ ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Male',,],
    family = 'Beta',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.95),
    prior = c(set_prior('normal(-0.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd'),
              set_prior('gamma(0.01, 0.01)', class='phi'),
              set_prior('lkj_corr_cholesky(1)', class='L')))),
  tar_target(m6, brm(
    formula = eigen ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) + (1|year),
    data = degree_df[ ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Female-Juvenile',,],
    family = 'Beta',
    chains=4,
    cores=4,
    prior = c(set_prior('normal(-0.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd'),
              set_prior('gamma(0.01, 0.01)', class='phi'),
              set_prior('lkj_corr_cholesky(1)', class='L')))),
  tar_target(m6b, brm(
    formula = eigen ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) + (1|year),
    data = degree_df[ ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Female-Juvenile' & yearSpan>10,,],
    family = 'Beta',
    chains=4,
    cores=4,
    prior = c(set_prior('normal(-0.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd'),
              set_prior('gamma(0.01, 0.01)', class='phi'),
              set_prior('lkj_corr_cholesky(1)', class='L')))),
  tar_target(m5biopsy, brm(
    formula = eigen | mi(eigen_SD) ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) +  (1|year),
    data = degree_df[Title %in% biopsy_titles & ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Male',,],
    family = 'Beta',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(-0.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd'),
              set_prior('gamma(0.01, 0.01)', class='phi'),
              set_prior('lkj_corr_cholesky(1)', class='L')))),
  tar_target(m6biopsy, brm(
    formula = eigen ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) + (1|year),
    data = degree_df[Title %in% biopsy_titles & ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Female-Juvenile',,],
    family = 'Beta',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(-0.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd'),
              set_prior('gamma(0.01, 0.01)', class='phi'),
              set_prior('lkj_corr_cholesky(1)', class='L')))),
  

  
  ## STRENGTH ##
  
  tar_target(m7, brm(
    formula = strength | mi(strength_SD) ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) + (1|year),
    data = degree_df[ ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Male',,],
    family = 'Gamma',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(0, 1)', class='Intercept'),
            set_prior('normal(0, 1)', class='b'),
            set_prior('exponential(1)', class='sd')))),
  tar_target(m8, brm(
    formula = strength | mi(strength_SD) ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) + (1|year),
    data = degree_df[ ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Female-Juvenile',,],
    family = 'Gamma',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.95),
    prior = c(set_prior('normal(0, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b'),
              set_prior('exponential(1)', class='sd')))),
  tar_target(m8b, brm(
    formula = strength | mi(strength_SD) ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) + (1|year),
    data = degree_df[ ageClass=='Adult' & numYears>1 & year>2000 & sex=='Female-Juvenile' & yearSpan>10,,],
    family = 'Gamma',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(0, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b'),
              set_prior('exponential(1)', class='sd')))),
  tar_target(m7biopsy, brm(
    formula = strength | mi(strength_SD) ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) + (1|year),
    data = degree_df[Title %in% biopsy_titles & ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Male',,],
    family = 'Gamma',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(0, 1)', class='Intercept'),
            set_prior('normal(0, 1)', class='b'),
            set_prior('exponential(1)', class='sd')))),
  tar_target(m8biopsy, brm(
    formula = strength | mi(strength_SD) ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) + (1|year),
    data = degree_df[Title %in% biopsy_titles & ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Female-Juvenile',,],
    family = 'Gamma',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(1, 0.5)', class='Intercept'),
              set_prior('normal(0, 1)', class='b'),
              set_prior('exponential(1)', class='sd')))),
  
  
  ## MEAN RELATIONSHIP STRENGTH ## 
  
  tar_target(m9, brm(
    formula = medianEdge ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) + (1|year),
    data = degree_df[ ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Male',,],
    family = 'Beta',
    chains=4,
    cores=4,
    prior = c(set_prior('normal(-0.5, 1)', class='Intercept'),
            set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
            set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
            set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
            set_prior('exponential(1)', class='sd'),
            set_prior('gamma(0.01, 0.01)', class='phi'),
            set_prior('lkj_corr_cholesky(1)', class='L')))),
  tar_target(m10, brm(
    formula = medianEdge ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) + (1|year),
    data = degree_df[ ageClass=='Adult' & numYears2>1  & year>2000 & sex=='Female-Juvenile',,],
    family = 'Beta',
    chains=4,
    cores=4,
    prior = c(set_prior('normal(-0.5, 1)', class='Intercept'),
            set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
            set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
            set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
            set_prior('exponential(1)', class='sd'),
            set_prior('gamma(0.01, 0.01)', class='phi'),
            set_prior('lkj_corr_cholesky(1)', class='L')))),
  tar_target(m10b, brm(
    formula = medianEdge ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) + (1|year),
    data = degree_df[ ageClass=='Adult' & numYears2>1  & year>2000 & sex=='Female-Juvenile' & yearSpan>10,,],
    family = 'Beta',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(-0.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd'),
              set_prior('gamma(0.01, 0.01)', class='phi'),
              set_prior('lkj_corr_cholesky(1)', class='L')))),
  tar_target(m9biopsy, brm(
    formula = medianEdge ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) + (1|year),
    data = degree_df[Title %in% biopsy_titles & ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Male',,],
    family = 'Beta',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(-0.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd'),
              set_prior('gamma(0.01, 0.01)', class='phi'),
              set_prior('lkj_corr_cholesky(1)', class='L')))),
  tar_target(m10biopsy, brm(
    formula = medianEdge ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) + (1|year),
    data = degree_df[Title %in% biopsy_titles & ageClass=='Adult' & numYears2>1  & year>2000 & sex=='Female-Juvenile',,],
    family = 'Beta',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(-0.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd'),
              set_prior('gamma(0.01, 0.01)', class='phi'),
              set_prior('lkj_corr_cholesky(1)', class='L')))),
  
  
  
  #######################################################################
  ## Figures
  
  tar_target(Figure1, save_figure('./Manuscript/Figures/Figure1.png',w=15,h=5,
                                  (within_between_plot(m3, 'degree', 'Number of partners','Male') | 
                                     within_between_plot(m4, 'degree', 'Number of partners','Female') |
                                     sex_coef_plot(m3, m4)) + plot_annotation(tag_levels = 'A',
                                                                              title = 'Number of social partners',
                                                                              theme = theme(plot.title = element_text(size = 18))))),
  
  
  tar_target(Figure2, save_figure('./Manuscript/Figures/Figure2.png',w=15,h=5,
                                  (within_between_plot(m5, 'eigen', 'Centrality','Male') | 
                                     within_between_plot(m6, 'eigen', 'Centrality','Female') |
                                     sex_coef_plot(m5, m6)) + plot_annotation(tag_levels = 'A',
                                                                              title = 'Eigenvectory centrality',
                                                                              theme = theme(plot.title = element_text(size = 18))))), 
  
  
  tar_target(Figure3, save_figure('./Manuscript/Figures/Figure3.png',w=15,h=5,
                                  (within_between_plot(m9, 'medianEdge', 'Mean bond strength','Male') | 
                                     within_between_plot(m10,  'medianEdge', 'Mean bond strength','Female') |
                                     sex_coef_plot(m9, m10)) + plot_annotation(tag_levels = 'A',
                                                                               title = 'Mean bond strength',
                                                                               theme = theme(plot.title = element_text(size = 18))))), 
  
  
  tar_target(Figure4, save_figure('./Manuscript/Figures/Figure4.png',w=5,h=5,
                                  pca_plot(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, biopsy_titles))), 
  
  
  tar_target(FigureS1, save_figure('./Manuscript/Figures/FigureS1.png',w=15,h=5,
                                  (within_between_plot(m1, 'gs', 'Group size','Male') | 
                                     within_between_plot(m2, 'gs', 'Group size','Female') |
                                     sex_coef_plot(m1, m2)) + plot_annotation(tag_levels = 'A',
                                                                              title = 'Group size',
                                                                              theme = theme(plot.title = element_text(size = 18))))), 
  
  
  tar_target(FigureS2, save_figure('./Manuscript/Figures/FigureS2.png',w=15,h=5,
                                  (within_between_plot(m7, 'strength', 'Strength','Male') | 
                                     within_between_plot(m8, 'strength', 'Strength','Female') |
                                     sex_coef_plot(m7, m8)) + plot_annotation(tag_levels = 'A',
                                                                              title = 'Social network strength',
                                                                              theme = theme(plot.title = element_text(size = 18))))), 
  
  
  #######################################################################
  ## Write supplement
  
  tar_quarto(
    supplement,
    file.path('Manuscript','Supplement_SocialAgeing.qmd')),  
  
  #######################################################################
  ## Write manuscript

  # tar_quarto(
  #   paper,
  #   file.path('Manuscript','MS_SocialAgeing.qmd')),
  
  tar_quarto(
    paper_iScience,
    file.path('Manuscript','MS_SocialAgeing_iScience.qmd'))


)

