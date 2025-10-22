
library(targets)

# Source
tar_source('R') # will do all in 'R' folder

# Seed
tar_option_set(seed = 1234)
options(future.globals.maxSize = 1024 * 1024 * 1024)  # Increase limit to 1 GB

# Variables
suppressMessages(set_cmdstan_path(path='C:/Users/sjfwa/AppData/Local/R/cmdstan-2.33.1')) # cmdstan path

# Targets
list(

  ###### What's group_df and why it's getting used where it is
  
  #######################################################################
  ## Read in data
  
  # Photographic and environmental data
  tar_target(raw_photoData, read_data_excel('./input/LV_SS_Primary_1988-2024.xlsx')),
  
  # Calculate Kappa statistic for classification of young animals
  tar_target(kappa, young_NBW_kappa(raw_photoData)),



  #######################################################################
  ## Define groups and social associations

  # process photo data
  tar_target(side_ungrouped, process_photoID(raw_photoData, 'Left', 'Gully')), # only left-sided photos from the Gully

  # Build nbw metadata
  tar_target(nbw_meta_raw, build_nbw_metadata(raw_photoData, 'left')), ######
  tar_target(nbw_meta, add_age_to_meta(nbw_meta_raw)), # Note, only including metadata for years after 2001

  # Defining groups
  tar_target(groups, carve_groups(side_ungrouped, nbw_meta)),
  tar_target(group_df_all, merge(groups, nbw_meta, by=c('Title','year'),all.x=TRUE)),
  tar_target(group_df_step1, group_df_all[!is.na(sex),,]), # cleave off individual-year combinations without sex information (including year < 2001)
  tar_target(group_df, group_df_step1[,numYears2:=length(unique(year)),by=Title]), # Only want individuals with at least 2 years in final dataframe to estimate change

  # Use tar_group to split dataframe by year
  tar_group_by(group_df_all_grouped, group_df_all, year),
  tar_target(year_key, year_key(group_df_all_grouped)), # track year-group IDs

  # Calculate social associations within each year-specific dataframe, then combine the results back together
  tar_target(group_associations, build_group_associations(group_df_all_grouped), pattern=map(group_df_all_grouped), iteration="list"),
  tar_target(group_associations_all, unlist_group_associations(group_associations)),
  tar_target(group_associations_all_aggregated, aggregate_group_associations(group_associations_all)), # aggregated version just useful for summaries



  #######################################################################
  ## Run social network models (estimate edge weights)

  # Run edge weight models using multi-annual network structure
  tar_target(brms_fit_group_multiAnnual, brm(
    formula = together ~ 1 + (1|dyad_annual),
    data = group_associations_all,
    family = bernoulli(),
    prior = c(
      set_prior("normal(-1.5, 1)", class = "Intercept"),
      set_prior("normal(0, 1)", class = "sd", group = "dyad_annual")), ####### took off coef=intercept here
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.9))),
  
  
  # Edge diagnostics # NOTE need to replace with edge_list_ma function which handles transformation to response scale better (see Bottlenose-SocialCausation project)
  # tar_target(edges_group_ma, edge_list_multiAnnual(brms_fit_group_multiAnnual, include_zeros=TRUE)),
  # tar_target(diagnostic_figure, save_figure('./Manuscript/Figures/Diagnostic.png',w=12,h=6,
  #                                               edge_diagnostics(edges_group_ma, group_associations_all_aggregated, multi=TRUE))),

  
  #######################################################################
  ## Prepare dataframes for modelling

  # Estimate nodal traits from annual edge weights
  tar_target(nodal_traits_multiAnnual, extract_trait_multiAnnual(brms_fit_group_multiAnnual)),

  # Combine estimates of annual nodal traits into dataframe
  tar_target(combined_traits_ma, combine_trait_dfs_ma(nodal_traits_multiAnnual, year_key)),

  # Finalize dataframe for modelling
  tar_target(model_df_ma_step1, merge(combined_traits_ma, nbw_meta, by=c('year','Title'))),
  tar_target(model_df_ma_step2, model_df_ma_step1[,numYears2:=.N,by=Title]), # Only want individuals with at least 2 years in final dataframe to estimate change
  tar_target(model_df_ma, add_degree_ma(side_ungrouped, group_associations_all, model_df_ma_step2, year_key)), ######

  # Extract biopsy details for models with genetically-sexed individuals only
  tar_target(biopsy_titles, biopsySexTitles(raw_photoData)),


  #######################################################################
  ## Model social traits as a function of minimum age

  ### Group size ###############################

  
  tar_target(ma_m1_minAge, brm(
    formula = gs ~  minimumAge + (1|Title) + (1|year),
    data = group_df[ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Male',,],
    family = 'Poisson',
    chains=4,
    cores=4,
    prior = c(set_prior('normal(1, 0.5)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='minimumAge'),
              set_prior('exponential(1)', class='sd')))),
  tar_target(ma_m1, brm(
    formula = gs ~  meanMinAge + deltaMinAge + (deltaMinAge|Title) +  (1|year),
    data = group_df[ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Male',,],
    family = 'Poisson',
    chains=4,
    cores=4,
    prior = c(set_prior('normal(1, 0.5)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('exponential(1)', class='sd'),
              set_prior('lkj_corr_cholesky(2)', class='L')))),
  tar_target(ma_m1_sd, brm(
    formula = gs ~  meanMinAge + minimumAge + (1|Title) +  (1|year),
    data = group_df[ ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Male',,],
    family = 'Poisson',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.95),
    prior = c(set_prior('normal(1, 0.5)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='minimumAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('exponential(1)', class='sd')))),
  tar_target(ma_m1_biopsy, brm(
    formula = gs ~  meanMinAge + deltaMinAge + (deltaMinAge|Title) +  (1|year),
    data = group_df[Title %in% biopsy_titles & ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Male',,],
    family = 'Poisson',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(0.5, 0.5)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('exponential(1)', class='sd'),
              set_prior('lkj_corr_cholesky(2)', class='L')))),
  tar_target(ma_m1_yearSpan, brm(
    formula = gs ~  meanMinAge + deltaMinAge + (deltaMinAge|Title) +  (1|year),
    data = group_df[ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Male' & yearSpan>10,,],
    family = 'Poisson',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(1, 0.5)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('exponential(1)', class='sd'),
              set_prior('lkj_corr_cholesky(2)', class='L')))),


  tar_target(ma_m2_minAge, brm(
    formula = gs ~  minimumAge + (1|Title) +  (1|year),
    data = group_df[ageClass=='Adult' & numYears2>1 & year>2000 & sex=='FemaleJ',,],
    family = 'Poisson',
    chains=4,
    cores=4,
    prior = c(set_prior('normal(1, 0.5)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='minimumAge'),
              set_prior('exponential(1)', class='sd')))),
  tar_target(ma_m2, brm(
    formula = gs ~  meanMinAge + deltaMinAge + (deltaMinAge|Title) +  (1|year),
    data = group_df[ageClass=='Adult' & numYears2>1 & year>2000 & sex=='FemaleJ',,],
    family = 'Poisson',
    chains=4,
    cores=4,
    prior = c(set_prior('normal(1, 0.5)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('exponential(1)', class='sd'),
              set_prior('lkj_corr_cholesky(2)', class='L')))),
  tar_target(ma_m2_sd, brm(
    formula = gs ~  meanMinAge + minimumAge + (1|Title) +  (1|year),
    data = group_df[ ageClass=='Adult' & numYears2>1 & year>2000 & sex=='FemaleJ',,],
    family = 'Poisson',
    chains=4,
    cores=4,
    prior = c(set_prior('normal(1, 0.5)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='minimumAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('exponential(1)', class='sd')))),
  tar_target(ma_m2_biopsy, brm(
    formula = gs ~  meanMinAge + deltaMinAge + (deltaMinAge|Title) +  (1|year),
    data = group_df[Title %in% biopsy_titles & ageClass=='Adult' & numYears2>1 & year>2000 & sex=='FemaleJ',,],
    family = 'Poisson',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(0.5, 0.5)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('exponential(1)', class='sd'),
              set_prior('lkj_corr_cholesky(2)', class='L')))),
  tar_target(ma_m2_yearSpan, brm(
    formula = gs ~  meanMinAge + deltaMinAge + (deltaMinAge|Title) +  (1|year),
    data = group_df[ageClass=='Adult' & numYears2>1 & year>2000 & sex=='FemaleJ' & yearSpan>10,,],
    family = 'Poisson',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(0.5, 0.5)', class='Intercept'),
              set_prior('normal(0, 0.5)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 0.5)', class='b', coef='meanMinAge'),
              set_prior('exponential(1)', class='sd'),
              set_prior('lkj_corr_cholesky(2)', class='L')))),


  ### Degree ###################################

  
  tar_target(ma_m3_minAge, brm(
    formula = degree ~ numSamplingPeriodsByYear + minimumAge + (1|Title) +  (1|year),
    data = model_df_ma[ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Male',,],
    family = 'Poisson',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(1.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='minimumAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd')))),
  tar_target(ma_m3, brm(
    formula = degree ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) +  (1|year),
    data = model_df_ma[ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Male',,],
    family = 'Poisson',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(1.5, 1)', class='Intercept'),
            set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
            set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
            set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
            set_prior('exponential(1)', class='sd'),
            set_prior('lkj_corr_cholesky(2)', class='L')))),
  tar_target(ma_m3_sd, brm(
    formula = degree ~ numSamplingPeriodsByYear + meanMinAge + minimumAge + (1|Title) +  (1|year),
    data = model_df_ma[ ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Male',,],
    family = 'Poisson',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(1.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='minimumAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd')))),
  tar_target(ma_m3_biopsy, brm(
    formula = degree ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) +  (1|year),
    data = model_df_ma[Title %in% biopsy_titles & ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Male',,],
    family = 'Poisson',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(1.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd'),
              set_prior('lkj_corr_cholesky(2)', class='L')))),
  tar_target(ma_m3_yearSpan, brm(
    formula = degree ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) +  (1|year),
    data = model_df_ma[ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Male' & yearSpan>10,,],
    family = 'Poisson',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(1.5, 1)', class='Intercept'),
              set_prior('normal(0, 0.5)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 0.5)', class='b', coef='meanMinAge'),
              set_prior('normal(0, 0.5)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd'),
              set_prior('lkj_corr_cholesky(2)', class='L')))),

  
  tar_target(ma_m4_minAge, brm(
    formula = degree ~ numSamplingPeriodsByYear + minimumAge + (1|Title) +  (1|year),
    data = model_df_ma[ageClass=='Adult' & numYears2>1 & year>2000 & sex=='FemaleJ',,],
    family = 'Poisson',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(1.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='minimumAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd')))),
  tar_target(ma_m4, brm(
    formula = degree ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) +  (1|year),
    data = model_df_ma[ageClass=='Adult' & numYears2>1 & year>2000 & sex=='FemaleJ',,],
    family = 'Poisson',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(1.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd'),
              set_prior('lkj_corr_cholesky(2)', class='L')))),
  tar_target(ma_m4_sd, brm(
    formula = degree ~ numSamplingPeriodsByYear + meanMinAge + minimumAge + (1|Title) +  (1|year),
    data = model_df_ma[ ageClass=='Adult' & numYears2>1 & year>2000 & sex=='FemaleJ',,],
    family = 'Poisson',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(1.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='minimumAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd')))),
  tar_target(ma_m4_biopsy, brm(
    formula = degree ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) +  (1|year),
    data = model_df_ma[Title %in% biopsy_titles & ageClass=='Adult' & numYears2>1 & year>2000 & sex=='FemaleJ',,],
    family = 'Poisson',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(1.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd'),
              set_prior('lkj_corr_cholesky(2)', class='L')))),  
  tar_target(ma_m4_yearSpan, brm(
    formula = degree ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) +  (1|year),
    data =  model_df_ma[ageClass=='Adult' & numYears2>1 & year>2000 & sex=='FemaleJ' & yearSpan>10,,],
    family = 'Poisson',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(1.5, 1)', class='Intercept'),
              set_prior('normal(0, 0.5)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 0.5)', class='b', coef='meanMinAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd'),
              set_prior('lkj_corr_cholesky(2)', class='L')))),


  ### Centrality ###############################
  
  # GAUSSIAN CHECK
  tar_target(tt_ma_m5_minAge, brm(
    formula = eigen | mi(eigen_SD) ~ numSamplingPeriodsByYear + minimumAge + (1|Title) +  (1|year),
    data = model_df_ma[ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Male',,],
    family = 'Gaussian',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(-0.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='minimumAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd')))),
  tar_target(tt_ma_m5, brm(
    formula = eigen | mi(eigen_SD) ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) +  (1|year),
    data = model_df_ma[ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Male',,],
    family = 'Gaussian',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(-0.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd'),
              set_prior('lkj_corr_cholesky(2)', class='L')))),

  tar_target(ma_m5_minAge, brm(
    formula = eigen | mi(eigen_SD) ~ numSamplingPeriodsByYear + minimumAge + (1|Title) +  (1|year),
    data = model_df_ma[ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Male',,],
    family = 'Beta',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(-0.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='minimumAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd'),
              set_prior('gamma(0.01, 0.01)', class='phi')))),
  tar_target(ma_m5, brm(
    formula = eigen | mi(eigen_SD) ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) +  (1|year),
    data = model_df_ma[ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Male',,],
    family = 'Beta',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(-0.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd'),
              set_prior('lkj_corr_cholesky(2)', class='L')))),
  tar_target(ma_m5_sd, brm(
    formula = eigen | mi(eigen_SD) ~ numSamplingPeriodsByYear + meanMinAge + minimumAge + (1|Title) +  (1|year),
    data = model_df_ma[ ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Male',,],
    family = 'Beta',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(-0.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='minimumAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd'),
              set_prior('gamma(0.01, 0.01)', class='phi')))),
  tar_target(ma_m5_biopsy, brm(
    formula = eigen | mi(eigen_SD) ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) +  (1|year),
    data = model_df_ma[Title %in% biopsy_titles & ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Male',,],
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
              set_prior('lkj_corr_cholesky(2)', class='L')))),  
  tar_target(ma_m5_yearSpan, brm(
    formula = eigen | mi(eigen_SD) ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) +  (1|year),
    data = model_df_ma[ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Male' & yearSpan>10,,],
    family = 'Beta',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(-0.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd'),
              # set_prior('gamma(0.01, 0.01)', class='phi'),
              set_prior('lkj_corr_cholesky(2)', class='L')))),

  
  tar_target(ma_m6_minAge, brm(
    formula = eigen | mi(eigen_SD) ~ numSamplingPeriodsByYear + minimumAge + (1|Title) +  (1|year),
    data = model_df_ma[ageClass=='Adult' & numYears2>1 & year>2000 & sex=='FemaleJ',,],
    family = 'Beta',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(-0.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='minimumAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd'),
              set_prior('gamma(0.01, 0.01)', class='phi')))),
  tar_target(ma_m6, brm(
    formula = eigen | mi(eigen_SD) ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) +  (1|year),
    data = model_df_ma[ageClass=='Adult' & numYears2>1 & year>2000 & sex=='FemaleJ',,],
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
              set_prior('lkj_corr_cholesky(2)', class='L')))),
  tar_target(ma_m6_sd, brm(
    formula = eigen | mi(eigen_SD) ~ numSamplingPeriodsByYear + meanMinAge + minimumAge + (1|Title) +  (1|year),
    data = model_df_ma[ ageClass=='Adult' & numYears2>1 & year>2000 & sex=='FemaleJ',,],
    family = 'Beta',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(-0.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='minimumAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd'),
              set_prior('gamma(0.01, 0.01)', class='phi')))),
  tar_target(ma_m6_biopsy, brm(
    formula = eigen | mi(eigen_SD) ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) +  (1|year),
    data = model_df_ma[Title %in% biopsy_titles & ageClass=='Adult' & numYears2>1 & year>2000 & sex=='FemaleJ',,],
    family = 'Beta',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(-0.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd'),
              # set_prior('gamma(0.01, 0.01)', class='phi'),
              set_prior('lkj_corr_cholesky(2)', class='L')))),
  tar_target(ma_m6_yearSpan, brm(
    formula = eigen | mi(eigen_SD) ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) +  (1|year),
    data = model_df_ma[ageClass=='Adult' & numYears2>1 & year>2000 & sex=='FemaleJ' & yearSpan>10,,],
    family = 'Beta',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(-0.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd'),
              # set_prior('gamma(0.01, 0.01)', class='phi'),
              set_prior('lkj_corr_cholesky(2)', class='L')))),


  ### Strength #################################

  
  tar_target(ma_m7_minAge, brm(
    formula = strength | mi(strength_SD) ~ numSamplingPeriodsByYear +  minimumAge + (1|Title) +  (1|year),
    data = model_df_ma[ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Male',,],
    family = 'Gamma',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(0, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='minimumAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd')))),
  tar_target(ma_m7, brm(
    formula = strength | mi(strength_SD) ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) +  (1|year),
    data = model_df_ma[ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Male',,],
    family = 'Gamma',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(0, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd'),
              set_prior('lkj_corr_cholesky(2)', class='L')))),
  tar_target(ma_m7_sd, brm(
    formula = strength | mi(strength_SD) ~ numSamplingPeriodsByYear + meanMinAge + minimumAge + (1|Title) +  (1|year),
    data = model_df_ma[ ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Male',,],
    family = 'Gamma',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(0, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='minimumAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd')))),
  tar_target(ma_m7_biopsy, brm(
    formula = strength | mi(strength_SD) ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) +  (1|year),
    data = model_df_ma[Title %in% biopsy_titles & ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Male',,],
    family = 'Gamma',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(0, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd'),
              set_prior('lkj_corr_cholesky(2)', class='L')))),
  tar_target(ma_m7_yearSpan, brm(
    formula = strength | mi(strength_SD) ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) +  (1|year),
    data = model_df_ma[ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Male' & yearSpan>10,,],
    family = 'Gamma',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(0, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd'),
              set_prior('lkj_corr_cholesky(2)', class='L')))),

  tar_target(ma_m8_minAge, brm(
    formula = strength | mi(strength_SD) ~ numSamplingPeriodsByYear + minimumAge + (1|Title) +  (1|year),
    data = model_df_ma[ageClass=='Adult' & numYears2>1 & year>2000 & sex=='FemaleJ',,],
    family = 'Gamma',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(0, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='minimumAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd')))),
  tar_target(ma_m8, brm(
    formula = strength | mi(strength_SD) ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) +  (1|year),
    data = model_df_ma[ageClass=='Adult' & numYears2>1 & year>2000 & sex=='FemaleJ',,],
    family = 'Gamma',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(0, 0.5)', class='Intercept'),
              set_prior('normal(0, 0.5)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 0.5)', class='b', coef='meanMinAge'),
              set_prior('normal(0, 0.5)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd'),
              set_prior('lkj_corr_cholesky(2)', class='L')))),
  tar_target(ma_m8_sd, brm(
    formula = strength | mi(strength_SD) ~ numSamplingPeriodsByYear + meanMinAge + minimumAge + (1|Title) +  (1|year),
    data = model_df_ma[ ageClass=='Adult' & numYears2>1 & year>2000 & sex=='FemaleJ',,],
    family = 'Gamma',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(0, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='minimumAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd')))),
  tar_target(ma_m8_biopsy, brm(
    formula = strength | mi(strength_SD) ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) +  (1|year),
    data = model_df_ma[Title %in% biopsy_titles & ageClass=='Adult' & numYears2>1 & year>2000 & sex=='FemaleJ',,],
    family = 'Gamma',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.999),
    prior = c(set_prior('normal(0, 0.5)', class='Intercept'),
              set_prior('normal(0, 0.5)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 0.5)', class='b', coef='meanMinAge'),
              set_prior('normal(0, 0.5)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd'),
              set_prior('lkj_corr_cholesky(2)', class='L')))),
  tar_target(ma_m8_yearSpan, brm(
    formula = strength | mi(strength_SD) ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|Title) +  (1|year),
    data = model_df_ma[ageClass=='Adult' & numYears2>1 & year>2000 & sex=='FemaleJ' & yearSpan>10,,],
    family = 'Gamma',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(0, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd'),
              set_prior('lkj_corr_cholesky(2)', class='L')))),


  ### Edge #####################################
 
  tar_target(ma_m9_minAge, brm(
    formula = medianEdge | mi(medianEdge_SD) ~  minimumAge + (1|Title) +  (1|year),
    data = model_df_ma[ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Male',,],
    family = 'Beta',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(0, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='minimumAge'),
              set_prior('exponential(1)', class='sd')))),
  tar_target(ma_m9, brm(
    formula = medianEdge | mi(medianEdge_SD) ~  meanMinAge + deltaMinAge + (deltaMinAge|Title) +  (1|year),
    data = model_df_ma[ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Male',,],
    family = 'Beta',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(-0.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('exponential(1)', class='sd'),
              set_prior('lkj_corr_cholesky(2)', class='L')))),
  tar_target(ma_m9_sd, brm(
    formula = medianEdge | mi(medianEdge_SD) ~  meanMinAge + minimumAge + (1|Title) +  (1|year),
    data = model_df_ma[ ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Male',,],
    family = 'Beta',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(-0.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='minimumAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('exponential(1)', class='sd'),
              set_prior('gamma(0.01, 0.01)', class='phi')))),
  tar_target(ma_m9_biopsy, brm(
    formula = medianEdge ~  meanMinAge + deltaMinAge + (deltaMinAge|Title) +  (1|year),
    data = model_df_ma[Title %in% biopsy_titles & ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Male',,],
    family = 'Beta',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(-0.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('exponential(1)', class='sd'),
              # set_prior('gamma(0.01, 0.01)', class='phi'),
              set_prior('lkj_corr_cholesky(2)', class='L')))),
  tar_target(ma_m9_yearSpan, brm(
    formula = medianEdge | mi(medianEdge_SD) ~  meanMinAge + deltaMinAge + (deltaMinAge|Title) +  (1|year),
    data = model_df_ma[ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Male' & yearSpan>10,,],
    family = 'Beta',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(-0.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('exponential(1)', class='sd'),
              # set_prior('gamma(0.01, 0.01)', class='phi'),
              set_prior('lkj_corr_cholesky(2)', class='L')))),

  

  tar_target(ma_m10_minAge, brm(
    formula = medianEdge | mi(medianEdge_SD) ~  minimumAge + (1|Title) +  (1|year),
    data = model_df_ma[ageClass=='Adult' & numYears2>1 & year>2000 & sex=='FemaleJ',,],
    family = 'Beta',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(0, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='minimumAge'),
              set_prior('exponential(1)', class='sd')))),    
  tar_target(ma_m10, brm(
    formula = medianEdge | mi(medianEdge_SD) ~  meanMinAge + deltaMinAge + (deltaMinAge|Title) +  (1|year),
    data = model_df_ma[ageClass=='Adult' & numYears2>1 & year>2000 & sex=='FemaleJ',,],
    family = 'Beta',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(-0.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('exponential(1)', class='sd'),
              set_prior('lkj_corr_cholesky(2)', class='L')))),
  tar_target(ma_m10_sd, brm(
    formula = medianEdge | mi(medianEdge_SD) ~  meanMinAge + minimumAge + (1|Title) +  (1|year),
    data = model_df_ma[ ageClass=='Adult' & numYears2>1 & year>2000 & sex=='FemaleJ',,],
    family = 'Beta',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(-0.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='minimumAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('exponential(1)', class='sd')))),
  tar_target(ma_m10_biopsy, brm(
    formula = medianEdge | mi(medianEdge_SD) ~  meanMinAge + deltaMinAge + (deltaMinAge|Title) +  (1|year),
    data = model_df_ma[Title %in% biopsy_titles & ageClass=='Adult' & numYears2>1 & year>2000 & sex=='FemaleJ',,],
    family = 'Beta',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(-0.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('exponential(1)', class='sd'),
              set_prior('lkj_corr_cholesky(2)', class='L')))),
  tar_target(ma_m10_yearSpan, brm(
    formula = medianEdge | mi(medianEdge_SD) ~  meanMinAge + deltaMinAge + (deltaMinAge|Title) +  (1|year),
    data = model_df_ma[ageClass=='Adult' & numYears2>1 & year>2000 & sex=='FemaleJ' & yearSpan>10,,],
    family = 'Beta',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(-0.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='deltaMinAge'),
              set_prior('normal(0, 1)', class='b', coef='meanMinAge'),
              set_prior('exponential(1)', class='sd'),
              set_prior('lkj_corr_cholesky(2)', class='L')))),


  
  #######################################################################
  ## Modelling non-linear trends in social behaviour with minimum age

  tar_target(m12_smooth2, brm(
    formula = gs ~ s(minimumAge, k=3, by=sex) + (1|Title) + (1|year),
    data = group_df[sex %in% c('Male','FemaleJ') & ageClass=='Adult' & year>2000,,],
    family = 'Poisson',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.95),
    prior = c(set_prior('normal(1, 1)', class='Intercept'),
              set_prior('exponential(1)', class='sd')))),
  
  tar_target(m34_smooth2, brm(
    formula = degree ~ numSamplingPeriodsByYear + s(minimumAge, k=3, by=sex) + (1|Title) + (1|year),
    data = model_df_ma[sex %in% c('Male','FemaleJ') &  ageClass=='Adult' & year>2000,,],
    family = 'Poisson',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.95),
    prior = c(set_prior('normal(0.5, 1)', class='Intercept'),
            set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
            set_prior('exponential(1)', class='sd')))),

  tar_target(m56_smooth2, brm(
    formula = eigen | mi(eigen_SD) ~ numSamplingPeriodsByYear + s(minimumAge, k=3, by=sex) + (1|Title) +  (1|year),
    data = model_df_ma[sex %in% c('Male','FemaleJ') &  ageClass=='Adult' & year>2000,,],
    family = 'Beta',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.95),
    prior = c(set_prior('normal(-0.5, 1)', class='Intercept'),
            set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
            set_prior('exponential(1)', class='sd')))),

  tar_target(m78_smooth2, brm(
    formula = strength | mi(strength_SD) ~ numSamplingPeriodsByYear + s(minimumAge, k=3, by=sex) + (1|Title) + (1|year),
    data = model_df_ma[sex %in% c('Male','FemaleJ') &  ageClass=='Adult' & year>2000,,],
    family = 'Gamma',
    chains=4,
    cores=4,
    control = list(adapt_delta = 0.99),
    prior = c(set_prior('normal(1.5, 1)', class='Intercept'),
              set_prior('normal(0, 1)', class='b', coef='numSamplingPeriodsByYear'),
              set_prior('exponential(1)', class='sd')))),

  tar_target(m910_smooth2, brm(
    formula = medianEdge | mi(medianEdge_SD) ~ s(minimumAge, k=3, by=sex) + (1|Title) + (1|year),
    data = model_df_ma[sex %in% c('Male','FemaleJ') &  ageClass=='Adult' & year>2000,,],
    family = 'Beta',
    chains=4,
    cores=4,
    prior = c(set_prior('normal(-1.5, 1)', class='Intercept'),
              set_prior('exponential(1)', class = 'sds'),
              set_prior('normal(0, 1)', class='b'),
              set_prior('exponential(1)', class='sd')))),

  

  #######################################################################
  ## Multi-response analyses

  
  # Simple total correlation between relationship strength and number partners
  tar_target(m_corr, brm(bf(scale(degree) ~ 1) + bf(scale(medianEdge) ~ 1) + set_rescor(TRUE),
    data   = model_df_ma,
    family = gaussian(),
    chains = 4, cores = 4, iter = 2000,
    prior = c(prior(lkj(2), class = "rescor")))),
  
  # Collapse group size into annual measure
  tar_target(group_df_with_means, group_df[,gs_annual:=mean(gs), by=c('Title','year')]),

  # Fit multi-response models
  tar_target(m_multi_m_ma, brm(formula=
                                 bf(gs_annual ~ meanMinAge + deltaMinAge + (deltaMinAge|p|Title) + (1|year), family='Gamma')+
                                 bf(medianEdge ~ meanMinAge + deltaMinAge + (deltaMinAge|p|Title) + (1|year), family='Beta')+
                                 bf(strength ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|p|Title) + (1|year), family='Gamma')+
                                 bf(eigen ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|p|Title) + (1|year), family='Beta')+
                                 bf(degree ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|p|Title) + (1|year), family='Poisson'),
                               data=multi_dataset_ma[ ageClass=='Adult' & numYears2>1 & year>2000 & sex=='Male',,],
                               iter=4000, chains=4, cores=4)),

  tar_target(m_multi_f_ma, brm(formula=
                                 bf(gs_annual ~ meanMinAge + deltaMinAge + (deltaMinAge|p|Title) + (1|year), family='Gamma')+
                                 bf(medianEdge ~ meanMinAge + deltaMinAge + (deltaMinAge|p|Title) + (1|year), family='Beta')+
                                 bf(strength ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|p|Title) + (1|year), family='Gamma')+
                                 bf(eigen ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|p|Title) + (1|year), family='Beta')+
                                 bf(degree ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|p|Title) + (1|year), family='Poisson'),
                               data=multi_dataset_ma[ ageClass=='Adult' & numYears2>1 & year>2000 & sex=='FemaleJ',,],
                               iter=4000, chains=4, cores=4)),

  tar_target(m_multi_both_ma, brm(formula=
                                    bf(gs_annual ~ meanMinAge + deltaMinAge + (deltaMinAge|p|Title) + (1|year), family='Gamma')+
                                    bf(medianEdge ~ meanMinAge + deltaMinAge + (deltaMinAge|p|Title) + (1|year), family='Beta')+
                                    bf(strength ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|p|Title) + (1|year), family='Gamma')+
                                    bf(eigen ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|p|Title) + (1|year), family='Beta')+
                                    bf(degree ~ numSamplingPeriodsByYear + meanMinAge + deltaMinAge + (deltaMinAge|p|Title) + (1|year), family='Poisson'),
                                  data=multi_dataset_ma[ ageClass=='Adult' & numYears2>1 & year>2000 & sex %in% c('Male','FemaleJ'),,],
                                  iter=4000, chains=4, cores=4)),
  
  tar_target(m_multi_both_ma_minimum, brm(formula=
                                    bf(gs_annual ~ minimumAge + (minimumAge|p|Title) + (1|year), family='Gamma')+
                                    bf(medianEdge ~ minimumAge + (minimumAge|p|Title) + (1|year), family='Beta')+
                                    bf(strength ~ numSamplingPeriodsByYear + minimumAge + (minimumAge|p|Title) + (1|year), family='Gamma')+
                                    bf(eigen ~ numSamplingPeriodsByYear + minimumAge + (minimumAge|p|Title) + (1|year), family='Beta')+
                                    bf(degree ~ numSamplingPeriodsByYear + minimumAge + (minimumAge|p|Title) + (1|year), family='Poisson'),
                                  data=multi_dataset_ma[ ageClass=='Adult' & numYears2>1 & year>2000 & sex %in% c('Male','FemaleJ'),,],
                                  iter=4000, chains=4, cores=4)),
  

  # Extract within-individual effects
  tar_target(within_effects_f_ma, extract_within_effects_sex_specific(multi_dataset_ma, m_multi_f_ma)),
  tar_target(within_effects_m_ma, extract_within_effects_sex_specific(multi_dataset_ma, m_multi_m_ma)),
  
  tar_target(within_effects_ma, extract_within_effects(multi_dataset_ma, m_multi_both_ma)),


  tar_target(multi_dataset_ma, merge(model_df_ma, unique(group_df_with_means[,c('Title','year','gs_annual'),]), all.x=TRUE)),


  # Version with minimum age 
  tar_target(within_effects_ma_minimum, extract_within_effects_minimum(multi_dataset_ma, m_multi_both_ma_minimum)),
  


  #######################################################################
  ## Figures

  tar_target(Figure1, save_figure('./Manuscript/Figures/Figure1.png',w=12,h=3.5,
                                      (smooth_plot(m12_smooth2, 'gs', 'Group size', FALSE, TRUE) |
                                         smooth_plot(m34_smooth2, 'degree', 'Number of social partners', FALSE, FALSE) |
                                         smooth_plot(m56_smooth2, 'eigen', 'Centrality', FALSE, FALSE) |
                                         smooth_plot(m78_smooth2, 'strength', 'Social strength', FALSE, FALSE) |
                                         smooth_plot(m910_smooth2, 'medianEdge', 'Mean relationship strength', FALSE, FALSE)))),

  tar_target(Figure2, save_figure('./Manuscript/Figures/Figure2.png',w=15,h=5,
                                  (within_between_plot(ma_m1, 'gs', 'Group size','Male') |
                                     within_between_plot(ma_m2, 'gs', 'Group size','Female') |
                                     (sex_coef_plot(ma_m1, ma_m2)) / var_plot(ma_m1, ma_m2)) +
                                    plot_annotation(tag_levels = 'A',
                                                    title = 'Group size',
                                                    theme = theme(plot.title = element_text(size = 18))))),

  tar_target(Figure3, save_figure('./Manuscript/Figures/Figure3.png',w=15,h=5,
                                  (within_between_plot(ma_m3, 'degree', 'Number of partners','Male') |
                                     within_between_plot(ma_m4, 'degree', 'Number of partners','Female') |
                                     (sex_coef_plot(ma_m3, ma_m4)) / var_plot(ma_m3, ma_m4)) +
                                    plot_annotation(tag_levels = 'A',
                                                    title = 'Number of social partners',
                                                    theme = theme(plot.title = element_text(size = 18))))),

  tar_target(Figure4, save_figure('./Manuscript/Figures/Figure4.png',w=15,h=5,
                                  (within_between_plot(ma_m7, 'strength', 'Strength','Male') |
                                     within_between_plot(ma_m8, 'strength', 'Strength','Female') |
                                     (sex_coef_plot(ma_m7, ma_m8)) / var_plot(ma_m7, ma_m8)) +
                                    plot_annotation(tag_levels = 'A',
                                                    title = 'Social network strength',
                                                    theme = theme(plot.title = element_text(size = 18))))),

  tar_target(Figure5, save_figure('./Manuscript/Figures/Figure5.png',w=15,h=5,
                                  (within_between_plot(ma_m5, 'eigen', 'Centrality','Male') |
                                     within_between_plot(ma_m6, 'eigen', 'Centrality','Female') |
                                     (sex_coef_plot(ma_m5, ma_m6)) / var_plot(ma_m5, ma_m6)) +
                                    plot_annotation(tag_levels = 'A',
                                                    title = 'Eigenvectory centrality',
                                                    theme = theme(plot.title = element_text(size = 18))))),

  tar_target(Figure6, save_figure('./Manuscript/Figures/Figure6.png',w=15,h=5,
                                  (within_between_plot(ma_m9, 'medianEdge', 'Mean relationship strength','Male') |
                                     within_between_plot(ma_m10,  'medianEdge', 'Mean bond strength','Female') |
                                     (sex_coef_plot(ma_m9, ma_m10)) / var_plot(ma_m9, ma_m10)) +
                                    plot_annotation(tag_levels = 'A',
                                                    title = 'Mean relationship strength',
                                                    theme = theme(plot.title = element_text(size = 18))))),


  #######################################################################
  ## Supplemental figures


  tar_target(FigureS1, save_figure('./Manuscript/Figures/FigureS1.png',w=11,h=8,
                                  plot_within_effects_sex_specific(within_effects_m_ma, within_effects_f_ma))),


  tar_target(FigureS2, save_figure('./Manuscript/Figures/FigureS2.png',w=11,h=8,
                                  plot_slope_correlations(m_multi_m_ma, m_multi_f_ma))),


  tar_target(FigureS3, save_figure('./Manuscript/Figures/FigureS3.png',w=12,h=6,
                                   (age_histogram(group_df, 'Male') | age_histogram(group_df, 'FemaleJ')))),



  #######################################################################
  ## Additional figures (unused)

  # Correlations of within-individual effects for deltaMinAge from model with both sexes combined
  tar_target(FigureX1, save_figure('./Manuscript/Figures/FigureX1.png',w=11,h=8,
                                   plot_within_effects(within_effects_ma))),
  
  # Correlations of within-individual effects for minimum age from model with both sexes combined
  tar_target(FigureX2, save_figure('./Manuscript/Figures/FigureX2.png',w=11,h=8,
                                   plot_within_effects(within_effects_ma_minimum))),

  # Smooth effects showing raw data
  tar_target(FigureX3, save_figure('./Manuscript/Figures/FigureX3.png', w=12, h=4,
                                   (smooth_plot(m12_smooth2, 'gs', 'Group size', TRUE, TRUE) |
                                      smooth_plot(m34_smooth2, 'degree', 'Number of social partners', TRUE, FALSE) |
                                      smooth_plot(m56_smooth2, 'eigen', 'Centrality', TRUE, FALSE) |
                                      smooth_plot(m78_smooth2, 'strength', 'Social strength', TRUE, FALSE) |
                                      smooth_plot(m910_smooth2, 'medianEdge', 'Mean relationship strength', TRUE, FALSE)))),

  # Plots of individual-specific slopes
  tar_target(FigureY1, save_figure('./Manuscript/Figures/Individual-slopes/slopes-m1.png', w=8, h=12, plot_within_slopes(ma_m1))),
  tar_target(FigureY2, save_figure('./Manuscript/Figures/Individual-slopes/slopes-m2.png', w=8, h=12, plot_within_slopes(ma_m2))),
  tar_target(FigureY3, save_figure('./Manuscript/Figures/Individual-slopes/slopes-m3.png', w=8, h=12, plot_within_slopes(ma_m3))),
  tar_target(FigureY4, save_figure('./Manuscript/Figures/Individual-slopes/slopes-m4.png', w=8, h=12, plot_within_slopes(ma_m4))),
  tar_target(FigureY5, save_figure('./Manuscript/Figures/Individual-slopes/slopes-m5.png', w=8, h=12, plot_within_slopes(ma_m5))),
  tar_target(FigureY6, save_figure('./Manuscript/Figures/Individual-slopes/slopes-m6.png', w=8, h=12, plot_within_slopes(ma_m6))),
  tar_target(FigureY7, save_figure('./Manuscript/Figures/Individual-slopes/slopes-m7.png', w=8, h=12, plot_within_slopes(ma_m7))),
  tar_target(FigureY8, save_figure('./Manuscript/Figures/Individual-slopes/slopes-m8.png', w=8, h=12, plot_within_slopes(ma_m8))),
  tar_target(FigureY9, save_figure('./Manuscript/Figures/Individual-slopes/slopes-m9.png', w=8, h=12, plot_within_slopes(ma_m9))),
  tar_target(FigureY10, save_figure('./Manuscript/Figures/Individual-slopes/slopes-m10.png', w=8, h=12, plot_within_slopes(ma_m10))),



  #######################################################################
  ## Write supplement

  tar_quarto(
    supplement,
    file.path('Manuscript','Supplement_SocialAgeing.qmd')),


  #######################################################################
  ## Write manuscript

  tar_quarto(
    paper_iScience,
    file.path('Manuscript','MS_SocialAgeing_iScience.qmd'))
  
  
)


















