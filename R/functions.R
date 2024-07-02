
# Custom functions to be called in '_targets.R'


###########################################################################
# Basic functions
###########################################################################


# Read data ---------------------------------------------------------------
read_data <- function(path) {
  d <- data.table(read.csv(path))
  return(d)
}




# Read data (excel) -------------------------------------------------------
read_data_excel <- function(path) {
  d <- data.table(read_xlsx(path))
  return(d)
}




# Save plot ---------------------------------------------------------------
save_figure <- function(path, w, h, call) {
  png(path, width=w, height=h, units='in', res=1000)
  print(call)
  dev.off()
}




# Save plot (transparent) -------------------------------------------------
save_figure_transparent <- function(path, w, h, call) {
  png(path, width=w, height=h, units='in', res=1000, bg = 'transparent')
  print(call)
  dev.off()
}




# In-line effect ----------------------------------------------------------
inline <- function(fit, var, p) {
  
  est <- summary(fit,prob=p)$fixed[paste(var),'Estimate']
  low <- summary(fit,prob=p)$fixed[paste(var),'l-90% CI']
  high <- summary(fit,prob=p)$fixed[paste(var),'u-90% CI']
  
  output <- paste(paste(format_number(est), ', ', sep = ''), 'CI ', paste(format_number(low), format_number(high), sep = ' -- '), sep = '')
  return(output)
  
}
# fit <- tar_read(m3)
# var <- 'numSamplingPeriodsByYear'
# p = 0.9




# Helper for inline effect ------------------------------------------------
format_number <- function(num) {
  if (abs(num) < 0.01) {
    
    # Convert to scientific notation with custom formatting
    formatted_number <- format(num, scientific = TRUE, trim = TRUE, digits = 2)
    # Extract mantissa and exponent
    parts <- strsplit(formatted_number, "e", fixed = TRUE)[[1]]
    # Print in the desired format
    out <- capture.output(cat(sprintf("%.2f x 10^%+d\n", as.numeric(parts[1]), as.numeric(parts[2]))))
    return(out)
    
  } else {
    return(round(num, 2))
  }
}




###########################################################################
# Project-specific functions
###########################################################################


# Coefficients table (brms) -----------------------------------------------
brms_table <- function(fit) {
  
  fixed <- data.table(as.data.frame(summary(fit, prob=0.9)$fixed)) 
  Title <- data.table(as.data.frame(summary(fit, prob=0.9)$random$Title)) 
  year <- data.table(as.data.frame(summary(fit, prob=0.9)$random$year)) 
  raw.table <- rbindlist(list(fixed,Title,year))
  
  stats.table <- cbind(c(row.names(summary(fit)$fixed),
                         row.names(summary(fit)$random$Title),
                         row.names(summary(fit)$random$year)), raw.table)
  
  #stats.table <- cbind(row.names(summary(fit)$fixed), stats.table)
  stats.table[,c('Rhat','Bulk_ESS','Tail_ESS'):=NULL,]
  names(stats.table) <- c("Term", "Estimate", "SE", "PI-Lower", "PI-Upper")
  
  nice_table(stats.table)
  
}
#fit <- tar_read(m5)




# Extract within-individual effects from each model -----------------------
model_to_slope <- function(fit, s) {

  # extract within-individual slopes from model
  c <- coef(fit)
  slopes <- data.frame(c$Title[,,'deltaMinAge'])
  slopes$Title <- rownames(slopes)
  slopes <- data.table(slopes)

  # add in response variable
  formula_str <- as.character(formula(fit))
  response_var <- strsplit(formula_str, " ~ ")[[1]][1]
  slopes[,var:=response_var,]
  slopes[,sex:=s,]
  
  return(slopes) #_#

}




# Create PCA plot ---------------------------------------------------------
pca_plot <- function(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, biopsy_titles) {
  
  df <- rbindlist(list(model_to_slope(m1, 'Male'),
                       model_to_slope(m3, 'Male'),
                       model_to_slope(m5, 'Male'),
                       model_to_slope(m7, 'Male'),
                       model_to_slope(m9, 'Male'),
                       
                       model_to_slope(m2, 'Female'),
                       model_to_slope(m4, 'Female'),
                       model_to_slope(m6, 'Female'),
                       model_to_slope(m8, 'Female'),
                       model_to_slope(m10, 'Female')))
  
  df[var=='eigen | mi(eigen_SD)', var:='eigen',]
  df[var=='strength | mi(strength_SD)', var:='strength',]
  
  # convert to wide format
  wide_df <- df[,c('Estimate','var','Title'),] %>%
    spread(key = var, value = Estimate)
  rownames(wide_df) <- wide_df$Title
  wide_df <- wide_df[2:6]
  
  wide_df_no_na <- wide_df[complete.cases(wide_df),]
  pca_res <- prcomp(wide_df_no_na, scale. = TRUE)
  
  wide_df_no_na$Title <- rownames(wide_df_no_na)
  wide_with_sex <- merge(wide_df_no_na, unique(df[,c('Title','sex'),]))
  
  wide_with_sex$biopsy <- (wide_with_sex$Title %in% biopsy_titles)
  wide_with_sex$biopsy <- ifelse(wide_with_sex$biopsy==TRUE,'Yes','No')
  
  with_genes <- wide_with_sex[wide_with_sex$biopsy=='Yes',]
  
  g <- autoplot(pca_res, data=wide_with_sex, colour='sex', shape='biopsy', size=1.7) + 
    scale_color_manual(values=c('#FDAE79','#4D1A7C')) +
    scale_shape_manual(values=c(16,3)) +
    annotate('text',x=0.1,y=0.19,label='Females', color='#FDAE79', size=5)+
    annotate('text',x=-0.25,y=-0.09,label='Males', color='#4D1A7C', size=5)+
    guides(color='none')+
    labs(shape='Genetically sexed')+
    theme_classic()+
    theme(legend.position=c(0.2,0.85))

  return(g)
  
}
# m1 <- tar_read(m1)
# m2 <- tar_read(m2)
# m3 <- tar_read(m3)
# m4 <- tar_read(m4)
# m5 <- tar_read(m5)
# m6 <- tar_read(m6)
# m7 <- tar_read(m7)
# m8 <- tar_read(m8)
# m9 <- tar_read(m9)
# m10 <- tar_read(m10)
# biopsy_titles <- tar_read(biopsy_titles)




# Create PCA table ---------------------------------------------------------
pca_table <- function(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, biopsy_titles) {
  
  df <- rbindlist(list(model_to_slope(m1, 'Male'),
                       model_to_slope(m3, 'Male'),
                       model_to_slope(m5, 'Male'),
                       model_to_slope(m7, 'Male'),
                       model_to_slope(m9, 'Male'),
                       
                       model_to_slope(m2, 'Female'),
                       model_to_slope(m4, 'Female'),
                       model_to_slope(m6, 'Female'),
                       model_to_slope(m8, 'Female'),
                       model_to_slope(m10, 'Female')))
  
  df[var=='eigen | mi(eigen_SD)', var:='eigen',]
  df[var=='strength | mi(strength_SD)', var:='strength',]
  
  # convert to wide format
  wide_df <- df[,c('Estimate','var','Title'),] %>%
    spread(key = var, value = Estimate)
  rownames(wide_df) <- wide_df$Title
  wide_df <- wide_df[2:6]
  
  wide_df_no_na <- wide_df[complete.cases(wide_df),]
  pca_res <- prcomp(wide_df_no_na, scale. = TRUE)
  
  wide_df_no_na$Title <- rownames(wide_df_no_na)
  wide_with_sex <- merge(wide_df_no_na, unique(df[,c('Title','sex'),]))
  
  b <- biopsy_titles
  
  wide_with_sex$biopsy <- (wide_with_sex$Title %in% b)
  wide_with_sex$biopsy <- ifelse(wide_with_sex$biopsy==TRUE,'Yes','No')
  
  with_genes <- wide_with_sex[wide_with_sex$biopsy=='Yes',]
  
  # make table for supplement
  d <- data.table(Measure=rownames(pca_res$rotation),pca_res$rotation)
  d[Measure=='degree',Measure:='Change in number of social partners']
  d[Measure=='eigen',Measure:='Change in network centrality']
  d[Measure=='gs',Measure:='Change in group size']
  d[Measure=='medianEdge',Measure:='Change in mean bond strength']
  d[Measure=='strength',Measure:='Change in network strength']
  
  nice_table(d)
  
}
# m1 <- tar_read(m1)
# m2 <- tar_read(m2)
# m3 <- tar_read(m3)
# m4 <- tar_read(m4)
# m5 <- tar_read(m5)
# m6 <- tar_read(m6)
# m7 <- tar_read(m7)
# m8 <- tar_read(m8)
# m9 <- tar_read(m9)
# m10 <- tar_read(m10)
# biopsy_titles <- tar_read(biopsy_titles)




# Compare calf/juv classifications ----------------------------------------
young_NBW_kappa <- function() {
  
  # Prep claire version
  
  claire <- data.table(read_xlsx('input/Calf_juv_updated_master.xlsx'))
  
  claire[, Claire_ageClass:=ifelse(grepl('claire-juv', `Keyword-export`),'Juv','none'),]
  claire[Claire_ageClass!='Juv', Claire_ageClass:=ifelse(grepl('claire-calf', `Keyword-export`),'Calf','none'),]
  claire[,.N,by='Claire_ageClass']
  
  
  claire[,year:=year(`Date Original`),]
  claire <- unique(claire[,c('Title','Claire_ageClass','year')])
  claire$Title <- as.character(claire$Title)
  unique(claire[,.N,by=c('year','Title')]$N) # expect 1
  # there is a 2 here now
  
  # 1990 28, 1994 390, 2016 6240
  
  # error log for Claire version
  claire <- claire[!(year==1990 & Title=='28' & Claire_ageClass=='none'),,] # assuming Juv
  claire <- claire[!(year==1994 & Title=='390' & Claire_ageClass=='none'),,] # assuming Juv
  claire <- claire[!(year==2016 & Title=='6240' & Claire_ageClass=='none'),,] # assuming Juv
  
  unique(claire[,.N,by=c('year','Title')]$N) # expect 1 now that errors have been fixed/assumedunique(claire[,.N,by=c('year','Title')]$N) # expect 1
  
  
  # Prep raw version
  
  original <- tar_read(raw_photoData)
  
  original[, ageClass:=ifelse(grepl('samJuv', Keyword.export),'Juv','none'),]
  original[ageClass!='Juv', ageClass:=ifelse(grepl('samCalf', Keyword.export),'Calf','none'),]
  original[,.N,by='ageClass']
  
  original[,year:=year(`Date.Original`),]
  original <- unique(original[,c('Title','ageClass','year')])
  original[,nByYear:=length(unique(ageClass)),by=c('Title','year')] # for checking
  unique(original[,.N,by=c('year','Title')]$N) # expect 1
  
  # error log for next year
  original <- original[!(year==2015 & Title=='5059' & ageClass=='none'),,]
  original <- original[!(year==2015 & Title=='6172' & ageClass=='none'),,]
  original <- original[!(year==2021 & Title=='6527' & ageClass=='none'),,]
  
  original <- original[Title!='see crops' & Title!='unk',,]
  claire <- claire[Title!='see crops' & Title!='unk',,]
  
  unique(original[,.N,by=c('year','Title')]$N) # expect 1
  
  
  compare <- merge(original[Title!='unk' & Title!='see crops',c('Title','ageClass','year'),], claire[Title!='unk' & Title!='see crops',,], by=c('Title','year'))
  
  
  return(kappam.fleiss(compare[,3:4], detail=TRUE))
  
}




# Process photo-ID observations -------------------------------------------
process_photoID <- function(lv_data, chosen_side, chosen_canyons){
  
  # check that Gully is keyworded appropriately in 2021 data
  
  # extract side information
  lv_data[,side:=ifelse(grepl("Left Side", Keyword.export),"Left","Right"),] # extract side from keyword list
  lv_data[,numSides:=length(unique(side)),by=Title] # identify IDs that are 2-sided
  
  # fix up date information
  lv_data[,dateTime:=as.POSIXct(Date.Original),] # reformat dateTimes
  lv_data[,day:=format(dateTime, format = "%Y-%m-%d")]
  lv_data[,year:=year(day),] # create a column of unique dates for use later
  lv_data[,numDaysByYear:=length(unique(day)),by=year]
  
  # remove observations with multiple or unknown whales
  lv_data <- lv_data[Title!="see crops",,] 
  lv_data <- lv_data[Title!="unk",,] 
  
  # add extra columns to indicate number of observations
  lv_data[,obsByYear:=.N,by=year] # count number of observations by year
  lv_data[,numDaysByTitle:=length(unique(day)),by=Title] # count number of days each individual was detected
  
  # add columns 
  lv_data[, canyon:=NA,]
  lv_data[, canyon:=ifelse(grepl('Gully', Keyword.export),'Gully',NA),]
  lv_data[, canyon:=ifelse(grepl('Shortland', Keyword.export),'Shortland',canyon),]
  lv_data[, canyon:=ifelse(grepl('Haldimand', Keyword.export),'Haldimand',canyon),]

  # subset by location and side  
  output <- lv_data[canyon%in%c(chosen_canyons),,]
  output <- output[side==chosen_side,,]
  
  return(output)
  
}
# lv_data <- tar_read(raw_photoData)
# chosen_side <- 'Left'
# chosen_canyons <- 'Gully'




# Carve identifications into groups ---------------------------------------
carve_groups <- function(photoData, meta) {
  
  photoData <- merge(photoData, meta[,c('Title','year','ageClass'),], by=c('Title','year'), all.x=TRUE) ## need to think on whether this is side-specific and so on? maybe other issues...
  #adults_only <- photoData[ageClass=='Adult',,]
  
  # sort data table
  setorder(photoData, dateTime)
  
  # Define time threshold (e.g., 2 hours)
  threshold <- 10 * 60  # 10 minutes in seconds
  
  # Carve groups based on time threshold
  photoData$t_group <- cut(photoData$dateTime, breaks = seq(min(photoData$dateTime), max(photoData$dateTime) + threshold, by = threshold), labels = FALSE)
  photoData[,myGroups:=.GRP,by=t_group]
  
  length(unique(photoData$myGroups)) # check number of timegroups

  for (i in 1:length(unique(photoData$myGroups))-1) {
    
    last_ID <- photoData[myGroups==i, max(dateTime)]
    next_ID <- photoData[myGroups==(i+1), min(dateTime)]
    
    if (as.numeric(difftime(next_ID, last_ID, units = 'mins')) < 5) {
      photoData[myGroups==i, myGroups:=i+1] 
    }
  }
  
  length(unique(photoData$myGroups)) # check how many timegroups remain 
  
  photoData[,groupSize:=as.numeric(length(unique(Title))),by=myGroups]
  hist(photoData$groupSize)
  
  photoData[,adultTitle:=ifelse(ageClass=='Adult',Title,NA),]
  #photoData[,groupSizeAdult:=as.numeric(length(unique(na.omit(adultTitle)))),by=myGroups]
  
  # subset to one row per group per individual
  group_ind_data <- unique(photoData[,c('Title','year', 'groupSize','myGroups')])
  
  group_ind_data[,gs:=groupSize,] # group size for each individual
  
  group_df <- unique(group_ind_data[,c('Title','myGroups','year','gs')])
  
  return(group_df)

}
# photoData <- tar_read(side_ungrouped)
# meta <- tar_read(nbw_meta) 




# Build binary association data -------------------------------------------
build_binary_association_data <- function(photoData) {
  
  # Step 1 - build dyadic dataframe (row = day-dyad combinations where at least one individual in dyad was observed)
  
  photoData <- data.table(photoData)
  
  # Note: currently hard-coding days as sampling periods
  dyads <- data.table(t(combn(photoData[,unique(Title),], 2))) # generate all dyads
  colnames(dyads) <- c('A', 'B') # rename ID columns
  df <- data.table(A = as.character(), B = as.character(), day = as.character(), A_observed = as.integer(), B_observed = as.integer())   # set up empty data frame for filling 
  
  for (d in unique(photoData$day)) { # create a new row for each observation day and each possible dyad in the dataset
    temp <- dyads[, day:=d, ]
    temp[, A_observed:=ifelse(A %in% photoData[day==d,unique(Title),],1,0)]
    temp[, B_observed:=ifelse(B %in% photoData[day==d,unique(Title),],1,0)]
    temp <- temp[A_observed==1 | B_observed==1,,] # prune out rows (day-dyad combinations) where neither member of dyad was observed
    df <- rbindlist(list(df, temp))
  }
  
  df[,index:=.I,] # add index for row number
  df[,bothPresent:=ifelse((A_observed==1 & B_observed==1),1,0),]
  
  df[,dyad:=paste(pmin(A,B),pmax(A,B),sep='_'),by=c('A','B')] # add useful dyad names
  
  
  # Step 2 - calculate associations
  
  df[bothPresent==1,together:=binary_association(photoData[day==day,,], A, B, 2), by='index']
  df[is.na(together),together:=0,] # assuming that presence on same day but not seen together is equivalent to not being present on same day
  
  
}
# side <- tar_read(side)
# photoData <- side[side$year==2016,,]




# Extract network traits --------------------------------------------------
extract_trait <- function(fit) {
  
  edge_samples <- data.frame(coef(fit, summary=FALSE)[[1]])
  colnames(edge_samples) <- gsub("X(\\d+_\\d+).Intercept", "\\1", colnames(edge_samples))
  
  N_draws <- 100
  
  # number of nodes (probably more efficient way to do this)
  nodes <- unique(c(sub("^(\\d+)_.*", "\\1", colnames(edge_samples)), sub("^\\d+_(\\d+)$", "\\1", colnames(edge_samples))))
  N_nodes <- length(nodes)
  
  # create empty matrix to store metric samples
  metric_samples <- matrix(0, N_draws, N_nodes)
  colnames(metric_samples) <- nodes
  
  metric_samples2 <- matrix(0, N_draws, N_nodes) 
  colnames(metric_samples2) <- nodes
  
  metric_samples3 <- matrix(0, N_draws, N_nodes) 
  colnames(metric_samples3) <- nodes

  for (i in 1:N_draws) {
    
    # pull draw of complete network (i.e., of each edge weight)
    draw_of_network <- as.data.frame(t(edge_samples[i,])) # probably more efficient way to do this (i.e., full matrix...)
    if (nrow(draw_of_network)==1) (rownames(draw_of_network) <- colnames(edge_samples)) # for 1991 with 2 photos... can simply delete this year too maybe
    colnames(draw_of_network) <- 'weight_sample_raw'
    
    # format IDs out of dyad label
    draw_of_network$A <- sub("^(\\d+)_.*", "\\1", rownames(draw_of_network)) # pull out ID A
    draw_of_network$B <- sub("^\\d+_(\\d+)$", "\\1", rownames(draw_of_network)) # pull out ID B
    
    
    ###########################################################
    ### Create main network with all Bayesian uncertainties ###
    
    # create edgelist
    net <- igraph::graph_from_edgelist(as.matrix(draw_of_network[,c('A','B'),]), directed=FALSE)
    
    # add in edgeweights
    igraph::E(net)$weight <- inv_logit(draw_of_network[,c('weight_sample_raw'),]) # inv-logit here because it's a bernoulli model # confirm makes sense to transform prior to strength calculation

    ######################################################################
    ### Create network where any dyads never seen together are given 0 ###
    
    # Need to pull actual non-zero dyads from data
    fitData <- data.table(fit$data)
    fitData[,nTogether:=sum(together),by=dyad]
    nonZero_dyads <- unique(fitData[nTogether>0,c('dyad')])$dyad
    
    NZ_draw_of_network <- data.table(copy(draw_of_network))
    NZ_draw_of_network[,dyad:=rownames(draw_of_network)]
    NZ_draw_of_network[,weight:=inv_logit(weight_sample_raw),] # inv-logit transformation
    # now set dyads that were never together as 0s
    NZ_draw_of_network[!(dyad %in% nonZero_dyads), weight:=NA,]
    
    NZ_net <- igraph::graph_from_edgelist(as.matrix(NZ_draw_of_network[,c('A','B'),]), directed=FALSE)
    #plot(NZ_net)
    
    igraph::E(NZ_net)$weight <- NZ_draw_of_network[,weight,] # inv-logit here because it's a binary model # confirm makes sense to transform prior to strength calculation
    NZ_net <- delete.edges(NZ_net, E(NZ_net)[is.na(E(NZ_net)$weight)])
    
    # For centrality
    vals <- igraph::eigen_centrality(net)$vector
    metric_samples[i, names(vals)] <- vals
    
    # For strength
    strengths <- igraph::strength(net)
    metric_samples2[i, names(strengths)] <- strengths
    
    # For mean non-0 edge weight
    for (name in names(strengths)) {
      
      # extract edge weights using igraph and cut below 5% threshold
      focal_weights <- data.table(dyad=incident(net,v=name), weight=(incident(net,v=name))$weight)
      #focal_weights_nonZero <- focal_weights[focal_weights>=0.05] # using 5% for Bayesian weight right now - could use real data (i.e., trim those never seen together or consider alternate threshold)
      focal_weights_nonZero <- data.table(dyad=incident(NZ_net,v=name), weight=(incident(NZ_net,v=name))$weight)
      
      # add mean from each draw of network to output matrix
      (metric_samples3[i, name] <- mean(focal_weights_nonZero$weight)) # ideally would look at edges never seen together and see range of values in each year

    }
    
  }
      
  ## FOR STANDARDIZATION ONLY ##
  # metric_samples <- standardize_matrix(metric_samples)
  # metric_samples2 <- standardize_matrix(metric_samples2)
  # metric_samples3 <- standardize_matrix(metric_samples3)
  ## FOR STANDARDIZATION ONLY ##

  # create more final output dataframe
  out <- data.table(ID = colnames(metric_samples)) 
  
  out$eigen <- apply(metric_samples, 2, mean) 
  out$eigen_SD <- apply(metric_samples, 2, sd) 
  
  out$strength <- apply(metric_samples2, 2, mean)  
  out$strength_SD <- apply(metric_samples2, 2, sd) 
  
  out$medianEdge <- apply(metric_samples3, 2, mean) # technically mean of mean edges
  out$medianEdge_SD <- apply(metric_samples3, 2, sd) 
 
  return(out)
  
}
# fit <- tar_read(brms_fit)[[2]]




# Standardize all values in matrix across entire matrix -------------------
standardize_matrix <- function(x) {
  
  return(x-mean(x))/mean(x)
  
}




# Extract numbers ---------------------------------------------------------
extract_numbers <- function(string) as.numeric(unlist(strsplit(string, "_")))




# Combine trait dataframes ------------------------------------------------
combine_trait_dfs <- function(nd, side) {
  
  combined <- rbindlist(nd, idcol=TRUE)
  combined[,year_group:=.GRP, by=.id]
  
  year_key <- data.frame(year_group = unique(combined$year_group), year = unique(side$year))
  combined <- merge(combined, year_key)
  
  combined[,Title:=ID,]
  combined[,ID:=NULL,]
  
  return(combined)
  
}
# nd <- tar_read(nodal_traits)
# side <- tar_read(side)




# Plot showing within and between effects ---------------------------------
within_between_plot <- function(fit, traitName, traitLabel, sex) {
  
  # data
  data <- data.table(fit$data)
  
  # # New data to predict for
  newdata <- expand.grid(deltaMinAge = seq(-10,10,1), # based on range in data
                         meanMinAge = mean(data$meanMinAge), 
                         Title = unique(fit$data$Title),
                         numSamplingPeriodsByYear = mean(data$numSamplingPeriodsByYear), # note this will give warning if model does not control for numSamplingPeriodsByYear, no problem
                         year = NA)
  newdata[[paste(traitName, '_SD', sep='')]] <- 0.001 # assign basically no measurement error - doesn't really seem to make a difference? maybe to CI?

  pred <- data.table(epred_draws(fit, newdata, re_formula=NULL)) 
  pred_means <- pred[, .(mean_epred=mean(.epred)),by=.(Title, deltaMinAge)] # extract means for each combination predictors (if not interested in plotting whole ribbon!)

  # next step: delimit within effects to span that they actually existed
  data[,floorMinAge:=min(meanMinAge+deltaMinAge),by=Title]
  data[,ceilingMinAge:=max(meanMinAge+deltaMinAge),by=Title]
  averageAges <- unique(data[,c('Title','meanMinAge', 'floorMinAge','ceilingMinAge'),])
  pred_means <- merge(pred_means, averageAges)
  pred_means[,reconstructed_minAge:=meanMinAge+deltaMinAge,]
  pred_means[,range:=((reconstructed_minAge >= floorMinAge) & (reconstructed_minAge <= ceilingMinAge)),] 

  newdata_btw <- expand.grid(deltaMinAge = 0, # because we're interested in between effects here?
                             meanMinAge = seq(2,36,1), 
                             Title = unique(fit$data$Title),
                             numSamplingPeriodsByYear = mean(data$numSamplingPeriodsByYear),
                             year = NA)
  newdata_btw[[paste(traitName, '_SD', sep='')]] <- 0.001 # assign basically no measurement error - doesn't really seem to make a difference? maybe to CI?
  pred_btw <- data.table(epred_draws(fit, newdata_btw, re_formula=NULL)) # NA for between effect? Changes magnitude but not shape necessarily

  pred_btw_means <- pred_btw[, .(mean_epred=mean(.epred)),by=.(meanMinAge)]
  
  # colors
  begin <- ifelse(sex=='Male', 0.2, 0.625)
  end <- ifelse(sex=='Male', 0.5, 0.925)
  # label
  plot_label <- ifelse(sex=='Male','M','F')
  label_col <- ifelse(sex=='Male', '#762F7C','#FCA797')
  
  p <- ggplot(pred_means[range==TRUE,,], aes(x=deltaMinAge+meanMinAge, y=mean_epred, group=Title, color=Title)) +

    # within effects
    geom_line(linewidth=1.25, alpha=0.75) +
    scale_color_viridis(option='A',discrete=TRUE,begin=begin, end=end)+
    geom_jitter(inherit.aes = FALSE, data=data,size=2,aes(x=deltaMinAge+meanMinAge,y=get(traitName)),color='black',alpha=0.1) +
    # between effect
    geom_line(inherit.aes = FALSE, data=pred_btw_means, aes(x=meanMinAge, y=mean_epred),linewidth=1, linetype='dashed', color='black')+
    labs(x='Minimum Age', y=paste(traitLabel))+
    annotate('text',x=Inf, y=Inf, size=8,label=plot_label,color=label_col, hjust=4,vjust=2)+
    theme_classic() +
    theme(strip.text=element_text(size=12),
          axis.text=element_text(size=12),
          axis.title=element_text(size=12),
          legend.position = "none") 
  p
  
  return(p)
  
}
# fit <- tar_read(m5)
# traitName <- 'eigen'
# traitLabel <- 'Centrality (std.)'
# sex='Female'




# Sex coefficient plot ----------------------------------------------------
sex_coef_plot <- function(fitM, fitF) {
  
  m <- data.table(mcmc_intervals(fitM, pars=c('b_meanMinAge', 'b_deltaMinAge', 'sd_Title__deltaMinAge'), point_est = 'median', prob=0.9, prob_outer = 0.9)$data) 
  f <- data.table(mcmc_intervals(fitF, pars=c('b_meanMinAge', 'b_deltaMinAge', 'sd_Title__deltaMinAge'), point_est = 'median', prob=0.9, prob_outer = 0.9)$data) 
  
  pars_df <- rbindlist(list(m,f))
  
  pars_df$model <- rep(c('Model M', 'Model F'), each = nrow(pars_df)/2)
  pars_df[parameter=='b_meanMinAge',parameter:='Between',]
  pars_df[parameter=='b_deltaMinAge',parameter:='Within',]
  pars_df[parameter=='sd_Title__deltaMinAge',parameter:='Variability',]
  
  pars_df$parameter <- factor(pars_df$parameter,levels = c('Variability', 'Within','Between'))

  # Adding an ID for each parameter for plotting
  #pars_df$param_id <- rep(1:(nrow(pars_df)/2), 2)
  
  #cols = c(viridis(1, begin=0.1,end=0.2, option='A'),viridis(1, begin=0.8,end=0.9, option='A')) # for flexible viridis approach
  cols <- c('#762F7C','#FCA797')
  
  # Plotting
  p <- ggplot(pars_df, aes(y = parameter, xmin = ll, xmax = hh, x = m, color=model)) +
    geom_errorbarh(aes(height = 0), position = position_dodge(width = 0.25), linewidth=1.75, alpha=0.75) +
    geom_point(position = position_dodge(width = 0.25), size=3.25) +
    geom_vline(xintercept = 0, linetype='dashed', linewidth=0.75)+
    scale_color_manual(values=c(cols[2],cols[1])) + 
    labs(x = 'Posterior', y = '') +
    #xlim(-0.1,0.1)+
    theme_classic()+theme(legend.position='none')+
    theme(axis.text = element_text(size=15))
  
  return(p)
  
}
# fitM <- tar_read(m5)
# fitF <- tar_read(m6)




# Pull out biopsy sex data ------------------------------------------------
# extract list of titles with biopsy sex information
biopsySexTitles <- function(d) {
  
  d[grepl("biopsy-Female", Keyword.export), biopsySex:='Y',]
  d[grepl("biopsy-Male", Keyword.export), biopsySex:='Y',]
  titles <- unique(d[Title!='unk' & biopsySex=='Y',Title,])
  
  return(titles)
  
  
}
#d <- tar_read(raw_photoData)




add_degree <- function(side, associations, model_df) {
  
  # loop through years
  for (i in 1:max(model_df$year_group)) {
    
    side <- data.table(side) 
    
    association_data <- associations[[i]] # pull out association target for year in question
    y <- association_data[1,year(day),] # pull out year from date data # Note would need to be changed for other sampling periods of course
    together <- association_data[together==1,,] # 'together' only
    
    # for each individual in model_df for year of question
    temp <- model_df[year==y,,]
    for (focal in unique(temp$Title)) {
      
      partner_df <- together[A==focal | B==focal,,]
      # add degree to updated model df
      model_df[year==y & Title==focal,degree:=max(0,length(unique(c(partner_df$A,partner_df$B)))-1),] # minus one so not to count focal individual, max 0 so that if empty, reverts to 0 not -1

      numDays <- side[year==y & Title==focal,length(unique(day)),]
      
      model_df[year==y & Title==focal,numSamplingPeriodsByYear:=numDays,]
      
    }
    
  }
  
  model_df[, standardized_numSamplingPeriodsByYear:=lapply(.SD, function(x) as.vector(scale(numSamplingPeriodsByYear))), by=year]

  return(model_df)
  
}
# side <- tar_read(side)
# associations <- tar_read(associations)
# model_df <- tar_read(model_df_temp)




# Build binary dataframe --------------------------------------------------
build_binary_df <- function(photoData) {
  
  # Note: currently hard-coding days as sampling periods
  dyads <- data.table(t(combn(photoData[,unique(Title),], 2))) # generate all dyads
  colnames(dyads) <- c('A', 'B') # rename ID columns
  df <- data.table(A = as.character(), B = as.character(), day = as.character(), A_observed = as.integer(), B_observed = as.integer())   # set up empty data frame for filling 
  
  for (d in unique(photoData$day)) { # create a new row for each observation day and each possible dyad in the dataset
    temp <- dyads[, day:=d, ]
    temp[, A_observed:=ifelse(A %in% photoData[day==d,unique(Title),],1,0)]
    temp[, B_observed:=ifelse(B %in% photoData[day==d,unique(Title),],1,0)]
    temp <- temp[A_observed==1 | B_observed==1,,] # prune out rows (day-dyad combinations) where neither member of dyad was observed
    df <- rbindlist(list(df, temp))
  }
  
  df[,index:=.I,] # add index for row number
  df[,bothPresent:=ifelse((A_observed==1 & B_observed==1),1,0),]
  
  df[,dyad:=paste(pmin(A,B),pmax(A,B),sep='_'),by=c('A','B')] # add useful dyad names
  
  return(df)
  
}
# photoData <- tar_read(side)




# Compute binary associations ---------------------------------------------
binary_association <- function(photos_sampling_period, ID1, ID2, time_limit) {
  
  # Pre-filter the data
  ID1_obs <- photos_sampling_period[(Title %in% ID1),,] # observations of ID 1 within sampling period
  ID2_obs <- photos_sampling_period[(Title %in% ID2),,] # observations of ID 2 within sampling period
  
  # Proceed only if there are any observations for ID1 or ID2 in the given sampling period (otherwise return 0 right away)
  if (nrow(ID1_obs) == 0 || nrow(ID2_obs) == 0) return(0)
  
  # Calculate the time difference between all combinations of observations
  time_diff <- outer(ID1_obs$dateTime, ID2_obs$dateTime, FUN = function(x, y) abs(as.numeric(difftime(x, y, units = 'mins'))))
  
  # Check if any time difference is within the time_limit
  if (any(time_diff <= time_limit)) return(1)
  
  return(0)
  
}




# Build NBW metadata ------------------------------------------------------
# takes raw LV output as input ('dt')
build_nbw_metadata <- function(dt, chosen_side) {
  
  # format dates
  dt[,dateTime:=as.POSIXct(Date.Original),] 
  dt[,day:=format(dateTime, format = "%Y-%m-%d")]
  dt[,year:=year(day),]
  
  # remove photographs that were processed as crops
  dt <- dt[Title!="see crops",,]
  dt <- dt[Title!="unk",,]
  
  # add left and right columns
  # dt=dt %>% mutate(side = ifelse(grepl("left", Keyword.export), "left", NA))
  # dt=dt %>% mutate(side = ifelse(grepl("right", Keyword.export), "right", NA))
  # 
  dt[grepl("Left", Keyword.export), side := "left"]
  dt[grepl("Right", Keyword.export), side := "right"]
  
  # subset to chosen side
  dt <- dt[side==chosen_side,,]
  
  # add columns with sex details
  dt=dt%>%mutate(sex = ifelse(grepl("FemaleJ", Keyword.export), "Female-Juvenile", NA))
  dt=dt%>%mutate(sex = ifelse(grepl("Male", Keyword.export), "Male", sex))
  dt$sex <- factor(dt$sex,levels=c("Female-Juvenile","Male"))

  # Juv and Calf classifications (using SFW's ratings)
  dt[,ageClass:=ifelse(grepl('samCalf', Keyword.export),'Calf','Adult'),by=c('Title','year')] # really more like "not young"
  dt[ageClass!='Calf',ageClass:=ifelse(grepl('samJuv', Keyword.export),'Juvenile',ageClass),by=c('Title','year')]
  dt[,young:=ageClass %in% c('Calf','Juvenile'),by=c('Title','year')]
  
  # add minimum age by year
  dt[,minYear:=min(year),by=Title]
  dt[,maxYear:=max(year),by=Title]
  dt[,yearSpan:=1+(maxYear-minYear),by=Title]
  dt[,catalogueAge:=year-minYear,] # Note: considering first year as "0"
  
  ## photo-ID error fix
  dt[Title==6527 & year==2021,ageClass:='Juvenile',]
  
  dt[, FirstYearAgeClass := unique(ageClass[minYear == year]), by = Title]
  
  dt[FirstYearAgeClass=='Calf',minimumAge:=catalogueAge,]
  dt[FirstYearAgeClass=='Juvenile',minimumAge:=catalogueAge+1,]
  dt[FirstYearAgeClass=='Adult',minimumAge:=catalogueAge+3,]
  
  # add day of year
  dt[,yday:=yday(day),]
  
  # diagnostic tests
  test_that("single sex classification for each individual", {
    expect_true(unique(dt[,length(unique(sex)),by=Title]$V1)==1)
  })
  
  test_that("single minimum age for each individual in each year", {
    expect_true(unique(dt[,length(unique(minimumAge)),by=c('Title', 'year'),]$V1)==1)
  })
  
  # sample one observation per individual per year and clean up dataset
  subsampled <- dt[dt[ , .I[sample(.N,1)] , by = c('Title','year')]$V1]
  meta <- subsampled[,c('Title', 'side', 'year', 'sex', 'minYear', 'maxYear', 'yearSpan', 'minimumAge', 'ageClass', 'yday'), ]
  
  return(meta)
  
}
# dt <- tar_read(raw_photoData)
# chosen_side = 'left'




# Add age info to metadata ------------------------------------------------
add_age_to_meta <- function(dt) {
  

  dt <- dt[year>=2000,,]
  dt <- dt[ageClass=='Adult',,] 
  
  # add number of years in long-term dataset
  dt[,numYears:=length(unique(year)),by=Title]
  
  # add mean minimum age for each individual
  dt[,meanMinAge:=mean(unique(minimumAge)),by=Title]
  # add delta age for each individual-year
  dt[,deltaMinAge:=minimumAge-meanMinAge,by=c('Title','year')]
  
  # sample one observation per individual per year and clean up dataset
  subsampled <- dt[dt[ , .I[sample(.N,1)] , by = c('Title','year')]$V1]
  meta_age <- subsampled[,c('Title', 'side', 'year', 'sex', 'numYears', 'minYear', 'maxYear', 'yearSpan', 'minimumAge', 'ageClass', 'meanMinAge', 'deltaMinAge', 'yday'), ]
  
  return(meta_age)
  
}
# dt <- tar_read(nbw_meta_raw)




# Coefficients table (cmdstan) --------------------------------------------
coef_table <- function(mcmc_object) {
  
  # create table
  precis_output <- precis(mcmc_object)
  stats.table <- data.table(data.frame(precis_output))
  stats.table <- cbind(rownames(precis_output), stats.table)
  #stats.table[,c('Rhat','Bulk_ESS','Tail_ESS'):=NULL,]
  names(stats.table) <- c("Parameter", "Mean", "SD", "CI-5.5", "CI-94.5", "R-hat", "ESS")
  
  return(nice_table(stats.table))
  
}




# Prune data to when both animals are known -------------------------------
both_catalogue_alive_assocData <- function(assocData, photoData) {
  
  assocData[,A_minYear:=photoData[Title==A,min(year),],by=A] # min years
  assocData[,B_minYear:=photoData[Title==B,min(year),],by=B]
  
  assocData[,A_maxYear:=photoData[Title==A,max(year),],by=A] # max years
  assocData[,B_maxYear:=photoData[Title==B,max(year),],by=B]
  
  assocData[,year:=year(day),] # add year column
  
  assocData[,bothAlive:=ifelse(all(year >= unique(A_minYear),
                                   year >= unique(B_minYear), 
                                   year <= unique(A_maxYear), 
                                   year <= unique(B_maxYear)),
                               TRUE,FALSE),
            by=c('A','B','year')]
  
  return(assocData)
  
}




print('Cleared functions')
