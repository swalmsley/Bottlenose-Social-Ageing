
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
inline <- function(fit, var, p, CI) {

  est <- summary(fit,prob=p)$fixed[paste(var),'Estimate']
  low <- summary(fit,prob=p)$fixed[paste(var),'l-90% CI']
  high <- summary(fit,prob=p)$fixed[paste(var),'u-90% CI']

  # report credible interval
  if (CI) output <- paste(paste(format_number(est), ', ', sep = ''), 'CI ', paste(format_number(low), format_number(high), sep = ' -- '), sep = '')

  # report probability of directional effect
  if (!CI) output <- paste(paste(format_number(est), ', ', sep = ''), p_dir <- paste('pd = ', round(as.numeric(p_direction(fit, parameters = paste(var)))*100,1), '%',sep=''), sep = '')

  return(output)

}
# fit <- tar_read(m3)
# var <- 'numSamplingPeriodsByYear'
# p = 0.9
# CI = FALSE




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



# Coefficients table (brms) -----------------------------------------------
brms_table_gam <- function(fit) {
  
  s <- summary(fit, prob = 0.9)
  
  fixed <- data.table(as.data.frame(s$fixed))
  Title <- data.table(as.data.frame(s$random$Title))
  year  <- data.table(as.data.frame(s$random$year))
  smooth <- if (!is.null(s$spec_pars)) data.table(as.data.frame(s$spec_pars)) else NULL
  
  raw.table <- rbindlist(list(fixed, Title, year, smooth), fill = TRUE)
  
  terms <- c(row.names(s$fixed),
             row.names(s$random$Title),
             row.names(s$random$year),
             if (!is.null(s$spec_pars)) row.names(s$spec_pars))
  
  stats.table <- cbind(terms, raw.table)
  setnames(stats.table, old = names(stats.table),
           new = c("Term", "Estimate", "SE", "PI-Lower", "PI-Upper", "Rhat", "Bulk_ESS", "Tail_ESS")[1:ncol(stats.table)])
  
  # drop diagnostics if present
  drop_cols <- intersect(c("Rhat","Bulk_ESS","Tail_ESS"), names(stats.table))
  if (length(drop_cols)) stats.table[, (drop_cols) := NULL]
  
  nice_table(stats.table)
}
#fit <- tar_read(m12_smooth2)





# Compare calf/juv classifications ----------------------------------------
young_NBW_kappa <- function(raw_photoData) {

  # Prep claire version

  claire <- data.table(read_xlsx('input/Calf_juv_updated_master.xlsx'))

  claire[, Claire_ageClass:=ifelse(grepl('claire-juv', `Keyword-export`),'Juv','none'),]
  claire[Claire_ageClass!='Juv', Claire_ageClass:=ifelse(grepl('claire-calf', `Keyword-export`),'Calf','none'),]
  claire[,.N,by='Claire_ageClass']


  claire[,year:=year(`Date Original`),]
  claire <- unique(claire[,c('Title','Claire_ageClass','year')])
  claire$Title <- as.character(claire$Title)
  unique(claire[,.N,by=c('year','Title')]$N) # expect 1

  # 1990 28, 1994 390, 2016 6240

  # error log for Claire version
  claire <- claire[!(year==1990 & Title=='28' & Claire_ageClass=='none'),,] # assuming Juv
  claire <- claire[!(year==1994 & Title=='390' & Claire_ageClass=='none'),,] # assuming Juv
  claire <- claire[!(year==2016 & Title=='6240' & Claire_ageClass=='none'),,] # assuming Juv

  unique(claire[,.N,by=c('year','Title')]$N) # expect 1 now that errors have been fixed/assumedunique(claire[,.N,by=c('year','Title')]$N) # expect 1
  ###### add these to update sheet

  # Prep raw version

  original <- raw_photoData

  original[, ageClass:=ifelse(grepl('ageJuv', Keyword.export),'Juv','none'),]
  original[ageClass!='Juv', ageClass:=ifelse(grepl('ageCalf', Keyword.export),'Calf','none'),]
  original[,.N,by='ageClass']

  original[,year:=year(`Date.Original`),]
  original <- unique(original[,c('Title','ageClass','year')])
  original[,nByYear:=length(unique(ageClass)),by=c('Title','year')] # for checking
  unique(original[,.N,by=c('year','Title')]$N) # expect 1

  # error log for next year
  original <- original[!(year==2015 & Title=='5059' & ageClass=='none'),,]
  original <- original[!(year==2015 & Title=='6172' & ageClass=='none'),,]
  original <- original[!(year==2021 & Title=='6527' & ageClass=='none'),,]
  ######
  
  original <- original[Title!='see crops' & Title!='unk',,]
  claire <- claire[Title!='see crops' & Title!='unk',,]

  unique(original[,.N,by=c('year','Title')]$N) # expect 1


  compare <- merge(original[Title!='unk' & Title!='see crops',c('Title','ageClass','year'),], claire[Title!='unk' & Title!='see crops',,], by=c('Title','year'))
  compare <- compare[year<=2023,,] # Claire only did classifications up to 2023 

  return(kappam.fleiss(compare[,3:4], detail=TRUE))

}
# raw_photoData <- tar_read(raw_photoData)



# Process photo-ID observations -------------------------------------------
process_photoID <- function(lv_data, chosen_side, chosen_canyons){

  # Manual edits for 2024
  # 1. Several photos with erroneous metadata, including times 
  photos_with_errors <- c('NBW_20240721_DSC_8272.NEF',
                          'NBW_20240721_DSC_8271.NEF',
                          'NBW_20240721_DSC_8267.NEF')
  lv_data <- lv_data[!(File.name %in% photos_with_errors),,]
  # 2. Fix to IDs 1317, 6018, 6402
  lv_data[Title==1317 & year(Date.Original)==2021, Title:=6018,]
  lv_data[Title==6402, Title:=6018, ]
  # all share unknown sex
  # 3. Fix to photograph with erroneous ID
  lv_data[File.name=='NBW_20240720_NB2_3650.JPG',Title:=6515,]

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

  ### Removing original version

  # # Define time threshold (e.g., 2 hours)
  # threshold <- 10 * 60  # 10 minutes in seconds
  #
  # # Carve groups based on time threshold
  # photoData$t_group <- cut(photoData$dateTime, breaks = seq(min(photoData$dateTime), max(photoData$dateTime) + threshold, by = threshold), labels = FALSE)
  # photoData[,myGroups:=.GRP,by=t_group]

  # Calculate the time difference in minutes
  photoData[, time_diff := c(0, diff(dateTime) / 60)]

  # Create a group identifier based on the threshold
  photoData[, myGroups := cumsum(time_diff > 10) + 1]

  # Summary statistics - group size and duration
  photoData[,groupSize:=as.numeric(length(unique(Title))),by=myGroups]
  photoData[, minTime:=min(dateTime), by=myGroups]
  photoData[, maxTime:=max(dateTime), by=myGroups]
  photoData[, groupDur:=as.numeric(maxTime-minTime, units='mins'), by=myGroups] # in minutes
  # extract unique
  u_groups <- unique(photoData[,c('myGroups','groupSize', 'groupDur')])
  # group size
  hist(u_groups$groupSize)
  range(u_groups$groupSize) # range 1-27
  median(u_groups$groupSize) # median 2
  mean(u_groups$groupSize) # mean 2.4
  # group duration
  hist(u_groups$groupDur)
  range(u_groups$groupDur) # range 0-74 mins
  median(u_groups$groupDur) # median 1 min
  mean(u_groups$groupDur) # mean 4.4 mins

  # subset to one row per group per individual
  group_ind_data <- unique(photoData[,c('Title','year', 'groupSize','myGroups')])

  group_ind_data[,gs:=groupSize,] # group size for each individual

  group_df <- unique(group_ind_data[,c('Title','myGroups','year','gs')])

  return(group_df)

}
# photoData <- tar_read(side_ungrouped)
# meta <- tar_read(nbw_meta)




# Function to create histogram of minimum age -----------------------------
age_histogram <- function(group_df, s) {

  col <- ifelse(s=='Male','#762F7C','#FDAE79')

  d <- unique(group_df[,c('Title','year','minimumAge','sex', 'ageClass', 'yearSpan', 'numYears')])
  g <- ggplot(d[sex==s & ageClass=='Adult' & numYears>1,,], aes(x=minimumAge)) +
    geom_histogram(binwidth = 2, fill=col) +
    labs(x='Minimum age', y='Number of annual observations')+
    theme_classic()+
    theme(strip.text=element_text(size=12),
          axis.text=element_text(size=12),
          axis.title=element_text(size=12),
          legend.position = "none")
  g
}
# group_df <- tar_read(group_df)
# s <- 'Male'




# Extract within-individual effects ---------------------------------------
# Note that this function uses "coef" and is intended to combine varying slopes with average trait-specific within-individual effects
extract_within_effects <- function(multi_dataset, m) {

  c <- coef(m, prob=c(0.95, 0.05))$Title

  gs <- data.frame(c[,,"gsannual_deltaMinAge"])
  colnames(gs) <- paste('gsannual_',colnames(gs),sep='')
  gs$Title <- rownames(gs)

  degree <- data.frame(c[,,"degree_deltaMinAge"])
  colnames(degree) <- paste('degree_',colnames(degree),sep='')
  degree$Title <- rownames(degree)

  strength <- data.frame(c[,,"strength_deltaMinAge"])
  colnames(strength) <- paste('strength_',colnames(strength),sep='')
  strength$Title <- rownames(strength)

  eigen <- data.frame(c[,,"eigen_deltaMinAge"])
  colnames(eigen) <- paste('eigen_',colnames(eigen),sep='')
  eigen$Title <- rownames(eigen)

  meanEdge <- data.frame(c[,,"medianEdge_deltaMinAge"])
  colnames(meanEdge) <- paste('meanEdge_',colnames(meanEdge),sep='')
  meanEdge$Title <- rownames(meanEdge)

  # combine effects for different traits
  combined <- Reduce(merge, list(gs, degree, strength, eigen, meanEdge))

  # incorporating sex information
  dt <- merge(combined, unique(multi_dataset[,c('Title','sex'),]),by='Title', all.x=TRUE, all.y=FALSE)

  return(dt)

}
# m <- tar_read(m_multi_both_ma)
# multi_dataset <- tar_read(multi_dataset_ma)


# Extract within-individual effects for minimum age ---------------------------------------
# Note that this function uses "coef" and is intended to combine varying slopes with average trait-specific within-individual effects
extract_within_effects_minimum <- function(multi_dataset, m) {
  
  c <- coef(m, prob=c(0.95, 0.05))$Title
  
  gs <- data.frame(c[,,"gsannual_minimumAge"])
  colnames(gs) <- paste('gsannual_',colnames(gs),sep='')
  gs$Title <- rownames(gs)
  
  degree <- data.frame(c[,,"degree_minimumAge"])
  colnames(degree) <- paste('degree_',colnames(degree),sep='')
  degree$Title <- rownames(degree)
  
  strength <- data.frame(c[,,"strength_minimumAge"])
  colnames(strength) <- paste('strength_',colnames(strength),sep='')
  strength$Title <- rownames(strength)
  
  eigen <- data.frame(c[,,"eigen_minimumAge"])
  colnames(eigen) <- paste('eigen_',colnames(eigen),sep='')
  eigen$Title <- rownames(eigen)
  
  meanEdge <- data.frame(c[,,"medianEdge_minimumAge"])
  colnames(meanEdge) <- paste('meanEdge_',colnames(meanEdge),sep='')
  meanEdge$Title <- rownames(meanEdge)
  
  # combine effects for different traits
  combined <- Reduce(merge, list(gs, degree, strength, eigen, meanEdge))
  
  # incorporating sex information
  dt <- merge(combined, unique(multi_dataset[,c('Title','sex'),]),by='Title', all.x=TRUE, all.y=FALSE)
  
  return(dt)
  
}
# m <- tar_read(m_multi_both_ma_minimum)
# multi_dataset <- tar_read(multi_dataset_ma)


# Extract within-individual effects ---------------------------------------
# Note that this function uses "coef" and is intended to combine varying slopes with average trait-specific within-individual effects
extract_within_effects_sex_specific <- function(multi_dataset, m) {
  
  c <- coef(m, prob=c(0.95, 0.05))$Title
  
  gs <- data.frame(c[,,"gsannual_deltaMinAge"])
  colnames(gs) <- paste('gsannual_',colnames(gs),sep='')
  gs$Title <- rownames(gs)
  
  degree <- data.frame(c[,,"degree_deltaMinAge"])
  colnames(degree) <- paste('degree_',colnames(degree),sep='')
  degree$Title <- rownames(degree)
  
  strength <- data.frame(c[,,"strength_deltaMinAge"])
  colnames(strength) <- paste('strength_',colnames(strength),sep='')
  strength$Title <- rownames(strength)
  
  eigen <- data.frame(c[,,"eigen_deltaMinAge"])
  colnames(eigen) <- paste('eigen_',colnames(eigen),sep='')
  eigen$Title <- rownames(eigen)
  
  meanEdge <- data.frame(c[,,"medianEdge_deltaMinAge"])
  colnames(meanEdge) <- paste('meanEdge_',colnames(meanEdge),sep='')
  meanEdge$Title <- rownames(meanEdge)
  
  # combine effects for different traits
  combined <- Reduce(merge, list(gs, degree, strength, eigen, meanEdge))
  
  # incorporating sex information
  dt <- merge(combined, unique(multi_dataset[,c('Title','sex'),]),by='Title', all.x=TRUE, all.y=FALSE)
  
  return(dt)
  
}
# multi_dataset <- tar_read(multi_dataset_ma)
# m <- tar_read(m_multi_f_ma)




# Plot within-individual effects ------------------------------------------
plot_within_effect <- function(w, traitX, traitY) {

  w$sex <- factor(w$sex, levels=c('FemaleJ','Male'))

  ggplot(w, aes(x=get(paste(traitX,'_Estimate',sep='')), y=get(paste(traitY, '_Estimate',sep='')), color=sex))+
    geom_point(alpha=0.75,size=2) +
    labs(x='',y='')+ ###### check these
    scale_color_manual(values = c('#FDAE79','#762F7C'))+
    theme_classic() +
    theme(legend.position='none')

}
# w <- tar_read(within_effects_ma)
# traitX <- 'degree'
# traitY <- 'meanEdge'




# Plot within-individual effects ------------------------------------------
plot_within_effect_sex_specific <- function(w_f, w_m, traitX, traitY) {
  
  w <- rbindlist(list(w_f, w_m))
  w$sex <- factor(w$sex, levels=c('FemaleJ','Male'))
  
  ggplot(w, aes(x=get(paste(traitX,'_Estimate',sep='')), y=get(paste(traitY, '_Estimate',sep='')), color=sex))+
    geom_point(alpha=0.4,size=2) +
    geom_vline(xintercept = 0, linetype='dashed')+ geom_hline(yintercept = 0, linetype='dashed')+
    labs(x='',y='')+ ###### check these
    scale_color_manual(values = c('#FDAE79','#762F7C'))+
    theme_classic() +
    theme(legend.position='none')
  
}
# w_f <- tar_read(within_effects_f_ma)
# w_m <- tar_read(within_effects_m_ma)
# traitX <- 'degree'
# traitY <- 'meanEdge'




# Create pairwise figure showing social ageing correlations ---------------
plot_within_effects <- function(w) {

  p1 <- plot_within_effect(w, 'degree', 'gsannual')
  p2 <- plot_within_effect(w, 'strength', 'gsannual')
  p3 <- plot_within_effect(w, 'eigen', 'gsannual')
  p4 <- plot_within_effect(w, 'meanEdge', 'gsannual')

  p5 <- plot_within_effect(w, 'strength', 'degree')
  p6 <- plot_within_effect(w, 'eigen', 'degree')
  p7 <- plot_within_effect(w, 'meanEdge', 'degree')

  p8 <- plot_within_effect(w, 'eigen', 'strength')
  p9 <- plot_within_effect(w, 'meanEdge', 'strength')

  p10 <- plot_within_effect(w, 'meanEdge', 'eigen')

  row_label_1 <- wrap_elements(panel = textGrob('Group size'))
  row_label_2 <- wrap_elements(panel = textGrob('N. Partners'))
  row_label_3 <- wrap_elements(panel = textGrob('Strength'))
  row_label_4 <- wrap_elements(panel = textGrob('Centrality'))

  col_label_1 <- wrap_elements(panel = textGrob('N. Partners'))
  col_label_2 <- wrap_elements(panel = textGrob('Strength'))
  col_label_3 <- wrap_elements(panel = textGrob('Centrality'))
  col_label_4 <- wrap_elements(panel = textGrob('Mean Relationship'))

  legend_plot <- ggplot(data.frame(x=c(1,1), y=c(0.5,-0.5), sex=c('f','m')), aes(x=x,y=y))+
    geom_point(aes(fill=sex), size=10, pch=22)+
    scale_fill_manual(values=c('#762F7C','#FDAE79'))+
    annotate('text',x=0,y=0.5,label='Male ageing correlations', size=2)+
    annotate('text',x=0,y=-0.5,label='Female ageing correlations', size=2)+
    xlim(-1,1.3)+
    ylim(-1,1)+
    theme(panel.background = element_blank(), axis.text=element_blank(), axis.ticks = element_blank(),
          axis.title=element_blank(), legend.position = 'none')

  p <- (col_label_1 | col_label_2 | col_label_3 | col_label_4 | plot_spacer()) /
    (p1 | p2 | p3 | p4 | row_label_1) /
    (legend_plot  | p5 | p6 | p7 | row_label_2) /
    (plot_spacer() | plot_spacer() | p8 | p9 | row_label_3) /
    (plot_spacer() | plot_spacer() | plot_spacer() | p10 | row_label_4) +
    plot_layout(widths=c(1,1,1,1,1))+
    plot_annotation(title='Correlated within-individual effects of age on social traits')

  p

  # return(p)

}
# w <- tar_read(within_effects_ma)
# w <- tar_read(within_effects_ma_minAge)




# Create pairwise figure showing social ageing correlations ---------------
plot_within_effects_sex_specific <- function(w_f, w_m) {
  
  p1 <- plot_within_effect_sex_specific(w_f, w_m, 'degree', 'gsannual')
  p2 <- plot_within_effect_sex_specific(w_f, w_m, 'strength', 'gsannual')
  p3 <- plot_within_effect_sex_specific(w_f, w_m, 'eigen', 'gsannual')
  p4 <- plot_within_effect_sex_specific(w_f, w_m, 'meanEdge', 'gsannual')
  
  p5 <- plot_within_effect_sex_specific(w_f, w_m, 'strength', 'degree')
  p6 <- plot_within_effect_sex_specific(w_f, w_m, 'eigen', 'degree')
  p7 <- plot_within_effect_sex_specific(w_f, w_m, 'meanEdge', 'degree')
  
  p8 <- plot_within_effect_sex_specific(w_f, w_m, 'eigen', 'strength')
  p9 <- plot_within_effect_sex_specific(w_f, w_m, 'meanEdge', 'strength')
  
  p10 <- plot_within_effect_sex_specific(w_f, w_m, 'meanEdge', 'eigen')
  
  row_label_1 <- wrap_elements(panel = textGrob('Group size'))
  row_label_2 <- wrap_elements(panel = textGrob('N. Partners'))
  row_label_3 <- wrap_elements(panel = textGrob('Strength'))
  row_label_4 <- wrap_elements(panel = textGrob('Centrality'))
  
  col_label_1 <- wrap_elements(panel = textGrob('N. Partners'))
  col_label_2 <- wrap_elements(panel = textGrob('Strength'))
  col_label_3 <- wrap_elements(panel = textGrob('Centrality'))
  col_label_4 <- wrap_elements(panel = textGrob('Mean Relationship'))
  
  legend_plot <- ggplot(data.frame(x=c(1,1), y=c(0.5,-0.5), sex=c('f','m')), aes(x=x,y=y))+
    geom_point(aes(fill=sex), size=10, pch=22)+
    scale_fill_manual(values=c('#762F7C','#FDAE79'))+
    annotate('text',x=0,y=0.5,label='Male ageing correlations', size=2)+
    annotate('text',x=0,y=-0.5,label='Female ageing correlations', size=2)+
    xlim(-1,1.3)+
    ylim(-1,1)+
    theme(panel.background = element_blank(), axis.text=element_blank(), axis.ticks = element_blank(),
          axis.title=element_blank(), legend.position = 'none')
  
  p <- (col_label_1 | col_label_2 | col_label_3 | col_label_4 | plot_spacer()) /
    (p1 | p2 | p3 | p4 | row_label_1) /
    (legend_plot  | p5 | p6 | p7 | row_label_2) /
    (plot_spacer() | plot_spacer() | p8 | p9 | row_label_3) /
    (plot_spacer() | plot_spacer() | plot_spacer() | p10 | row_label_4) +
    plot_layout(widths=c(1,1,1,1,1))+
    plot_annotation(title='Correlated within-individual effects of age on social traits')
  
  p
  
  # return(p)
  
}
# w_f <- tar_read(within_effects_f_ma)
# w_m <- tar_read(within_effects_m_ma)

# w_f <- tar_read(within_effects_f_ma_minAge)
# w_nm <- tar_read(within_effects_m_ma_minAge)




# Extract network traits --------------------------------------------------
extract_trait_multiAnnual <- function(fit) {
  
  # extract all edge samples
  all_edge_samples <- data.frame(coef(fit, summary=FALSE)[[1]])
  colnames(all_edge_samples) <- str_extract(colnames(all_edge_samples), "\\d+_\\d+\\.\\d+")
  
  # extract year group and proper dyad names
  d <- data.table(fit$data)
  d[,year_group:=str_extract(dyad_annual, "(?<=-)[0-9]+$"), by=.I]
  d[,dyad:=str_extract(dyad_annual, "^\\d+_\\d+"), by=.I]  
  
  # year groups to loop through
  years <- unique(d$year_group)
  years <- years[!years %in% c(4, 18)]
  
  # loop through years
  for (y in years) {
    
    edge_samples <- all_edge_samples[,(str_extract(colnames(all_edge_samples), "(?<=\\.)\\d+")==y)]
    colnames(edge_samples) <- str_remove(colnames(edge_samples), "\\.\\d+")
    
    N_draws <- 1000
    
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
      fitData <- d[year_group==y,,]
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
      NZ_net <- delete_edges(NZ_net, E(NZ_net)[is.na(E(NZ_net)$weight)])
      
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
    
    # create more final output dataframe
    out <- data.table(ID = colnames(metric_samples))
    
    out$eigen <- apply(metric_samples, 2, mean)
    out$eigen_SD <- apply(metric_samples, 2, sd)
    
    out$strength <- apply(metric_samples2, 2, mean)
    out$strength_SD <- apply(metric_samples2, 2, sd)
    
    out$medianEdge <- apply(metric_samples3, 2, mean) # technically mean of mean edges
    out$medianEdge_SD <- apply(metric_samples3, 2, sd)
    
    ## FOR STANDARDIZATION ONLY ##
    # metric_samples <- standardize_matrix(metric_samples)
    # metric_samples2 <- standardize_matrix(metric_samples2)
    # metric_samples3 <- standardize_matrix(metric_samples3)
    ## FOR STANDARDIZATION ONLY ##
    
    out[, year_group:=y, ]
    
    if (exists('combined')) (combined <- rbindlist(list(combined, out)))
    if (!exists('combined')) (combined <- out)
    
  }
  
  return(combined)
  
}
# fit <- tar_read(brms_fit_group_multiAnnual)



# Extract network traits --------------------------------------------------
extract_trait_multiAnnual_imputation <- function(fit, N_draws) {
  
  # extract all edge samples
  all_edge_samples <- data.frame(coef(fit, summary=FALSE)[[1]])
  colnames(all_edge_samples) <- str_extract(colnames(all_edge_samples), "\\d+_\\d+\\.\\d+")
  
  # extract year group and proper dyad names
  d <- data.table(fit$data)
  d[,year_group:=str_extract(dyad_annual, "(?<=-)[0-9]+$"), by=.I]
  d[,dyad:=str_extract(dyad_annual, "^\\d+_\\d+"), by=.I]  
  
  # Randomly select draws to use in advance
  chosen_draws <- sample(1:nrow(all_edge_samples), size = N_draws, replace = FALSE)    
  
  # year groups to loop through
  years <- unique(d$year_group)
  years <- years[!years %in% c(4, 18)]
  
  # Initialize result
  result <- list()
  
  # loop through draws
  for (draw in chosen_draws) {
    
    # set output to NULL
    output <- NULL
    
    # loop through years
    for (y in years) {
      
      # extract relevant edge weights for given year
      annual_edge_samples <- all_edge_samples[,(str_extract(colnames(all_edge_samples), "(?<=\\.)\\d+")==y)]
      colnames(annual_edge_samples) <- str_remove(colnames(annual_edge_samples), "\\.\\d+")
      
      # pull draw of complete network (i.e., of each edge weight)
      draw_of_network <- as.data.frame(t(annual_edge_samples[draw,])) # probably more efficient way to do this (i.e., full matrix...)
      if (nrow(draw_of_network)==1) (rownames(draw_of_network) <- colnames(annual_edge_samples)) # for 1991 with 2 photos... can simply delete this year too maybe
      colnames(draw_of_network) <- 'weight_sample_raw'
      
      # format IDs out of dyad label
      draw_of_network$A <- sub("^(\\d+)_.*", "\\1", rownames(draw_of_network)) # pull out ID A
      draw_of_network$B <- sub("^\\d+_(\\d+)$", "\\1", rownames(draw_of_network)) # pull out ID B
      
      ###########################################################
      ### Create main network with all Bayesian uncertainties ###
      
      # create network
      net <- igraph::graph_from_edgelist(as.matrix(draw_of_network[,c('A','B'),]), directed=FALSE)
      # add in edgeweights
      igraph::E(net)$weight <- inv_logit(draw_of_network[,c('weight_sample_raw'),]) # inv-logit here because it's a bernoulli model # confirm makes sense to transform prior to strength calculation
      
      # create version of network without 0s
      fitData <- d[year_group==y,,]
      fitData[,nTogether:=sum(together),by=dyad]
      nonZero_dyads <- unique(fitData[nTogether>0,c('dyad')])$dyad
      
      NZ_draw_of_network <- data.table(copy(draw_of_network))
      NZ_draw_of_network[,dyad:=rownames(draw_of_network)]
      NZ_draw_of_network[,weight:=inv_logit(weight_sample_raw),] # inv-logit transformation
      # now set dyads that were never together as 0s
      NZ_draw_of_network[!(dyad %in% nonZero_dyads), weight:=NA,]
      
      NZ_net <- igraph::graph_from_edgelist(as.matrix(NZ_draw_of_network[,c('A','B'),]), directed=FALSE)
      #plot(NZ_net)
      
      igraph::E(NZ_net)$weight <- NZ_draw_of_network[,weight,] ###### inv-logit here because it's a binary model # confirm makes sense to transform prior to strength calculation
      NZ_net <- delete_edges(NZ_net, E(NZ_net)[is.na(E(NZ_net)$weight)])
      
      # plot network if desired
      # plot(net)
      
      ## Add traits
      
      # Strength
      strength_df <- as.data.frame(igraph::strength(net))
      strength_df$Title <- rownames(strength_df)
      strength_df <- data.table(strength_df)
      colnames(strength_df) <- c('strength', 'Title')
      
      # Centrality
      centrality_df <- as.data.frame(igraph::eigen_centrality(net)$vector)
      centrality_df$Title <- rownames(centrality_df)
      centrality_df <- data.table(centrality_df)
      colnames(centrality_df) <- c('eigen', 'Title')
      centrality_df[eigen==1,eigen:=eigen-0.0001,]
      centrality_df[eigen==0,eigen:=eigen+0.0001,]
      
      
      # Mean relationship 
      # For mean non-0 edge weight
      meanRel_df <- data.table(Title=strength_df$Title)
      
      for (name in strength_df$Title) {

        # extract edge weights using igraph and cut below 5% threshold
        focal_weights <- data.table(dyad=incident(net,v=name), weight=(incident(net,v=name))$weight)
        #focal_weights_nonZero <- focal_weights[focal_weights>=0.05] # using 5% for Bayesian weight right now - could use real data (i.e., trim those never seen together or consider alternate threshold)
        focal_weights_nonZero <- data.table(dyad=incident(NZ_net,v=name), weight=(incident(NZ_net,v=name))$weight)
        
        # save out mean edge weight by ID
        meanRel_df[Title==name,meanEdge:=mean(focal_weights_nonZero$weight),]
        meanRel_df[is.nan(meanEdge),meanEdge:=NA,]
      
      }

      # format nodal measures
      nodals <- merge(strength_df, centrality_df, by='Title')
      nodals <- merge(nodals, meanRel_df, by='Title')
      nodals$year_group <- y
    
      # Save results from each year to draw-specific output
      output <- rbindlist(list(output, nodals)) # if output is null at start, will simply give trait_df
      
    }
    
    # Save results from single draw to main results
    result <- c(result, list(output[,ID:=Title,]))
    # result <- c(result, list(as.data.frame(output)))
    
    
  }
  
  return(result)
  
}
# fit <- tar_read(brms_fit_group_multiAnnual)
# N_draws <- 3





combine_trait_dfs_ma <- function(nd, key) {
  
  key[,year_group:=as.character(year_group),]
  combined <- nd[,Title:=ID,]
  combined <- merge(combined, key, by='year_group')
  combined[,ID:=NULL,]
  
  return(combined)
  
}
# nd <- tar_read(nodal_traits_multiAnnual)
# combined_traits <- tar_read(year_key)




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
  label_col <- ifelse(sex=='Male', '#762F7C','#FDAE79')

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
# fit <- tar_read(ma_m5)
# traitName <- 'eigen'
# traitLabel <- 'Centrality (std.)'
# sex='Male'



# Sex coefficient plot ----------------------------------------------------
sex_coef_plot <- function(fitM, fitF) {

  m <- data.table(mcmc_intervals(fitM, pars=c('b_meanMinAge', 'b_deltaMinAge'), point_est = 'median', prob=0.9, prob_outer = 0.9)$data)
  f <- data.table(mcmc_intervals(fitF, pars=c('b_meanMinAge', 'b_deltaMinAge'), point_est = 'median', prob=0.9, prob_outer = 0.9)$data)

  pars_df <- rbindlist(list(m,f))

  pars_df$model <- rep(c('Model M', 'Model F'), each = nrow(pars_df)/2)
  pars_df[parameter=='b_meanMinAge',parameter:='Between',]
  pars_df[parameter=='b_deltaMinAge',parameter:='Within',]

  pars_df$parameter <- factor(pars_df$parameter,levels = c('Within','Between'))

  # Adding an ID for each parameter for plotting
  #pars_df$param_id <- rep(1:(nrow(pars_df)/2), 2)

  #cols = c(viridis(1, begin=0.1,end=0.2, option='A'),viridis(1, begin=0.8,end=0.9, option='A')) # for flexible viridis approach
  cols <- c('#762F7C','#FDAE79')

  # Plotting
  p <- ggplot(pars_df, aes(y = parameter, xmin = ll, xmax = hh, x = m, color=model)) +
    geom_errorbarh(aes(height = 0), position = position_dodge(width = 0.25), linewidth=1.75, alpha=0.75) +
    geom_point(position = position_dodge(width = 0.25), size=3.25) +
    geom_vline(xintercept = 0, linetype='dashed', linewidth=0.75)+
    scale_color_manual(values=c(cols[2],cols[1])) +
    labs(x = 'Posterior', y = '') +
    #xlim(-0.1,0.1)+
    theme_classic()+theme(legend.position='none')+
    theme(axis.text = element_text(size=12))

  return(p)

}
# fitM <- tar_read(m5)
# fitF <- tar_read(m6)




# Variability plot --------------------------------------------------------
var_plot <- function(fitM, fitF) {

  # set up colors
  cols <- c('#FDAE79','#762F7C')

  ## density for variability
  m <- data.table(as_draws_df(fitM)$sd_Title__deltaMinAge)
  m[,sex:='m',]
  f <- data.table(as_draws_df(fitF)$sd_Title__deltaMinAge)
  f[,sex:='f',]

  d <- rbindlist(list(m,f))
  d$sex <- factor(d$sex, levels=c('f','m'))

  p <- ggplot(d, aes(x=V1, group=sex, fill=sex, color=sex))+
    geom_density(alpha=0.75)+
    scale_color_manual(values=cols)+
    scale_fill_manual(values=cols)+
    labs(x='Variability of within-individual slopes', y='Density')+
    theme_classic()+
    theme(axis.text = element_text(size=15),
          legend.position='none')

  return(p)

}
# fitM <- tar_read(m5)
# fitF <- tar_read(m6)




# plot correlations in social ageing effects
plot_slope_correlation <- function(parName, male_m, female_m) {

  # set up colors
  cols <- c('#FDAE79','#762F7C')

  d_m <- as_draws_df(male_m)
  d_f <- as_draws_df(female_m)

  m_effect <- d_m[,parName]
  m_effect$sex <- 'm'

  f_effect <- d_f[,parName]
  f_effect$sex <- 'f'

  effects <- rbindlist(list(m_effect, f_effect))
  effects$sex <- factor(effects$sex, levels=c('f','m'))

  p <- ggplot(effects, aes(x=get(parName), group=sex, fill=sex, color=sex))+
    geom_density(alpha=0.75)+
    scale_color_manual(values=cols)+
    scale_fill_manual(values=cols)+
    labs(x=paste(parName))+
    labs(x='Correlation', y='Density')+
    theme_classic()+
    geom_vline(xintercept=0, linetype='dashed')+
    xlim(-1,1)+
    theme(axis.title=element_text(size=7),axis.text=element_text(size=7), legend.position='none')

  p

}
# parName <- 'cor_Title__eigen_deltaMinAge__degree_deltaMinAge'
# male_m <- tar_read(m_multi_m_ma)
# female_m <- tar_read(m_multi_f_ma)




# Create pairwise figure showing social ageing correlations ---------------
plot_slope_correlations <- function(male_m, female_m) {

  p1 <- plot_slope_correlation('cor_Title__gsannual_deltaMinAge__degree_deltaMinAge', male_m, female_m)
  p2 <- plot_slope_correlation('cor_Title__gsannual_deltaMinAge__strength_deltaMinAge', male_m, female_m)
  p3 <- plot_slope_correlation('cor_Title__gsannual_deltaMinAge__eigen_deltaMinAge', male_m, female_m)
  p4 <- plot_slope_correlation('cor_Title__gsannual_deltaMinAge__medianEdge_deltaMinAge', male_m, female_m)

  p5 <- plot_slope_correlation('cor_Title__strength_deltaMinAge__degree_deltaMinAge', male_m, female_m)
  p6 <- plot_slope_correlation('cor_Title__eigen_deltaMinAge__degree_deltaMinAge', male_m, female_m)
  p7 <- plot_slope_correlation('cor_Title__medianEdge_deltaMinAge__degree_deltaMinAge', male_m, female_m)

  p8 <- plot_slope_correlation('cor_Title__strength_deltaMinAge__eigen_deltaMinAge', male_m, female_m)
  p9 <- plot_slope_correlation('cor_Title__medianEdge_deltaMinAge__strength_deltaMinAge', male_m, female_m)

  p10 <- plot_slope_correlation('cor_Title__medianEdge_deltaMinAge__eigen_deltaMinAge', male_m, female_m)


  row_label_1 <- wrap_elements(panel = textGrob('Group size'))
  row_label_2 <- wrap_elements(panel = textGrob('N. Partners'))
  row_label_3 <- wrap_elements(panel = textGrob('Strength'))
  row_label_4 <- wrap_elements(panel = textGrob('Centrality'))

  col_label_1 <- wrap_elements(panel = textGrob('N. Partners'))
  col_label_2 <- wrap_elements(panel = textGrob('Strength'))
  col_label_3 <- wrap_elements(panel = textGrob('Centrality'))
  col_label_4 <- wrap_elements(panel = textGrob('Mean Relationship'))


  legend_plot <- ggplot(data.frame(x=c(1,1), y=c(0.5,-0.5), sex=c('f','m')), aes(x=x,y=y))+
    geom_point(aes(fill=sex), size=10, pch=22)+
    scale_fill_manual(values=c('#762F7C','#FDAE79'))+
    annotate('text',x=0,y=0.5,label='Male ageing correlations', size=2)+
    annotate('text',x=0,y=-0.5,label='Female ageing correlations', size=2)+
    xlim(-1,1.3)+
    ylim(-1,1)+
    theme(panel.background = element_blank(), axis.text=element_blank(), axis.ticks = element_blank(),
          axis.title=element_blank(), legend.position = 'none')
  legend_plot

  p <- (col_label_1 | col_label_2 | col_label_3 | col_label_4 | plot_spacer()) /
    (p1 | p2 | p3 | p4 | row_label_1) /
    (legend_plot  | p5 | p6 | p7 | row_label_2) /
    (plot_spacer() | plot_spacer() | p8 | p9 | row_label_3) /
    (plot_spacer() | plot_spacer() | plot_spacer() | p10 | row_label_4) +
    plot_layout(widths=c(1,1,1,1,1))+
    plot_annotation(title='Correlated within-individual effects of age on social traits')

  return(p)

}
# male_m <- tar_read(m_multi_m_ma)
# female_m <- tar_read(m_multi_f_ma)




# Pull out biopsy sex data ------------------------------------------------
# extract list of titles with biopsy sex information
biopsySexTitles <- function(d) {

  d[grepl("biopsy-Female", Keyword.export), biopsySex:='Y',]
  d[grepl("biopsy-Male", Keyword.export), biopsySex:='Y',]
  titles <- unique(d[Title!='unk' & biopsySex=='Y',Title,])

  return(titles)


}
#d <- tar_read(raw_photoData)



# Plot GAMs ---------------------------------------------------------------
smooth_plot <- function(fit, traitName, traitLabel, show_data, add_legend) {
  
  # colors
  cols <- c('Male' = '#762F7C', 'FemaleJ' = '#db874f')
  
  # data
  data <- data.table(fit$data)
  
  # # New data to predict for
  newdata <- expand.grid(minimumAge = seq(0,37,1), # based on range in data
                         Title = NA,
                         sex = unique(data$sex),
                         numSamplingPeriodsByYear = mean(data$numSamplingPeriodsByYear), # note this will give warning if model does not control for numSamplingPeriodsByYear, no problem
                         year = NA)
  newdata[[paste(traitName, '_SD', sep='')]] <- 0.001 # assign basically no measurement error - doesn't really seem to make a difference? maybe to CI?
  
  pred <- data.table(epred_draws(fit, newdata, re_formula=NA))
  pred_means <- pred[, .(mean_epred=mean(.epred)),by=c('sex', 'minimumAge')] # extract means for each combination predictors (if not interested in plotting whole ribbon!)
  
  pred$sex <- factor(pred$sex, levels=c('Male', 'FemaleJ'))
  
  g <- ggplot(pred, aes(x=minimumAge, y=.epred, color=sex, fill=sex)) +
    
    stat_lineribbon(alpha=0.1, .width=seq(0.1, 0.9, by=0.1)) +
    stat_lineribbon(alpha=0.8, .width=c(0.0)) +
    
    scale_color_manual(values=cols)+
    scale_fill_manual(values=cols)+
    
    labs(x='Minimum age', y=paste(traitLabel))+
    
    theme_classic() + theme(legend.position='none')
  
  if (show_data==TRUE) (g <- g + geom_jitter(inherit.aes = FALSE, data=data,size=1,aes(x=minimumAge,y= get(traitName), color=sex),alpha=0.05))
  
  if (add_legend==TRUE) (g <- g + add_phylopic(name='Hyperoodon ampullatus', x=5, y=6.5, width = 12, alpha=0.85, fill = '#762F7C') + 
                           add_phylopic(name='Hyperoodon ampullatus', x=5, y=6, width = (12*0.85), alpha=0.85, fill= '#db874f') + 
                           annotate('text',x=17,y=6.5,label='Males', color='#762F7C', size=3)+
                           annotate('text',x=17,y=6,label='Females', color='#db874f', size=3))

  g
  
  return(g)
                             
    
}
# fit <- tar_read(m12_smooth2)
# traitName <- 'gs'
# traitLabel <- 'Group size'
# show_data <- FALSE
# add_legend <- TRUE

# fit <- tar_read(m12_smooth2_freek)
# fit <- tar_read(m12_smooth2_no2024)
# fit <- tar_read(m12_smooth2_2024)
# fit <- tar_read(m12_smooth2_2023)


# Add degree to dataframe -------------------------------------------------
###### This needs to be checked over carefully
add_degree_ma <- function(side, associations, model_df, year_key) {
  
  associations <- merge(associations, year_key, by='year_group')
  
  for (y in unique(model_df$year)) {
    
    association_data <- associations[year==y,,] # pull out association target for year in question
    together <- association_data[together==1,,] # 'together' only
    
    temp <- model_df[year==y,,]
    
    for (focal in unique(temp$Title)) {
      
      partner_df <- together[A==focal | B==focal,,]
      # add degree to updated model df
      
      ###### important: because group sampling period here, can have multiple rows per dyad per day, so unique() is necessary...
      model_df[year==y & Title==focal,degree:=max(0,length(unique(c(partner_df$A,partner_df$B)))-1),] # minus one so not to count focal individual, max 0 so that if empty, reverts to 0 not -1
      
      numDays <- side[year==y & Title==focal,length(unique(day)),] 
      
      numGroups <- length(unique(c(association_data[A_present & A==focal, unique(group_ID), ], association_data[B_present & B==focal, unique(group_ID), ])))
      
      model_df[year==y & Title==focal,numDaysByYear:=numDays,]
      model_df[year==y & Title==focal,numSamplingPeriodsByYear:=numGroups,]
      
    }
    
  }
  
  model_df[, standardized_numSamplingPeriodsByYear:=lapply(.SD, function(x) as.vector(scale(numSamplingPeriodsByYear))), by=year]
  
  return(model_df)
  
}
# side <- tar_read(side_ungrouped)
# associations <- tar_read(group_associations_all)
# model_df <- tar_read(model_df_ma_step2)
# year_key <- tar_read(year_key)



# # Build NBW metadata ------------------------------------------------------
# # takes raw LV output as input ('dt')
# build_nbw_metadata <- function(dt, chosen_side) {
# 
#   # format dates
#   dt[,dateTime:=as.POSIXct(Date.Original),]
#   dt[,day:=format(dateTime, format = "%Y-%m-%d")]
#   dt[,year:=year(day),]
# 
#   # remove photographs that were processed as crops
#   dt <- dt[Title!="see crops",,]
#   dt <- dt[Title!="unk",,]
# 
#   # add left and right columns
#   # dt=dt %>% mutate(side = ifelse(grepl("left", Keyword.export), "left", NA))
#   # dt=dt %>% mutate(side = ifelse(grepl("right", Keyword.export), "right", NA))
#   #
#   dt[grepl("Left", Keyword.export), side := "left"]
#   dt[grepl("Right", Keyword.export), side := "right"]
# 
#   # subset to chosen side
#   dt <- dt[side==chosen_side,,]
# 
#   # add columns with sex details
#   dt=dt%>%mutate(sex = ifelse(grepl("FemaleJ", Keyword.export), "Female-Juvenile", NA))
#   dt=dt%>%mutate(sex = ifelse(grepl("Male", Keyword.export), "Male", sex))
#   dt$sex <- factor(dt$sex,levels=c("Female-Juvenile","Male"))
# 
#   # Juv and Calf classifications (using SFW's ratings)
#   dt[,ageClass:=ifelse(grepl('samCalf', Keyword.export),'Calf','Adult'),by=c('Title','year')] # really more like "not young"
#   dt[ageClass!='Calf',ageClass:=ifelse(grepl('samJuv', Keyword.export),'Juvenile',ageClass),by=c('Title','year')]
#   dt[,young:=ageClass %in% c('Calf','Juvenile'),by=c('Title','year')]
# 
#   # add minimum age by year
#   dt[,minYear:=min(year),by=Title]
#   dt[,maxYear:=max(year),by=Title]
#   dt[,yearSpan:=1+(maxYear-minYear),by=Title]
#   dt[,catalogueAge:=year-minYear,] # Note: considering first year as "0"
# 
#   ## photo-ID error fix
#   dt[Title==6527 & year==2021,ageClass:='Juvenile',]
# 
#   dt[, FirstYearAgeClass := unique(ageClass[minYear == year]), by = Title]
# 
#   dt[FirstYearAgeClass=='Calf',minimumAge:=catalogueAge,]
#   dt[FirstYearAgeClass=='Juvenile',minimumAge:=catalogueAge+1,]
#   dt[FirstYearAgeClass=='Adult',minimumAge:=catalogueAge+3,]
# 
#   # add day of year
#   dt[,yday:=yday(day),]
# 
#   # diagnostic tests
#   test_that("single sex classification for each individual", {
#     expect_true(unique(dt[,length(unique(sex)),by=Title]$V1)==1)
#   })
# 
#   test_that("single minimum age for each individual in each year", {
#     expect_true(unique(dt[,length(unique(minimumAge)),by=c('Title', 'year'),]$V1)==1)
#   })
# 
#   # sample one observation per individual per year and clean up dataset
#   subsampled <- dt[dt[ , .I[sample(.N,1)] , by = c('Title','year')]$V1]
#   meta <- subsampled[,c('Title', 'side', 'year', 'sex', 'minYear', 'maxYear', 'yearSpan', 'minimumAge', 'ageClass', 'yday'), ]
# 
#   return(meta)
# 
# }
# # dt <- tar_read(raw_photoData)
# # chosen_side = 'left'



# Build NBW metadata ------------------------------------------------------
# takes raw LV output as input ('dt')
build_nbw_metadata <- function(dt, chosen_side) {
  
  # Edits for 2024
  dt[Title==1317 & year(Date.Original)==2021, Title:=6018,]
  dt[Title==6402, Title:=6018, ]
  
  # # Manual edits for 2024
  # # 1. Several photos with erroneous metadata, including times 
  # photos_with_errors <- c('NBW_20240721_DSC_8272.NEF',
  #                         'NBW_20240721_DSC_8271.NEF',
  #                         'NBW_20240721_DSC_8267.NEF')
  # dt <- dt[!(File.name %in% photos_with_errors),,]
  # # 2. Fix to IDs 1317, 6018, 6402
  # dt[Title==1317 & year(Date.Original)==2021, Title:=6018,]
  # dt[Title==6402, Title:=6018, ]
  # # 3. Fix to photograph with erroneous ID
  # dt[File.name=='NBW_20240720_NB2_3650.JPG',Title:=6515,]
  
  # format dates
  dt[,dateTime:=as.POSIXct(Date.Original),]
  dt[,day:=format(dateTime, format = "%Y-%m-%d")]
  dt[,year:=year(day),]
  
  # remove photographs that were processed as crops
  dt <- dt[Title!="see crops",,]
  dt <- dt[Title!="unk",,]
  
  # DO RESIDENCY BEFORE restricting by side 
  dt[,residency:=ifelse(length(unique(year))==1, 'Transient', 'Resident'),by=Title]
  dt[residency=='Transient' & year %in% c(1988, 2024),residency:='Unknown',] # Changed to 2024, now most recent year of data
  
  # add left and right columns
  dt[grepl("Left", Keyword.export), side := "left"]
  dt[grepl("Right", Keyword.export), side := "right"]
  
  # subset to chosen side
  dt <- dt[side==chosen_side,,]
  
  # reliability
  dt[,reliability:='Not reliable',]
  dt[grepl('Notch|Back Indent|Nick', Keyword.export), reliability:='Reliable',by=c('Title','year')] # Changed to include Nick in 2025
  
  # add columns with sex details
  dt[, sex:=ifelse(grepl("Unknown", Keyword.export), "Unknown", NA), ]
  dt[, sex:=ifelse(grepl("FemaleJ", Keyword.export), "FemaleJ", sex), ]
  dt[, sex:=ifelse(grepl("Male", Keyword.export), "Male", sex), ]
  # old version below
  # dt=dt%>%mutate(sex = ifelse(grepl("FemaleJ", Keyword.export), "Female-Juvenile", NA))
  # dt=dt%>%mutate(sex = ifelse(grepl("Male", Keyword.export), "Male", sex))
  # dt$sex <- factor(dt$sex,levels=c("Female-Juvenile","Male"))
  
  # Step 1: Assign 'Calf' or 'Adult' based on 'Keyword.export'
  dt[, ageClass := ifelse(grepl('ageCalf', Keyword.export), 'Calf', 'Adult')]
  
  # Step 2: Further refine 'ageClass' for non-'Calf' rows to include 'Juvenile'
  dt[ageClass == 'Adult' & grepl('ageJuv', Keyword.export), ageClass := 'Juvenile']
  
  # Step 3: Create 'young' column indicating if the individual is either 'Calf' or 'Juvenile'
  dt[, young := ageClass %in% c('Calf', 'Juvenile')]
  
  # # add minimum age by year
  dt[,minYear:=min(year),by=Title]
  dt[,maxYear:=max(year),by=Title]
  dt[,yearSpan:=1+(maxYear-minYear),by=Title]
  dt[,catalogueAge:=year-minYear,] # Note: considering first year as "0"
  
  dt[, FirstYearAgeClass := unique(ageClass[minYear == year]), by = Title]
  
  dt[FirstYearAgeClass=='Calf',minimumAge:=catalogueAge,]
  dt[FirstYearAgeClass=='Juvenile',minimumAge:=catalogueAge+1,]
  dt[FirstYearAgeClass=='Adult',minimumAge:=catalogueAge+3,]
  
  dt[,minAge_in_2030:=(2030-year)+minimumAge,]
  
  # add day of year
  dt[,yday:=yday(day),]
  
  # diagnostic tests
  test_that("single sex classification for each individual", {
    expect_true(unique(dt[,length(unique(sex)),by=Title]$V1)==1)
  })
  
  # test_that("single minimum age for each individual in each year", {
  #   expect_true(unique(dt[,length(unique(minimumAge)),by=c('Title', 'year'),]$V1)==1)
  # })
  #
  # sample one observation per individual per year and clean up dataset
  subsampled <- dt[dt[ , .I[sample(.N,1)] , by = c('Title','year')]$V1]
  meta <- subsampled[,c('Title', 'side', 'year', 'sex', 'minYear', 'maxYear', 'yearSpan', 'ageClass', 'yday', 'minimumAge', 'minAge_in_2030', 'residency', 'reliability'), ]
  
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




# Calculate associations based on group observations ----------------------
build_group_associations <- function(g) {

  # convert to data.table for easier processing
  g <- data.table(g)

  # unique individuals
  individuals <- g[,unique(Title),]

  # generate all possible dyads
  possible_dyads <- data.table(t(combn(individuals, 2)))
  colnames(possible_dyads) <- c('A','B')

  # get unique group IDs
  group_ids <- unique(g$myGroups)
  
  # clear output
  output <- NULL

  # loop through and build dataframe of associations
  for (i in unique(group_ids)) {

    # individuals present in given group
    seen_individuals <- g[myGroups==i,unique(Title),]

    # possible dyads FOR each group (for computational efficiency, only including dyads where at least one individual is present)
    local_dyads <- possible_dyads[(A %in% seen_individuals) | (B %in% seen_individuals),, ]

    # add group ID
    local_dyads[,group_ID:=i,]

    # check whether A and B are present and calculate associations ('together')
    local_dyads[,A_present:=ifelse(A %in% seen_individuals,1,0),by=.I]
    local_dyads[,B_present:=ifelse(B %in% seen_individuals,1,0),by=.I]
    local_dyads[,together:=ifelse(A_present==1 & B_present==1,1,0),by=.I]

    # add dyad name
    local_dyads[,dyad:=paste(pmin(A,B),pmax(A,B),sep='_'),by=c('A','B')]

    # build up output iteratively
    if (exists('output')) (output <- rbindlist(list(output, local_dyads)))
    if (!(exists('output'))) (output <- local_dyads)

  }

  return(output)
  

  ###### test this carefully

}
# g <- tar_read(group_df_all_grouped)
# g <- g[g$tar_group==1,]




##################################################
# Network model with partial pooling across years
# [But year-specific dyads are kept separate]
unlist_group_associations <- function(group_associations) {
  
  # loop through list to extract and combine annual associations
  for (i in 1:length(group_associations)) {
    
    temp <- group_associations[[i]]
    temp$year_group <- i 
    
    if (exists('output')) output <- rbindlist(list(output, temp))
    if (!exists('output')) output <- temp
    
  }
  
  # add year-specific dyad IDs
  output[,dyad_annual:=paste(dyad, year_group, sep='-')]
  
  return(output)
  
  
}
# group_associations <- tar_read(group_associations)




# Create key to link year_groups with years  ------------------------------
year_key <- function(g) {
  
  key <- data.table(unique(g[,c('tar_group','year'),]))
  colnames(key) <- c('year_group','year')
  key <- key[order(key$year)]
  
}
# g <- tar_read(group_df_all_grouped)




# Combine group associations for summarizing ------------------------------
aggregate_group_associations <- function(g) {
  
  g[, sum_together:=sum(together), by=c('year_group','dyad')]
  g[, sum_opps:=.N, by=c('year_group','dyad')]
  
  aggregated <- unique(g[, c('A','B','sum_together','sum_opps','dyad', 'dyad_annual','year_group')])
  aggregated[,together:=sum_together,]
  
  return(aggregated)
  
}
# g <- tar_read(group_associations_all)




# Plot credible intervals of individual-specific slopes -------------------
plot_within_slopes <- function(fit) {
  
  # Extract coefficients
  ranef_slopes <- coef(fit, prob = c(0.05, 0.95))$Title  # Replace 'Title' with actual grouping factor
  
  # Extract the relevant coefficient for deltaMinAge
  delta_slopes <- ranef_slopes[,, 'deltaMinAge']  # Adjust if using a nested structure
  
  # Convert to dataframe
  df <- data.frame(
    Title = rownames(delta_slopes),  # Assuming individuals/groups as row names
    Estimate = delta_slopes[, "Estimate"],
    CI_low = delta_slopes[, "Q5"], 
    CI_high = delta_slopes[, "Q95"]
  )
  
  g <- ggplot(df, aes(x = reorder(Title, Estimate), y = Estimate)) +
    geom_pointrange(aes(ymin = CI_low, ymax = CI_high, color = (CI_low > 0 | CI_high < 0)), 
                    size = 0.2, linewidth=0.75) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth=0.1) +
    coord_flip() +
    scale_color_manual(values = c("grey", "red"), labels = c("Includes 0", "Excludes 0")) +
    labs(x = "ID", y = "Random Slope for deltaMinAge", 
         title = "Credible Intervals of Random Slopes for deltaMinAge",
         color = "CI Relationship to Zero") +
    theme_minimal()
  
  return(g)
  
}
# fit <- tar_read(ma_m2)



edge_diagnostics <- function(edgelist, associations, multi) {
  
  edgelist$year_group <- as.numeric(edgelist$year_group)
  combined <- merge(edgelist, associations, by=c('dyad', 'year_group'))
  
  combined[,never_together:=ifelse(together==0,'never together','together at least once'),]
  
  g <- ggplot(combined, aes(x=sum_opps, y=edge, group=never_together, color=never_together)) + 
    facet_wrap(~never_together)+
    geom_linerange(aes(ymin = lower_edge, ymax = upper_edge), linewidth = 1, color = 'yellow') +  # Add vertical error bars for CIs
    geom_jitter(alpha=0.1, color='grey30') + 
    geom_smooth(linewidth=2, se = FALSE) + 
    scale_color_manual(values=c('turquoise', 'navy'))+
    xlim(0,25)+
    theme_classic()
  
  return(g)
  
}
# edgelist <- tar_read(edges_group_ma)
# associations <- tar_read(group_associations_all_aggregated)
# multi <- TRUE




# Pull out edges
edge_list_multiAnnual <- function(fit, include_zeros) {
  
  # identify non-zero dyads
  m_data <- data.table(fit$data)
  m_data[,dyad:=dyad_annual,]
  m_data[,nAssociations:=sum(together),by=dyad_annual]
  nz_dyads <- m_data[nAssociations>=1,unique(dyad_annual),]
  
  # for multi-annual only
  nz_dyads <- str_extract(nz_dyads, "^\\d+_\\d+")
  
  # extract edge weights
  edges <- data.frame(coef(fit))
  edges$dyad <- rownames(edges)
  edges <- data.table(edges)
  edges[,edge:=inv_logit(dyad_annual.Estimate.Intercept),] ###### is this correct for the zero-inflated model as well?
  edges[,lower_edge:=inv_logit(dyad_annual.Q2.5.Intercept),] ###### is this correct for the zero-inflated model as well?
  edges[,upper_edge:=inv_logit(dyad_annual.Q97.5.Intercept),] ###### is this correct for the zero-inflated model as well?
  
  ###### think on whether I need to manage transformation of error some other way?
  
  edges <- merge(edges, unique(m_data[,c('dyad', 'nAssociations'),]), by='dyad')
  
  # for multi-annual
  edges[,year_group:=str_extract(dyad, "(?<=-)[0-9]+$"), by=.I]
  edges[,dyad:=str_extract(dyad, "^\\d+_\\d+"), by=.I]
  
  ###### so sd of transformed samples quite different from est. error -- need to figure out what's going on here and which we want to use
  # ask on brms forum or ask Jordan Hart directly
  
  # exclude non-zero dyads if necessary
  if (!include_zeros) (edges <- edges[dyad %in% nz_dyads,,])
  
  return(edges)
  
  
}
# fit <- tar_read(brms_fit_group_multiAnnual)
# include_zeros <- TRUE



print('Cleared functions')
