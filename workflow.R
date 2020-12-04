# Next steps
# 5. Describe two patient-cases using breakDown lib:
#    - One patient with just one visit to the hospital
#    - Another patient with the maximum number of visits in the testing dataset
# 6. Write conclusions

# R options
options(warn = -1, scipen = 999)
if(!require(pacman)){install.packages('pacman')}
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, vroom, psych,
                                caret, caretEnsemble, corrplot,
                                ranger, fastcluster, sparklyr,
                                fastDummies, ape, Boruta,
                                future, future.apply, furrr,
                                lsr, RColorBrewer, DT, skimr,
                                naniar, ggrepel, DALEX,
                                breakDown, pdp, png))

## --------------------------------------------------- ##
## Data obtention
## --------------------------------------------------- ##

if(!file.exists(paste0(getwd(),'/dataset_diabetes/diabetic_data.csv'))){
  # Download data
  url <- 'https://archive.ics.uci.edu/ml/machine-learning-databases/00296/dataset_diabetes.zip'
  download.file(url, destfile = paste0(getwd(),'/dataset_diabetes.zip'))
  # Unzip files
  unzip(paste0(getwd(),'/dataset_diabetes.zip'))
  file.remove(paste0(getwd(),'/dataset_diabetes.zip'))
  # Load the data
  tbl <- vroom::vroom(paste0(getwd(),'/dataset_diabetes/diabetic_data.csv'), delim = ',')
} else {
  # Load the data
  tbl <- vroom::vroom(paste0(getwd(),'/dataset_diabetes/diabetic_data.csv'), delim = ',')
}
# Replace '?' character by NA's
tbl <- tbl %>% dplyr::na_if(y = '?')

cat(paste0('Dataset dimensions\n'))
print(dim(tbl))

## --------------------------------------------------- ##
## Data pre-processing
## --------------------------------------------------- ##
# Identify and remove variables without or with low variance
zvar <- tbl %>% caret::nzv()
cat(paste0('Removed features due to the low or null variance\n'))
print(names(tbl)[zvar])
tbl  <- tbl[,-zvar]; rm(zvar)

# Replace the codes by the right categories to the *admission* variables
# Load ID identifiers
idi  <- vroom::vroom(paste0(getwd(),'/dataset_diabetes/IDs_mapping.csv'), delim = ',')
mtch <- which(is.na(idi$admission_type_id))
# admision_type_id
idi1 <- idi[1:(mtch[1]-1),]
# discharge_disposition_id
idi2 <- idi[(mtch[1]+1):(mtch[2]-1),]
names(idi2) <- as.character(idi2[1,])
idi2 <- idi2[-1,]
# admission_source_id
idi3 <- idi[(mtch[2]+1):nrow(idi),]
names(idi3) <- as.character(idi3[1,])
idi3 <- idi3[-1,]; rm(mtch)

tbl$admission_type_id <- factor(tbl$admission_type_id)
levels(tbl$admission_type_id) <- idi1$description
levels(tbl$admission_type_id)[levels(tbl$admission_type_id) %in% c('NULL','Not Available')] <- NA
# Create a new category which comprises the levels: 'Newborn','Trauma Center','Not Mapped'
levels(tbl$admission_type_id)[levels(tbl$admission_type_id) %in% c('Newborn','Trauma Center','Not Mapped')] <- 'Other'

tbl$discharge_disposition_id <- factor(tbl$discharge_disposition_id)
levels(tbl$discharge_disposition_id) <- idi2$description[match(levels(tbl$discharge_disposition_id), idi2$discharge_disposition_id)]
levels(tbl$discharge_disposition_id)[levels(tbl$discharge_disposition_id) %in% c('NULL','Not Available')] <- NA
# Create new categories: Discharged, Expired, Hospice, and Admitted
levels(tbl$discharge_disposition_id)[grep(pattern = '[dD][iI][sS][cC][hH][aA][rR][gG][eE][dD]', x = levels(tbl$discharge_disposition_id))] <- 'Discharged'
levels(tbl$discharge_disposition_id)[grep(pattern = '[eE][xX][pP][iI][rR][eE][dD]', x = levels(tbl$discharge_disposition_id))] <- 'Expired'
levels(tbl$discharge_disposition_id)[grep(pattern = '[hH][oO][sS][pP][iI][cC][eE]', x = levels(tbl$discharge_disposition_id))] <- 'Hospice'
levels(tbl$discharge_disposition_id)[c(2,3,5)] <- 'Admitted'

tbl$admission_source_id <- factor(tbl$admission_source_id)
levels(tbl$admission_source_id) <- idi3$description[match(levels(tbl$admission_source_id), idi3$admission_source_id)]
levels(tbl$admission_source_id)[levels(tbl$admission_source_id) %in% c('NULL','Not Available')] <- NA
# Create new categories: Transfer, Referral, and Other
levels(tbl$admission_source_id)[grep(pattern = '[tT][rR][aA][nN][sS][fF][eE][rR]', x = levels(tbl$admission_source_id))] <- 'Transfer'
levels(tbl$admission_source_id)[grep(pattern = '[rR][eE][fF][eE][rR][rR][aA][lL]', x = levels(tbl$admission_source_id))] <- 'Referral'
levels(tbl$admission_source_id)[4:8] <- 'Other'
rm(idi, idi1, idi2, idi3)

## Omit Unknown category in gender
# It only has 3 registers
cat(paste0('Omit Unknown category in gender due to low variation\n'))
print(table(tbl$gender))
tbl$gender[grep(pattern = 'Unknown', x = tbl$gender)] <- NA

## Identify variables with more than 25% of missing data
msg <- tbl %>%
  apply(X = ., MARGIN = 2, function(x){sum(is.na(x))/nrow(.)}) %>%
  sort(decreasing = T) %>%
  .[. > 0.25]
cat(paste0('Removed features due to high proportion of missing data\n'))
print(msg)

naniar::vis_miss(tbl, warn_large_data = F)

## Remove variables with more than 25% of missing data
tbl <- tbl[,-which(names(tbl) %in% names(msg))]; rm(msg)

## Create a new feature: How many times the same patient have visited the hospital?
# Select the ones with the maximum number of days interned in the hospital
visits <- tbl$patient_nbr %>% table %>% sort(decreasing = T) %>% base::as.data.frame()
names(visits)[1] <- 'patient_nbr'
visits$patient_nbr <- visits$patient_nbr %>% as.character() %>% as.numeric()
unq_vs <- visits %>% dplyr::filter(Freq == 1)
visits <- visits %>% dplyr::filter(Freq > 1)

tbl_unq <- tbl %>% dplyr::filter(patient_nbr %in% unq_vs$patient_nbr)
tbl_unq$number_visits <- 1

tbl_dup <- 1:nrow(visits) %>%
  purrr::map(.f = function(i){
    df <- tbl %>%
      dplyr::filter(patient_nbr == visits$patient_nbr[i]) %>%
      .[which.max(.$time_in_hospital)[1],]
    df$number_visits <- visits$Freq[i]
    return(df)
  }) %>%
  dplyr::bind_rows()

tbl <- rbind(tbl_dup, tbl_unq)
# Number of unique patients
cat(paste0('Unique number of patients after removing duplicated information\n'))
print(nrow(tbl))
rm(tbl_dup, tbl_unq, unq_vs, visits)

## Transform character to factors
tbl[sapply(tbl, is.character)] <- lapply(tbl[sapply(tbl, is.character)], as.factor)

## Create a new feature: numerical age using the midpoint of each age interval
tbl$age_num <- tbl$age %>%
  gsub('\\[', '', .) %>%
  gsub('\\)', '', .) %>%
  strsplit(., split = '-') %>%
  purrr::map(.f = function(int){
    x <- as.numeric(int)
    return(mean(x))
  }) %>%
  unlist

# The following features: diag_1, diag_2, and diag_3 make
# reference of the three initial detected diagnoses. They
# have more than 700 categories. Additionally, the feature
# 'number_diagnoses' captures the sum of all of detected
# diagnoses. For that reason, diag_1,diag_2, and diag_3
# were excluded of the analysis
tbl$diag_1 <- tbl$diag_2 <- tbl$diag_3 <- NULL

# Discard expired and hospiced people from analysis
tbl$discharge_disposition_id <- tbl$discharge_disposition_id %>% as.character
tbl <- tbl %>%
  dplyr::filter(!(discharge_disposition_id %in% c('Expired','Hospice')))
tbl$discharge_disposition_id <- tbl$discharge_disposition_id %>% factor()

## Response variable
# Considering that the readmission before and after 30 days
# are both bad situations for a patient. The analysis will
# be carried out using a binary classification, merging the
# admissions before and after 30 days in one category
tbl$readmitted <- ifelse(test = tbl$readmitted %in% c('<30','>30'),
                         yes  = 1,
                         no   = 0)
# Data proportion of response category
cat(paste0('Response variable proportion\n'))
print(table(tbl$readmitted)/nrow(tbl))

# Remove ID variables
tbl$encounter_id <- NULL
tbl$patient_nbr  <- NULL

# Transform categorical variables to dummies 
tbl_num <- tbl %>% fastDummies::dummy_cols(ignore_na = T)
names(tbl_num)
tbl_num <- tbl_num %>%
  dplyr::select(time_in_hospital:number_diagnoses,
                number_visits,age_num,
                race_AfricanAmerican:race_Other,
                gender_Female,gender_Male,
                admission_type_id_Emergency:admission_type_id_Other,
                discharge_disposition_id_Admitted:`discharge_disposition_id_Not Mapped`,
                admission_source_id_Referral:admission_source_id_Other,
                `A1Cresult_>7`:diabetesMed_Yes,
                readmitted)

# Based on the proportion of complete observations which is 82%
# the final decision is not impute missing data
cat('Proportion of complete observations\n')
print(nrow(tbl_num[complete.cases(tbl_num),])/nrow(tbl_num))

tbl_num_full <- tbl_num %>% tidyr::drop_na()

# Data proportion of response category
cat(paste0('Response variable proportion\n'))
print(table(tbl_num_full$readmitted)/nrow(tbl_num_full))

## Identify and remove categories without or with low variance
zvar <- tbl_num_full %>% caret::nzv()
cat('Dummy categories without or with low variance\n')
print(names(tbl_num_full)[zvar])
tbl_num_full <- tbl_num_full[,-zvar]; rm(zvar)

cat(paste0('Dataset dimensions\n'))
print(dim(tbl_num_full))

out_dir <- paste0(getwd(),'/processed_data')
if(!dir.exists(out_dir)){dir.create(out_dir, recursive = T)}

if(!file.exists(paste0(out_dir,'/diabetic_data_processed_1.csv'))){
  vroom::vroom_write(x = tbl, paste0(out_dir,'/diabetic_data_processed_1.csv'), delim = ',')
}
if(!file.exists(paste0(out_dir,'/diabetic_data_processed_2.csv'))){
  vroom::vroom_write(x = tbl_num, paste0(out_dir,'/diabetic_data_processed_2.csv'), delim = ',')
}
if(!file.exists(paste0(out_dir,'/diabetic_data_processed_2_complete.csv'))){
  vroom::vroom_write(x = tbl_num_full, paste0(out_dir,'/diabetic_data_processed_2_complete.csv'), delim = ',')
}
rm(out_dir)

## --------------------------------------------------- ##
## Exploratory Data Analysis
## --------------------------------------------------- ##

# This function brings a complete description of the dataset
# for both categorical and numerical variables
skimr::skim(tbl)

## --------------------------------------------------- ##
## Correlation Analysis
## --------------------------------------------------- ##

tbl_num_full <- vroom::vroom('./processed_data/diabetic_data_processed_2_complete.csv', delim = ',')

# This function approximates the correlation calculation for
# mixed variables (numerical and categorical). The following
# cases are considered:
# - Both variables numeric: non-parametric Spearman correlation
# - One variable categorical, the other one numerical: ANOVA
# test getting the R2 and applying the square root
# - Both variables categorical: Cramer's V test to measure
# association between two nominal variables
# Taken from https://gist.github.com/talegari/b514dbbc651c25e2075d88f31d48057b
source(paste0(getwd(),'/scripts/cor2_modified.R'))

m <- cor2(df = tbl_num_full)
corrplot::corrplot(m,
                   type   = "upper",
                   method = "square",
                   order  = "hclust",
                   col    = brewer.pal(n = 8, name = "RdBu"),
                   tl.cex	= .5)

# PCA analysis
if(FALSE){
  sc      <- sparklyr::spark_connect(master = 'local')
  tbl_spk <- sparklyr::copy_to(sc, tbl_num_full, overwrite = T)
  pca_mdl <- sparklyr::ml_pca(tbl_spk)
  print(pca_mdl)
  cc <- pca_mdl$pc %>% base::as.data.frame()
  plot(cc[,1:2], xlim = c(-1, 1), ylim = c(-1, 1))
  
  # Taken from: https://www.researchgate.net/deref/http%3A%2F%2Fdx.doi.org%2F10.13140%2FRG.2.2.33337.26720
  circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
    r = diameter / 2
    tt <- seq(0,2*pi,length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
  }
  circ <- circleFun(c(0,0),2,npoints = 500)
  vars.p <- ggplot() +
    geom_path(data = circ, aes(x, y), lty = 2, color = "grey", alpha = 0.7) +
    geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
    geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
    geom_segment(data = cc, aes(x = 0, xend = PC1, y = 0, yend = PC2), arrow = arrow(length = unit(0.025, "npc"), type = "open"), lwd = 0.5) +
    geom_text_repel(data = cc[(cc$PC1 < -0.1)|(cc$PC1 > 0.1)|(cc$PC2 < -0.1)|(cc$PC2 > 0.1),], aes(x = PC1 * 1.1, y = PC2 * 1.1, label = rownames(cc[(cc$PC1 < -0.1)|(cc$PC1 > 0.1)|(cc$PC2 < -0.1)|(cc$PC2 > 0.1),]), check_overlap = T, size = 1)) +
    xlab("PC 1") +
    ylab("PC 2") +
    coord_equal() +
    labs(size = NULL) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(fill = "transparent"),
          legend.position = "none")
  vars.p
}

## --------------------------------------------------- ##
## Feature selection
## --------------------------------------------------- ##

out_dir <- paste0(getwd(),'/results')
if(!dir.exists(out_dir)){dir.create(out_dir, recursive = T)}
if(!file.exists(paste0(out_dir,'/feature_selection_simulation.csv'))){
  
  # Generate a local cluster of processors
  future::plan(multiprocess, workers = future::availableCores()-1)
  # Define a seed to do the results replicable
  set.seed(1)
  # Generate 20 random seeds from a uniform distribution
  seeds <- round(runif(20) * 1000, 0)
  # Run Boruta algorithm of feature selection using 20
  # random subsets of 2000 patients information and
  # save final decision: features CONFIRMED, TENTATIVE,
  # and REJECTED
  selected_fts <- seeds %>%
    furrr::future_map(.f = function(seed){
      # Fix each seed
      set.seed(seed)
      # Obtain a random sample without replacement of 2000 observations
      smp <- sample(x = 1:nrow(tbl_num_full), size = 2000, replace = F)
      # Run Boruta algorithm of feature selection over this random sample
      fts <- Boruta::Boruta(formula = readmitted ~ ., data = tbl_num_full[smp,])
      # Get the final decision
      res <- fts$finalDecision
      return(res)
    })
  
  fts <- selected_fts %>% as.data.frame()
  fts$Feature <- rownames(fts)
  rownames(fts) <- 1:nrow(fts)
  colnames(fts)[1:20] <- paste0('Run_',1:20)
  
  write.csv(fts, paste0(out_dir,'/feature_selection_simulation.csv'), row.names = F)
} else {
  fts <- read.csv(paste0(out_dir,'/feature_selection_simulation.csv'))
}

fts %>%
  tidyr::pivot_longer(cols = Run_1:Run_20) %>%
  dplyr::group_by(Feature, value) %>%
  dplyr::summarise(Count = n()) %>%
  ggplot2::ggplot(aes(x = Feature, y = Count, fill = factor(value, levels = c('Confirmed','Tentative','Rejected')))) +
  ggplot2::scale_fill_brewer(palette = 'Set1', direction = -1) +
  ggplot2::geom_bar(stat = 'identity') +
  ggplot2::coord_flip() +
  ggplot2::labs(fill = 'Decision')

fnl_fst <- fts %>%
  tidyr::pivot_longer(cols = Run_1:Run_20) %>%
  dplyr::group_by(Feature, value) %>%
  dplyr::summarise(Count = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Feature) %>%
  dplyr::summarise(Decision = value[which.max(Count)]) %>%
  dplyr::filter(Decision %in% c('Confirmed','Tentative')) %>%
  dplyr::pull(Feature)
fnl_fst[1] <- gsub('\`','',fnl_fst[1])
cat('Selected features from Boruta algorithm\n')
print(fnl_fst)

tbl_num_full_sel <- tbl_num_full[,c(fnl_fst,'readmitted')]

## --------------------------------------------------- ##
## Classification analysis
## --------------------------------------------------- ##

if(!file.exists('./results/model_accuracies.csv')){
  
  # Spark connection
  if(!exists('sc')){
    sc <- sparklyr::spark_connect(master = 'local', version = '3.0.1')
  }
  tbl_num_full_sel_spk <- sparklyr::copy_to(sc, tbl_num_full_sel, overwrite = T)
  
  # Define partitions of the dataset: 70% for training, 30% for testing
  partitions <- tbl_num_full_sel_spk %>%
    sparklyr::sdf_partition(training = 0.7, test = 0.3, seed = 1099)
  
  # Fit a logistic regression model
  fit_log <- partitions$training %>% ml_logistic_regression(readmitted ~ .)
  # Evaluate one model manually
  # prd_log <- ml_predict(fit_log, partitions$test)
  # ml_binary_classification_evaluator(x = prd_log, label_col = 'readmitted', raw_prediction_col = 'prediction')
  
  # Fit a neural network
  fit_ann <- partitions$training %>% ml_multilayer_perceptron_classifier(readmitted ~ ., layers = c(13, 26, 26, 2))
  
  # Fit a decision tree model
  fit_dct <- partitions$training %>% ml_decision_tree(readmitted ~ ., type = 'classification')
  
  # Fit a random forest model
  fit_rdf <- partitions$training %>% ml_random_forest(readmitted ~ ., type = 'classification', num_trees = 100)
  
  # Fit a gradient boosting trees model
  fit_gbt <- partitions$training %>% ml_gradient_boosted_trees(readmitted ~ ., type = 'classification')
  
  metrics <- data.frame(Model = c('Logistic regression', 'Neural network','Decision tree', 'Random forest', 'Gradient boosting trees'),
                        Training = c(ml_evaluate(fit_log, partitions$training)$accuracy(),
                                     ml_evaluate(fit_ann, partitions$training)$Accuracy,
                                     ml_evaluate(fit_dct, partitions$training)$Accuracy,
                                     ml_evaluate(fit_rdf, partitions$training)$Accuracy,
                                     ml_evaluate(fit_gbt, partitions$training)$Accuracy),
                        Testing  = c(ml_evaluate(fit_log, partitions$test)$accuracy(),
                                     ml_evaluate(fit_ann, partitions$test)$Accuracy,
                                     ml_evaluate(fit_dct, partitions$test)$Accuracy,
                                     ml_evaluate(fit_rdf, partitions$test)$Accuracy,
                                     ml_evaluate(fit_gbt, partitions$test)$Accuracy))
  write.csv(metrics, './results/model_accuracies.csv', row.names = F)
} else {
  metrics <- read.csv('./results/model_accuracies.csv')
}

metrics %>%
  tidyr::pivot_longer(cols = 2:3, names_to = 'Partition', values_to = 'Accuracy') %>%
  ggplot2::ggplot(aes(x = reorder(Model, Accuracy), y = Accuracy, fill = Partition)) +
  ggplot2::geom_bar(stat = "identity", position = position_dodge()) +
  ggplot2::scale_fill_brewer(palette = "Paired") +
  ggplot2::coord_flip() +
  ggplot2::theme_bw() +
  ggplot2::xlab('')

cat('Based on the accuracy metric the best model is the Gradient boosting trees\n')

if(FALSE){
  # Interpretation of the results
  test_set <- partitions$test %>% dplyr::collect()
  custom_predict <- function(model, new_data){
    
    new_data_spark <- copy_to(sc, new_data, name="spark_temp1b3c4e6", overwrite = T)
    spark_tbl_scored <- ml_predict(model, new_data_spark)
    res <- as.numeric(as.data.frame(spark_tbl_scored %>% select(prediction) %>% collect())$prediction)
    # dplyr::db_drop_table(con = sc, table = "spark_temp1b3c4e6")
    
    return(res)
  }
  
  # Explain
  xpl_spark_gbt <- DALEX::explain(model            = fit_gbt,
                                  data             = test_set %>% select(-readmitted),
                                  y                = test_set$readmitted,
                                  predict_function = custom_predict,
                                  label            = 'Gradient Boosting Trees model',
                                  type             = 'classification')
  
}

if(!file.exists('./results/gbt_variable_importance.png')){
  # Variable importance
  vim_spark_gbt <- DALEX::variable_importance(xpl_spark_gbt)
  plot(vim_spark_gbt)
} else {
  vimp <- './results/gbt_variable_importance.png'
  knitr::include_graphics(vimp)
}

if(!file.exists('./results/pdp_gbt.png')){
  # Predictor-response relationship
  pdp_gbt  <- DALEX::variable_effect_partial_dependency(xpl_spark_gbt,
                                                        variables = c("number_visits",
                                                                      "number_inpatient",
                                                                      "number_diagnoses",
                                                                      "time_in_hospital",
                                                                      "number_emergency",
                                                                      "age_num"))
  plot(pdp_gbt)
} else {
  pdp_m <- './results/pdp_gbt.png'
  knitr::include_graphics(pdp_m)
}

if(!(file.exists('./results/prediction_decomposition_obs1.png') &
     file.exists('./results/prediction_decomposition_obs2.png'))){
  # Prediction understanding
  row_00001 <- test_set[1,]
  row_06560 <- test_set[6560,]
  pd_spark_gbt1 <- breakDown::break_down(xpl_spark_gbt, new_observation = row_00001)
  pd_spark_gbt2 <- breakDown::break_down(xpl_spark_gbt, new_observation = row_06560)
  plot(pd_spark_gbt1)
  plot(pd_spark_gbt2)
} else {
  bd1 <- './results/prediction_decomposition_obs1.png'
  knitr::include_graphics(bd1)
  bd2 <- './results/prediction_decomposition_obs2.png'
  knitr::include_graphics(bd2)
}

if(exists('sc')){
  spark_disconnect(sc) 
}
