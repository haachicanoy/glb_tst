## Steps to follow

# 8. Correlation analysis among all variables:
# 9. Having a good database I have to:
# 	- Define how to deal with missing data: imputing is the best option
#   - Define how to balance the unbalanced class: downsampling, upsampling, mixed approach
#   - Identify and check outliers
#   - Classification model: 
#   - Clustering analysis: k-means

# R options
options(warn = -1, scipen = 999)
if(!require(pacman)){install.packages('pacman')}
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, vroom, psych,
                                caret, caretEnsemble, corrplot,
                                ranger, fastcluster, sparklyr,
                                fastDummies, ape, Boruta,
                                future, future.apply, furrr,
                                cor2))

## --------------------------------------------------- ##
## Obtain the data
## --------------------------------------------------- ##
if(!file.exists(paste0(getwd(),'/dataset_diabetes/diabetic_data.csv'))){
  ## Download data
  url <- 'https://archive.ics.uci.edu/ml/machine-learning-databases/00296/dataset_diabetes.zip'
  download.file(url, destfile = paste0(getwd(),'/dataset_diabetes.zip'))
  ## Unzip files
  unzip(paste0(getwd(),'/dataset_diabetes.zip'))
  file.remove(paste0(getwd(),'/dataset_diabetes.zip'))
  ## Load the data
  tbl <- vroom::vroom(paste0(getwd(),'/dataset_diabetes/diabetic_data.csv'), delim = ',')
} else {
  ## Load the data
  tbl <- vroom::vroom(paste0(getwd(),'/dataset_diabetes/diabetic_data.csv'), delim = ',')
}
## Replace '?' character by NA's
tbl <- tbl %>% dplyr::na_if(y = '?')

cat(paste0('The dataset is conformed by: ',dim(tbl), ' records vs features\n'))

## --------------------------------------------------- ##
## Data pre-processing
## --------------------------------------------------- ##
## Identify and remove variables without or with low variance
zvar <- tbl %>% caret::nzv()
cat(paste0('The following features were removed: ',names(tbl)[zvar], ' due to the low or insignificant variance\n'))
tbl  <- tbl[,-zvar]; rm(zvar)

## Replace the codes by the right categories to the *admission* variables
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
cat(paste0('Omit Unknown category in gender. The distribution of this feature is as follows\n'))
table(tbl$gender)
tbl$gender[grep(pattern = 'Unknown', x = tbl$gender)] <- NA

## Identify variables with more than 25% of missing data
msg <- tbl %>%
  apply(X = ., MARGIN = 2, function(x){sum(is.na(x))/nrow(.)}) %>%
  sort(decreasing = T) %>%
  .[. > 0.25]
cat(paste0('The following features were removed: ',names(msg), ' due to high amount of missing data\n'))

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
cat(paste0('The number of unique patients after removing duplicated information is: ',nrow(tbl), '\n'))
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
cat(paste0('The data proportion in the response variable is as follows: \n'))
table(tbl$readmitted)/nrow(tbl)

# Remove ID variables
tbl$encounter_id <- NULL
tbl$patient_nbr  <- NULL

# Quick summary of the data
summary(tbl)

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

# Proportion of complete observations
nrow(tbl_num[complete.cases(tbl_num),])/nrow(tbl_num)

tbl_num_full <- tbl_num %>% tidyr::drop_na()

# Data proportion of response category
cat(paste0('The data proportion in the response variable, leaving just complete data is as follows: \n'))
table(tbl_num_full$readmitted)/nrow(tbl_num_full)

## Identify and remove categories without or with low variance
zvar <- tbl_num_full %>% caret::nzv()
tbl_num_full <- tbl_num_full[,-zvar]; rm(zvar)

out_dir <- paste0(getwd(),'/processed_data')
if(!dir.exists(out_dir)){dir.create(out_dir, recursive = T)}

vroom::vroom_write(x = tbl, paste0(out_dir,'/diabetic_data_processed_1.csv'), delim = ',')
vroom::vroom_write(x = tbl_num, paste0(out_dir,'/diabetic_data_processed_2.csv'), delim = ',')
vroom::vroom_write(x = tbl_num_full, paste0(out_dir,'/diabetic_data_processed_2_complete.csv'), delim = ',')

## --------------------------------------------------- ##
## Exploratory Data Analysis
## --------------------------------------------------- ##

cor2


tbl_num_full_sp <- tbl_num_full[,-ncol(tbl_num_full)]
tbl_num_less_30 <- tbl_num_full[tbl_num_full$readmitted == '<30',]
tbl_num_less_30$readmitted <- NULL
tbl_num_less_30 <- tbl_num_less_30[,-caret::nzv(tbl_num_less_30)]
tbl_num_more_30 <- tbl_num_full[tbl_num_full$readmitted == '>30',]
tbl_num_more_30$readmitted <- NULL
tbl_num_more_30 <- tbl_num_more_30[,-caret::nzv(tbl_num_more_30)]
tbl_num_no_rdms <- tbl_num_full[tbl_num_full$readmitted == 'NO',]
tbl_num_no_rdms$readmitted <- NULL
tbl_num_no_rdms <- tbl_num_no_rdms[,-caret::nzv(tbl_num_no_rdms)]

cmat_less30 <- tbl_num_less_30 %>%
  dplyr::select(names(tbl_num_more_30)) %>%
  cor(method = 'spearman')
cmat_more30 <- tbl_num_more_30 %>%
  cor(method = 'spearman')
cmat_nordsm <- tbl_num_no_rdms %>%
  dplyr::select(names(tbl_num_more_30)) %>%
  cor(method = 'spearman')

mantel.test(cmat_less30, cmat_more30, nperm = 999, graph = T)
mantel.test(cmat_less30, cmat_nordsm, nperm = 999, graph = T)
mantel.test(cmat_more30, cmat_nordsm, nperm = 999, graph = T)

cmat_less30 %>%
  corrplot()

## --------------------------------------------------- ##
## Feature selection
## --------------------------------------------------- ##

future::plan(multiprocess, workers = future::availableCores()-1)

set.seed(1)
seeds <- round(runif(20) * 1000, 0)

selected_fts <- seeds %>%
  furrr::future_map(.f = function(seed){
    set.seed(seed)
    smp <- sample(x = 1:nrow(tbl_num_full), size = 2000, replace = F)
    fts <- Boruta::Boruta(formula = readmitted ~ ., data = tbl_num_full[smp,])
    res <- fts$finalDecision
    return(res)
  })

# set.seed(1)
# seeds <- round(runif(20) * 1000, 0)
# selected_fts <- seeds %>%
#   purrr::map(.f = function(seed){
#     set.seed(seed)
#     smp <- sample(x = 1:nrow(tbl_num_full_bin), size = 2000, replace = F)
#     fts <- Boruta::Boruta(formula = readmitted ~ ., data = tbl_num_full_bin[smp,])
#     res <- fts$finalDecision
#     return(res)
#   })

fts <- selected_fts %>% as.data.frame()
fts$Feature <- rownames(fts)
rownames(fts) <- 1:nrow(fts)
colnames(fts)[1:20] <- paste0('Run_',1:20)
fts %>%
  tidyr::pivot_longer(cols = Run_1:Run_20) %>%
  dplyr::group_by(Feature, value) %>%
  dplyr::summarise(Count = n()) %>%
  dplyr::filter(!(value == 'Rejected' & Count == 20)) %>%
  ggplot2::ggplot(aes(x = reorder(Feature, -Count), y = Count, fill = value)) +
  ggplot2::geom_bar(stat = 'identity') +
  ggplot2::coord_flip()

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

tbl_num_sel_bin <- tbl_num_full_bin[,c(fnl_fst,'readmitted')]

## Spark connection

sc <- sparklyr::spark_connect(master = 'local')

tbl_spk    <- sparklyr::copy_to(sc, tbl_num_full_bin, overwrite = T)
partitions <- tbl_spk %>%
  sparklyr::sdf_partition(training = 0.7, test = 0.3, seed = 1099)

tbl_fnl_spk <- sparklyr::copy_to(sc, tbl_num_sel_bin, overwrite = T)
partitions2 <- tbl_fnl_spk %>%
  sparklyr::sdf_partition(training = 0.7, test = 0.3, seed = 1099)

# Perform PCA
pca_model <- sparklyr::ml_pca(partitions$training)
print(pca_model)
pca_model$model$pc %>% View()

tbl_less_30_spk <- sparklyr::copy_to(sc, tbl_num_less_30, overwrite = T)
tbl_more_30_spk <- sparklyr::copy_to(sc, tbl_num_more_30, overwrite = T)
tbl_no_rdms_spk <- sparklyr::copy_to(sc, tbl_num_no_rdms, overwrite = T)

pca_less_30 <- sparklyr::ml_pca(tbl_less_30_spk)
pca_less_30$explained_variance %>% head(4)
r_less_30 <- pca_less_30$pc %>% as.data.frame %>% .[,1:5]
r_less_30$variable <- rownames(r_less_30)
r_less_30$readmitted <- '<30'

pca_more_30 <- sparklyr::ml_pca(tbl_more_30_spk)
pca_more_30$explained_variance %>% head(4)
r_more_30 <- pca_more_30$pc %>% as.data.frame %>% .[,1:5]
r_more_30$variable <- rownames(r_more_30)
r_more_30$readmitted <- '>30'

pca_no_rdms <- sparklyr::ml_pca(tbl_no_rdms_spk)
pca_no_rdms$explained_variance %>% head(4)
r_no_rdms <- pca_no_rdms$pc %>% as.data.frame %>% .[,1:5]
r_no_rdms$variable <- rownames(r_no_rdms)
r_no_rdms$readmitted <- 'NO'

r_pca <- rbind(r_less_30, r_more_30, r_no_rdms)

r_pca %>%
  ggplot(aes(x = reorder(variable, -PC2), y = PC2)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~readmitted)

# Fit a linear model to the training dataset
fit <- partitions$training %>%
  ml_logistic_regression(readmitted ~ .)

fit

pred <- ml_predict(fit, partitions$test)
ml_binary_classification_eval(pred)

# Fit a linear model to the training dataset selected features
fit2 <- partitions2$training %>%
  dplyr::select(-number_visits) %>%
  ml_logistic_regression(readmitted ~ .)

fit2$coefficients

pred2 <- ml_predict(fit2, partitions2$test)
ml_binary_classification_eval(pred2)

# Fit a GBM model to the training dataset selected features
gbm_fit <- partitions2$training %>%
  sparklyr::ml_gradient_boosted_trees(readmitted ~ .)

pred3 <- ml_predict(gbm_fit, partitions2$test)
ml_binary_classification_eval(x = pred3, label_col = 'readmitted', prediction_col = "prediction")

ml_binary_classification_evaluator(x = pred2, label_col = 'readmitted', raw_prediction_col = 'prediction')
ml_binary_classification_evaluator(x = pred3, label_col = 'readmitted', raw_prediction_col = 'prediction')


ml_multiclass_classification_evaluator(pred)


pat_no <- tbl_cmp %>%
  dplyr::filter(readmitted == 'NO')
pat_no$readmitted <- NULL
dst_no <- pat_no %>% daisy(x = ., metric="gower")

pat_less_30 <- tbl_cmp %>%
  dplyr::filter(readmitted == '<30')
pat_less_30$readmitted <- NULL


pat_more_30 <- tbl_cmp %>%
  dplyr::filter(readmitted == '>30')
pat_more_30$readmitted <- NULL


fastcluster::hclust()


spark_disconnect(sc)

# tst <- tbl[complete.cases(tbl),] %>% base::as.data.frame()
# # cor2(df = tst)
# 
# # Numeric variables
# colnames(tbl[sapply(tbl, is.numeric)])
# tbl[sapply(tbl, is.numeric)] %>%
#   cor(use = 'pairwise.complete.obs', method = 'spearman') %>%
#   corrplot::corrplot.mixed()
# 
# # Categorical variables
# cat_vars <- tbl[sapply(tbl, is.factor)]
# cat_vars <- cat_vars[complete.cases(cat_vars),]
# cat_vars <- cat_vars %>% as.data.frame()
# p.chisq <- matrix(0, nrow = ncol(cat_vars), ncol = ncol(cat_vars), byrow = T)
# for(i in 1:ncol(cat_vars)){
#   for(j in 1:ncol(cat_vars)){
#     p.chisq[i,j] <- round(chisq.test(cat_vars[,i],cat_vars[,j])$p.value,3)
#   }
# }; rm(i); rm(j)
# 
# diag(p.chisq) <- 1
# colnames(p.chisq) <- colnames(cat_vars)
# rownames(p.chisq) <- colnames(cat_vars)
# 
# ord <- p.chisq %>%
#   corrplot::corrMatOrder(order = 'AOE')
# 
# p.chisq2 <- p.chisq[ord,ord]
# corrplot::corrplot.mixed(p.chisq2)
# 
# # Select just complete data
# tbl_cmp <- tbl[complete.cases(tbl),]
# ids     <- tbl_cmp[,1:2]
# tbl_cmp <- tbl_cmp[,-c(1:2)]
# 
# # modelr::crossv_kfold()
# train_control<- caret::trainControl(method = "cv",
#                                     number = 10,
#                                     savePredictions = TRUE)
# 
# library(doParallel)
# cl <- makePSOCKcluster(3)
# registerDoParallel(cl)
# 
# model <- caret::train(readmitted ~ ., data=tbl_cmp, trControl=train_control, method = "gbm")
# 
# stopCluster(cl)
# 
# set.seed(1235)
# inTrain <- caret::createDataPartition(y = tbl_cmp$readmitted, p = 0.7, list = FALSE)
# training <- tbl_cmp[inTrain,]
# testing  <- tbl_cmp[-inTrain,]
# 
# control_prmt <- caret::trainControl(method          = "LGOCV",
#                                     p               = 0.7,
#                                     number          = 10,
#                                     savePredictions = "final",
#                                     verboseIter     = TRUE)
# 
# model_list <- caretEnsemble::caretList(
#     readmitted ~ .,
#     data       = training,
#     trControl  = control_prmt,
#     tuneList   = list(ranger = caretModelSpec(method = "ranger", importance = "impurity")),
#     methodList = c("svmRadial", "gbm", "cforest", "hdda", "xgbTree", "xgbLinear")
#   )
```
