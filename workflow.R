## Steps to follow

# 7. Remove ID variables from analysis:
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
suppressMessages(pacman::p_load(tidyverse, vroom, psych, caret, caretEnsemble, corrplot, ranger, fastcluster, sparklyr))

## Obtain the data
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

# psych::describe(tbl)

## Data pre-processing

## Identify and remove variables without variance
zvar <- tbl %>% caret::nzv()
tbl  <- tbl[,-zvar]

## Replace the codes by the right categories to the *admision* variables
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
idi3 <- idi3[-1,]

tbl$admission_type_id <- factor(tbl$admission_type_id)
levels(tbl$admission_type_id) <- idi1$description
levels(tbl$admission_type_id)[levels(tbl$admission_type_id) %in% c('NULL','Not Available')] <- NA
# Create a new category which comprises the levels: 'Newborn','Trauma Center','Not Mapped'
levels(tbl$admission_type_id)[levels(tbl$admission_type_id) %in% c('Newborn','Trauma Center','Not Mapped')] <- 'Other'

tbl$discharge_disposition_id <- factor(tbl$discharge_disposition_id)
levels(tbl$discharge_disposition_id) <- idi2$description[match(levels(tbl$discharge_disposition_id), idi2$discharge_disposition_id)]
levels(tbl$discharge_disposition_id)[levels(tbl$discharge_disposition_id) %in% c('NULL','Not Available')] <- NA
grep(pattern = 'discharged ', x = levels(tbl$discharge_disposition_id))

tbl$admission_source_id <- factor(tbl$admission_source_id)
levels(tbl$admission_source_id) <- idi3$description[match(levels(tbl$admission_source_id), idi3$admission_source_id)]
levels(tbl$admission_source_id)[levels(tbl$admission_source_id) %in% c('NULL','Not Available')] <- NA
rm(idi, idi1, idi2, idi3)

## Identify variables with more than 15% of missing data
msg <- tbl %>%
  apply(X = ., MARGIN = 2, function(x){sum(is.na(x))/nrow(.)}) %>%
  sort(decreasing = T) %>%
  .[. > 0.15]
print(msg)

# Remove variables with more than 15% of missing data
tbl <- tbl[,-which(names(tbl) %in% names(msg))]

# Check the confusion matrix
# table(tbl$race, tbl$readmitted) %>%
#   base::data.frame() %>%
#   ggplot(aes(x = Var1, y = Var2, fill = Freq)) +
#   geom_tile() +
#   coord_equal() +
#   scale_fill_gradientn(name = "", colors = terrain.colors(10)) +
#   scale_x_discrete(name = "") +
#   scale_y_discrete(name = "")

# How many times the same patient visited the hospital
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
rm(tbl_dup, tbl_unq, unq_vs, visits)

## Transform character to factors
tbl[sapply(tbl, is.character)] <- lapply(tbl[sapply(tbl, is.character)], as.factor)

## Create numerical age variables using the midpoint of each age interval
tbl$age_num <- tbl$age %>%
  gsub('\\[', '', .) %>%
  gsub('\\)', '', .) %>%
  strsplit(., split = '-') %>%
  purrr::map(.f = function(int){
    x <- as.numeric(int)
    return(mean(x))
  }) %>%
  unlist

# Response variable
tbl$readmitted %in% c('<30','>30')

## Exploratory Data Analysis

summary(tbl)
tbl$encounter_id <- NULL
tbl$patient_nbr  <- NULL

tbl$gender[grep(pattern = 'Unknown', x = tbl$gender)] <- NA

tbl$discharge_disposition_id %>% table() %>% sort(decreasing = T) %>% View()
tbl$admission_source_id %>% table() %>% sort(decreasing = T)
tbl$diag_1 %>% unique %>% length()
tbl$diag_2 %>% unique %>% length()
tbl$diag_3 %>% unique %>% length()

## Correlation analysis

tbl_cmp <- tbl[complete.cases(tbl),]

sc <- sparklyr::spark_connect(master = 'local')

tbl_spk <- sparklyr::copy_to(sc, tbl_cmp)
partitions <- tbl_spk %>%
  sparklyr::sdf_partition(training = 0.7, test = 0.3, seed = 1099)

# fit a linear model to the training dataset
fit <- partitions$training %>%
  ml_random_forest(readmitted ~ ., type = 'classification')
fit

pred <- ml_predict(fit, partitions$test)

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
