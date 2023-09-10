#install packages

install.packages('tidyverse')
install.packages('httr')
install.packages('jsonlite')
install.packages('highcharter')
install.packages('dplyr')
install.packages('lubridate')

#open packages
library(tidyverse)
library(httr)
library(jsonlite)
library(highcharter)
library(dplyr)
library(lubridate)

#Load dataframes
df1 <- read_tsv("prad_tcga_pan_can_atlas_2018_clinical_data.tsv")
df2 <- read_tsv("prad_tcga_pub_clinical_data.tsv")
df3 <- read_csv("tcga_prad_from_xml.csv")

#Select key variables from dataframes
table1 <- df1 %>%
  select(2,4,7,15,16,21,31,40,44,45)

table2 <- df2 %>%
  select(2,4,11,22,25,40,41,42,49,50,54,70,74,78,88,89,92)

table3 <- df3 %>% 
  select(14,16,20,22,45,48,57,58,64,65,67,71:89)

# Join tables
join1 <- left_join(table1, table2, by = "Patient ID")
finaltable <- right_join(join1, table3, by = "Patient ID")

# Rename columns
finaltable <- finaltable %>% rename(
  patient_id = `Patient ID`, diagnosis_age = `Diagnosis Age`,
  disease_free_status = `Disease Free Status`)

finaltable <- finaltable %>% rename(
  erg_status = `Erg status`, etv1_status = `ETV1 Status`,
  etv4_status = `ETV4 Status`, foxa1_mutation = `FOXA1 mutation`,
  spop_mutation = `SPOP Mutation`, pten_mutation = `PTEN Mutation`
)

# Filter values TRUE for radiotherapy and hormone therapy
filtered_data <- finaltable %>%
  filter(radiotherapy == TRUE, hormone_therapy == TRUE)

#Select key variables from filtered data
fd1 <- filtered_data %>%
  select(1,2,3,6,7,8,11,12,15,23,24,26,29,31,32,33,35,37,46,54,55)

# Clean variable names
install.packages('Janitor')
library(janitor)

Clean_Data <- fd1 %>%
  clean_names()

# Edit column values
Clean_Data$bone_scan_results[
  Clean_Data$bone_scan_results == 'Normal (no evidence of prostate cancer) [cM0]'
] <- 'Normal'

Clean_Data$bone_scan_results[
  Clean_Data$bone_scan_results == 'Abnormal (not related to prostate cancer)'
] <- 'Abnormal'

# Convert outcome variable to factor with levels "No" and "Yes"
pfs_factor <- factor(Clean_Data$pfs_status, levels = c(1, 0), labels = c("Yes", "No"))
pfs <- pfs_factor


#Exploratory Data Analysis
install.packages('visdat')
library(visdat)
install.packages('ggplot2')
library(ggplot2)

# visualise final dataset
vis_dat(Clean_Data)
vis_dat(Clean_Data, sort_type = FALSE)

# explore missing data
vis_miss(Clean_Data)
vis_miss(Clean_Data, sort_miss = TRUE)

# Create bar plot
ggplot(Clean_Data, aes(x=race_list)) +
  geom_bar()

Clean_Data %>% count(race_list)

ggplot(Clean_Data, aes(x=diagnosis_age)) +
  geom_histogram(binwidth = 5)

Clean_Data %>% 
  ggplot(aes(diagnosis_age, aneuploidy_score, color = pfs)) +
  geom_point() +
  facet_wrap(~race_list)

Clean_Data %>% count(diagnosis_age)

Clean_Data %>% 
  ggplot(aes(diagnosis_age, reviewed_gleason_sum, color = pfs)) +
  geom_point() +
  facet_wrap(~race_list)

# Select Numeric Variables
numeric_data <- select_if(Clean_Data, is.numeric)

# Compute the correlation matrix for numeric variables
cor_matrix <- cor(numeric_data)

# Print the correlation matrix
print(cor_matrix)

#convert categorical variables to factor and relevel outcome variable
data_new <- 
  Clean_Data %>% 
  select(
    patient_id:biochemical_recurrence,
    pfs_status
  ) %>% 
  mutate(
    across(where(is_character), as_factor),
    pfs_status = factor(pfs_status, levels = c("0", "1"))
  )


# Set the random number stream using `set.seed()` so that the results can be reproduced later
set.seed(123)

# Save the split information for a 75/25 split of the data
prad_split <- initial_split(data_new, prop = 0.75, strata = pfs_status)


# Generate resulting dataset from the intial split
prad_train <- training(prad_split)
prad_test  <- testing(prad_split)

# Bootstrap resampling on the training set
set.seed(234)
prad_boot <- bootstraps(prad_train, times = 10)

# Logistic regression model specification
glm_model <- logistic_reg() %>%
  set_engine("glm") %>% 
  set_mode("classification")

#create recipe
recipe_object <- 
  recipe(pfs_status ~ ., data = prad_train) %>% 
  update_role(patient_id, new_role = "label") %>% 
  step_normalize(all_numeric_predictors()) %>% 
  step_corr(all_numeric_predictors()) %>%
  step_impute_knn(all_predictors()) %>% 
  step_dummy(all_nominal_predictors()) %>% 
  step_smote(pfs_status)

#recipe post variable importance
recipe_object <- 
  recipe(pfs_status ~ patient_id + residual_tumor +
           number_of_lymphnodes_positive_by_he +
           spop_mutation + fraction_genome_altered + 
           tp53_mutation + biochemical_recurrence, data = prad_train) %>% 
  update_role(patient_id, new_role = "label") %>% 
  step_normalize(all_numeric_predictors()) %>% 
  step_corr(all_numeric_predictors()) %>%
  step_impute_knn(all_predictors()) %>% 
  step_dummy(all_nominal_predictors()) %>% 
  step_smote(pfs_status)

# Set up workflow
prad_wf <- workflow() %>%
  add_model(glm_model) %>% 
  add_recipe(recipe_object)

prad_wf

#fit model to bootstrap resamples
glm_rs <- 
  prad_wf %>%
  fit_resamples(
    resamples = prad_boot, 
    control = control_resamples(save_pred = TRUE, verbose = TRUE),
    metrics = metric_set(accuracy, sens, spec, roc_auc)
  )

#Evaluate logistic regression model

collect_metrics(glm_rs)

collect_predictions(glm_rs)

# confusion matrix
glm_rs %>% 
  collect_predictions() %>% 
  conf_mat(
    truth = pfs_status,
    estimate = .pred_class
  ) %>%
  autoplot(type = "heatmap")

#roc_curve
glm_rs %>%
  collect_predictions() %>%
  roc_curve(pfs_status, .pred_0) %>%
  autoplot()

# variable importance lr --------------------------------------------------

# get the prepped data set for assessing variable importance
vip_data <- 
  recipe_object %>% 
  prep() %>% 
  bake(new_data = NULL) %>% 
  select(!patient_id)

# fit model to the prepped training data NOT resamples
glm_fit <- 
  glm_model %>% 
  fit(pfs_status ~ ., data = vip_data)

# get the coefficients and model stats
tidy(glm_fit)

#filter out intercept
tidyintercept <- tidy(glm_fit) %>%
  filter(term != "(Intercept)")

print(tidyintercept)

library(vip)

# plot variable importance
glm_fit %>% 
  vip(geom = "point")

install.packages("randomForest")
library(randomForest)
install.packages("ranger")
library(ranger)

# Random Forest model specification
rf_model <- 
  rand_forest() %>%
  set_engine("ranger") %>% 
  set_mode("classification")

recipe_object <- 
  recipe(pfs_status ~ ., data = prad_train) %>% 
  update_role(patient_id, new_role = "label") %>% 
  step_normalize(all_numeric_predictors()) %>% 
  step_corr(all_numeric_predictors()) %>%
  step_impute_knn(all_predictors()) %>% 
  step_dummy(all_nominal_predictors()) %>% 
  step_smote(pfs_status)

#post variable importance
recipe_object <- 
  recipe(pfs_status ~ patient_id + number_of_lymphnodes_positive_by_he +
           ar_score + biochemical_recurrence +
           spop_mutation + erg_status +
           person_neoplasm_cancer_status, data = prad_train) %>% 
  update_role(patient_id, new_role = "label") %>% 
  step_normalize(all_numeric_predictors()) %>% 
  step_corr(all_numeric_predictors()) %>%
  step_impute_knn(all_predictors()) %>% 
  step_dummy(all_nominal_predictors()) %>% 
  step_smote(pfs_status) 

#rf model workflow
prad_wf <- workflow() %>%
  add_model(rf_model) %>% 
  add_recipe(recipe_object)

prad_wf

rf_rs <- prad_wf %>%
  update_model(rf_model) %>%
  fit_resamples(
    resamples = prad_boot,
    control = control_resamples(save_pred = TRUE, verbose = TRUE),
    metrics = metric_set(accuracy, sens, spec, roc_auc)
  )

#Evaluate Random Forest model

collect_metrics(rf_rs)

collect_predictions(rf_rs)

# confusion matrix
rf_rs %>% 
  collect_predictions() %>% 
  conf_mat(
    truth = pfs_status,
    estimate = .pred_class
  ) %>%
  autoplot(type="heatmap")

#roc_curve
rf_rs %>%
  collect_predictions() %>%
  roc_curve(pfs_status, .pred_0) %>%
  autoplot()

# variable importance rf --------------------------------------------------
install.packages("vip")
library(vip)

# get the prepped data set for assessing variable importance
vip_data <- 
  recipe_object %>% 
  prep() %>% 
  bake(new_data = NULL) %>% 
  select(!patient_id)

rf_model %>% 
  set_engine("ranger", importance = "permutation") %>% 
  fit(pfs_status ~ ., data = vip_data) %>% 
  vip(geom = "point")

install.packages("xgboost")
library(xgboost)

# xgboost model specification
xgb_model <- 
  boost_tree() %>%
  set_engine("xgboost") %>% 
  set_mode("classification")


recipe_object <- 
  recipe(pfs_status ~ ., data = prad_train) %>% 
  update_role(patient_id, new_role = "label") %>% 
  step_normalize(all_numeric_predictors()) %>% 
  step_corr(all_numeric_predictors()) %>%
  step_impute_knn(all_predictors()) %>% 
  step_dummy(all_nominal_predictors()) %>% 
  step_smote(pfs_status)

#post variable-importance
recipe_object <- 
  recipe(pfs_status ~ patient_id + number_of_lymphnodes_positive_by_he +
           ar_score + mutation_count +
           tp53_mutation, data = prad_train) %>% 
  update_role(patient_id, new_role = "label") %>% 
  step_normalize(all_numeric_predictors()) %>% 
  step_corr(all_numeric_predictors()) %>%
  step_impute_knn(all_predictors()) %>% 
  step_dummy(all_nominal_predictors()) %>% 
  step_smote(pfs_status) 

#xgb worklow
prad_wf <- workflow() %>%
  add_model(xgb_model) %>% 
  add_recipe(recipe_object)

prad_wf

xgb_rs <- prad_wf %>%
  update_model(xgb_model) %>%
  fit_resamples(
    resamples = prad_boot, 
    control = control_resamples(save_pred = TRUE, verbose = TRUE),
    metrics = metric_set(accuracy, sens, spec, roc_auc)
  )

collect_metrics(xgb_rs)

collect_predictions(xgb_rs)

# confusion matrix
xgb_rs %>% 
  collect_predictions() %>% 
  conf_mat(
    truth = pfs_status,
    estimate = .pred_class
  ) %>%
  autoplot(type="heatmap")

#roc curve
xgb_rs %>%
  collect_predictions() %>%
  roc_curve(pfs_status, .pred_0) %>%
  autoplot()

# variable importance xgb --------------------------------------------------

# get the prepped data set for assessing variable importance
vip_data <- 
  recipe_object %>% 
  prep() %>% 
  bake(new_data = NULL) %>% 
  select(!patient_id)

xgb_model %>% 
  set_engine("xgboost") %>% 
  fit(pfs_status ~ ., data = vip_data) %>% 
  vip(geom = "point")
