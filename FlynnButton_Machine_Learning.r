# This script was developed as part of a data science course project.
# It uses 2021 Canadian Census microdata to predict low-income status (LIMBM)
# among adults in Nova Scotia. The file includes data cleaning, feature
# engineering, and training/evaluation of three classifiers: Random Forest,
# XGBoost, and Elastic Net logistic regression. Hyperparameters are tuned
# via cross-validation using the tidymodels framework, with performance
# compared using confusion matrices and ROC curves.
#
# Author: Flynn Button (Sole Author)
# ============================================================



rm(list = ls())

library(randomForest)
library(tidymodels)
library(dplyr)
library(ranger) 
library(tibble)
library(ggplot2) 
library(yardstick)
library(doParallel)

data <- read.csv("/Users/Flynn/Desktop/BSE/Data Science/Data/cen21_ind_98m0001x_part_rec21/data_donnees_2021_ind_v2.csv")



# Select only observations for Nova Scotia
my_data <- data[data$PR == 46,]

# Drop weighting variables and province
weight_index <- c(128:144, 101)
my_data <- my_data[-weight_index]



#Restricting to adults:
my_data <- my_data[my_data$AGEGRP > 6,]
my_data <- my_data[,-weight_index]

my_data[my_data == 88888888 | my_data == 8888 |my_data == 88] <- NA


obs_with_na <- which(rowSums(is.na(my_data)) > 0)


# Step 2: Identify columns with >10% NA in those rows
cols_with_na <- which(colSums(is.na(my_data[obs_with_na, ])) / length(obs_with_na) > 0.1)
cols_with_na  # indices of columns

names(my_data)[cols_with_na]

# Step 3: Check counts of NA per column
colSums(is.na(my_data[obs_with_na, cols_with_na]))

high_na_obs <- which(rowSums(is.na(my_data)) > 5)

my_data <- my_data[-high_na_obs,]


#88888888 is used for when an observation is unavailable, is NA. Only 395 observations with NA values
#99999999 is used when the question is not applicable, in questions about 



poverty_rate <- sum(my_data$LoMBM_2018==2)/nrow(my_data)

poverty_na <- which(my_data$LoMBM_2018==8 & my_data$LoMBM_2018==9)

class(my_data$LoMBM_2018)

hist(my_data$LoMBM_2018)




my_data$YRIM <- as.numeric(as.character(my_data$YRIM))

# Original binning
my_data$YRIM_bin <- case_when(
  my_data$YRIM == 9999                              ~ 0,   # Non-immigrant
  my_data$YRIM >= 1995 & my_data$YRIM <= 1999      ~ 10,
  my_data$YRIM >= 2000 & my_data$YRIM <= 2004      ~ 11,
  my_data$YRIM >= 2005 & my_data$YRIM <= 2009      ~ 12,
  my_data$YRIM >= 2010 & my_data$YRIM <= 2014      ~ 13,
  my_data$YRIM >= 2015 & my_data$YRIM <= 2019      ~ 14,
  my_data$YRIM >= 2020 & my_data$YRIM <= 2021      ~ 15,
  TRUE                                              ~ my_data$YRIM   # bins 1–9
)

# Group bins with descriptive labels
my_data$YRIM_bin_grouped <- case_when(
  my_data$YRIM_bin == 0               ~ "Non-Immigrant",
  my_data$YRIM_bin %in% c(1,2)       ~ "Early Census Cohort 1",
  my_data$YRIM_bin %in% c(3,4)       ~ "Early Census Cohorts 2",
  my_data$YRIM_bin %in% c(5,6)       ~ "Cohorts 3",
  my_data$YRIM_bin %in% c(7,8)       ~ "Cohorts 4",
  my_data$YRIM_bin %in% c(9,10)      ~ "Cohorts 5",
  my_data$YRIM_bin %in% c(11,12)     ~ "Cohorts 6",
  my_data$YRIM_bin %in% c(13,14,15)  ~ "Recent Immigrants 2010–2021",
  TRUE                               ~ "Other"
)

# Convert to factor
my_data$YRIM_bin_grouped <- factor(my_data$YRIM_bin_grouped)

# Check the distribution
table(my_data$YRIM_bin_grouped)


my_data$YRIM_bin_grouped <- factor(my_data$YRIM_bin_grouped)

# Check counts
table(my_data$YRIM_bin_grouped)



library(dplyr)

my_data %>%
  group_by(YRIM_bin_grouped) %>%
  summarise(count = n()) %>%
  arrange(YRIM_bin_grouped)

my_data$YRIM[is.na(my_data$YRIM_bin) ==1]


#------Income Variable Indices
income_variables <- c(
  "CFInc", "CFInc_AT", "CHDBN", "COVID_ERB", "CQPPB",
  "CapGn", "ChldC", "EFDecile", "EFInc", "EFInc_AT",
  "EICBN", "EfDIMBM_2018", "EmpIn", "GovtI", "GTRfs",
  "HHInc", "HHInc_AT", "HHMRKINC", "IncTax", "Invst",
  "LICO_BT", "LICO_AT", "LoLIMA", "LoLIMB", 
  "MrkInc", "OASGI", "OtInc", "Retir", "SempI",
  "TotInc", "TotInc_AT", "Wages", "PPSORT", "YRIM_bin"
)

housing_variables<- c("HCORENEED_IND","NOS") 


columns_to_remove <- which(names(my_data) %in% c(income_variables, housing_variables))

my_data <- my_data[,-columns_to_remove]

my_data <- my_data %>% 
  select(-YRIM)

obs_with_na <- which(rowSums(is.na(my_data)) > 0)
my_data <- my_data[-obs_with_na,]

sum(is.na(my_data))


# Recode LoMBM_2018 to 0/1
my_data$LoMBM <- ifelse(my_data$LoMBM_2018 == 2, 1, 0)

# Convert to factor for classification
my_data$LoMBM <- factor(my_data$LoMBM, levels = c(0, 1))


# Check counts
table(my_data$LoMBM)


my_data <- my_data %>% 
  select(-LoMBM_2018)





which(rowSums(is.na(my_data))>0)

#-----------      Initial Plots


lom_counts <- my_data %>%
  count(LoMBM) %>%
  mutate(percent = n / sum(n) * 100,
         label = paste0(LoMBM, " (", round(percent, 1), "%)"))

# Pie chart
ggplot(lom_counts, aes(x = "", y = percent, fill = factor(LoMBM))) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5)) +
  scale_fill_brewer(palette = "Set2", name = "Individual counted in Poverty") +
  theme_void() +
  labs(title = "Distribution of Poverty Status (LoMBM)")


#Plot LoMBM by YRIM_bin
lom_prop <- my_data %>%
  group_by(YRIM_bin_grouped) %>%
  summarise(prop_1 = mean(LoMBM == 1), .groups = "drop")  # proportion of 1s


ggplot(lom_prop, aes(x = YRIM_bin_grouped, y = prop_1, fill = YRIM_bin_grouped)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Proportion of Poverty (LoMBM = 1) by Year of Immigration Bin",
    x = "Year of Immigration Bin",
    y = "Percentage of LoMBM = 1",
    fill = "YRIM_bin"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set2")

library(tidymodels)


set.seed(123)

#######     Split Data: Test/Train
split <- initial_split(my_data, prop = 0.7, strata = LoMBM)
train_data <- training(split)
test_data  <- testing(split)

#For entire data set
table(my_data$LoMBM)
prop.table(table(my_data$LoMBM))

# For training set
table(train_data$LoMBM)            
prop.table(table(train_data$LoMBM))  

# For test set
table(test_data$LoMBM)
prop.table(table(test_data$LoMBM))


rf_recipe <- recipe(LoMBM ~ ., data = train_data) %>%
  step_zv() %>%                        # remove zero-variance predictors
  step_impute_median(all_numeric_predictors()) %>%
  step_impute_mode(all_nominal_predictors())


rf_model <- rand_forest(
  trees = 100,        # smaller number of trees for tuning
  mtry  = tune(),
  min_n = tune()
) %>%
  set_engine("ranger", importance = "impurity", num.threads = parallel::detectCores()) %>%
  set_mode("classification")

rf_wf <- workflow() %>%
  add_recipe(rf_recipe) %>%
  add_model(rf_model)


#Hyperparameter tuning###

rf_grid <- grid_random(
  mtry(range = c(2, ncol(train_data) - 1)),
  min_n(range = c(2, 10)),
  size = 10
)

cv_splits <- vfold_cv(train_data, v = 3, strata = LoMBM)

rf_tune <- tune_grid(
  rf_wf,
  resamples = cv_splits,
  grid = rf_grid,
  metrics = metric_set(accuracy, roc_auc)
)

# Best hyperparameters
best_params <- rf_tune %>%
  collect_metrics() %>%
  filter(.metric == "accuracy") %>%
  arrange(desc(mean)) %>%
  slice_head(n = 1)


rf_final <- rand_forest(
  trees = 500,         # higher number of trees for final fit
  mtry  = best_params$mtry,
  min_n = best_params$min_n
) %>%
  set_engine("ranger", importance = "impurity", num.threads = parallel::detectCores()) %>%
  set_mode("classification")

rf_workflow_final <- workflow() %>%
  add_recipe(rf_recipe) %>%
  add_model(rf_final)

rf_fit_final <- fit(rf_workflow_final, data = train_data)

#     RF.   Prediction
rf_test_preds <- predict(rf_fit_final, test_data, type = "prob") %>%
  bind_cols(test_data %>% select(LoMBM)) %>%
  mutate(
    .pred_class = factor(ifelse(.pred_1 > 0.5, 1, 0), levels = levels(LoMBM))
  )


conf_mat_rf <- conf_mat(rf_test_preds, truth = LoMBM, estimate = .pred_class)
conf_mat_rf


prop.table(conf_mat_rf$table, margin = 2)


#           ROC Curve

roc_data <- roc_curve(rf_test_preds, truth = LoMBM, .pred_1)
autoplot(roc_data)


#       Importance Plot

rf_fit_final$fit$fit$fit %>%  # extract ranger model
  ranger::importance() %>%
  enframe(name = "Variable", value = "Importance") %>%
  arrange(desc(Importance)) %>%
  ggplot(aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Random Forest Variable Importance")













################ Boosting Model ###########


num_neg <- sum(train_data$LoMBM == 0)
num_pos <- sum(train_data$LoMBM == 1)
scale_weight <- num_neg / num_pos
cat("Scale_pos_weight:", scale_weight, "\n")


xgb_recipe <- recipe(LoMBM ~ ., data = train_data) %>%
  step_zv() %>%
  step_impute_median(all_numeric_predictors()) %>%
  step_impute_mode(all_nominal_predictors()) %>%
  step_dummy(all_nominal_predictors())


xgb_model <- boost_tree(
  trees = 200,          
  tree_depth = tune(),
  learn_rate = tune(),
  loss_reduction = tune(),
  sample_size = tune(),
  stop_iter = 10
) %>%
  set_engine("xgboost", scale_pos_weight = scale_weight) %>% 
  set_mode("classification")



xgb_wf <- workflow() %>%
  add_recipe(xgb_recipe) %>%
  add_model(xgb_model)


cv_splits <- vfold_cv(train_data, v = 5, strata = LoMBM)


xgb_grid <- grid_random(
  tree_depth(range = c(3, 10)),
  learn_rate(range = c(0.01, 0.2)),
  loss_reduction(range = c(0, 5)),
  sample_size = sample_prop(range = c(0.6, 1)),
  size = 20
)


xgb_tune <- tune_grid(
  xgb_wf,
  resamples = cv_splits,
  grid = xgb_grid,
  metrics = metric_set(roc_auc, accuracy)
)


best_params <- select_best(xgb_tune, metric = "roc_auc")
xgb_final <- finalize_workflow(xgb_wf, best_params)


xgb_fit_final <- fit(xgb_final, data = train_data)

xgb_fit_final$fit$fit$fit$params



xgb_preds <- predict(xgb_fit_final, test_data, type = "prob") %>%
  bind_cols(test_data %>% select(LoMBM)) %>%
  mutate(
    .pred_class = factor(ifelse(.pred_1 > 0.5, 1, 0), levels = c(0, 1))
  )


conf_mat_xgb <- conf_mat(xgb_preds, truth = LoMBM, estimate = .pred_class)
conf_mat_xgb


library(ggplot2)

# Elastic Net------------


enet_recipe <- recipe(LoMBM ~ ., data = train_data) %>%
  step_zv() %>%                           # remove zero-variance predictors
  step_impute_median(all_numeric_predictors()) %>%
  step_impute_mode(all_nominal_predictors()) %>%
  step_dummy(all_nominal_predictors())    # convert categorical to dummy variables


enet_model <- logistic_reg(
  penalty = tune(),   # lambda
  mixture = tune()    # alpha (0 = ridge, 1 = lasso)
) %>%
  set_engine("glmnet") %>%
  set_mode("classification")


enet_wf <- workflow() %>%
  add_recipe(enet_recipe) %>%
  add_model(enet_model)


# 5. Hyperparameter tuning

enet_grid <- grid_regular(
  penalty(range = c(-4, 0)),   # log10 scale for lambda
  mixture(range = c(0,1)),     # alpha 0 to 1
  levels = 5
)

cv_splits <- vfold_cv(train_data, v = 3, strata = LoMBM)


enet_tune <- tune_grid(
  enet_wf,
  resamples = cv_splits,
  grid = enet_grid,
  metrics = metric_set(roc_auc, accuracy)
)



best_params <- enet_tune %>% select_best(metric = "roc_auc")
enet_final <- finalize_workflow(enet_wf, best_params)


enet_fit <- fit(enet_final, data = train_data)

# Making prediction on testing data


# Probabilities
enet_probs <- predict(enet_fit, test_data, type = "prob") %>%
  bind_cols(test_data %>% select(LoMBM))

# Predicted classes
enet_preds <- predict(enet_fit, test_data, type = "class") %>%
  bind_cols(test_data %>% select(LoMBM))



# Confusion matrix

conf_mat_enet <- conf_mat(enet_preds, truth = LoMBM, estimate = .pred_class)

# ROC Curve
roc_curve(enet_probs, truth = LoMBM, .pred_1) %>% autoplot()


library(dplyr)
library(ggplot2)
library(purrr)

# Put confusion matrices in a named list
cm_list <- list(
  RandomForest = conf_mat_rf,
  XGBoost      = conf_mat_xgb,
  ElasticNet   = conf_mat_enet
)

# Tidy all confusion matrices
cm_tidy <- imap_dfr(cm_list, ~{
  df <- as.data.frame(.x$table)  # convert table to data frame
  # rename the columns explicitly
  colnames(df) <- c("Prediction", "Reference", "Freq")
  
  df %>%
    mutate(
      Model = .y,
      Reference = as.factor(Reference),
      Prediction = as.factor(Prediction)
    )
})

# Convert counts to proportions per Reference
cm_tidy <- cm_tidy %>%
  group_by(Model, Reference) %>%
  mutate(Proportion = Freq / sum(Freq)) %>%
  ungroup()

# Plot heatmaps
ggplot(cm_tidy, aes(x = Prediction, y = Reference, fill = Proportion)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Proportion, 2)), size = 5) +
  facet_wrap(~Model) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  labs(
    title = "Confusion Matrix Comparison",
    x = "Predicted",
    y = "Actual",
    fill = "Proportion"
  )


library(vip)
library(patchwork)   


rf_engine <- rf_fit_final$fit$fit$fit

rf_plot <- vip(
  rf_engine,
  num_features = 10,
  aesthetics = list(fill = "steelblue")
) +
  theme_minimal(base_size = 13) +
  labs(title = "Random Forest Variable Importance") +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title.y = element_blank()
  )


xgb_engine <- xgb_fit_final$fit$fit$fit

xgb_plot <- vip(
  xgb_engine,
  num_features = 10,
  aesthetics = list(fill = "darkorange")
) +
  theme_minimal(base_size = 13) +
  labs(title = "XGBoost Variable Importance") +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title.y = element_blank()
  )


rf_plot + xgb_plot
