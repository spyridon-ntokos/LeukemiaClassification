## Spyridon Ntokos
## CYO Project - Leukemia Classification
## HarvardX: PH125.9x - Capstone Project
## https://github.com/spyridon-ntokos

###############################################
# Create train & validation (indepedent) sets #
###############################################

if(!require(tidyverse)) install.packages("tidyverse", repos = "http://cran.us.r-project.org")
if(!require(caret)) install.packages("caret", repos = "http://cran.us.r-project.org")
if(!require(ggplot2)) install.packages("ggplot2", repos = "http://cran.us.r-project.org")
if(!require(naniar)) install.packages("naniar", repos = "http://cran.us.r-project.org")
if(!require(matrixStats)) install.packages("matrixStats", repos = "http://cran.us.r-project.org")
if(!require(gridExtra)) install.packages("gridExtra", repos = "http://cran.us.r-project.org")


# Gene expression dataset (Golub et al.):
# https://www.kaggle.com/crawford/gene-expression/data

dl <- tempfile()
download.file("https://storage.googleapis.com/kaggle-data-sets/1868/3249/bundle/archive.zip?GoogleAccessId=web-data@kaggle-161607.iam.gserviceaccount.com&Expires=1590574751&Signature=pNKbnTCWdotyK7OrAVw%2FRlxbTyU2%2BbtseYhX3yd0ahsPrCubhMpDU1U5GtHevYPq%2FaO0jTZzQ0dcXTBOYtR%2BrnUrM3Gcg3N%2FIrWJvw8EP%2FThEllTskhx6MsmeynQ%2BR8cyt6vWWR9S6%2Fm73OBLs7E1xBo6wpurikEP%2Fq%2FHx9Q5ugdIMtPO57Kw6N4pkzDBWMyPmKBlOXj3vtRTyvHgORFYdcdhhFcGcRq06CjEImO%2FbPn1jt2C%2BB5ksmnV2ICcbTIQ9kCn6ayiY%2B3Eai1tCyLTasS18E2diMEZ95qnY5vQxw5lDXJwq3G%2F4HQ1o90%2FDmtngJ9YI1hidl%2Fk3Xb%2BiyWYg%3D%3D&response-content-disposition=attachment%3B+filename%3Dgene-expression.zip", dl)

# Extract .csv files
actual <- read.csv(unzip(dl, "actual.csv"), stringsAsFactors = FALSE)
train_set <- read.csv(unzip(dl, "data_set_ALL_AML_train.csv"), stringsAsFactors = FALSE)
validation <- read.csv(unzip(dl, "data_set_ALL_AML_independent.csv"), stringsAsFactors = FALSE)

rm(dl)

####################
# Data preparation #
####################

## Data understanding
head(actual)
head(train_set)
head(validation)

dim(as.matrix(actual))
dim(as.matrix(train_set))
dim(as.matrix(validation))

## Reshaping the Data
original_train <- train_set

# Divide 'actual' table 
actual_train <- actual[1:38, ]
actual_validation <- actual[39:72, ]

# Removing irrelevant 'call' columns
call_columns <- grep(pattern = "call", x = names(train_set))
train_set <- train_set[,-call_columns]

# Store 'Gene Access Name' column as new column headers
gene_id <- as.vector(t(train_set[,2]))

# Keep only numeric data
train_set <- train_set[, - c(1,2)]

# Transpose table and name new columns according to 'gene_id'
# and rows according to patient ID
train_set <- as.data.frame(t(train_set))
colnames(train_set) <- gene_id
train_set <- cbind(patient = row.names(train_set),
                   cancer_type = actual_train$cancer,
                   train_set)

# Drop 'X' from patient IDs
train_set <- train_set %>%
  separate(patient, c(NA, "patient"), sep = "X")

# Any missing values?
anyNA(train_set)
try(gg_miss_upset(train_set))
train_set[, 'Var.1650']
gene_id[1649]

# Remove empty column & confirm
train_set <- train_set[, - c(1650)]
vis_miss(train_set) +
  theme(axis.text.x=element_blank()) +
  ggtitle("Missing data in train set?")

# Confirm changes in table
head(train_set) %>% as_tibble()

# Repeat for validation set:

# First, make sure genes in validation set are also in training set
validation <- validation %>% 
  semi_join(original_train, by = "Gene.Accession.Number")

# Removing irrelevant 'call' columns
call_columns <- grep(pattern = "call", x = names(validation))
validation <- validation[,-call_columns]

# Store 'Gene Access Name' column as new column headers
gene_id <- as.vector(t(validation[,2]))

# Keep only numeric data
validation <- validation[, - c(1,2)]

# Transpose table and name new columns according to 'gene_id'
# and rows according to patient ID
validation <- as.data.frame(t(validation))
colnames(validation) <- gene_id
validation <- cbind(patient = row.names(validation), 
                    cancer_type = actual_validation$cancer,
                    validation)

# Drop 'X' from patient IDs
validation <- validation %>%
  separate(patient, c(NA, "patient"), sep = "X")

# Any missing values?
vis_miss(validation) +
  theme(axis.text.x = element_blank()) +
  ggtitle("Missing data in validation set?")

# Confirm changes in table
head(validation) %>% as_tibble()

# Tidy workspace, remove anything unnecessary
rm(gene_id, call_columns, original_train)

#############################
# Exploratory Data Analysis #
#############################

## PCA analysis (quick)

# Matrix transformation
x <- as.matrix(train_set[, - c(1,2)])
x_centered <- sweep(x, 2, colMeans(x), FUN = "-")
x_scaled <- sweep(x_centered, 2, colSds(x), FUN = "/")

# Proportion of variance
pca <- prcomp(x_scaled)
summary(pca)

var_explained <- cumsum(pca$sdev^2/sum(pca$sdev^2))
plot(var_explained)

# Detect non-overlapping PCs
data.frame(type = actual_train$cancer, pca$x) %>%
  gather(key = "PC", value = "value", -type) %>%
  ggplot(aes(PC, value, fill = type)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45))

# Plot non-overlapping PCs against each other
data.frame(pca$x[,c(1,3)], type = actual_train$cancer) %>%
  ggplot(aes(PC1, PC3, color = type)) +
  geom_point()

# Remove anything unnecessary
rm(var_explained, x, x_centered, x_scaled, pca)

## Top 10 genes expressed by leukemia type
train_set %>%
  select(-patient) %>%
  group_by(cancer_type) %>%
  summarise_all(list(mean)) %>%
  gather("gene", "mean", - cancer_type) %>%
  group_by(cancer_type) %>%
  top_n(10) %>%
  ggplot(aes(gene, fill = mean)) +
  geom_bar() +
  facet_wrap(~cancer_type) +
  coord_flip() +
  ggtitle("Top 10 Genes Expressed in AML and ALL patients") +
  labs(x = "Gene Accession Number",
       caption = "Based on mean expression value per leukemia type") +
  scale_fill_continuous("Mean expression", low = "steelblue3", high = "indianred4") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

## Top 10 highest gene expressions...
# ...in ALL patients
all_top_genes <- train_set %>%
  filter(cancer_type == "ALL") %>%
  gather("gene", "value", - c(patient, cancer_type)) %>% 
  top_n(10, value) %>%
  ggplot(aes(patient, value)) +
  geom_point(aes(col = gene)) +
  xlab("ALL-patient") +
  ylab("Max. expressed value") +
  scale_color_brewer("Gene Accession Numner", palette = "Dark2") +
  ggtitle("Top 10 highest gene expressions")
  
# ...in AML patients
aml_top_genes <- train_set %>%
  filter(cancer_type == "AML") %>%
  gather("gene", "value", - c(patient, cancer_type)) %>% 
  top_n(10, value) %>%
  ggplot(aes(patient, value)) +
  geom_point(aes(col = gene)) +
  xlab("AML-patient") +
  ylab("Max. expressed value") +
  scale_color_brewer("Gene Accession Numner", palette = "Dark2")

grid.arrange(all_top_genes, aml_top_genes, ncol = 1)

# Remove anything unnecessary 
rm(all_top_genes, aml_top_genes)

##################
# Model Training #
##################

if(!require(naivebayes)) install.packages("naivebayes")
if(!require(kernlab)) install.packages("kernlab")
if(!require(RSNNS)) install.packages("RSNNS")
if(!require(deepnet)) install.packages("deepnet")
if(!require(caTools)) install.packages("caTools")
if(!require(gbm)) install.packages("gbm")
if(!require(C50)) install.packages("C50")
if(!require(plyr)) install.packages("plyr")

## Model preparation

train_canc <- train_set[, !(colnames(train_set) == "patient")]

# First partition the data into train and test sets
set.seed(1998, sample.kind = "Rounding")
test_index <- createDataPartition(train_canc$cancer_type,
                                  p = 0.15, list = FALSE)
test_canc <- train_canc[test_index, ] %>%
  mutate(cancer_type = as.factor(cancer_type))
train_canc <- train_canc[-test_index, ] %>%
  mutate(cancer_type = as.factor(cancer_type))

rm(test_index)

## Setup metric and control
control <- trainControl(method = "repeatedcv", 
                        number = 10, repeats = 3)
metric <- "Accuracy"

## Setup models

# Naive Bayes
set.seed(1998, sample.kind = "Rounding")
fit.naiveBayes <-caret::train(cancer_type ~ .,
                              data = train_canc,
                              method = "naive_bayes",
                              trControl = control,
                              metric = metric)
preds.naiveBayes <- predict(fit.naiveBayes, test_canc)
caret::confusionMatrix(preds.naiveBayes, test_canc$cancer_type)$overall

# k-nearest neighbors
set.seed(1998, sample.kind = "Rounding")
fit.knn <- caret::train(cancer_type ~ .,
                        data = train_canc,
                        method = "knn",
                        trControl = control,
                        metric = metric,
                        tuneGrid = data.frame(k = seq(5, 31, 2)))
fit.knn$bestTune
preds.knn <- predict(fit.knn, test_canc)
caret::confusionMatrix(preds.knn, test_canc$cancer_type)$overall

# Support Vector Machines
set.seed(1998, sample.kind = "Rounding")
fit.svm <- caret::train(cancer_type ~ .,
                        data = train_canc,
                        method = "svmRadial",
                        trControl = control,
                        metric = metric)
preds.svm <- predict(fit.svm, test_canc)
caret::confusionMatrix(preds.svm, test_canc$cancer_type)$overall

# Multilayer Perceptrons
set.seed(1998, sample.kind = "Rounding")
fit.mlp <- caret::train(cancer_type ~ .,
                        data = train_canc,
                        method = "mlp",
                        trControl = control,
                        metric = metric)
preds.mlp <- predict(fit.mlp, test_canc)
caret::confusionMatrix(preds.mlp, test_canc$cancer_type)$overall

# Deep Neural Network
set.seed(1998, sample.kind = "Rounding")
fit.dnn <- caret::train(cancer_type ~ .,
                        data = train_canc,
                        method = "dnn",
                        trControl = control,
                        metric = metric)
preds.dnn <- predict(fit.dnn, test_canc)
caret::confusionMatrix(preds.dnn, test_canc$cancer_type)$overall

# Generalized Linear Model
set.seed(1998, sample.kind = "Rounding")
fit.glm <- caret::train(cancer_type ~ .,
                        data = train_canc,
                        method = "glm",
                        trControl = control,
                        metric = metric)
preds.glm <- predict(fit.glm, test_canc)
caret::confusionMatrix(preds.glm, test_canc$cancer_type)$overall

# Boosted Logistic Regression
set.seed(1998, sample.kind = "Rounding")
fit.lb <- caret::train(cancer_type ~ .,
                       data = train_canc,
                       method = "LogitBoost",
                       trControl = control,
                       metric = metric)
preds.lb <- predict(fit.lb, test_canc)
caret::confusionMatrix(preds.lb, test_canc$cancer_type)$overall

# Stochastic Gradient Boosting
set.seed(1998, sample.kind = "Rounding")

gbm_grid <- expand.grid(
  n.trees = 5,
  interaction.depth = 3,
  shrinkage = 0.03,
  n.minobsinnode = 2)

fit.gbm <- caret::train(cancer_type ~ .,
                        data = train_canc,
                        method = "gbm",
                        trControl = control,
                        metric = metric,
                        tuneGrid = gbm_grid,
                        verbose = FALSE)

preds.gbm <- predict(fit.gbm, test_canc)
caret::confusionMatrix(preds.gbm, test_canc$cancer_type)$overall

# C5.0 (Decision Tree algorithm)
set.seed(1998, sample.kind = "Rounding")
fit.c5 <- caret::train(cancer_type ~ .,
                       data = train_canc,
                       method = "C5.0",
                       trControl = control,
                       metric = metric,
                       verbose = FALSE)
preds.c5 <- predict(fit.c5, test_canc)
caret::confusionMatrix(preds.c5, test_canc$cancer_type)$overall

plot(fit.c5, main = "Decision Tree")

# Random Forest
set.seed(1998, sample.kind = "Rounding")
fit.rf <- caret::train(cancer_type ~ .,
                       data = train_canc,
                       method = "rf",
                       trControl = control,
                       metric = metric,
                       verbose = FALSE)
preds.rf <- predict(fit.rf, test_canc)
caret::confusionMatrix(preds.rf, test_canc$cancer_type)$overall

plot(fit.rf$finalModel, main = "Random Forest")

# K-means Clustering
predict_kmeans <- function(x, k) {
  centers <- k$centers    # extract cluster centers
  
  # calculate distance to cluster centers
  distances <- sapply(1:nrow(x), function(i){
    apply(centers, 1, function(y) dist(rbind(x[i,], y)))
  })
  max.col(-t(distances))  # select cluster with min distance to center
}

train_canc_m <- as.matrix(train_canc[, -1])
test_canc_m <- as.matrix(test_canc[, -1])

set.seed(1998, sample.kind = "Rounding")
k <- kmeans(train_canc_m, centers = 2)
kmeans_preds <- ifelse(predict_kmeans(test_canc_m, k) == 1, "ALL", "AML")

caret::confusionMatrix(as.factor(kmeans_preds), as.factor(test_canc$cancer_type))$overall

rm(train_canc_m, test_canc_m, k, kmeans_preds)

# Model results
results <- resamples(list(naiveBayes = fit.naiveBayes,
                          knn = fit.knn,
                          svm = fit.svm,
                          mlp = fit.mlp,
                          dnn = fit.dnn,
                          glm = fit.glm,
                          logitBoost = fit.lb,
                          gbm = fit.gbm,
                          C5.0 = fit.c5,
                          randomForest = fit.rf))
dotplot(results, main = "Model Accuracy Results")

rm(fit.gbm, fit.dnn, fit.glm,
   preds.gbm, preds.dnn, preds.glm)

# Create an ensemble with 7 best-performing models
ensemble <- cbind(knn = preds.knn == "ALL",
                  naiveBayes = preds.naiveBayes == "ALL",
                  randomForest = preds.rf == "ALL",
                  logitBoost = preds.lb == "ALL",
                  C5.0 = preds.c5 == "ALL",
                  svm = preds.svm == "ALL",
                  mlp = preds.mlp == "ALL")
preds.ensemble <- ifelse(rowMeans(ensemble) > 0.5, "ALL", "AML")
caret::confusionMatrix(as.factor(preds.ensemble), 
                       test_canc$cancer_type)$overall

rm(train_canc, test_canc)

# Apply final models to "unknown" data in validation set

preds.knn <- predict(fit.knn, validation)
results <- data.frame(Method = "k-Nearest Neighbors",
                      Accuracy = caret::confusionMatrix(preds.knn, validation$cancer_type)$overall["Accuracy"])

preds.naiveBayes <- predict(fit.naiveBayes, validation)
results <- dplyr::bind_rows(results,
                            data.frame(Method = "Naive Bayes",
                                       Accuracy = caret::confusionMatrix(preds.naiveBayes, validation$cancer_type)$overall["Accuracy"]))

preds.rf <- predict(fit.rf, validation)
results <- dplyr::bind_rows(results,
                            data.frame(Method = "Random Forest",
                                       Accuracy = caret::confusionMatrix(preds.rf, validation$cancer_type)$overall["Accuracy"]))

preds.lb <- predict(fit.lb, validation)
results <- dplyr::bind_rows(results,
                            data.frame(Method = "Boosted Logistic Regression",
                                       Accuracy = caret::confusionMatrix(preds.lb, validation$cancer_type)$overall["Accuracy"]))

preds.c5 <- predict(fit.c5, validation)
results <- dplyr::bind_rows(results,
                            data.frame(Method = "C5.0 (Decision Tree Algorithm)",
                                       Accuracy = caret::confusionMatrix(preds.c5, validation$cancer_type)$overall["Accuracy"]))

preds.svm <- predict(fit.svm, validation)
results <- dplyr::bind_rows(results,
                            data.frame(Method = "Support Vector Machines",
                                       Accuracy = caret::confusionMatrix(preds.svm, validation$cancer_type)$overall["Accuracy"]))

preds.mlp <- predict(fit.mlp, validation)
results <- dplyr::bind_rows(results,
                            data.frame(Method = "Multilayer Perceptons",
                                       Accuracy = caret::confusionMatrix(preds.mlp, validation$cancer_type)$overall["Accuracy"]))

ensemble <- cbind(knn = preds.knn == "ALL",
                  naiveBayes = preds.naiveBayes == "ALL",
                  randomForest = preds.rf == "ALL",
                  logitBoost = preds.lb == "ALL",
                  C5.0 = preds.c5 == "ALL",
                  svm = preds.svm == "ALL",
                  mlp = preds.mlp == "ALL")
preds.ensemble <- ifelse(rowMeans(ensemble) > 0.5, "ALL", "AML")
results <- dplyr::bind_rows(results,
                            data.frame(Method = "Ensemble",
                                       Accuracy = caret::confusionMatrix(as.factor(preds.ensemble), validation$cancer_type)$overall["Accuracy"]))

# Visualise results & compare with actual data
results %>% knitr::kable()

data.frame(PatientID = actual_validation$patient,
           Predicted = preds.ensemble,
           Actual = actual_validation$cancer)
