library(keras)
library(tfdatasets)
library(readr)
library(tidyverse)
library(rsample)
library(tensorflow)

setwd("~/bio/Article-PeerJ-codes")
datasetLUSCLit <- read_csv("LUSCGenes.csv")
datasetLUSC <- read_csv("GENESLUSCSELEC.csv")
validacao<-read_csv("validacaoLUSC.csv")

datasetLUSCLit <- subset(datasetLUSCLit, select= c(TP53,NF1,ARID1A,RB1,CDKN2A,PIK3CA,NFE2L2,PTEN,KMT2D,FAT1,NOTCH1,KDM6A,HRAS,vital_status) )
datasetLUSC32 <- subset(datasetLUSC, select= c(PKHD1,CNTNAP5,HCN1,ERICH3,RELN,DNAH8,PKHD1L1,DNAH5,KMT2D,FAM135B,SYNE1,SI,CDH10,PAPPA2,DAMTS12,RYR3,MUC17,PCDH15,PCLO,COL11A1,NAV3,SPTA1,FLG,XIRP2,ZFHX4,USH2A,LRP1B,RYR2,CSMD3,MUC16,TP53,TTN,vital_status) )
datasetLUSC15 <-  subset(datasetLUSC, select= c(TP53,MUC16,LRP1B,MUC17,CDH10,FAM135B,DAMTS12,PKHD1,HCN1,RYR2,SYNE1,KMT2D,PAPPA2,SI,CSMD3,vital_status) )
validacaoLUSC32<- subset(validacao, select= c(PKHD1,CNTNAP5,HCN1,ERICH3,RELN,DNAH8,PKHD1L1,DNAH5,KMT2D,FAM135B,SYNE1,SI,CDH10,PAPPA2,DAMTS12,RYR3,MUC17,PCDH15,PCLO,COL11A1,NAV3,SPTA1,FLG,XIRP2,ZFHX4,USH2A,LRP1B,RYR2,CSMD3,MUC16,TP53,TTN,vital_status) ) 
validacaoLUSC15<- subset(validacao, select= c(TP53,MUC16,LRP1B,MUC17,CDH10,FAM135B,DAMTS12,PKHD1,HCN1,RYR2,SYNE1,KMT2D,PAPPA2,SI,CSMD3,vital_status) ) 

#datasetLUSCLit
# first we split between training and testing sets
split <- initial_split(datasetLUSCLit, prop = 4/5)
train <- training(split)
write_csv(train, file = "lusc_train.csv")
test <- testing(split)
write_csv(test, file = "lusc_test.csv")


TRAIN_DATA_URL<- "lusc_train.csv"
TEST_DATA_URL<- "lusc_test.csv"

train_file_path <- get_file("lusc_train.csv",TRAIN_DATA_URL )
test_file_path <- get_file("lusc_test.csv", TEST_DATA_URL )


train_dataset <- make_csv_dataset(
  train_file_path, 
  field_delim = ",",
  batch_size = 5, 
  num_epochs = 1,
)

#test_dataset <- train_dataset <- make_csv_dataset(
test_dataset <- make_csv_dataset(
  test_file_path, 
  field_delim = ",",
  batch_size = 5, 
  num_epochs = 1
)


train_dataset %>% 
  reticulate::as_iterator() %>% 
  reticulate::iter_next() %>% 
  reticulate::py_to_r()

spec <- feature_spec(train_dataset, vital_status ~ .) %>% 
  step_numeric_column(all_numeric(), normalizer_fn = scaler_standard()) %>% 
  step_categorical_column_with_vocabulary_list(all_nominal()) %>% 
  step_indicator_column(all_nominal())

spec <- fit(spec)
layer <- layer_dense_features(feature_columns = dense_features(spec))
train_dataset %>% 
  reticulate::as_iterator() %>% 
  reticulate::iter_next() %>% 
  layer()


model <- keras_model_sequential() %>% 
  layer_dense_features(feature_columns = dense_features(spec)) %>%
  layer_dropout(rate = 0.2) %>%
  layer_dense(units = 50, activation = "relu") %>%
  layer_dropout(rate = 0.2) %>%
  layer_dense(units = 30, activation = "relu") %>% 
  layer_dropout(rate = 0.2) %>%
  layer_dense(units = 20, activation = "relu") %>% 
  layer_dropout(rate = 0.2) %>%
  layer_dense(units = 1, activation = "sigmoid")

model %>% compile(
  loss = "binary_crossentropy",
  optimizer = "adam",
  metrics = "accuracy"
)

history <- model %>% 
  fit(
    train_dataset %>% dataset_use_spec(spec) %>% dataset_shuffle(500),
    epochs = 100,
    validation_data = test_dataset %>% dataset_use_spec(spec),
    verbose = 2,
  )
summary(model)

model %>% evaluate(test_dataset %>% dataset_use_spec(spec), verbose = 0)

batch <- test_dataset %>% 
  reticulate::as_iterator() %>% 
  reticulate::iter_next() %>% 
  reticulate::py_to_r()
predict(model, batch)

plot(history)

test_predictions <- model %>% predict(test, select=c(vital_status))
test_predictions[ , 1]

real <- subset(test, select= c(vital_status)) 
#real[,1]
#View(pred_dataset )
##prediction 
pred_dataset = round(test_predictions)

#Display table with all predicted and actual values
final_statusLUSCLit=cbind (real, pred_dataset )
colnames(final_statusLUSCLit) = c("Real","Predicao")
final_statusLUSCLit
final_statusLUSCLit$score <- test_predictions[, 1]

library(caret)
cfm=caret::confusionMatrix(table(final_statusLUSCLit$Predicao, final_statusLUSCLit$Real))
print(cfm)

#datasetLUSC32
# first we split between training and testing sets
split <- initial_split(datasetLUSC32, prop = 4/5)
train <- training(split)
write_csv(train, file = "LUSC32_train.csv")
test <- testing(split)
write_csv(test, file = "LUSC32_test.csv")


TRAIN_DATA_URL<- "LUSC32_train.csv"
TEST_DATA_URL<- "LUSC32_test.csv"

train_file_path <- get_file("LUSC32_train.csv",TRAIN_DATA_URL )
test_file_path <- get_file("LUSC32_test.csv",TEST_DATA_URL )


train_dataset <- make_csv_dataset(
  train_file_path, 
  field_delim = ",",
  batch_size = 6, 
  num_epochs = 1,
)

#test_dataset <- train_dataset <- make_csv_dataset(
test_dataset <- make_csv_dataset(
  test_file_path, 
  field_delim = ",",
  batch_size = 6, 
  num_epochs = 1
)


train_dataset %>% 
  reticulate::as_iterator() %>% 
  reticulate::iter_next() %>% 
  reticulate::py_to_r()

spec <- feature_spec(train_dataset, vital_status ~ .)

spec <- feature_spec(train_dataset, vital_status ~ .) %>% 
  step_numeric_column(all_numeric(), normalizer_fn = scaler_standard()) %>% 
  step_categorical_column_with_vocabulary_list(all_nominal()) %>% 
  step_indicator_column(all_nominal())

spec <- fit(spec)

layer <- layer_dense_features(feature_columns = dense_features(spec))
train_dataset %>% 
  reticulate::as_iterator() %>% 
  reticulate::iter_next() %>% 
  layer()


model <- keras_model_sequential() %>% 
  layer_dense_features(feature_columns = dense_features(spec)) %>%
  layer_dropout(rate = 0.2) %>%
  layer_dense(units = 100, activation = "relu") %>%
  layer_dropout(rate = 0.2) %>%
  layer_dense(units = 50, activation = "relu") %>% 
  layer_dropout(rate = 0.2) %>%
  layer_dense(units = 20, activation = "relu") %>% 
  layer_dropout(rate = 0.2) %>%
  layer_dense(units = 1, activation = "sigmoid")

model %>% compile(
  loss = "binary_crossentropy",
  optimizer = "adam",
  metrics = "accuracy"
)

history <- model %>% 
  fit(
    train_dataset %>% dataset_use_spec(spec) %>% dataset_shuffle(500),
    epochs = 100,
    validation_data = test_dataset %>% dataset_use_spec(spec),
    verbose = 2,
  )
summary(model)

model %>% evaluate(test_dataset %>% dataset_use_spec(spec), verbose = 0)

batch <- test_dataset %>% 
  reticulate::as_iterator() %>% 
  reticulate::iter_next() %>% 
  reticulate::py_to_r()

predict(model, batch)

plot(history)

test_predictions <- model %>% predict(test, select=c(vital_status))
#test_predictions[ , 1]

real <- subset(test, select = c(vital_status)) 
#real[,1]
#View(pred_dataset )
##prediction 
pred_dataset <- round(test_predictions)
final_statusLUSC32  <- as.data.frame(matrix(nrow = dim(real)[1], ncol = 3), )
names(final_statusLUSC32) <- c("Real","Predicao", "score")
final_statusLUSC32$Real <- real$vital_status
final_statusLUSC32$Predicao <- pred_dataset
final_statusLUSC32$score <- test_predictions[, 1]  

library(caret)
cfm <- caret::confusionMatrix(table(final_statusLUSC32$Predicao, final_statusLUSC32$Real))
print(cfm)

#datasetLUSC32 with Validation LUSC-KR
# first we split between training and testing sets
write_csv(datasetLUSC32, file = "lusc_train.csv")
write_csv(validacaoLUSC32, file = "lusc_test.csv")

TRAIN_DATA_URL<- "lusc_train.csv"
TEST_DATA_URL<- "lusc_test.csv"

train_file_path <- get_file("lusc_train.csv",TRAIN_DATA_URL )
test_file_path <- get_file("lusc_test.csv",TEST_DATA_URL )


train_dataset <- make_csv_dataset(
  train_file_path, 
  field_delim = ",",
  batch_size = 6, 
  num_epochs = 1,
)

#test_dataset <- train_dataset <- make_csv_dataset(
test_dataset <- make_csv_dataset(
  test_file_path, 
  field_delim = ",",
  batch_size = 6, 
  num_epochs = 1
)


train_dataset %>% 
  reticulate::as_iterator() %>% 
  reticulate::iter_next() %>% 
  reticulate::py_to_r()

spec <- feature_spec(train_dataset, vital_status ~ .)

spec <- feature_spec(train_dataset, vital_status ~ .) %>% 
  step_numeric_column(all_numeric(), normalizer_fn = scaler_standard()) %>% 
  step_categorical_column_with_vocabulary_list(all_nominal()) %>% 
  step_indicator_column(all_nominal())

spec <- fit(spec)
layer <- layer_dense_features(feature_columns = dense_features(spec))
train_dataset %>% 
  reticulate::as_iterator() %>% 
  reticulate::iter_next() %>% 
  layer()


model <- keras_model_sequential() %>% 
  layer_dense_features(feature_columns = dense_features(spec)) %>%
  layer_dropout(rate = 0.2) %>%
  layer_dense(units = 100, activation = "relu") %>%
  layer_dropout(rate = 0.2) %>%
  layer_dense(units = 50, activation = "relu") %>% 
  layer_dropout(rate = 0.2) %>%
  layer_dense(units = 20, activation = "relu") %>% 
  layer_dropout(rate = 0.2) %>%
  layer_dense(units = 1, activation = "sigmoid")

model %>% compile(
  loss = "binary_crossentropy",
  optimizer = "adam",
  metrics = "accuracy"
)

history <- model %>% 
  fit(
    train_dataset %>% dataset_use_spec(spec) %>% dataset_shuffle(500),
    epochs = 100,
    validation_data = test_dataset %>% dataset_use_spec(spec),
    verbose = 2,
  )
summary(model)

model %>% evaluate(test_dataset %>% dataset_use_spec(spec), verbose = 0)

batch <- test_dataset %>% 
  reticulate::as_iterator() %>% 
  reticulate::iter_next() %>% 
  reticulate::py_to_r()
predict(model, batch)

plot(history)

test_predictions <- model %>% predict(validacaoLUSC32, select=c(vital_status))
test_predictions[ , 1]

real <- subset(validacaoLUSC32, select= c(vital_status)) 
#real[,1]
#View(pred_dataset )
##prediction 
pred_dataset = round(test_predictions)

final_statusVal32  <- as.data.frame(matrix(nrow = dim(real)[1], ncol = 2), )
names(final_statusVal32) <- c("Real","Predicao")
final_statusVal32$Real <- real$vital_status
final_statusVal32$Predicao <- pred_dataset
final_statusVal32$score <- test_predictions[, 1]
library(caret)
cfm=caret::confusionMatrix(table(final_statusVal32$Predicao, final_statusVal32$Real))
print(cfm)

#datasetLUSC15
# first we split between training and testing sets
split <- initial_split(datasetLUSC15, prop = 4/5)
train <- training(split)
write_csv(train, file = "lusc_train.csv")
test <- testing(split)
write_csv(test, file = "lusc_test.csv")


TRAIN_DATA_URL<- "lusc_train.csv"
TEST_DATA_URL<- "lusc_test.csv"

train_file_path <- get_file("lusc_train.csv",TRAIN_DATA_URL )
test_file_path <- get_file("lusc_test.csv",TEST_DATA_URL )


train_dataset <- make_csv_dataset(
  train_file_path, 
  field_delim = ",",
  batch_size = 6, 
  num_epochs = 1,
)

#test_dataset <- train_dataset <- make_csv_dataset(
test_dataset <- make_csv_dataset(
  test_file_path, 
  field_delim = ",",
  batch_size = 6, 
  num_epochs = 1
)


train_dataset %>% 
  reticulate::as_iterator() %>% 
  reticulate::iter_next() %>% 
  reticulate::py_to_r()

spec <- feature_spec(train_dataset, vital_status ~ .)

spec <- feature_spec(train_dataset, vital_status ~ .) %>% 
  step_numeric_column(all_numeric(), normalizer_fn = scaler_standard()) %>% 
  step_categorical_column_with_vocabulary_list(all_nominal()) %>% 
  step_indicator_column(all_nominal())

spec <- fit(spec)
layer <- layer_dense_features(feature_columns = dense_features(spec))
train_dataset %>% 
  reticulate::as_iterator() %>% 
  reticulate::iter_next() %>% 
  layer()


model <- keras_model_sequential() %>% 
  layer_dense_features(feature_columns = dense_features(spec)) %>%
  layer_dense(units = 100, activation = "relu") %>%
  layer_dropout(rate = 0.2) %>%
  layer_dense(units = 50, activation = "relu") %>% 
  layer_dropout(rate = 0.2) %>%
  layer_dense(units = 20, activation = "relu") %>% 
  layer_dropout(rate = 0.2) %>%
  layer_dense(units = 1, activation = "sigmoid")

model %>% compile(
  loss = "binary_crossentropy",
  optimizer = "adam",
  metrics = "accuracy"
)

history <- model %>% 
  fit(
    train_dataset %>% dataset_use_spec(spec) %>% dataset_shuffle(500),
    epochs = 100,
    validation_data = test_dataset %>% dataset_use_spec(spec),
    verbose = 2,
  )
summary(model)

model %>% evaluate(test_dataset %>% dataset_use_spec(spec), verbose = 0)

batch <- test_dataset %>% 
  reticulate::as_iterator() %>% 
  reticulate::iter_next() %>% 
  reticulate::py_to_r()
predict(model, batch)

plot(history)

test_predictions <- model %>% predict(test, select=c(vital_status))
test_predictions[ , 1]

real <- subset(test, select= c(vital_status)) 
#real[,1]
#View(pred_dataset )
##prediction 
pred_dataset = round(test_predictions)

final_statusLUSC15  <- as.data.frame(matrix(nrow = dim(real)[1], ncol = 2), )
names(final_statusLUSC15) <- c("Real","Predicao")
final_statusLUSC15$Real <- real$vital_status
final_statusLUSC15$Predicao <- pred_dataset[, 1]
final_statusLUSC15$score <- test_predictions[, 1]
library(caret)
cfm=caret::confusionMatrix(table(final_statusLUSC15$Predicao, final_statusLUSC15$Real))
print(cfm)

#datasetLUSC15 with Validation LUSC-KR
# first we split between training and testing sets
write_csv(datasetLUSC15, file = "lusc_train.csv")
write_csv(validacaoLUSC15, file = "lusc_test.csv")

TRAIN_DATA_URL<- "lusc_train.csv"
TEST_DATA_URL<- "lusc_test.csv"

train_file_path <- get_file("lusc_train.csv",TRAIN_DATA_URL )
test_file_path <- get_file("lusc_test.csv",TEST_DATA_URL )


train_dataset <- make_csv_dataset(
  train_file_path, 
  field_delim = ",",
  batch_size = 6, 
  num_epochs = 1,
)

#test_dataset <- train_dataset <- make_csv_dataset(
test_dataset <- make_csv_dataset(
  test_file_path, 
  field_delim = ",",
  batch_size = 6, 
  num_epochs = 1
)


train_dataset %>% 
  reticulate::as_iterator() %>% 
  reticulate::iter_next() %>% 
  reticulate::py_to_r()

#spec <- feature_spec(train_dataset, vital_status ~ .)

spec <- feature_spec(train_dataset, vital_status ~ .) %>% 
  step_numeric_column(all_numeric(), normalizer_fn = scaler_standard()) %>% 
  step_categorical_column_with_vocabulary_list(all_nominal()) %>% 
  step_indicator_column(all_nominal())

spec <- fit(spec)
layer <- layer_dense_features(feature_columns = dense_features(spec))
train_dataset %>% 
  reticulate::as_iterator() %>% 
  reticulate::iter_next() %>% 
  layer()

model <- keras_model_sequential() %>% 
  layer_dense_features(feature_columns = dense_features(spec)) %>%
  layer_dense(units = 100, activation = "relu") %>%
  layer_dropout(rate = 0.2) %>%
  layer_dense(units = 50, activation = "relu") %>% 
  layer_dropout(rate = 0.2) %>%
  layer_dense(units = 20, activation = "relu") %>% 
  layer_dropout(rate = 0.2) %>%
  layer_dense(units = 1, activation = "sigmoid")

model %>% compile(
  loss = "binary_crossentropy",
  optimizer = "adam",
  metrics = "accuracy"
)

history <- model %>% 
  fit(
    train_dataset %>% dataset_use_spec(spec) %>% dataset_shuffle(500),
    epochs = 100,
    validation_data = test_dataset %>% dataset_use_spec(spec),
    verbose = 2,
  )
summary(model)

model %>% evaluate(test_dataset %>% dataset_use_spec(spec), verbose = 0)

batch <- test_dataset %>% 
  reticulate::as_iterator() %>% 
  reticulate::iter_next() %>% 
  reticulate::py_to_r()
predict(model, batch)

plot(history)

test_predictions <- model %>% predict(validacaoLUSC15, select=c(vital_status))
test_predictions[ , 1]

real <- subset(validacaoLUSC15, select= c(vital_status)) 
#real[,1]
#View(pred_dataset )
##prediction 
pred_dataset = round(test_predictions)

final_statusVal15  <- as.data.frame(matrix(nrow = dim(real)[1], ncol = 2), )
names(final_statusVal15) <- c("Real","Predicao")
final_statusVal15$Real <- real$vital_status
final_statusVal15$Predicao <- pred_dataset
final_statusVal15$score <- test_predictions[, 1]
library(caret)
cfm=caret::confusionMatrix(table(final_statusVal15$Predicao, final_statusVal15$Real))
print(cfm)

#ROC curve Figure 4B
#library(pROC)

library(ROCit)

save(final_statusLUSCLit, final_statusLUSC32, final_statusLUSC15, final_statusVal32, final_statusVal15, file = "metrics_noleaked.RData", compress = T)

rocit_Lit <- rocit(score = final_statusLUSCLit$score, 
                   class = final_statusLUSCLit$Real, 
                   method = "emp")

rocit_LUSC32 <- rocit(score = final_statusLUSC32$score, 
                   class = final_statusLUSC32$Real, 
                   method = "emp")

rocit_LUSC15 <- rocit(score = final_statusLUSC15$score, 
                      class = final_statusLUSC15$Real, 
                      method = "emp")

rocit_LUSCVal32 <- rocit(score = final_statusVal32$score, 
                      class = final_statusVal32$Real, 
                      method = "emp")

rocit_LUSCVal15 <- rocit(score = final_statusVal15$score, 
                         class = final_statusVal15$Real, 
                         method = "emp")

## Plot ROC curve
pdf("fig/AUC_noleaked.pdf", width = 9.2, height = 5.5)
plot(rocit_Lit, col = c("red", 1), legend = FALSE, YIndex = FALSE)
#lines(rocit_LUSC32$TPR~rocit_LUSC32$FPR,  col = "blue", lwd = 2)
#lines(rocit_LUSC15$TPR~rocit_LUSC15$FPR, col = "darkgreen", lwd = 2)
lines(rocit_LUSCVal32$TPR~rocit_LUSCVal32$FPR, col = "darkorange", lwd = 2)
lines(rocit_LUSCVal15$TPR~rocit_LUSCVal15$FPR, col = "purple", lwd = 2)

legend("bottomright", legend=c(paste0("LUSC Literature genes (AUC=", round(rocit_Lit$AUC, digits = 3), ")"),
                               #paste0("LUSC 32 genes (AUC=", round(rocit_LUSC32$AUC, digits = 3), ")"),
                               #paste0("LUSC 15 genes (AUC=", round(rocit_LUSC15$AUC, digits = 3), ")"),
                               paste0("LUSC-KR 32 genes (AUC=", round(rocit_LUSCVal32$AUC, digits = 3), ")"),
                               paste0("LUSC-KR 15 genes (AUC=", round(rocit_LUSCVal15$AUC, digits = 3), ")")),
       col=c("red",
             #"blue",
             #"darkgreen",
             "darkorange",
             "purple"), 
       lwd=2)
dev.off()
