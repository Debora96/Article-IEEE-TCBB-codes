library(keras)
library(tfdatasets)
library(readr)
library(tidyverse)
library(rsample)
library(tensorflow)


datasetLUSCLit <- read_csv("../data/LUSCGenes.csv")
datasetLUSC <- read_csv("../data/GENESLUSCSELEC.csv")
validacao<-read_csv("../data/validacaoLUSC.csv")

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


TRAIN_DATA_URL<- "../data/lusc_train.csv"
TEST_DATA_URL<- "../data/lusc_test.csv"

train_file_path <- get_file("lusc_train.csv",TRAIN_DATA_URL )
test_file_path <- get_file("lusc_test.csv",TEST_DATA_URL )


train_dataset <- make_csv_dataset(
  train_file_path, 
  field_delim = ",",
  batch_size = 5, 
  num_epochs = 1,
)

test_dataset <- train_dataset <- make_csv_dataset(
  test_file_path, 
  field_delim = ",",
  batch_size = 5, 
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
real[,1]
View(pred_dataset )
##prediction 
pred_dataset = round(test_predictions)

#Display table with all predicted and actual values
final_statusLUSCLit=cbind (real, pred_dataset )
colnames(final_statusLUSCLit) = c("Real","Predição")

library(caret)
cfm=caret::confusionMatrix(table(final_statusLUSCLit$Predição, final_statusLUSCLit$Real))
print(cfm)

#datasetLUSC32
# first we split between training and testing sets
split <- initial_split(datasetLUSC32, prop = 4/5)
train <- training(split)
write_csv(train, file = "lusc_train.csv")
test <- testing(split)
write_csv(test, file = "lusc_test.csv")


TRAIN_DATA_URL<- "../data/lusc_train.csv"
TEST_DATA_URL<- "../data/lusc_test.csv"

train_file_path <- get_file("lusc_train.csv",TRAIN_DATA_URL )
test_file_path <- get_file("lusc_test.csv",TEST_DATA_URL )


train_dataset <- make_csv_dataset(
  train_file_path, 
  field_delim = ",",
  batch_size = 6, 
  num_epochs = 1,
)

test_dataset <- train_dataset <- make_csv_dataset(
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
test_predictions[ , 1]

real <- subset(test, select= c(vital_status)) 
real[,1]
View(pred_dataset )
##prediction 
pred_dataset = round(test_predictions)

final_statusLUSC32=cbind (real, pred_dataset )
(final_statusLUSC32) = c("Real","Predição")

library(caret)
cfm=caret::confusionMatrix(table(final_statusLUSC32$Predição, final_statusLUSC32$Real))
print(cfm)

#datasetLUSC32 with Validation LUSC-KR
# first we split between training and testing sets
write_csv(datasetLUSC32, file = "lusc_train.csv")
write_csv(validacaoLUSC32, file = "lusc_test.csv")

TRAIN_DATA_URL<- "../data/lusc_train.csv"
TEST_DATA_URL<- "../data/lusc_test.csv"

train_file_path <- get_file("lusc_train.csv",TRAIN_DATA_URL )
test_file_path <- get_file("lusc_test.csv",TEST_DATA_URL )


train_dataset <- make_csv_dataset(
  train_file_path, 
  field_delim = ",",
  batch_size = 6, 
  num_epochs = 1,
)

test_dataset <- train_dataset <- make_csv_dataset(
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
real[,1]
View(pred_dataset )
##prediction 
pred_dataset = round(test_predictions)

final_statusVal32=cbind (real, pred_dataset )
colnames(final_statusVal32) = c("Real","Predição")

library(caret)
cfm=caret::confusionMatrix(table(final_statusVal32$Predição, final_statusVal32$Real))
print(cfm)

#datasetLUSC15
# first we split between training and testing sets
split <- initial_split(datasetLUSC15, prop = 4/5)
train <- training(split)
write_csv(train, file = "lusc_train.csv")
test <- testing(split)
write_csv(test, file = "lusc_test.csv")


TRAIN_DATA_URL<- "../data/lusc_train.csv"
TEST_DATA_URL<- "../data/lusc_test.csv"

train_file_path <- get_file("lusc_train.csv",TRAIN_DATA_URL )
test_file_path <- get_file("lusc_test.csv",TEST_DATA_URL )


train_dataset <- make_csv_dataset(
  train_file_path, 
  field_delim = ",",
  batch_size = 6, 
  num_epochs = 1,
)

test_dataset <- train_dataset <- make_csv_dataset(
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
real[,1]
View(pred_dataset )
##prediction 
pred_dataset = round(test_predictions)

final_statusLUSC15=cbind (real, pred_dataset )
colnames(final_statusLUSC15) = c("Real","Predição")

library(caret)
cfm=caret::confusionMatrix(table(final_statusLUSC15$Predição, final_statusLUSC15$Real))
print(cfm)

#datasetLUSC15 with Validation LUSC-KR
# first we split between training and testing sets
write_csv(datasetLUSC15, file = "lusc_train.csv")
write_csv(validacaoLUSC15, file = "lusc_test.csv")

TRAIN_DATA_URL<- "../data/lusc_train.csv"
TEST_DATA_URL<- "../data/lusc_test.csv"

train_file_path <- get_file("lusc_train.csv",TRAIN_DATA_URL )
test_file_path <- get_file("lusc_test.csv",TEST_DATA_URL )


train_dataset <- make_csv_dataset(
  train_file_path, 
  field_delim = ",",
  batch_size = 6, 
  num_epochs = 1,
)

test_dataset <- train_dataset <- make_csv_dataset(
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

test_predictions <- model %>% predict(validacaoLUSC15, select=c(vital_status))
test_predictions[ , 1]

real <- subset(validacaoLUSC15, select= c(vital_status)) 
real[,1]
View(pred_dataset )
##prediction 
pred_dataset = round(test_predictions)

final_statusVal15=cbind (real, pred_dataset )
colnames(final_statusVal15) = c("Real","Predição")

library(caret)
cfm=caret::confusionMatrix(table(final_statusVal15$Predição, final_statusVal15$Real))
print(cfm)

#ROC curve Figure 4B
library(pROC)

rocobj1 <- plot.roc(final_statusLUSCLit$Real, final_statusLUSCLit$Predição, percent=TRUE, col="#000000" ,ci=TRUE)
rocobj2 <- lines.roc(final_statusLUSC32$Real, final_statusLUSC32$Predição, percent=TRUE, col="#008600" ,ci=TRUE)
rocobj3 <- lines.roc(final_statusLUSC15$Real, final_statusLUSC15$Predição, percent=TRUE, col="#0000FF" ,ci=TRUE)
rocobj4 <- lines.roc(final_statusVal32$Real, final_statusVal32$Predição, percent=TRUE, col="#800000" ,ci=TRUE)
rocobj5 <- lines.roc(final_statusVal15$Real, final_statusVal15$Predição, percent=TRUE, col="#FFFF00" ,ci=TRUE)

legend("bottomright", legend=c("LUSC lit","LUSC  32","LUSC 15","LUSC-KR 32","LUSC-KR 15"), 
       col=c("#000000","#008600","#0000FF","#800000","#FFFF00" ), lwd=2)

