##----------------------------------------------------------------------------#
##data loading

#loading library
library(DNAshapeR)
if (!requireNamespace("Biostrings", quietly = TRUE)) {
    BiocManager::install("Biostrings")
}
library(Biostrings)
library(caret)

#reading in Mad, Max, Myc seq:
madSequence <- readDNAStringSet("/Users/aly/Downloads/QBIO481-master/gcPBM/Mad.txt.fa")
maxSequence <- readDNAStringSet("/Users/aly/Downloads/QBIO481-master/gcPBM/Max.txt.fa")
mycSequence <- readDNAStringSet("/Users/aly/Downloads/QBIO481-master/gcPBM/Myc.txt.fa")

#1mer
mad1Mer <- oligonucleotideFrequency(madSequence, width=1)
max1Mer <- oligonucleotideFrequency(maxSequence, width=1)
myc1Mer <- oligonucleotideFrequency(mycSequence, width=1)

#shape mad:
temporaryFile <- tempfile(fileext = ".fa")
writeXStringSet(madSequence, filepath = temporaryFile)

madShapeList <- tryCatch({
  getShape(temporaryFile)
}, error = function(e) {
  message("Error encountered: ", conditionMessage(e))
  NULL
})

file.remove(temporaryFile)

#shape max:
temporaryFile <- tempfile(fileext = ".fa")
writeXStringSet(maxSequence, filepath = temporaryFile)

maxShapeList <- tryCatch({
  getShape(temporaryFile)
}, error = function(e) {
  message("Error encountered: ", conditionMessage(e))
  NULL
})

file.remove(temporaryFile)

#shape myc:
temporaryFile <- tempfile(fileext = ".fa")
writeXStringSet(mycSequence, filepath = temporaryFile)

mycShapeList <- tryCatch({
  getShape(temporaryFile)
}, error = function(e) {
  message("Error encountered: ", conditionMessage(e))
  NULL
})

file.remove(temporaryFile)

#1mer+shape for mad
madShapeMatrix <- do.call(cbind, madShapeList)
mad1MerShape <- cbind(as.matrix(mad1Mer), madShapeMatrix)
#max
maxShapeMatrix <- do.call(cbind, maxShapeList)
max1MerShape <- cbind(as.matrix(max1Mer), maxShapeMatrix)
#myc
mycShapeMatrix <- do.call(cbind, mycShapeList)
myc1MerShape <- cbind(as.matrix(myc1Mer), mycShapeMatrix)

##----------------------------------------------------------------------------#
##caret, QUESTION 3

#reading in data+intensities
#mad
madData <- read.table("/Users/aly/Downloads/QBIO481-master/gcPBM/Mad.txt", sep = "\t", header = FALSE)
madSeq <- madData$V1
madIntensity <- madData$V2
#max
maxData <- read.table("/Users/aly/Downloads/QBIO481-master/gcPBM/Max.txt", sep = "\t", header = FALSE)
maxSeq <- maxData$V1
maxIntensity <- maxData$V2
#myc
mycData <- read.table("/Users/aly/Downloads/QBIO481-master/gcPBM/Myc.txt", sep = "\t", header = FALSE)
mycSeq <- mycData$V1
mycIntensity <- mycData$V2


#train ctrl 10fold crossval 
trainCtrl <- trainControl(method = "cv", number = 10)
set.seed(0)
#mad
modelMad1Mer <- train(mad1Mer, madIntensity, method = "ridge",
                      trControl = trainCtrl, metric = "Rsquared")
#max
modelMax1Mer <- train(max1Mer, maxIntensity, method = "ridge",
                      trControl = trainCtrl, metric = "Rsquared")
#myc
modelMyc1Mer <- train(myc1Mer, mycIntensity, method = "ridge",
                      trControl = trainCtrl, metric = "Rsquared")
#print
print(modelMad1Mer)
print(modelMax1Mer)
print(modelMyc1Mer)

#for printing out Rsquared values for 1-mer + shape, i tried for like 3 hours
#i couldn't get it to work :( the attempt is shown below

#fixing mad1MerShape
mad1MerShape <- mad1MerShape[, -c(5, 6, 39, 40, 41)]
mad1MerShape <- na.omit(mad1MerShape)
#modeling 1-mer+shape
modelMad1MerShape <- train(mad1MerShape, madIntensity, method = "rf",
                           trControl = trainCtrl, metric = "Rsquared")

#fixing max1MerShape
max1MerShape <- max1MerShape[, -c(5, 6, 39, 40, 41)]
max1MerShape <- na.omit(max1MerShape)
#modeling 1-mer+shape
modelMax1MerShape <- train(max1MerShape, maxIntensity, method = "rf",
                           trControl = trainCtrl, metric = "Rsquared")

#fixing myc1MerShape
myc1MerShape <- myc1MerShape[, -c(5, 6, 39, 40, 41)]
myc1MerShape <- na.omit(myc1MerShape)
#modeling 1-mer+shape
modelMyc1MerShape <- train(myc1MerShape, mycIntensity, method = "rf",
                           trControl = trainCtrl, metric = "Rsquared")
#print
print(modelMad1MerShape)
print(modelMax1MerShape)
print(modelMyc1MerShape)

##----------------------------------------------------------------------------#
##plotting, QUESTION 4

#i do not have 1-mer + shape data unfortunately, so this will just be code for a plot
madPlotData <- data.frame(
  model = c('1-mer', '1-mer', '1-mer', '1-mer+shape', '1-mer+shape', '1-mer+shape'),
  ind = c('RMSE', 'Rsquared', 'MAE', 'RMSE', 'Rsquared', 'MAE'),
  val = c(0.2100877, 0.04016394, 0.1701872, NA, NA, NA)
)

maxPlotData <- data.frame(
  model = c('1-mer', '1-mer', '1-mer', '1-mer+shape', '1-mer+shape', '1-mer+shape'),
  ind = c('RMSE', 'Rsquared', 'MAE', 'RMSE', 'Rsquared', 'MAE'),
  val = c(0.7966562, 0.01251410, 0.6468446, NA, NA, NA)
)

mycPlotData <- data.frame(
  model = c('1-mer', '1-mer', '1-mer', '1-mer+shape', '1-mer+shape', '1-mer+shape'),
  ind = c('RMSE', 'Rsquared', 'MAE', 'RMSE', 'Rsquared', 'MAE'),
  val = c(0.7905888, 0.006532783, 0.6396503, NA, NA, NA)
)

#plotting
ggplot(madPlotData, aes(x = ind, y = val, fill = model))
ggplot(maxPlotData, aes(x = ind, y = val, fill = model))
ggplot(mycPlotData, aes(x = ind, y = val, fill = model))

##----------------------------------------------------------------------------#
##downloading, QUESTION 5

#bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

#DNAshapeR
BiocManager::install("DNAshapeR")

#carat
install.packages("caret", dependencies = TRUE)

#bound.fa and unbound.fa obtained

##----------------------------------------------------------------------------#
##in vivo analysis, QUESTION 6

library(DNAshapeR)

#loading in data
boundSequence <- readDNAStringSet("/Users/aly/Downloads/QBIO481-master/CTCF/bound.fa")
unboundSequence <- readDNAStringSet("/Users/aly/Downloads/QBIO481-master/CTCF/unbound.fa")

#list intermediate
boundList <- as.list(boundSequence)
unboundList <- as.list(unboundSequence)

#shape list
boundShapeList <- lapply(boundList, getShape)
unboundShapeList <- lapply(unboundList, getShape)

#plotting
#bound 
boundEnsembleMGW <- ensembleMatrix(lapply(boundShapeList, function(x) x$MGW))
plotShape(boundEnsembleMGW, plot.type=1, ylab="MGW")
#unbound
unboundEnsembleMGW <- ensembleMatrix(lapply(unboundShapeList, function(x) x$MGW))
plotShape(unboundEnsembleMGW, plot.type=1, ylab="MGW")

##i tried to get MGW, couldn't because i cannot figure out shaping for the life of me
##it has been like 5 hours(?) i am struggling :(

##----------------------------------------------------------------------------#
##prediction models for in vitro data, question 7

library(pROC)

#log regression
madOut <- mad1Mer$target
pModelMad1Mer <- glm(madOut ~ ., data = mad1Mer, family = binomial)

maxOut <- max1Mer$target
pModelMax1Mer <- glm(maxOut ~ ., data = max1Mer, family = binomial)

mycOut <- myc1Mer$target
pModelMyc1Mer <- glm(mycOut ~ ., data = myc1Mer, family = binomial)

#combined shape data + logr
combinedShapes <- data.frame(mad1MerShape, max1MerShape, myc1MerShape)
modelAll1MerShape <- glm(madOut ~ ., data = combinedShapes, family = binomial)

#predicting probabilities
predictMad1Mer <-  predict(pModelMad1Mer, newdata = mad1Mer, type = "response")
predictMax1Mer <-  predict(pModelMax1Mer, newdata = max1Mer, type = "response")
predictMyc1Mer <-  predict(pModelMyc1Mer, newdata = myc1Mer, type = "response")

#predicting shape
predictShape <-  predict(modelAll1MerShape, newdata = combinedShapes, type = "response")

#ROC curve:
rocMad1Mer <- roc(madOut, predictMad1Mer)
rocMax1Mer <- roc(maxOut, predictMax1Mer)
rocMyc1Mer <- roc(mycOut, predictMyc1Mer)
rocAllShape <- roc(madOut, predictShape)

#plot
plot(rocMad1Mer, col = "red", main = "ROC Curves")
plot(rocMax1Mer, col = "green", main = "ROC Curves")
plot(rocMyx1Mer, col = "yellow", main = "ROC Curves")
lines(rocAllShape, col = "blue")

#auc
aucMad <- auc(rocMad1Mer)
aucMax <- auc(rocMax1Mer)
aucMyc <- auc(rocMyc1Mer)
aucShape <- auc(rocAllShape)

#print
print(aucMad)
print(aucMax)
print(aucMyc)
print(aucShape)