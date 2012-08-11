rm(list=ls())
load("fdpa.report.Rdata")
f <- function(x)return(signif(as.double(x), digits = 5))
load("readpan20120114 md5.cc4a88a33be3296333121e9c4d524a8f")
load("fixAlleleData")
source("internals.R")
allele.setNames <- alleles <- names(fixAlleleData)
auc <- double(14)
K1 <- readpan20120114$K1

peptideList <- alleleList <- rValues <- cvInds <- c()
for(al in alleles){
    peptideList <- c(peptideList, fixAlleleData[[al]]$peptides)
    alleleList <- c(alleleList, rep(al, length(fixAlleleData[[al]]$peptides)))
    rValues <- c(rValues, fixAlleleData[[al]]$rValues)
    cvInds <- c(cvInds, fixAlleleData[[al]]$cvInds)
}
peptideStrings <- unique(peptideList)
names(peptideStrings) <- peptideStrings
names(auc) <- allele.setNames
weights <- rmse <- auc
for(nm in allele.setNames){
    index = (alleleList == nm)
    auc[nm] <- getAUC.guoxin(predicted = fdpa.report$predict[index],
                   realValue = rValues[index],
                   threshold = 1-(log(500)/log(50000)) # 0.425625189808507
               )
    rmse[nm] <- getRMSE(predicted = fdpa.report$predict[index],
                   realValue = rValues[index])
    weights[nm] <- sum(index)
    write(paste(nm, weights[nm], f(rmse[nm]), f(auc[nm]), sep = " & "), file="")
}
print("average")
print("auc")
print(mean(auc))
print("rmse")
print(mean(rmse))
print("weighted average")
weights <- weights/sum(weights)
print("auc")
print(sum(auc * weights))
print("rmse")
print(sum(rmse * weights))
