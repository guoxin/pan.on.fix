load("readpan20120114 md5.cc4a88a33be3296333121e9c4d524a8f")
load("fixAlleleData")
source("internals.R")
beta.set <- 0.02 * (1:5)
lambda.set <- exp(seq(-16,-10))
alleles <- names(fixAlleleData)
alleleStrings <- readpan20120114$allele.set[alleles]

K1 <- readpan20120114$K1

peptideList <- c()
alleleList  <- c()
rValues     <- c()
cvInds      <- c()

for(al in alleles){
    peptideList <- c(peptideList, fixAlleleData[[al]]$peptides)
    alleleList <- c(alleleList, rep(al, length(fixAlleleData[[al]]$peptides)))
    rValues <- c(rValues, fixAlleleData[[al]]$rValues)
    cvInds <- c(cvInds, fixAlleleData[[al]]$cvInds)
}
peptideStrings <- unique(peptideList)
names(peptideStrings) <- peptideStrings

fdpa.report <- LooPanReport(
    peptideStrings = peptideStrings,
    alleleStrings = alleleStrings,
    peptideList = peptideList,
    alleleList = alleleList,
    K1 = readpan20120114$K1,
    rValues = rValues,
    cvInds = cvInds,
    beta.set = beta.set,
    lambda.set = lambda.set,
    threshold = (1 - log(500)/log(50000))
)
save(fdpa.report, file = "fdpa.report.Rdata")
