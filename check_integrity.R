#Script en R para limpipar datos:

cnames = c("CHR","SNP","BP","A1","A2","FREQ","BETA","SE","P")
args <- commandArgs(trailingOnly=TRUE)
gwas <- args[1]
out <- args[2]
xPx <-as.numeric(args[3])
yPx <- as.numeric(args[4])
factor_division <- as.numeric(args[5])


library(data.table)

data <- fread(gwas, data.table=FALSE);

checkColumns <- all(cnames %in% colnames(data));
if(!checkColumns) {
  missingCols <- which(!(cnames %in% colnames(data)))
  message <- paste(cnames[missingCols], collapse=" and ")
  stop(paste("Error: Columns", message, "are missing from your data." ))
}

data <- data[complete.cases(data),]
if(nrow(data)<10) {
  stop("Many missing data or very few SNPs in data file")
}

checkFreq <- all(is.numeric(data$FREQ))
if(!checkFreq) {
  stop("FREQ column should contain only numbers")
}

checkFreqBounds <- all(data$FREQ >0 & data$FREQ <1)
if(!checkFreqBounds) {
  stop("FREQ column should contain only numbers between 0 and 1 exclusive")
}

checkBeta <- all(is.numeric(data$BETA))
if(!checkBeta) {
  stop("BETA column should contain only numbers")
}

checkSE <- all(is.numeric(data$SE))
if(!checkSE) {
  stop("SE column should contain only numbers")
}

checkP <- all(is.numeric(data$P))
if(!checkP) {
  data$P <- as.numeric(data$P)
  data$P[data$P>=0 & data$P<1e-300] <-1e-300
  # stop("P-value column should contain only numbers")
}

checkPBounds <- all(data$P <=1 & data$P >0)
if(!checkPBounds) {
  stop("P-value column should contain only numbers between 0 and 1 inclusive")
}
print(head(data))

data<-data[order(data$P),]

write.table(data[,cnames], file=out, quote=F, row.names=F, col.names=T, sep="\t")

inicio_scrp <- as.numeric(Sys.time())


chrm <- c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15',
          '16','17','18','19','20','21','22','X','Y')

#posc <- as.matrix(c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,
#          145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,
#          90338345,83257441,80373285,58617616,64444167,46709983,50818468,156040895,57227415))


posc <- as.matrix(c(0,248956422,242193529,198295559,190214555,181538259,170805979,159345973,
          145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,
          90338345,83257441,80373285,58617616,64444167,46709983,50818468,156040895))

rownames(posc)<-chrm

data$manhattanPos <- data$BP + posc[as.character(data$CHR),1]

xLen <- as.integer(xPx/factor_division)
yLen <- as.integer(yPx/factor_division)

vMax <- -log10(min(data$P))
vMin <- -log10(max(data$P))
bpPerPixel_X = sum(posc[,1])/xLen
valuesPerPixel_y  = (vMax-vMin)/yLen

data$mX <- round(data$manhattanPos/bpPerPixel_X)
data$mY <- round(-log10(data$P)/valuesPerPixel_y)

d2<-data[!duplicated(data[,c("mX","mY")]),]

cnames <- c(cnames,'AbslPost')

write.table(d2, file=paste0(out,".manhattan"), quote=F, row.names=F, col.names=T, sep="\t")

print(as.numeric(Sys.time())-inicio_scrp)
