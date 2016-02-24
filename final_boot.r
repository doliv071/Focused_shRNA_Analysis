### set the working directory to where your raw data is
setwd("~/Documents/eCAT/Focused shRNA Library Project/Boot_Focused_Analysis/")

### get some packages
library(MASS); library(reshape2)
source('~/R/Useful_Functions/zipper.r')
source('PenaRProgs.R')
source('~/R/Useful_Functions/upper_quartile_norm.r')

### read in raw data and fix format
Emily.Output <- read.csv(file = 'raw_data.csv',header=T)
eo.trim <- Emily.Output
### make names and get rid of name column
row.names(eo.trim) <- eo.trim[,1]; eo.trim <- eo.trim[,-1]; eo.trim <- as.matrix(eo.trim)
### perform upper quartile normalization 
eo.trim <- eo.trim+1
eo.trim <- upperQuartileNorm(eo.trim)

### make a universal unique Gene ID list, shRNA list, and unique shRNAs list
geneids <- read.table(file = 'uniq_gene_names.txt', header = F, stringsAsFactors = F)
uniq_ids <- unique(geneids)

### full list of gene names
name.unstack <- NULL
for(i in 1:nrow(uniq_ids)){
    reps <- NULL
    if (uniq_ids[i,] == "RPL36AP37"){
        reps <- rep(uniq_ids[i,], 5)
    }
    else if (uniq_ids[i,] == "ATP5J2-PTCD1"){
        reps <- rep(uniq_ids[i,], 2)
    }
    else if (uniq_ids[i,] == "DAZ2"){
        reps <- rep(uniq_ids[i,], 3)
    }
    else if (uniq_ids[i,] == "SYNJ2BP-COX16"){
        reps <- rep(uniq_ids[i,], 4)
    }
    else if (uniq_ids[i,] == "LOC100506455"){
        reps <- rep(uniq_ids[i,], 2)
    }
    else if (uniq_ids[i,] == "LOC100653008"){
        reps <- rep(uniq_ids[i,], 3)
    }
    else{
        reps <- rep(uniq_ids[i,], 6)
    }
    name.unstack <- rbind(name.unstack, as.matrix(reps))
}

### split dataset into individual cancer cell lines
data_MDA <- eo.trim[,1:12]; data_BJ <- eo.trim[,13:23]; 
data_HT1080 <- eo.trim[,24:34]; data_HCT116 <- eo.trim[,35:46];
data_PC3 <- eo.trim[,47:ncol(eo.trim)]

### get rid of 6th day.
data_MDA_good <- data_MDA[,1:10];
data_BJ_good <- data_BJ[,1:10]; 
data_HT1080_good <- data_HT1080[,1:9]; 
data_HT1080_good <- cbind(data_HT1080_good, as.matrix(rowMeans(data_HT1080[,10:11])))
colnames(data_HT1080_good) <- c("HT1080.1A", "HT1080.1B", "HT1080.2A", "HT1080.2B", "HT1080.3A",
                                "HT1080.3B", "HT1080.4A", "HT1080.4B", "HT1080.5A", "HT1080.5B")
data_HCT116_good <- data_HCT116[,1:10]; 
data_PC3_good <- data_PC3[,1:10]

datasets <- list(data_PC3_good, data_BJ_good, data_MDA_good, data_HT1080_good, data_HCT116_good)
### stack bio-rep names
name.stack <- NULL
for(i in 1:nrow(uniq_ids)){
    reps <- NULL
    if (uniq_ids[i,] == "RPL36AP37"){
        reps <- rep(uniq_ids[i,], 10)
    }
    else if (uniq_ids[i,] == "ATP5J2-PTCD1"){
        reps <- rep(uniq_ids[i,], 4)
    }
    else if (uniq_ids[i,] == "DAZ2"){
        reps <- rep(uniq_ids[i,], 6)
    }
    else if (uniq_ids[i,] == "SYNJ2BP-COX16"){
        reps <- rep(uniq_ids[i,], 8)
    }
    else if (uniq_ids[i,] == "LOC100506455"){
        reps <- rep(uniq_ids[i,], 4)
    }
    else if (uniq_ids[i,] == "LOC100653008"){
        reps <- rep(uniq_ids[i,], 6)
    }
    else{
        reps <- rep(uniq_ids[i,], 12)
    }
    name.stack <- rbind(name.stack, as.matrix(reps))
}

### combine biological replicates
### first, split all cell lines into bio-rep1 and bio-rep2
MDA.bio1 <- NULL; MDA.bio2 <- NULL; BJ.bio1 <- NULL; BJ.bio2 <- NULL;
HT.bio1 <- NULL; HT.bio2 <- NULL; HCT.bio1 <- NULL; HCT.bio2 <- NULL;
PC3.bio1 <- NULL; PC3.bio2 <- NULL
for (i in 1:(ncol(data_MDA_good)/2)){
    MDA.bio1 <- cbind(MDA.bio1, data_MDA_good[,(2*i)-1])
    MDA.bio2 <- cbind(MDA.bio2, data_MDA_good[,2*i])
    BJ.bio1  <- cbind(BJ.bio1, data_BJ_good[,(2*i)-1])
    BJ.bio2  <- cbind(BJ.bio2, data_BJ_good[,2*i])
    HT.bio1  <- cbind(HT.bio1, data_HT1080_good[,(2*i)-1])
    HT.bio2  <- cbind(HT.bio2, data_HT1080_good[,2*i])
    HCT.bio1 <- cbind(HCT.bio1, data_HCT116_good[,(2*i)-1])
    HCT.bio2 <- cbind(HCT.bio2, data_HCT116_good[,2*i])
    PC3.bio1 <- cbind(PC3.bio1, data_PC3_good[,(2*i)-1])
    PC3.bio2 <- cbind(PC3.bio2, data_PC3_good[,2*i])
}
### restack bioreps.
MDA.stacked <- zipper(MDA.bio1, MDA.bio2, along=1); 
rownames(MDA.stacked) <- make.unique(rownames(MDA.stacked));
MDA.stacked <- as.data.frame(MDA.stacked);
MDA.stacked <- cbind(name.stack, MDA.stacked); 

BJ.stacked  <- zipper(BJ.bio1,  BJ.bio2,  along=1)
rownames(BJ.stacked) <- make.unique(rownames(BJ.stacked));
BJ.stacked <- as.data.frame(BJ.stacked);
BJ.stacked <- cbind(name.stack, BJ.stacked); 

HT.stacked  <- zipper(HT.bio1,  HT.bio2,  along=1)
rownames(HT.stacked) <- make.unique(rownames(HT.stacked));
HT.stacked <- as.data.frame(HT.stacked);
HT.stacked <- cbind(name.stack, HT.stacked); 

HCT.stacked <- zipper(HCT.bio1, HCT.bio2, along=1)
rownames(HCT.stacked) <- make.unique(rownames(HCT.stacked));
HCT.stacked <- as.data.frame(HCT.stacked);
HCT.stacked <- cbind(name.stack, HCT.stacked); 

PC3.stacked <- zipper(PC3.bio1, PC3.bio2, along=1)
rownames(PC3.stacked) <- make.unique(rownames(PC3.stacked));
PC3.stacked <- as.data.frame(PC3.stacked);
PC3.stacked <- cbind(name.stack, PC3.stacked); 


gendat <- function(cancer = cancer, normal = normal, gene = "METTL18", BootRep = 10000) {
    dat1 <- melt(cancer[which(cancer$name.stack == gene),])[,2:3]
    dat2 <- melt(normal[which(normal$name.stack == gene),])[,2:3]
    colnames(dat1) <- c("x", "y")
    colnames(dat2) <- c("x", "y")
    dat1$x <- as.numeric(dat1$x)
    dat2$x <- as.numeric(dat2$x)
    BootPValue(d1 = dat1, d2 = dat2, BootRep = BootRep)
}


PC3.out <- list()
for (i in 1:nrow(uniq_ids)) {
    PC3.out[[i]] <- gendat(cancer = PC3.stacked, normal = BJ.stacked, 
                           gene = uniq_ids[,1][i], BootRep = 10000)
}

MDA.out <- list()
for (i in 1:nrow(uniq_ids)) {
    MDA.out[[i]] <- gendat(cancer = MDA.stacked, normal = BJ.stacked, 
                           gene = uniq_ids[,1][i], BootRep = 10000)
}

HT.out <- list()
for (i in 1:nrow(uniq_ids)) {
    HT.out[[i]] <- gendat(cancer = HT.stacked, normal = BJ.stacked, 
                          gene = uniq_ids[,1][i], BootRep = 10000)
}

HCT.out <- list()
for (i in 1:nrow(uniq_ids)) {
    HCT.out[[i]] <- gendat(cancer = HCT.stacked, normal = BJ.stacked, 
                           gene = uniq_ids[,1][i], BootRep = 10000)
}
names(PC3.out) <- uniq_ids[,1]; names(MDA.out) <- uniq_ids[,1];
names(HT.out)  <- uniq_ids[,1]; names(HCT.out) <- uniq_ids[,1];

# lots of plots
Boot.out <- list(PC3.out, MDA.out, HT.out, HCT.out)
names(Boot.out) <- c("PC3_Boot_Results", "MDA_Boot_Results", "HT1080_Boot_Results", "HCT116_Boot_Results")
for (j in 1:length(Boot.out)){
    pdf(file = paste("~/Desktop/", names(Boot.out)[j], "_plots.pdf", sep = ""))
    for(i in 1:length(Boot.out[[j]])) {
        plot(as.data.frame(Boot.out[[j]][[i]][7])[,1], 
             as.data.frame(Boot.out[[j]][[i]][7])[,2], pch = 0, 
             main = paste(names(Boot.out[[1]])[i]),
             xlab = "Passages", ylab = "Normalized shRNA Counts")
        abline(coef = c(Boot.out[[j]][[i]][3], Boot.out[[j]][[i]][4]), col = "blue", lwd = 2)
        points(as.data.frame(Boot.out[[j]][[i]][8])[,1], as.data.frame(Boot.out[[j]][[i]][8])[,2], pch = 1)
        abline(coef = c(Boot.out[[j]][[i]][5], Boot.out[[j]][[i]][6]), col = "red", lwd = 2)
        legend("topright", legend=c("Cancer Counts", "Cancer Slope", "Normal Counts", "Normal Slope", NA, 
                                    paste("Pval = ", as.character(Boot.out[[j]][[i]][2]), sep ="")), 
               pch=c(0,NA,1,NA,NA,NA), lty = c(NA,1,NA,1,NA,NA), col = c("black", "blue", "black", "red", NA, NA),
               cex = 0.75, bg = "white")
    }
    dev.off()
}

for (j in 1:length(Boot.out)){
    final.out <- NULL
    for (i in 1:length(Boot.out[[j]])) {
        if (is.null(final.out)){
            final.out <- cbind(names(Boot.out[[j]])[i], as.data.frame(Boot.out[[j]][[i]][1:6]))
        }
        else {
            final.out <- rbind(final.out, cbind(names(Boot.out[[j]])[i], as.data.frame(Boot.out[[j]][[i]][1:6])))
        }
    }
    write.table(final.out, file = paste("~/Desktop/", names(Boot.out)[j], ".txt", sep = ""), row.names = F)
}


# final.out <- cbind(PC3.out[order(PC3.out$name),], cbind(MDA.out[order(MDA.out$name),2], cbind(HCT.out[order(HCT.out$name),2], HT.out[order(HT.out$name),2])))
# colnames(final.out) <- c("Gene", "PC3.Pvalue","MDA.Pvalue","HCT116.Pvalue","HT1080.Pvalue")
# write.table(final.out, file = "/home/david/Desktop/Focused_Analysis/final.boot.txt", sep = "\t", row.names = F)

################################################################################
####################### End of bootstrapped LRT analysis #######################
################################################################################

################################################################################
########################## Start edgeR analysis ################################
################################################################################
setwd("~/Documents/eCAT/Focused shRNA Library Project/Boot_Focused_Analysis/")
library ("edgeR")
Emily.Output <- read.csv(file = 'raw_data.csv',header=T)
eo.trim <- Emily.Output
### make names and get rid of name column
row.names(eo.trim) <- eo.trim[,1]
eo.trim <- eo.trim[,-1]
eo.trim <- as.matrix(eo.trim)
# remove pass 6
eo.trim <- eo.trim[,c(1:10, 13:22, 24:31, 33:34, 35:44, 47:56)]


Focused <- DGEList(eo.trim)
Focused$samples$group <- c(rep("MDA", 10), rep("BJ", 10), rep("HT", 10), 
                           rep("HCT", 10), rep("PC3", 10))
Focused$samples$time <- rep(1:5, each = 2, times = 5)
Focused$samples$bio.rep <- rep(1:2, times = 25)
sel <- rowSums(cpm(Focused$counts) > 0.5) >= 3
Focused <- Focused[ sel, ]

par(mfrow = c(2, 1))
barplot(colSums(Focused$counts), las = 2, main = "Counts per Sample", cex.names = 0.5, 
        cex.axis = 0.8, ylim = c(0, max(Focused$samples$lib.size)*1.2))

barplot(c(rowSums(Focused$counts[,1:10]), rowSums(Focused$counts[,11:20]), 
          rowSums(Focused$counts[,21:30]), rowSums(Focused$counts[,31:40]),
          rowSums(Focused$counts[,41:50])), las = 2, main = "Counts per shRNA", 
        cex.names = 0.5, cex.axis = 0.8, names.arg = c(rep("MDA", 1269), rep("BJ", 1269),
                                                   rep("HT", 1269), rep("HCT", 1269),
                                                   rep("PC3", 1269)))
par(mfrow = c(1,1), pty = "s")
plotMDS(Focused, labels = rownames(Focused$samples))

Focused <- calcNormFactors(Focused, method = "TMM")
design <- model.matrix(~Focused$samples$group*Focused$samples$time)
Focusedglm = estimateDisp(Focused, design)
plotBCV(Focusedglm, main = "Focused screen: BCV Plot")
fit = glmFit(Focusedglm, design)
HCT.lrt <- glmLRT(fit, coef = 7)
HT.lrt  <- glmLRT(fit, coef = 8) 
MDA.lrt <- glmLRT(fit, coef = 9)
PC3.lrt <- glmLRT(fit, coef = 10)

thresh = 0.05
top.HCT = topTags(HCT.lrt, n = Inf)
top.HCT.ids = rownames(top.HCT$table[top.HCT$table$FDR < thresh,]) 
plotSmear(HCT.lrt, de.tags = top.HCT.ids, pch = 20, cex = 0.6, main = "Focused screen: Slope vs logCPM")
abline(h = c(-1, 0, 1), col = c("dodgerblue", "yellow", "dodgerblue"), lty = 2)

genesymbols <- name.unstack[sel,]
genesymbollist = list()

for (i in unique(genesymbols)) genesymbollist[[i]] = which(genesymbols == i)


roast.HCT = mroast(Focusedglm, index = genesymbollist, design, contrast = 7,  nrot = 10000)
roast.HT  = mroast(Focusedglm, index = genesymbollist, design, contrast = 8,  nrot = 10000)
roast.MDA = mroast(Focusedglm, index = genesymbollist, design, contrast = 9,  nrot = 10000)
roast.PC3 = mroast(Focusedglm, index = genesymbollist, design, contrast = 10, nrot = 10000)

venn.list <- list(rownames(roast.HCT[which(roast.HCT$FDR.Mixed < 0.05 & roast.HCT$Direction == "Down"),]), 
                  rownames(roast.HT[which(roast.HT$FDR.Mixed   < 0.05 & roast.HT$Direction == "Down"),]),
                  rownames(roast.MDA[which(roast.MDA$FDR.Mixed < 0.05 & roast.MDA$Direction == "Down"),]),
                  rownames(roast.PC3[which(roast.PC3$FDR.Mixed < 0.05 & roast.PC3$Direction == "Down"),]))
library("gplots")
venn.out <- venn(venn.list)


pdf(file = "~/Desktop/HCT_edgeR.pdf")
par(mfrow = c(3,1))
for(i in 1:length(unique(genesymbols))){
    barcodeplot(HCT.lrt$table$logFC, index = genesymbollist[[i]], 
                main = paste("Barcodeplot for ", names(genesymbollist[i])),
                labels = c("Positive logFC", "Negative logFC"))
}
dev.off()


################################################################################
####################### End of edgeR shRNA-seq analysis ########################
################################################################################

################################################################################
##################### Start Multiple T-test analysis ###########################
################################################################################
setwd("~/Documents/eCAT/Focused shRNA Library Project/Boot_Focused_Analysis/")

### get some packages
library(MASS); library(reshape2)
source('~/R/Useful_Functions/zipper.r')
source('PenaRProgs.R')
source('~/R/Useful_Functions/upper_quartile_norm.r')
#install if necessary
source('http://bioconductor.org/biocLite.R')
biocLite('preprocessCore')
#load package
library(preprocessCore)


### read in raw data and fix format
Emily.Output <- read.csv(file = 'raw_data.csv',header=T)
eo.trim <- Emily.Output
### make names and get rid of name column
row.names(eo.trim) <- eo.trim[,1]
eo.trim <- eo.trim[,-1]
eo.trim <- as.matrix(eo.trim)
# remove pass 6
eo.trim <- eo.trim[,c(1:10, 13:22, 24:31, 33:34, 35:44, 47:56)]
rn.eo <- rownames(eo.trim)
cn.eo <- colnames(eo.trim)
eo.trim <- normalize.quantiles(eo.trim)
rownames(eo.trim) <- rn.eo
colnames(eo.trim) <- cn.eo
eo.trim <- as.data.frame(eo.trim)

temp <- rownames(eo.trim)
for (i in 1:9){
    temp[i] <- substr(temp[i], 2, nchar(temp[i])-1)
}
for (i in 10:99){
    temp[i] <- substr(temp[i], 3, nchar(temp[i])-2)
}
for (i in 100:999){
    temp[i] <- substr(temp[i], 4, nchar(temp[i])-3)
}
for(i in 1000:length(temp)){
    temp[i] <- substr(temp[i], 5,nchar(temp[i])-4)
}
temp <- as.factor(temp)

eo.trim.2 <- cbind(temp, as.data.frame(eo.trim))
colnames(eo.trim.2)[1] <- "GeneID"

# function does multiple T-test...
multi.t <- function( test.data, cell.line, gene, GeneID ) {
    all.tests <- NULL
    cancer <- test.data[which(test.data$GeneID == gene), grep(cell.line, colnames(test.data))]
    normal <- test.data[which(test.data$GeneID == gene), grep("BJ", colnames(test.data))]
    for (i in 1:(ncol(cancer)/2)) {
        stat  <- t.test(stack(cancer[ , (i*2-1):(i*2)])$values, 
                        stack(normal[ , (i*2-1):(i*2)])$values)$statistic
        p.val <- t.test(stack(cancer[ , (i*2-1):(i*2)])$values, 
                        stack(normal[ , (i*2-1):(i*2)])$values)$p.value
        if (is.null(all.tests)) {
            all.tests <- data.frame(Gene = gene, t.stat = stat, P.value = p.val)
        }
        else {
            all.tests <- rbind(all.tests, data.frame(Gene = gene, t.stat = stat, P.value = p.val))
        }
    }
    
    return(all.tests)
}

all.genes <- list()
Cells <- list("MDA", "HT", "HCT", "PC")
Genes <- list()
for (i in 1:length(Cells)){
    temp.j <- list()
    for(j in 1:length(unique(eo.trim.2$GeneID))){
        temp.j[[j]] <- multi.t(eo.trim.2, Cells[[i]], as.character(unique(eo.trim.2$GeneID)[j]), GeneID)
    }
    all.genes[[i]] <- temp.j
}

MB231.t  <- data.frame(matrix(unlist(all.genes[[1]]), nrow = 1075, byrow = T))
colnames(MB231.t) <- c("Gene", "t.stat", "P.value")
HT1080.t <- data.frame(matrix(unlist(all.genes[[2]]), nrow = 215, byrow = T))
colnames(HT1080.t) <- c("Gene", "t.stat", "P.value")
HCT116.t <- data.frame(matrix(unlist(all.genes[[3]]), nrow = 215, byrow = T))
colnames(HCT116.t) <- c("Gene", "t.stat", "P.value")
PC3.t    <- data.frame(matrix(unlist(all.genes[[4]]), nrow = 215, byrow = T))
colnames(PC3.t) <- c("Gene", "t.stat", "P.value")

head(MB231.t[order(MB231.t$P.value),])
head(HCT116.t[order(HCT116.t$P.value),])
head(HT1080.t[order(HT1080.t$P.value),])
head(PC3.t[order(PC3.t$P.value),])


pdf(file = "MDA_t-test.pdf")
for(i in 1:length(unique(eo.trim.2$GeneID))){
    plot(all.genes[[1]][[i]][,2], main = as.character(unique(eo.trim.2$GeneID)[i]),
         ylim = c(-3,3))
    abline(h = 0)
    if(any(all.genes[[1]][[i]][,3] < 0.05)) {
        points(which(all.genes[[1]][[i]][,3] < 0.05),
               all.genes[[1]][[i]][ which(all.genes[[1]][[i]][,3] < 0.05), 2], 
               pch = 20, col = "red")
    }
}
dev.off()

pdf(file = "HT1080_t-test.pdf")
for(i in 1:length(unique(eo.trim.2$GeneID))){
    plot(all.genes[[2]][[i]][,2], main = as.character(unique(eo.trim.2$GeneID)[i]),
         ylim = c(-3,3))
    abline(h = 0)
    if(any(all.genes[[2]][[i]][,3] < 0.05)) {
        points(which(all.genes[[2]][[i]][,3] < 0.05),
               all.genes[[2]][[i]][ which(all.genes[[2]][[i]][,3] < 0.05), 2], 
               pch = 20, col = "red")
    }
}
dev.off()

pdf(file = "HCT116_t-test.pdf")
for(i in 1:length(unique(eo.trim.2$GeneID))){
    plot(all.genes[[3]][[i]][,2], main = as.character(unique(eo.trim.2$GeneID)[i]),
         ylim = c(-3,3))
    abline(h = 0)
    if(any(all.genes[[3]][[i]][,3] < 0.05)) {
        points(which(all.genes[[3]][[i]][,3] < 0.05),
               all.genes[[3]][[i]][ which(all.genes[[3]][[i]][,3] < 0.05), 2], 
               pch = 20, col = "red")
    }
}
dev.off()

pdf(file = "PC3_t-test.pdf")
for(i in 1:length(unique(eo.trim.2$GeneID))){
    plot(all.genes[[4]][[i]][,2], main = as.character(unique(eo.trim.2$GeneID)[i]),
         ylim = c(-3,3))
    abline(h = 0)
    if(any(all.genes[[4]][[i]][,3] < 0.05)) {
        points(which(all.genes[[4]][[i]][,3] < 0.05),
               all.genes[[4]][[i]][ which(all.genes[[4]][[i]][,3] < 0.05), 2], 
               pch = 20, col = "red")
    }
}
dev.off()

################################################################################
####################### End of Multiple T-test analysis ########################
################################################################################

################################################################################
################### Start Bootstrap LRT redo analysis ##########################
################################################################################
setwd("~/Documents/eCAT/Focused shRNA Library Project/Boot_Focused_Analysis/")
load(".RData")
### get some packages
library(MASS); library(reshape2)
source('~/R/Useful_Functions/zipper.r')
source('PenaRProgs.R')
source('~/R/Useful_Functions/upper_quartile_norm.r')
library(preprocessCore)


### read in raw data and fix format
Emily.Output <- read.csv(file = 'raw_data.csv',header=T)
eo.trim <- Emily.Output
### make names and get rid of name column
row.names(eo.trim) <- eo.trim[,1]
eo.trim <- eo.trim[,-1]
eo.trim <- as.matrix(eo.trim)
# remove pass 6
eo.trim <- eo.trim[,c(1:10, 13:22, 24:31, 33:34, 35:44, 47:56)]
rn.eo <- rownames(eo.trim)
cn.eo <- colnames(eo.trim)
eo.trim <- normalize.quantiles(eo.trim)+15
eo.trim <- sqrt(eo.trim)
# eo.trim <- log2(eo.trim)
rownames(eo.trim) <- rn.eo
colnames(eo.trim) <- cn.eo
eo.trim <- as.data.frame(eo.trim)

temp <- rownames(eo.trim)
for (i in 1:9){
    temp[i] <- substr(temp[i], 2, nchar(temp[i])-1)
}
for (i in 10:99){
    temp[i] <- substr(temp[i], 3, nchar(temp[i])-2)
}
for (i in 100:999){
    temp[i] <- substr(temp[i], 4, nchar(temp[i])-3)
}
for(i in 1000:length(temp)){
    temp[i] <- substr(temp[i], 5,nchar(temp[i])-4)
}
temp <- as.factor(temp)

eo.trim.2 <- cbind(temp, as.data.frame(eo.trim))
colnames(eo.trim.2)[1] <- "GeneID"

### make a universal unique Gene ID list, shRNA list, and unique shRNAs list
geneids <- read.table(file = 'uniq_gene_names.txt', header = F, stringsAsFactors = F)
uniq_ids <- unique(geneids)


gendat.2 <- function(data = data, cancer = "cancer", normal = "normal", gene = "METTL18", BootRep = 10000) {
    dat1 <- melt(data[which(data$GeneID == gene), grep(cancer, colnames(data))])
    dat2 <- melt(data[which(data$GeneID == gene), grep(normal, colnames(data))])
    colnames(dat1) <- c("x", "y")
    colnames(dat2) <- c("x", "y")
    dat1$x <- as.numeric(dat1$x)
    dat2$x <- as.numeric(dat2$x)
    Boot.Result <- BootPValue(d1 = dat1, d2 = dat2, BootRep = BootRep)
    return(c(Gene = gene, Cancer.Slope = Boot.Result$cancerSlope, 
           Normal.Slope = Boot.Result$normalSlope, LRT = Boot.Result$obsval,
           P.value = Boot.Result$bpvalue ))
}

# function does bootstrapped analysis but with better normalization...

all.genes <- list()
Cells <- list("MDA", "HT", "HCT", "PC")
for (i in 1:length(Cells)){
    temp.j <- list()
    for(j in 1:dim(uniq_ids)[1]){
        temp.j[[j]] <- gendat.2(data = eo.trim.2, cancer = Cells[i], normal = "BJ", 
                                gene = uniq_ids[j,], BootRep = 10000)
    }
    all.genes[[i]] <- temp.j
}

MB231.lrt  <- data.frame(matrix(unlist(all.genes[[1]]), nrow = 215, byrow = T))
colnames(MB231.lrt) <- c("Gene", "Cancer.Slope", "Normal.Slope", "LRT", "P.value")
HT1080.lrt <- data.frame(matrix(unlist(all.genes[[2]]), nrow = 215, byrow = T))
colnames(HT1080.lrt) <- c("Gene", "Cancer.Slope", "Normal.Slope", "LRT", "P.value")
HCT116.lrt <- data.frame(matrix(unlist(all.genes[[3]]), nrow = 215, byrow = T))
colnames(HCT116.lrt) <- c("Gene", "Cancer.Slope", "Normal.Slope", "LRT", "P.value")
PC3.lrt    <- data.frame(matrix(unlist(all.genes[[4]]), nrow = 215, byrow = T))
colnames(PC3.lrt) <- c("Gene", "Cancer.Slope", "Normal.Slope", "LRT", "P.value")

log.MB231.lrt  <- data.frame(matrix(unlist(log2.genes[[1]]), nrow = 215, byrow = T))
colnames(log.MB231.lrt) <- c("Gene", "Cancer.Slope", "Normal.Slope", "LRT", "P.value")
log.HT1080.lrt <- data.frame(matrix(unlist(log2.genes[[2]]), nrow = 215, byrow = T))
colnames(log.HT1080.lrt) <- c("Gene", "Cancer.Slope", "Normal.Slope", "LRT", "P.value")
log.HCT116.lrt <- data.frame(matrix(unlist(log2.genes[[3]]), nrow = 215, byrow = T))
colnames(log.HCT116.lrt) <- c("Gene", "Cancer.Slope", "Normal.Slope", "LRT", "P.value")
log.PC3.lrt    <- data.frame(matrix(unlist(log2.genes[[4]]), nrow = 215, byrow = T))
colnames(log.PC3.lrt) <- c("Gene", "Cancer.Slope", "Normal.Slope", "LRT", "P.value")

write.table(x = log.MB231.lrt, file = "log2_trans_result_MB231.txt")
write.table(x = log.HT1080.lrt, file = "log2_trans_result_HT1080.txt")
write.table(x = log.HCT116.lrt, file = "log2_trans_result_HCT116.txt")
write.table(x = log.PC3.lrt, file = "log2_trans_result_PC3.txt")

write.table(x = MB231.lrt, file = "sqrt_trans_result_MB231.txt")
write.table(x = HT1080.lrt, file = "sqrt_trans_result_HT1080.txt")
write.table(x = HCT116.lrt, file = "sqrt_trans_result_HCT116.txt")
write.table(x = PC3.lrt, file = "sqrt_trans_result_PC3.txt")


### exploring the mean/var relationship.
# temp.mean <- NULL
# temp.var <- NULL
# for(i in 1:dim(uniq_ids)[1]){
#     temp.mean[i] <- mean(stack(eo.trim.2[which(eo.trim.2$GeneID == uniq_ids[i,]), grep("MDA", colnames(eo.trim.2))])$values)
#     temp.var[i]  <- var( stack(eo.trim.2[which(eo.trim.2$GeneID == uniq_ids[i,]), grep("MDA", colnames(eo.trim.2))])$values)
# }
# temp.2 <- cbind(temp.mean, temp.var)
# plot(temp.2)



