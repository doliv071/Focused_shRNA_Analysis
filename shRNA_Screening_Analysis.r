#  ____  _     _                                            _    _   _          
# |  _ \(_) __| |   _   _  ___  _   _    _ __ ___  __ _  __| |  | |_| |__   ___ 
# | | | | |/ _` |  | | | |/ _ \| | | |  | '__/ _ \/ _` |/ _` |  | __| '_ \ / _ \
# | |_| | | (_| |  | |_| | (_) | |_| |  | | |  __/ (_| | (_| |  | |_| | | |  __/
# |____/|_|\__,_|   \__, |\___/ \__,_|  |_|  \___|\__,_|\__,_|   \__|_| |_|\___|
#                  |___/                                                       
#                     _                  ___ 
#  _ __ ___  __ _  __| |_ __ ___   ___  |__ \
# | '__/ _ \/ _` |/ _` | '_ ` _ \ / _ \   / /
# | | |  __/ (_| | (_| | | | | | |  __/  |_| 
# |_|  \___|\__,_|\__,_|_| |_| |_|\___|  (_) 

shRNA_Screening_Analysis <- function(dat1 = "control.csv", dat2 = "treat.csv", 
                                     genes = "targets.txt", BootRep = 10000, thresh = 10, 
                                     passages = 5, bio.reps = 2, adjust = "BH") {
    stopifnot(file.exists(genes))
    gene.vector <- read.table(file = genes, stringsAsFactors = F)
    gene.vector <- gene.vector[,1]
    return.data <- NULL
    for(i in 1:length(gene.vector)){
        return.data <- rbind(return.data, shRNA_Screening_Function2(dat1 = dat1, dat2 = dat2, gene = gene.vector[i],  
                                                      passages = passages, thresh = thresh, 
                                                      bio.reps = bio.reps, BootRep = BootRep))
    }
    return.data <- as.data.frame(return.data)
    return.data$AdjustedPvalue <- p.adjust(return.data$Pvalue, method = adjust)
    rownames(return.data) <- gene.vector
    return(return.data)
}

##########################################################################################

shRNA_Screening_Function2 <- function(dat1 = dat1, dat2 = dat2, gene = gene, BootRep = BootRep, 
                                      passages = passages, thresh = thresh, bio.reps = bio.reps){
    stopifnot(is.character(gene), file.exists(dat1), file.exists(dat2), is.numeric(BootRep))
    norm.control.data <- shRNA_Screening_Function1(dat = dat1, thresh = thresh, passages = passages, bio.reps = bio.reps)
    norm.treated.data <- shRNA_Screening_Function1(dat = dat2, thresh = thresh, passages = passages, bio.reps = bio.reps)
    stopifnot(!any(norm.control.data$Var1 != norm.treated.data$Var1))
    norm.control.data.gene <- norm.control.data[grep(paste0("*",gene,"*"), norm.control.data$Var1),3:4]
    norm.control.data.gene <- norm.control.data.gene[which(complete.cases(norm.control.data.gene)),]
    norm.treated.data.gene <- norm.treated.data[grep(paste0("*",gene,"*"), norm.treated.data$Var1),3:4]
    norm.treated.data.gene <- norm.treated.data.gene[which(complete.cases(norm.treated.data.gene)),]
    BootPValue(d1 = norm.treated.data.gene, d2 = norm.control.data.gene, BootRep = BootRep)
}

##########################################################################################

shRNA_Screening_Function1 <- function(dat = "data.csv", thresh = thresh, passages = passages, bio.reps = bio.reps) {
    require(MASS); require(reshape2); require(preprocessCore)
    stopifnot(file.exists(dat), is.numeric(thresh), is.numeric(passages), is.numeric(bio.reps))
    dat <- as.matrix(read.csv(dat, row.names = 1))
    stopifnot(is.matrix(dat), !is.null(rownames(dat)))
    shRNAs <- rownames(dat); CN <- colnames(dat)
    dat[which(dat < thresh)] <- NA
    dat.2 <- normalize.quantiles(dat)
    rownames(dat.2) <- shRNAs
    colnames(dat.2) <- CN
    dat.3 <- melt(dat.2)
    dat.3$x <- rep(1:passages, each = bio.reps)
    colnames(dat.3)[3] <- "y"
    bc <- boxcox(dat.3$y~dat.3$x, plotit = F)
    dat.3$y <- dat.3$y^bc$x[which.max(bc$y)]
    return(dat.3)
}

##########################################################################################

BootPValue <- function(d1=dat1, d2=dat2, BootRep=10000) {
    
    out = LikRatio(d1,d2)
    obsNegLogLR = out$NegLogLR
    
    n1 = out$UnRestrict$n1
    n2 = out$UnRestrict$n2
    
    b10 = out$UnRestrict$b10
    b11 = out$UnRestrict$b11
    b20 = out$UnRestrict$b20
    b21 = out$UnRestrict$b21
    
    b10res = out$Restrict$b10
    b20res = out$Restrict$b20
    b11res = out$Restrict$b11
    b21res = out$Restrict$b21
    
    BootVals = rep(0,BootRep)
    resid1 = d1$y-(b10+b11*d1$x)
    resid2 = d2$y-(b20+b21*d2$x)
    allresids = c(resid1,resid2)
    
    for (m in 1:BootRep) {
        resstar1    = sample(allresids,n1)
        resstar2    = sample(allresids,n2)
        y1star      = b10res + b11res*d1$x + resstar1
        y2star      = b20res + b21res*d2$x + resstar2
        dat1star    = data.frame(x = d1$x, y = y1star)
        dat2star    = data.frame(x = d2$x,y = y2star)
        outstar     = LikRatio(dat1star, dat2star)
        BootVals[m] = outstar$NegLogLR
    }
    bootPvalue = sum((BootVals >= obsNegLogLR))/BootRep
    return(cbind(ControlSlope = b21, TreatedSlope = b11, LRTstatistics = obsNegLogLR, Pvalue = bootPvalue))
}

##########################################################################################

LikRatio <- function(d1 = dat1, d2 = dat2) {
    out1 = LREstimates(d1)
    out2 = LREstimates(d2)

    n1    = out1$n
    n2    = out2$n
    x1bar = out1$barx
    y1bar = out1$bary
    x2bar = out2$barx
    y2bar = out2$bary
    b10   = out1$b0
    b11   = out1$b1
    b20   = out2$b0
    b21   = out2$b1
    sig2hat1 = out1$sig2hat
    sig2hat2 = out2$sig2hat
    sig2hat  = (n1/(n1+n2))*sig2hat1 + (n2/(n1+n2))*sig2hat2

    b10res = b10
    b11res = b11
    b20res = b20
    b21res = b21
    sig2hatres = sig2hat
    rat1   = out1$SXX/(out1$SXX + out2$SXX)
    rat2   = 1 - rat1
    
    if((b11 < 0) && (b21 > 0) && (b21 < abs(b11))) {
        b1dhat = rat1*b11 - rat2*b21
        b11res = b1dhat
        b21res = -b1dhat
    }
    
    if((b11 < 0) && (b21 < 0) && (b21 > b11)) {
        b1dhat = rat1*b11 + rat2*b21
        b11res = b1dhat
        b21res = b1dhat
    }
    
    b10res = y1bar - b11res*x1bar
    b20res = y2bar - b21res*x2bar
    sig2hatres = (sum((d1$y - (b10res + b11res*d1$x))^2) + sum((d2$y - (b20res + b21res*d2$x))^2))/(n1 + n2)
    NegLogLR = (n1 + n2)*(log(sig2hatres) - log(sig2hat))
    UnRestrict = list(n1 = n1, n2 = n2, b10 = b10, b11 = b11, b20 = b20, 
                      b21 = b21, sig2hat = sig2hat)
    Restrict=list(b10res = b10res, b11res = b11res, b20res = b20res, 
                  b21res = b21res, sig2hatres = sig2hatres)
    
    return(list(UnRestrict=UnRestrict,Restrict=Restrict,NegLogLR=NegLogLR))
}

##########################################################################################

LREstimates <- function(dat) {
    x    = dat$x
    y    = dat$y
    n    = length(x)
    barx = mean(x)
    bary = mean(y)
    SXX  = (n - 1)*var(x)
    SYY  = (n - 1)*var(y)
    SXY  = (n - 1)*cov(x, y)
    b1   = SXY/SXX
    b0   = bary - b1*barx
    sig2hat = sum((y - (b0 + b1*x))^2)/n
    return(list(n = n, barx = barx, bary = bary, SXX = SXX, SYY = SYY, SXY = SXY,
                b1 = b1, b0 = b0, sig2hat = sig2hat))
}
