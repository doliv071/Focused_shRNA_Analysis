####Programs for David
####Prepared by Piaomu and Edsel

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

LikRatio <- function(d1 = dat1, d2 = dat2) {
    out1 = LREstimates(d1)
    out2 = LREstimates(d2)
    ###UnRestricted Estimates
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
    ###Restricted Estimates
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
    ###Negative of Log-Likelihood Ratio
    NegLogLR = (n1 + n2)*(log(sig2hatres) - log(sig2hat))
    UnRestrict = list(n1 = n1, n2 = n2, b10 = b10, b11 = b11, b20 = b20, 
                      b21 = b21, sig2hat = sig2hat)
    Restrict=list(b10res = b10res, b11res = b11res, b20res = b20res, 
                  b21res = b21res, sig2hatres = sig2hatres)

    return(list(UnRestrict=UnRestrict,Restrict=Restrict,NegLogLR=NegLogLR))
}

###
###Generates two data sets for testing purposes
###
# DatGener <- function(n1=100,b10=2,b11=-2,n2=100,b20=-3,b21=4,sig2=.5,plotok=F)
# {
# x1 = rnorm(n1)
# err1=rnorm(n1,0,sig2)
# y1=b10+b11*x1+err1
# dat1=data.frame(x=x1,y=y1)
# 
# x2 = rnorm(n2)
# err2=rnorm(n2,0,sig2)
# y2=b20+b21*x2+err2
# dat2=data.frame(x=x2,y=y2)
# 
# if(plotok) {
# close.screen(all.screens=T)
# split.screen(c(1,2))
# screen(1)
# plot(x1,y1)
# screen(2)
# plot(x2,y2)
# }
# 
# return(list(dat1=dat1,dat2=dat2))
# }

###
### Sampling Distribution of -2*log(LR)
###

SampDist <- function(M=10000,n1=100,b10=2,b11=-2,n2=100,b20=-3,b21=4,sig2=.5,plotok=F) {
    NegLRStats = rep(0, M)
    for(m in 1:M) {
        datout=DatGener(n1=n1,b10=b10,b11=b11,n2=n2,b20=b20,b21=b21,sig2=sig2,plotok=plotok)
        outLR=LikRatio(datout$dat1,datout$dat2)
        NegLRStats[m] = outLR$NegLogLR
    }

    close.screen(all.screens=T)
    hist(NegLRStats,main="Histogram of Negative of Log of LR")

    return(NegLRStats)
}

####
#### P-Value via Bootstrap
####
BootPValue <- function(d1=dat1,d2=dat2,BootRep=10000) {

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
# 
#     close.screen(all.screens = T)
# 
#     hist(BootVals, col = "gray", breaks = 100)
#     abline(v = obsNegLogLR, col="red")
#
    bootPvalue = sum((BootVals >= obsNegLogLR))/BootRep

    return(list(obsval = obsNegLogLR, bpvalue = bootPvalue, cancerIntercept = b10, 
                cancerSlope = b11, normalIntercept = b20, normalSlope = b21, cancerData = d1, normalData = d2))
}
