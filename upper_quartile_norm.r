# upper quartile normalization
upperQuartileNorm <- function(dat) {
    stopifnot(is.matrix(dat))
    # divide each shRNA count by the upper quartile (75th percentile)
    upperQ <- apply(eo.trim, 2, FUN = fivenum)[4,]
    dat.norm <- t(t(dat)/upperQ)
    dat.norm.scale <- round(dat.norm * mean(temp), 0)
    
}