##'combining similarity matrix to similarity score
##'
##'Functions for combining similarity matrix to similarity score
##'
##'
##'@param SimScores similarity matrix
##'@param combine combine method
##'@return similarity value
##'@export
##'@author Guangchuang Yu \url{http://ygc.name}
combineScores <- function(SimScores, combine) {

    if (length(combine) == 0) {  #if not define combine
        return(round(SimScores, digits=3))
    }

    ## if combine was defined...
    if(!sum(!is.na(SimScores))) return (NA)

    if (is.vector(SimScores) || nrow(SimScores)==1 || ncol(SimScores)==1) {
        if (combine == "avg") {
            return(round(mean(SimScores, na.rm=TRUE), digits=3))
        } else {
            return (round(max(SimScores, na.rm=TRUE), digits=3))
        }
    }


    row.na.idx <- apply(SimScores, 1, function(i) all(is.na(i)))
    if (any(row.na.idx)) {
        SimScores <- SimScores[-which(row.na.idx), ]
    }

    col.na.idx <- apply(SimScores, 2, function(i) all(is.na(i)))
    if (any(col.na.idx)) {
        SimScores <- SimScores[ , -which(col.na.idx)]
    }

    if (combine        == "avg") {
        result   <- mean(SimScores, na.rm=TRUE)
    } else if (combine == "max") {
        result   <- max(SimScores, na.rm=TRUE)
    } else if (combine == "rcmax") {
        rowScore <- mean(apply(SimScores, 1, max, na.rm=TRUE))
        colScore <- mean(apply(SimScores, 2, max, na.rm=TRUE))
        result   <- max(rowScore, colScore)
    } else if (combine == "rcmax.avg" || combine == "BMA") {
        result   <- sum( apply(SimScores, 1, max, na.rm=TRUE),
                        apply(SimScores, 2, max, na.rm=TRUE)
                        ) / sum(dim(SimScores))
    }

    return (round(result, digits=3))
}
