soft <- function(x,d){
    return(sign(x)*pmax(0, abs(x)-d))
}

mean_na <- function(vec){
    return(mean(vec[!is.na(vec)]))
}

safesvd <- function(x){
    i <- 1
    out <- try(svd(x), silent=TRUE)
    while(i<10 && inherits(out, "try-error")){
        out <- try(svd(x), silent=TRUE)
        i <- i+1
    }
    if(inherits(out, "try-error")) out <- svd(matrix(rnorm(nrow(x)*ncol(x)), ncol=ncol(x)))
    return(out)
}

l2n <- function(vec){
    a <- sqrt(sum(vec^2))
    if(a==0) a <- .05
    return(a)
}

# return the value such that ||u||_1 = c_u
BinarySearch <- function(argu,sumabs){
    if(l2n(argu)==0 || sum(abs(argu/l2n(argu)))<=sumabs) return(0)
    lam1 <- 0
    lam2 <- max(abs(argu))-1e-5
    iter <- 1
    while(iter < 150){
        su <- soft(argu,(lam1+lam2)/2)
        if(sum(abs(su/l2n(su)))<sumabs){
            lam2 <- (lam1+lam2)/2
        } else {
            lam1 <- (lam1+lam2)/2
        }
        if((lam2-lam1)<1e-6) return((lam1+lam2)/2)
        iter <- iter+1
    }
    warning("Didn't quite converge")
    return((lam1+lam2)/2)
}


CheckPMDV <-  function(v,w,K){
    if(!is.null(v) && is.matrix(v) && ncol(v)>=K){
        v <- matrix(v[,1:K], ncol=K)
    } else if(ncol(w)>nrow(w)){
        v <- matrix(t(w)%*%(safesvd(w%*%t(w))$v[,1:K]),ncol=K)
        if(sum(is.na(v))>0) v <- matrix(safesvd(w)$v[,1:K], ncol=K)
        v <- sweep(v,2,apply(v, 2, l2n), "/")
        if(sum(is.na(v))>0) stop("some are NA")
    } else if (ncol(w)<=nrow(w)){
        v <- matrix(safesvd(t(w)%*%w)$v[,1:K],ncol=K)
    }
    return(v)
}

CheckPMDU <-  function(u,w,K){
    if(!is.null(u) && is.matrix(u) && ncol(u)>=K){
        u <- matrix(u[,1:K], ncol=K)
    } else {
        u <- matrix(safesvd(w%*%t(w))$u[,1:K],ncol=K)
    }
    return(u)
}




get_determinant <- function(V){
    
    R <- length(V) # list of 3 modalities
    value <- c()
    for (i in 1:(R-1)){
        v1 <- V[[i]] # G*K
        for (j in (i+1):R){
            v2 <- V[[j]] # G*K
            r <- cor(v1, v2) # K*K
            r = abs(r)
            tmp <- 2*mean(diag(r)) - mean(r)
            value <- c(value, tmp)
        }
    }
    return(value)
    
}



get_determinant1 <- function(V){
    
    R <- length(V) # list of 3 modalities
    value <- c()
    for (i in 1:(R-1)){
        v1 <- V[[i]] # G*K
        for (j in (i+1):R){
            v2 <- V[[j]] # G*K
            r <- cor(v1, v2) # K*K
            value <- c(value, det(r))
        }
    }
    return(value)
    
}



# - same effect size
# - coupling parameter
# - different effect size, check u vector.
# - sparsity parameters between multiview and singleview model: same sparsity parameter
# - single view sparsity parameters to choose coupling parameter. sparsity parameter first and then coupling parameter.
# - put those parameters in the results.
# - draw gamma - determinant.
# - pmd for single view then run MVE.





