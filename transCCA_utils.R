
#' Calculate Pseudo Inverse Square Root of a Matrix
#'
#' This function computes the pseudo inverse of the square root of a given square matrix.
#' It is based on eigenvalue decomposition, and a tolerance level is used to determine
#' the number of eigenvalues to include.
#'
#' @param Sigma A square matrix for which the pseudo inverse square root is required.
#' @param tol Tolerance level to decide the number of eigenvalues to include.
#' @return The pseudo inverse square root of the matrix.
#' @examples
#' Sigma <- matrix(c(4, 1, 1, 3), nrow = 2)
#' pseudo_inverse_root(Sigma)
#' @export
pseudo_inverse_root <- function(Sigma, tol = 0.999){
    
    if (nrow(Sigma) != ncol(Sigma)) stop("Please provide a square matrix!")
    
    eigen_Sigma <- eigen(Sigma)
    L <- which(cumsum(eigen_Sigma$values) / sum(eigen_Sigma$values) > tol)[1]
    tmp_inverse_root <- eigen_Sigma$vectors[,1:L] %*% 
        diag(1/sqrt(eigen_Sigma$values[1:L])) %*% 
        t(eigen_Sigma$vectors[,1:L])
    return(tmp_inverse_root)
    
}



get_pvalue <- function(null_scores, res_transCCA, K = NULL, method = "q2"){
    
    if (is.list(res_transCCA$res$v)){
        
        R <- length(res_transCCA$res$v)
        pv <- list()
        for (r in 1:R){
            q2 <- res_transCCA$q2_score[[r]]
            null <- null_scores$q2_score[[r]][,1:K]
            pv[[r]] <- rowMeans(t(null) > q2)
        }
        
    } else {
        
        q2 <- res_transCCA$q2_score
        null <- null_scores$q2_score[[1]][,1:K]
        pv <- rowMeans(t(null) > q2)
        pv <- list(pv)
        
    }
    return(pv)
    
}


get_null_scores <- function(num, K, sumabs,
                            method = "diag",
                            R, N, P, 
                            MAF,
                            Sigma_EE_R){
    
    ss_score <- vector(mode='list', length=R)
    q2_score <- vector(mode='list', length=R)
    for (i in 1:num){
        null.data <- Generate_NullData_multi(R = R, N = N, P = P, 
                                             MAF = MAF,
                                             Sigma_EE_R = Sigma_EE_R)
        
        data.list <- lapply(1:R, function(r){ 
            datals <- null.data[[r]]
            list("G" = datals$G_scaled, "E" = datals$E_trans_scaled) 
        })
        if (method != "diag"){
            data.list <- lapply(data.list, function(dl){ 
                list("G" = dl$G, "E" = dl$E, "Sigma_EE" = cova(dl$E)) 
            })
        }
        res_transCCA <- transCCA(data.list, K = K, sumabs = sumabs)
        if (R == 1){
            ss_score[[1]] <- rbind(ss_score[[1]], res_transCCA$ss_score)
            q2_score[[1]] <- rbind(q2_score[[1]], res_transCCA$q2_score)
        } else {
            for (r in 1:R){
                ss_score[[r]] <- rbind(ss_score[[r]], res_transCCA$ss_score[[r]])
                q2_score[[r]] <- rbind(q2_score[[r]], res_transCCA$q2_score[[r]])
            }
        }
        
        if ( i%%1==0 ){
            message(paste("Finish", i, "out for", num))
        }
        
    }
    scores <- list("ss_score" = ss_score,
                   "q2_score" = q2_score)
    return(scores)
}


sample_MAF_density <- function(d, cdf) {
    u <- runif(1)
    sampled_indices <- which.min(abs(cdf - u))
    sampled_values <- d$x[sampled_indices]
    if (sampled_values<=0.01 ||sampled_values>=0.5)
        return(sample_MAF_density(d,cdf))
    return(sampled_values)
}

# --- null simulation
Generate_NullData_one <- function(n, # sample size
                                  p, # number of SNPs
                                  MAF,
                                  Sigma_EE){
    # - initial generated data
    data <- list()
    
    ########### - step 1. generate independent 50% SNPs with MAF in [0.1, 0.4] and 50% SNPs with MAF in [0.01, 0.1]
    MAF=as.numeric(MAF[!is.na(MAF)])
    d <- density(MAF)
    cdf <- cumsum(d$y) / sum(d$y)
    G <- lapply(1:p, function(i) rbinom(n, 2,sample_MAF_density(d,cdf)))
    G <- do.call(cbind, G)
    data$G_scaled <- apply(G, 2, scale) # standardize 
    
    # step 2. generate independent trans-gene expressions 
    E_trans <- MASS::mvrnorm(n=n, mu=rep(0,ncol(Sigma_EE)), Sigma = Sigma_EE)
    data$E_trans_scaled <- E_trans
    
    return(data)
}


Generate_NullData_multi <- function(R, # number of views
                                    N, # sample size
                                    P, # number of SNPs
                                    MAF,
                                    Sigma_EE_R){
    
    data.list <- list()
    
    for (i.r in 1:R){
        
        # - check if same N
        if (length(N) == 1){
            n <- N
        } else if (length(N) == R){
            n <- N[i.r]
        } else {
            stop("Please provide only one sample size if all views consistent, otherwise, provide R sample size!")
        }
        # - check if same P
        if (length(P) == 1){
            p <- P
        } else if (length(P) == R){
            p <- P[i.r]
        } else {
            stop("Please provide only one P if all views consistent, otherwise, provide R Ps!")
        }
        # - check if same Sigma_EE
        if (is.matrix(Sigma_EE_R)){
            Sigma_EE <- Sigma_EE_R
        } else {
            if (length(Sigma_EE_R) == 1){
                Sigma_EE <- Sigma_EE_R
            } else if (length(Sigma_EE_R) == R){
                Sigma_EE <- Sigma_EE_R[[i.r]]
            } else {
                stop("Please provide only one Sigma_EE_R if all views consistent, otherwise, provide R Sigma_EE_R!")
            }
        }
        
        # single view data generation
        data.list[[i.r]] <- Generate_NullData_one(n=n, p=p, MAF = MAF, Sigma_EE = Sigma_EE)
        
    }
    
    return(data.list)
}




