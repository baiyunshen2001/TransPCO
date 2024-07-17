transCCA_mv <- function(W,
                        K = NULL, select.style = "svd", 
                        niter = 20, V = NULL,
                        select.method = "max", var_exp = 2/3, # how to select K
                        null_thres = 0, # null filtering, need to check
                        sumabs = NULL, # set up subabs = 0.4 in PMD
                        sumabsu = NULL, sumabsv = NULL,
                        upos=FALSE, uneg=FALSE, vpos=FALSE, vneg=FALSE,
                        gamma = 0,
                        verbose = FALSE
){
    
    R <- length(W)
    # ------- upos, uneg, vpos, vneg
    # upos: constrain the elements of u to be positive.
    # uneg: constrain the elements of u to be negative.
    # vpos: constrain the elements of v to be positive.
    # vneg: constrain the elements of v to be negative.
    if(upos&&uneg || vpos&&vneg) stop("Cannot contrain elements to be both positive and negative!")
    
    # -------- sumabs, sumabsu, sumabsv
    # sumabs: a measure of sparsity for u and v vectors, between 0 and 1. When sumbas is specified, and sumabsu and sumabsv are NULL, 
    # then sumabsu = sqrt(p)*sumabs abd sumabsv = sqrt(g)*sumabs. 
    # If sumabs is specified, then sumabsu and sumabsv should be NULL. Or if sumabsu and sumabsv are specified, then sumabs should be NULL.
    if (is.null(sumabsu)){
        if (!is.null(sumabs)){
            sumabsu <- sqrt(nrow(W[[1]]))*sumabs
        } else {
            stop("Provide sumabs or sumabsu for the sparsity of u!")
        }
    } else {
        if (!is.null(sumabs)){
            warning("Both sumabs and sumabsu are provided, sumabs was ignored.")
        }
        if ( (sumabsu<1) | (sumabsu>sqrt(nrow(W[[1]]))) ){ stop("sumabsu must be between 1 and sqrt(p).")}
    }
    if (is.null(sumabsv)){
        if (!is.null(sumabs)){
            sumabsv <- sqrt(ncol(W[[1]]))*sumabs
            sumabsV <- rep(sumabsv, 3)
        } else {
            stop("Provide sumabs or sumabsu for the sparsity of u!")
        }
    } else {
        if (!is.null(sumabs)){
            warning("Both sumabs and sumabsu are provided, sumabs was ignored.")
        }
        if (length(sumabsv) == 1){
            warning("You only provide one sumabs, that is, all modalities will have the same sparsity for V.")
        } else {
            if (length(sumabsv) != R){
                stop("If you provide multiple sumabsv, you should provide the same number of modalities.")
            }
        }
        sumabsV <- c()
        for (i in 1:R){
            if ( (sumabsv[i]<1) | (sumabsv[i]>sqrt(ncol(W[[1]]))) ){ 
                stop(paste0("sumabsv must be between 1 and sqrt(g) for modality ", i))
            }
            sumabsV[i] <- sumabsv[i]
        }
        
    }
    
    # ------- V
    # V: the first right singular vectors of data. V is used as the initial value for the iterative PMD algorithm.
    # If W is large, then this step can be time-consuming; therefore, if PMD is to be run multiple times, then V should be computed once and saved.
    if (!is.null(V)){
        if (!is.list(V)){
            stop("V should be a list of first right singular vectors of the list of W with the same length.")
        }
        else {
            if (length(W)!=length(V)){stop("V should have the same length of W.")}
        }
    } else {
        V <- vector(mode = "list", length = length(W))
    }
    V <- lapply(1:R, function(i){
        w <- W[[i]]
        v <- V[[i]]
        return(CheckPMDV(v, w, K))
    })
    
    P <- nrow(W[[1]])
    G <- ncol(W[[1]])
    Ds <- vector(mode = "list", length = R)
    Us <- matrix(0, nrow = P, ncol = K)
    Vs <- vector(mode = "list", length = R)
    W_use <- W
    keep_update <- 1:R
    MVE <- list()
    
    for(k in 1:K){
        
        old.V <- lapply(keep_update, function(i) rnorm(G))
        V_k <- lapply(V[keep_update], function(v) v[,k,drop=FALSE])
        diff.V <- mean(sapply(keep_update, function(i) sum(abs(old.V[[i]]-V_k[[i]]))))
        check <- diff.V
        
        for (iter in 1:niter){
            
            if (diff.V > 1e-7){
                old.V <- V_k
                # - updating u
                argu <- lapply(keep_update, function(i){
                    w <- W_use[[i]]
                    v <- V_k[[i]]
                    return(w%*%v)
                })
                argu <- do.call(cbind, argu)
                argu <- rowmeans(argu)
                if(upos) argu <- pmax(argu,0)
                if(uneg) argu <- pmin(argu,0)
                lamu <- BinarySearch(argu,sumabsu)
                su <- soft(argu, lamu)
                u <- matrix(su/l2n(su), ncol=1)
                # - done updating u
                
                # - updating v
                for (i in keep_update){
                    pos <- setdiff(keep_update, i)
                    tmp <- lapply(pos, function(it){ gamma*V_k[[it]] })
                    tmp <- do.call(cbind, tmp)
                    coupling_v <- rowSums(tmp) # - take sum and smaller gammas
                    # coupling_v <- rowMeans(tmp)
                    # coupling_v <- 0
                    # for (ii in keep_update){
                    #     v1 <- V_k[[ii]]
                    #     for (jj in setdiff(keep_update, ii)){
                    #         v2 <- V_k[[jj]]
                    #         if (coupling == "mean"){
                    #             coupling_v <- coupling_v + gamma*mean(v1*v2) # - take mean and keep the same gammas
                    #         } else if (coupling == "sum"){
                    #             coupling_v <- coupling_v + gamma*sum(v1*v2) # - take sum and smaller gammas
                    #         }
                    #     }
                    # } 
                    argv <- t(u) %*% W_use[[i]] + coupling_v
                    if(vpos) argv <- pmax(argv,0)
                    if(vneg) argv <- pmin(argv,0)
                    lamv <- BinarySearch(argv, sumabsV[i])
                    sv <- soft(argv, lamv)
                    v <- matrix(sv/l2n(sv), ncol=1)
                    V_k[[i]] <- v
                }
                # - done updating V
                diff.V <- mean(sapply(keep_update, function(i) sum(abs(old.V[[i]]-V[[i]][,k]))))
                check <- c(check, diff.V)
            }
            
        }
        
        
        for (i in keep_update){
            Ds[[i]][k] <- as.numeric(t(u) %*% (W_use[[i]] %*% V_k[[i]]))
            Us[,k] <- u
            Vs[[i]] <- cbind(Vs[[i]], V_k[[i]])
            W_use[[i]] <- W_use[[i]] - Ds[[i]][k]*u%*%t(V_k[[i]])
        }
        
        # # - check the total explained of W by the first k componenets
        # # - if W has been explained after 80%, we will stop 
        # # Based on RRMSE in flashr paper
        # U_k <- Us[,c(1:k),drop=FALSE]
        # W_k_hat <- lapply(keep_update, function(ii){
        #     D_k <- matrix(0, nrow = k, ncol = k)
        #     diag(D_k) <- Ds[[ii]][1:k]
        #     return(U_k %*% D_k %*% t(Vs[[ii]]))
        # })
        # # - we need to consider the sparsity.
        # MVE[[k]] <- lapply(keep_update, function(ii){
        #     sum((W_k_hat[[ii]])^2) / sum(W[[ii]]^2) 
        # })
        
    }
    
    results <- list("u" = Us,
                    "v" = Vs,
                    "d" = Ds)
    return(results)
    
}





