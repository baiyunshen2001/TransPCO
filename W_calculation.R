library(Rfast)
W_calculation <- function(data.list,
                          tol = 0.999,
                          sig.cut = 1e-4){
    
    if (!is.list(data.list)){
        stop("Input data.list must be the list containing data for all views!")
    }
    
    R <- length(data.list)
    W <- vector(mode = "list", length = R)
    for (i in 1:R){
        
        data.tmp <- data.list[[i]]
        Sigma_GE = data.tmp$Sigma_GE
        Z = data.tmp$Z
        N = data.tmp$N
        MAF = data.tmp$MAF 
        sdY = data.tmp$sdY
        Beta = data.tmp$Beta
        seBeta = data.tmp$seBeta
        Pvalue = data.tmp$Pvalue 
        G = data.tmp$G
        E = data.tmp$E
        log2FC = data.tmp$log2FC
        N1 = data.tmp$N1
        N0 = data.tmp$N0
        Sigma_GG = data.tmp$Sigma_GG
        Sigma_EE = data.tmp$Sigma_EE
        
        W[[i]] <- W_calculation_single(Sigma_GE = Sigma_GE, Z = Z, N = N, MAF = MAF, sdY = sdY,
                                       Beta = Beta, seBeta = seBeta, Pvalue = Pvalue, G = G, E = E,
                                       log2FC = log2FC, N1 = N1, N0 = N0,
                                       Sigma_GG = Sigma_GG, Sigma_EE = Sigma_EE,
                                       tol = tol, sig.cut = sig.cut)
        
    }
    
    return(W)
}


W_calculation_single <- function(Sigma_GE = NULL, 
                                 Z = NULL, N = NULL, 
                                 MAF = NULL, sdY = NULL,
                                 Beta = NULL, seBeta = NULL, Pvalue = NULL, 
                                 G = NULL, E = NULL,
                                 log2FC = NULL, N1 = NULL, N0 = NULL,
                                 Sigma_GG = NULL,
                                 Sigma_EE = NULL,
                                 tol = 0.999,
                                 sig.cut = 1e-4){
    
    if (!is.null(Sigma_GE) && 
        is.null(G) && is.null(E) &&
        is.null(Z) && is.null(N) && 
        is.null(MAF) && is.null(sdY) &&
        is.null(Beta) && is.null(seBeta) && is.null(Pvalue) &&
        is.null(log2FC) && is.null(N1) && is.null(N0)){ 
        
        Sigma_GE <- Sigma_GE 
        
        if (is.null(Sigma_GG) && is.null(Sigma_EE)){
            
            W <- Sigma_GE
            
        } else if (is.null(Sigma_GG) && !is.null(Sigma_EE)) {
            
            Sigma_EE_inv_root <- pseudo_inverse_root(Sigma_EE, tol = tol)
            W <- Sigma_GE %*% Sigma_EE_inv_root
            
        } else if (is.null(Sigma_EE) && !is.null(Sigma_GG)){
            
            Sigma_GG_inv_root <- pseudo_inverse_root(Sigma_GG, tol = tol)
            W <- Sigma_GG_inv_root %*% Sigma_GE
            
        } else {
            
            Sigma_GG_inv_root <- pseudo_inverse_root(Sigma_GG, tol = tol)
            Sigma_EE_inv_root <- pseudo_inverse_root(Sigma_EE, tol = tol)
            W <- Sigma_GG_inv_root %*% Sigma_GE %*% Sigma_EE_inv_root
            
        }
        
    } else if (!is.null(G) && !is.null(E) &&
               is.null(Z) && is.null(N) && 
               is.null(Sigma_GE) && 
               is.null(MAF) && is.null(sdY) &&
               is.null(Beta) && is.null(seBeta) && is.null(Pvalue) &&
               is.null(log2FC) && is.null(N1) && is.null(N0)){
        
        n <- nrow(G)
        
        if (is.null(Sigma_GG) && is.null(Sigma_EE)){
            
            W <- t(G) %*% E / n
            
        } else if (is.null(Sigma_GG) && !is.null(Sigma_EE)) {
            
            if ( all(Sigma_EE == "insample") ){ Sigma_EE <- cova(E) } 
            Sigma_EE_inv_root <- pseudo_inverse_root(Sigma_EE, tol = tol)
            E_tilde <- E %*% Sigma_EE_inv_root
            W <- t(G) %*% E_tilde / n
            
        } else if (is.null(Sigma_EE) && !is.null(Sigma_GG)){
            
            if ( all(Sigma_GG == "insample") ){ Sigma_GG <- cova(G) } 
            Sigma_GG_inv_root <- pseudo_inverse_root(Sigma_GG, tol = tol)
            G_tilde <- G %*% Sigma_GG_inv_root
            W <- t(G_tilde) %*% E / n
            
        } else {
            
            Sigma_EE <- cova(E)
            Sigma_EE_inv_root <- pseudo_inverse_root(Sigma_EE, tol = tol)
            E_tilde <- E %*% Sigma_EE_inv_root
            Sigma_GG <- cova(G)
            Sigma_GG_inv_root <- pseudo_inverse_root(Sigma_GG, tol = tol)
            G_tilde <- G %*% Sigma_GG_inv_root
            W <- t(G_tilde) %*% E_tilde / nrow(G)
            
        }
        
    } else if (!is.null(Z) && !is.null(N) && 
               is.null(Sigma_GE) && 
               is.null(MAF) && is.null(sdY) &&
               is.null(Beta) && is.null(seBeta) && is.null(Pvalue) &&
               is.null(G) && is.null(E) &&
               is.null(log2FC) && is.null(N1) && is.null(N0)){
        
        Sigma_GE <- Z *(1 / sqrt((N - 2 + Z^2)))
        
        if (is.null(Sigma_GG) && is.null(Sigma_EE)){
            
            W <- Sigma_GE
            
        } else if (is.null(Sigma_GG) && !is.null(Sigma_EE)) {
            
            Sigma_EE_inv_root <- pseudo_inverse_root(Sigma_EE, tol = tol)
            W <- Sigma_GE %*% Sigma_EE_inv_root
            
        } else if (is.null(Sigma_EE) && !is.null(Sigma_GG)){
            
            Sigma_GG_inv_root <- pseudo_inverse_root(Sigma_GG, tol = tol)
            W <- Sigma_GG_inv_root %*% Sigma_GE
            
        } else {
            
            Sigma_GG_inv_root <- pseudo_inverse_root(Sigma_GG, tol = tol)
            Sigma_EE_inv_root <- pseudo_inverse_root(Sigma_EE, tol = tol)
            W <- Sigma_GG_inv_root %*% Sigma_GE %*% Sigma_EE_inv_root
            
        }
        
    } else if (!is.null(Z) && !is.null(N) && 
               is.null(Sigma_GE) && 
               !is.null(MAF) && is.null(sdY) &&
               is.null(Beta) && is.null(seBeta) && is.null(Pvalue) &&
               is.null(G) && is.null(E) &&
               is.null(log2FC) && is.null(N1) && is.null(N0)){
        
        D_G <- diag(sqrt(2 * MAF * (1-MAF)))
        if (!is.null(sdY)) {D_E_inv <- diag(1/sdY)
        } else {D_E_inv <- diag(ncol(Z))}
        one_p_g <- matrix(1, nrow = nrow(Z), ncol = ncol(Z))
        Sigma_GE <- (D_G %*% one_p_g %*% D_E_inv) * (Z *(1 / sqrt((N - 2 + Z^2))))
        
        if (is.null(Sigma_GG) && is.null(Sigma_EE)){
            
            W <- Sigma_GE
            
        } else if (is.null(Sigma_GG) && !is.null(Sigma_EE)) {
            
            Sigma_EE_inv_root <- pseudo_inverse_root(Sigma_EE, tol = tol)
            W <- Sigma_GE %*% Sigma_EE_inv_root
            
        } else if (is.null(Sigma_EE) && !is.null(Sigma_GG)){
            
            Sigma_GG_inv_root <- pseudo_inverse_root(Sigma_GG, tol = tol)
            W <- Sigma_GG_inv_root %*% Sigma_GE
            
        } else {
            
            Sigma_GG_inv_root <- pseudo_inverse_root(Sigma_GG, tol = tol)
            Sigma_EE_inv_root <- pseudo_inverse_root(Sigma_EE, tol = tol)
            W <- Sigma_GG_inv_root %*% Sigma_GE %*% Sigma_EE_inv_root
            
        }
        
    } else if (!is.null(Beta) && 
               is.null(N) && is.null(MAF) && is.null(sdY) &&
               is.null(Z) && is.null(seBeta) && is.null(Pvalue) &&
               is.null(Sigma_GE) && 
               is.null(G) && is.null(E) &&
               is.null(log2FC) && is.null(N1) && is.null(N0)){
        
        Sigma_GE <- Beta 
        
        if (is.null(Sigma_GG) && is.null(Sigma_EE)){
            
            W <- Sigma_GE
            
        } else if (is.null(Sigma_GG) && !is.null(Sigma_EE)) {
            
            Sigma_EE_inv_root <- pseudo_inverse_root(Sigma_EE, tol = tol)
            W <- Sigma_GE %*% Sigma_EE_inv_root
            
        } else if (is.null(Sigma_EE) && !is.null(Sigma_GG)){
            
            Sigma_GG_inv_root <- pseudo_inverse_root(Sigma_GG, tol = tol)
            W <- Sigma_GG_inv_root %*% Sigma_GE
            
        } else {
            
            Sigma_GG_inv_root <- pseudo_inverse_root(Sigma_GG, tol = tol)
            Sigma_EE_inv_root <- pseudo_inverse_root(Sigma_EE, tol = tol)
            W <- Sigma_GG_inv_root %*% Sigma_GE %*% Sigma_EE_inv_root
            
        }
        
    } else if (!is.null(Beta) && !is.null(MAF) &&
               is.null(N) && 
               is.null(Z) && is.null(seBeta) && is.null(Pvalue) &&
               is.null(Sigma_GE) && 
               is.null(G) && is.null(E) &&
               is.null(log2FC) && is.null(N1) && is.null(N0)){
        
        D_G <- diag(sqrt(2 * MAF * (1-MAF)))
        if (!is.null(sdY)) {D_E_inv <- diag(1/sdY)
        } else {D_E_inv <- diag(ncol(Beta))}
        one_p_g <- matrix(1, nrow = nrow(Beta), ncol = ncol(Beta))
        Sigma_GE <- (D_G %*% one_p_g %*% D_E_inv) * Beta 
        
        if (is.null(Sigma_GG) && is.null(Sigma_EE)){
            
            W <- Sigma_GE
            
        } else if (is.null(Sigma_GG) && !is.null(Sigma_EE)) {
            
            Sigma_EE_inv_root <- pseudo_inverse_root(Sigma_EE, tol = tol)
            W <- Sigma_GE %*% Sigma_EE_inv_root
            
        } else if (is.null(Sigma_EE) && !is.null(Sigma_GG)){
            
            Sigma_GG_inv_root <- pseudo_inverse_root(Sigma_GG, tol = tol)
            W <- Sigma_GG_inv_root %*% Sigma_GE
            
        } else {
            
            Sigma_GG_inv_root <- pseudo_inverse_root(Sigma_GG, tol = tol)
            Sigma_EE_inv_root <- pseudo_inverse_root(Sigma_EE, tol = tol)
            W <- Sigma_GG_inv_root %*% Sigma_GE %*% Sigma_EE_inv_root
            
        }
        
    } else if (!is.null(Beta) && !is.null(N) && 
               (!is.null(Z) | !is.null(Pvalue) | !is.null(seBeta)) && 
               is.null(MAF) &&
               is.null(Sigma_GE) && 
               is.null(G) && is.null(E) &&
               is.null(log2FC) && is.null(N1) && is.null(N0)){
        
        if (!is.null(Pvalue)){
            Z <- sign(Beta) * sqrt(qchisq(Pvalue, 1, lower.tail = FALSE))
            seBeta <- Beta / Z
        } else if (!is.null(Z)){
            seBeta <- Beta / Z
            Pvalue <- pchisq(Z^2, 1, lower.tail = FALSE)
        } else if(!is.null(seBeta)){
            Z <- Beta / seBeta
            Pvalue <- pchisq(Z^2, 1, lower.tail = FALSE)
        }
        if (!is.null(sdY)) {
            D_E_inv <- diag(1/sdY)
            S_E <- t(replicate(nrow(Beta), sdY^2))
        } else {
            D_E_inv <- diag(ncol(Beta))
            S_E <- matrix(1, nrow = nrow(Beta), ncol = ncol(Beta))
        }
        v_tmp <- sqrt(S_E / (Beta^2 + (N-2)*seBeta^2))
        
        # add-on if sigma^2 NAN, we set as 0.
        v_tmp[is.nan(v_tmp)] <- 0
        
        sig.num <- colSums(Pvalue < sig.cut)
        pos <- which(sig.num != 0)
        if (length(pos) != 0){
            v_tmp <- v_tmp[, pos, drop = FALSE]
            D_G <- diag(rowMeans(v_tmp))
        } else {
            g.min <- apply(Pvalue, 2, min)
            n.min <- ceiling(0.2 * length(g.min))
            pos <- rank(g.min)[1:n.min]
            v_tmp <- v_tmp[, pos, drop = FALSE]
            D_G <- diag(rowMeans(v_tmp))
        }
        one_p_g <- matrix(1, nrow = nrow(Beta), ncol = ncol(Beta))
        Sigma_GE <- (D_G %*% one_p_g %*% D_E_inv) * Beta 
        
        if (is.null(Sigma_GG) && is.null(Sigma_EE)){
            
            W <- Sigma_GE
            
        } else if (is.null(Sigma_GG) && !is.null(Sigma_EE)) {
            
            Sigma_EE_inv_root <- pseudo_inverse_root(Sigma_EE, tol = tol)
            W <- Sigma_GE %*% Sigma_EE_inv_root
            
        } else if (is.null(Sigma_EE) && !is.null(Sigma_GG)){
            
            Sigma_GG_inv_root <- pseudo_inverse_root(Sigma_GG, tol = tol)
            W <- Sigma_GG_inv_root %*% Sigma_GE
            
        } else {
            
            Sigma_GG_inv_root <- pseudo_inverse_root(Sigma_GG, tol = tol)
            Sigma_EE_inv_root <- pseudo_inverse_root(Sigma_EE, tol = tol)
            W <- Sigma_GG_inv_root %*% Sigma_GE %*% Sigma_EE_inv_root
            
        }
        
    } else if (!is.null(log2FC) && !is.null(N1) && !is.null(N0) &&
               is.null(Sigma_GE) && 
               is.null(G) && is.null(E) &&
               is.null(Z) && is.null(N) && 
               is.null(MAF) && is.null(sdY) &&
               is.null(Beta) && is.null(seBeta) && is.null(Pvalue)){
        
        D_G <- diag(sqrt(N1*N0) / (N1+N0))
        Sigma_GE <- D_G %*% log2FC 
        
        if (is.null(Sigma_GG) && is.null(Sigma_EE)){
            
            W <- Sigma_GE
            
        } else if (is.null(Sigma_GG) && !is.null(Sigma_EE)) {
            
            Sigma_EE_inv_root <- pseudo_inverse_root(Sigma_EE, tol = tol)
            W <- Sigma_GE %*% Sigma_EE_inv_root
            
        } else if (is.null(Sigma_EE) && !is.null(Sigma_GG)){
            
            Sigma_GG_inv_root <- pseudo_inverse_root(Sigma_GG, tol = tol)
            W <- Sigma_GG_inv_root %*% Sigma_GE
            
        } else {
            
            Sigma_GG_inv_root <- pseudo_inverse_root(Sigma_GG, tol = tol)
            Sigma_EE_inv_root <- pseudo_inverse_root(Sigma_EE, tol = tol)
            W <- Sigma_GG_inv_root %*% Sigma_GE %*% Sigma_EE_inv_root
            
        }
    }
    
    return(W)
}

