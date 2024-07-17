transCCA_sv <- function(W,
                        K = NULL, select.style = "svd", 
                        select.method = "max", var_exp = 2/3, # how to select K
                        null_thres = 0, # null filtering, need to check
                        sumabs = NULL, # set up subabs = 0.4 in PMD
                        verbose = FALSE, ... # ... is for the parameters of PMD function
                        
){
    
    if (select.style == "svd"){
        
        s_value = svd(W)$d
        
    } else if (select.style == "eigen"){
        
        if (nrow(W) >= ncol(W)){
            WW <- t(W) %*% W
        } else {
            WW <- W %*% t(W)
        }
        
        s_value <- eigen(WW)$values
    }
    
    if (is.null(K)){
        
        if (select.method == "max"){
            
            d1 <- -diff(s_value)
            K <- which.max(d1)
            
        } else if (select.method == "var_explain"){
            
            K <- which(cumsum(s_value) / sum(s_value) > var_exp )[1]
            
        }
        
    } 
    
    if (max(abs(W)) < null_thres){
        
        u_variants <- NULL
        v_genes <- NULL
        d <- NULL
        
    } else {
        
        # cross-validation for choosing the best sparsity parameters
        if (is.null(sumabs)){
            min.abs <- max(1/sqrt(nrow(W)), 1/sqrt(ncol(W)))+0.01
            cv.out <- PMD.cv(W, type="standard", center = FALSE, 
                             sumabss=seq(min.abs, 0.8, len=20), trace = FALSE)
            sumabs <- cv.out$bestsumabs
        }
        
        # run PMD based on the best sparsity parameters
        if (!verbose){
            suppressWarnings(
                invisible(capture.output(
                    res <- PMD(x = W, K = K, sumabs = sumabs,
                               type = "standard", center = FALSE, ...),
                    type = c("output", "message")
                ))
            )
        } else {
            res <- PMD(x = W, K = K, sumabs = sumabs, 
                       type = "standard", center = FALSE, ...)
        }
        u_variants <- res$u
        v_genes <- res$v
        d <- res$d
        
    }
    
    results <- list("u" = u_variants,
                    "v" = v_genes,
                    "d" = d,
                    # "W" = W,
                    "s_values" = s_value,
                    "sumabs_best" = sumabs)
    return(results)
    
}

