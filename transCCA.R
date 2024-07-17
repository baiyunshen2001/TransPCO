library(PMA)
transCCA <- function(data.list, 
                     tol = 0.999,
                     sig.cut = 1e-4, # W calculation
                     K = NULL, select.style = "svd", 
                     select.method = "max", var_exp = 2/3, # how to select K
                     niter = 20, V = NULL,
                     null_thres = 0, # null filtering, need to check
                     sumabs = NULL, # set up subabs = 0.4 in PMD
                     sumabsu = NULL, sumabsv = NULL,
                     upos=FALSE, uneg=FALSE, vpos=FALSE, vneg=FALSE,
                     gamma = NULL,
                     gammas = NULL,
                     coupling = "mean",
                     verbose = FALSE, ...
){
    
    
    if (!is.list(data.list)){
        stop("Input data.list must be the list containing all views!")
    }
    
    W <- W_calculation(data.list, tol = tol, sig.cut = sig.cut)
    R <- length(W)
    
    # - check if the number of SNPs are the same in all views
    pp <- sapply(W, nrow)
    if (length(unique(pp)) != 1){
        stop("Please provide the same number of SNP in each data.list!")
    }
    
    if (length(W) == 1){
        
        W <- W[[1]]
        res <- transCCA_sv(W, K = K, select.style = select.style,
                           select.method = select.method, var_exp = var_exp,
                           null_thres = null_thres, sumabs = sumabs, verbose = verbose, ...)
        sumabs = res$sumabs_best
        sumabsu = sumabs*sqrt(nrow(W))
        sumabsv = sumabs*sqrt(ncol(W))
        ss_tmp <- diag(t(res$u) %*% W %*% res$v)
        q2_tmp <- ss_tmp^2 / sqrt(colSums((t(W)%*%res$u)^2)) / sqrt(colSums((W%*%res$v)^2))

        ll <- list("res" = res,
                   "sumabs" = sumabs,
                   "sumabsu" = sumabsu,
                   "sumabsv" = sumabsv,
                   "ss_score" = ss_tmp,
                   "q2_score" = q2_tmp,
                   "u_id" = rownames(W),
                   "v_id" = colnames(W))
        
    } else {
        
        # - if different genes, we need to impute 0 into the matrix
        gg <- sapply(W, ncol)
        if (length(unique(gg)) != 1){
            g.names <- lapply(W, colnames)
            if (any(sapply(g.names, is.null))){
                stop("Please provide the gene names for your data.list, since you have different number of genes!")
            }
            gnames <- unique(unlist(g.names))
            W <- lapply(W, function(w) {
                w_final <- matrix(0, nrow = nrow(w), ncol = length(gnames))
                pos <- match(colnames(w), gnames)
                w_final[,pos] <- w
                return(w_final)
            })
        } else {
            g.names <- lapply(W, colnames)
            gnames <- unique(unlist(g.names))
        }
        
        # if we need cv for sumabsu and sumabsv, we need write a new function "transCCA_mv_cv" for tuning the best combinations of sumabsu and sumabsv
        # we can also run "PMD.cv" three times, then, for sumabsu = mean(sumabsu from R modalities).
        # cross-validation for choosing the best sparsity parameters
        P <- nrow(W[[1]])
        G <- ncol(W[[1]])
        if (is.null(sumabs)){
            min.abs <- max(1/sqrt(P), 1/sqrt(G))+0.01
            sumabs_tmp <- c()
            for (i in 1:length(W)){
                cv.out <- PMD.cv(W[[i]], type="standard", center = FALSE, 
                                 sumabss=seq(min.abs, 0.8, len=20), trace = FALSE)
                sumabs_tmp[i] <- cv.out$bestsumabs
            }
            sumabsv <- lapply(1:length(W), function(i) sumabs_tmp[i]*sqrt(G))
            sumabsu <- mean(sumabs_tmp)*sqrt(P)
        }
        
        #############################################
        # test: for different sparsity parameters, loss_gamma change a lot.
        # need to fix
        # current: keep the same sparsity, save gamma=0 results, and gamma=gamma_best without choosing from 0.
        #############################################
        
        if (is.null(gamma)){
            
            if (!is.null(gammas)){
                gammas <- gammas
            } else {
                # warning("Do not provide both coupling parameter gamma and gammas. ",
                #                 "Will choose the best gamma from a default sequence gammas.")
                gammas = c(0, 10^seq(-7,-1), 10^seq(0,2))
            }
            
        } else {
            gammas <- gamma
        }
        
        loss_gamma <- c()
        res <- list()
        flag = 1
        for (gamma in gammas){
            tmp <- transCCA_mv(W, K = K, niter = niter, V = V, sumabs = sumabs, sumabsu = sumabsu, sumabsv = sumabsv,
                               upos=upos, uneg=uneg, vpos=vpos, vneg=vneg, gamma = gamma)
            res[[flag]] <- tmp
            flag = flag+1
            V_tmp <- tmp$v
            res_det <- get_determinant(V_tmp)
            loss_gamma <- c(loss_gamma, sum(res_det))
        }
        
        pos <- which.max(loss_gamma)
        diff.max <- (loss_gamma[pos] - loss_gamma[1:(pos-1)])/length(W)
        pos1 <- which(diff.max < 0.01)
        if (length(pos1) == 0){
            gamma_best <- gammas[pos]
            res_gammabest <- res[[pos]]
        } else {
            gamma_best <- gammas[pos1[1]]
            res_gammabest <- res[[pos[1]]]
        }

        
        # pos <- which(loss_gamma>0.9*R)
        # if (length(pos) >= 1){
        #     gamma_best <- gammas[pos[1]]
        #     res_gammabest <- res[[pos[1]]]
        # } else {
        #     pos <- which.max(loss_gamma)
        #     diff.max <- (loss_gamma[pos] - loss_gamma[1:(pos-1)])/R
        #     pos1 <- which(diff.max < 0.01)
        #     if (length(pos1) == 0){
        #         gamma_best <- gammas[pos]
        #         res_gammabest <- res[[pos]]
        #     } else {
        #         gamma_best <- gammas[pos1[1]]
        #         res_gammabest <- res[[pos[1]]]
        #     }
        # }
        
        ss_tmp <- lapply(1:length(W), function(r){
            u <- res_gammabest$u
            v <- res_gammabest$v[[r]]
            diag(t(u) %*% W[[r]] %*% v)
        })
        
        q2_tmp <- lapply(1:length(W), function(r){
            u <- res_gammabest$u
            v <- res_gammabest$v[[r]]
            ww <- W[[r]]
            ss_tmp[[r]]^2 / sqrt(colSums((t(ww)%*%u)^2)) / sqrt(colSums((ww%*%v)^2))
        })
        
        ll <- list("res_gammabest" = res_gammabest,
                   "sumabs" = sumabs,
                   "sumabsu" = sumabsu,
                   "sumabsv" = sumabsv,
                   "gamma" = gamma_best,
                   "loss_gamma" = loss_gamma,
                   "ss_score" = ss_tmp,
                   "q2_score" = q2_tmp,
                   "u_id" = rownames(W[[1]]),
                   "v_id" = gnames)
        
    }
    
    return(ll)
}
