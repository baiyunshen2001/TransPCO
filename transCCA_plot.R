get_pvalue_q2 <- function(q2, null_q2){
    pv=rowMeans(t(null_q2) > q2)
    return(pv)
}

library(ggplot2)
library(dplyr)


plot_function_boxplot=function(q2, null_q2, alpha = 0.005){
    order_col=order(colMeans(null_q2),decreasing = TRUE)
    order_null_q2=null_q2[,order_col]
    order_q2=q2[order_col]
    p_value=get_pvalue_q2(order_q2,order_null_q2)
    colmax=apply(rbind(order_null_q2,order_q2),2,max)
    original_par <- par()
    par(mar = c(7, 4.1, 4.1, 2.1))
    boxplot(order_null_q2,cex=1.5,lwd=1.5,col="grey", pch=19,ylim=c(0,1.1),xlab=paste("Top",ncol(null_q2),"components"),ylab="Sparse canonical correlation values")
    points(x = 1:ncol(order_null_q2), y =order_q2,pch=24 ,cex=1.5,col="red",bg="red")
    points(x = which(p_value<=alpha), y = colmax[which(p_value<=0.05)]+0.05,pch=8 ,cex=1.5,col="black")
    legend(x="topright",pch=17,col="red",legend = "Observed q^2",xpd = TRUE,horiz = TRUE)
    par(mar=original_par$mar)
}


# plot_function_barplot=function(q2, null_q2, alpha = 0.005){
#     
#     order_col=order(colMeans(null_q2),decreasing = TRUE)
#     order_null_q2=null_q2[,order_col]
#     order_q2=q2[order_col]
#     p_value=get_pvalue_q2(order_q2,order_null_q2)
#     colmax=apply(rbind(colMeans(order_null_q2),order_q2),2,max)
#     colerror=apply(order_null_q2,2,sd)
#     bar_mean=barplot(height=colMeans(order_null_q2),border = "grey",
#                      ylim = c(0,1.1),
#                      xlab=paste("Top",ncol(null_q2),"components"),
#                      ylab="Sparse canonical correlation values")
#     points(x = as.vector(bar_mean), y =order_q2,pch=24 ,cex=1.5,col="red",bg="red")
#     points(x = bar_mean[which(p_value<=alpha)], y = colmax[which(p_value<=0.05)]+0.05,pch=8 ,cex=1.5,col="black")
#     arrows(x0 = bar_mean, y0 = colMeans(order_null_q2) - 1.96*colerror, 
#            x1 = bar_mean, y1 = colMeans(order_null_q2) + 1.96*colerror, 
#            angle = 90, code = 3, length = 0.1,lwd=3.0)
#     legend(x="topright",pch=c(17,15),col=c("red","grey"),
#            legend = c("Observed q^2","Null q^2 distribution"))
#     
# }

plot_function_barplot <- function(q2, null_q2, alpha = 0.005) {
    # order_col <- order(colMeans(null_q2), decreasing = TRUE)
    # order_null_q2 <- null_q2[, order_col]
    # order_q2 <- q2[order_col]
    order_null_q2 <- null_q2
    order_q2 <- q2
    p_value <- get_pvalue_q2(order_q2, order_null_q2)  # Assuming this function is defined correctly
    
    # Data preparation for ggplot2 with specific aesthetics for legend
    bar_data <- data.frame(
        Component = factor(1:ncol(order_null_q2), labels = paste0("V", 1:ncol(order_null_q2))),
        Mean = colMeans(order_null_q2),
        SD = apply(order_null_q2, 2, sd),
        Observed = order_q2,
        P_value = p_value,
        Type = rep("Null q2 distribution", ncol(order_null_q2)),  # Constant for fill
        ObsType = rep("Observed q2", ncol(order_null_q2))  # Constant for color
    )
    
    # Adding confidence intervals
    bar_data$Lower <- bar_data$Mean - 1.96 * bar_data$SD
    bar_data$Upper <- bar_data$Mean + 1.96 * bar_data$SD
    
    # Plot using ggplot2
    p <- ggplot(bar_data, aes(x = Component)) +
        geom_bar(aes(y = Mean, fill = Type), stat = "identity") +
        geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
        geom_point(aes(y = Observed, color = ObsType), size = 12, shape = 17) +  # Triangle for observed points
        scale_fill_manual("Legend", values = c("Null q2 distribution" = "grey")) +
        scale_color_manual("Legend", values = c("Observed q2" = "red")) +
        geom_text(data = bar_data[bar_data$P_value <= alpha, ], aes(y = 0, label = "*"), vjust = 1.5, size = 15, color = "blue") +
        labs(x = NULL, y = "Sparse canonical correlation values", title = NULL) +
        theme_minimal() +
        theme(axis.title.x = element_text(size = 20, face = "bold"),
              axis.title.y = element_text(size = 24, face = "bold"),
              axis.text.x = element_text(angle = 0, hjust = 1, size = 20),
              axis.text.y = element_text(size = 20),
              plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
              legend.position = "topright") 
    
    print(p)
}

