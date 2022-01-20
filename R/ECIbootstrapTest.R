###########
# ECI

#' @export
bootstrap.pval <- function(data){
  B = 1000
  pE <- mean(data)
  N = sum(data < pE)
  M = sum(data < 0 )
  alpha1 = pnorm(-2*qnorm(N/B) + qnorm(M/B))
  if(pE < 0) alpha1 = 1 - alpha1
  alpha = 2*alpha1
  return(alpha)
}

#' @export
bootstrap.qval <- function(pval){
  sortP <- sort(pval)
  q <- c()
  N <- length(pval)
  for(i in 1:N){
    q[i] <- sortP[i]*N/i
  }

  for(i in 1:N){
    q[i] <- min(q[(i):N])
  }
  names(q) <- names(sortP)
  return(q)
}

#' ECI bootstrap test
#'
#' This function performs a molecular equivalence test based on the ECI statistic using a bootstrapt test.
#'
#' @param gene_list1 A matrix of differential gene expression values for the first study
#' @param gene_list2 A matrix of differential gene expression values for the second study
#' @param geneExpr1 A gene expression matrix with genes as rows and subjects as columns.
#' @param geneExpr2 A gene expression matrix with genes as rows and subjects as columns.
#' @param targets1 A data.frame of subject "ID" and "Type" as columns. "Type" has to be "normal" or "tumor".
#' @param targets2 A data.frame of subject "ID" and "Type" as columns. "Type" has to be "normal" or "tumor".
#' @param filter A boolean parameter indicating if the returned values should by filtered by genes expressing
#' being differentially expressed in at least one of the studies.
#' @return A matrix of ECI values, confidence interval, p value, and q value.
#' @export
ECIbootstrapTest <- function(gene_list1,gene_list2, geneExpr1,geneExpr2,targets1,targets2, filter = TRUE){
  #get beta values and pvalues from both studies
  beta1 = gene_list1$log2FC
  names(beta1) = rownames(gene_list1)
  pval1 = gene_list1$pval
  names(pval1) = rownames(gene_list1)
  beta2 = gene_list2$log2FC
  names(beta2) = rownames(gene_list2)
  pval2 = gene_list2$pval
  names(pval2) = rownames(gene_list2)

  #get eci for all genes
  eci <- ECEA::getECI(beta1,beta2,pval1,pval2)

  n <- 1000
  len <- dim(gene_list1)[1]
  #id <- seq(1, len)

  group1 = which(targets1[,2] == "tumor")
  group2 = which(targets1[,2] == "normal")
  group3 = which(targets2[,2] == "tumor")
  group4 = which(targets2[,2] == "normal")

  bootstrap <- matrix(NA, ncol = n, nrow = len)
  for(i in 1:n){
    #if (i %% 10 == 0){
    #  print(i)
    #}
    #Bootstrap sampling
    new1 = c(sample(group1,length(group1),replace = TRUE),sample(group2,length(group2),replace = TRUE))
    targets1new = targets1[new1,]
    geneExpr1new = geneExpr1[,new1]
    new2 = c(sample(group3,length(group3),replace = TRUE),sample(group4,length(group4),replace = TRUE))
    targets2new = targets2[new2,]
    geneExpr2new = geneExpr2[,new2]
    #diff gene expression
    gene_list_B1 = diffExpr(geneExpr1new,targets1new)
    gene_list_B2 = diffExpr(geneExpr2new,targets2new)
    #ECI
    beta_B1 = gene_list_B1$log2FC
    names(beta_B1) = rownames(gene_list_B1)
    pval_B1 = gene_list_B1$pval
    names(pval_B1) = rownames(gene_list_B1)
    beta_B2 = gene_list_B2$log2FC
    names(beta_B2) = rownames(gene_list_B2)
    pval_B2 = gene_list_B2$pval
    names(pval_B2) = rownames(gene_list_B2)

    #get eci for all genes
    bootstrap[,i] <- ECEA::getECI(beta_B1,beta_B2,pval_B1,pval_B2)
  }

  #confidence intervals
  CI <- matrix(NA,ncol = 2, nrow = len)
  pval <- c()
  for(i in 1:len){
    CI[i,] <- coxed::bca(bootstrap[i,], conf.level = 0.95)
    pval[i] <- bootstrap.pval(bootstrap[i,])
  }

  result <- data.frame(ECI = eci, CI.LL = CI[,1], CI.UL = CI[,2], p_value = pval)
  rownames(result) <- rownames(gene_list1)

  # output only those ECI where at least one of the genes has abs(log2FC) > 1 and pval < 0.05
  if(filter){
    result <- result[(abs(gene_list1$log2FC) > 1 & gene_list1$pval < 0.05) | (abs(gene_list2$log2FC) > 1 & gene_list2$pval < 0.05),]
  }
  p <- result$p_value
  names(p) <- rownames(result)
  qvalue <- bootstrap.qval(p)
  qvalue <- qvalue[names(p)]
  result$q_value <- qvalue
  return(result)
}
