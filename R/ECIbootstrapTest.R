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
#' @param data1 A matrix with attributes (e.g. genes) as rows and subjects as columns.
#' @param data2 A matrix with attributes (e.g. genes) as rows and subjects as columns.
#' @param targets1 A data.frame of subject "ID" and "Type" as columns. "Type" has to be "normal" or "tumor".
#' @param targets2 A data.frame of subject "ID" and "Type" as columns. "Type" has to be "normal" or "tumor".
#' @param alpha Significance level for the confidence interval.
#' @param analysisFunc Function to analyze the data. Should return "ES" (e.g. log2FC) and "pval" as columns of a matrix
#' @return The function returns a list with the returned values for the raw data analysis (e.g. differential gene expression analysis) for both studies and a matrix with the results for the ECI bootstrap test with ECI values, the Bca confidence interval, p-pvalue, and q-value of each ECI value. 
#' @export
#' @examples 
#' library(truncnorm)
#' # generate sample data
#' sample1 <- data.frame(ID = seq(1,10),Type = rep(c("normal","tumor"), each = 5))
#' sample2 <- data.frame(ID = seq(1,10),Type = rep(c("normal","tumor"), each = 5))
#' 
#' data1 <- matrix(NA,nrow = 5,ncol = 10)
#' data2 <- matrix(NA,nrow = 5,ncol = 10)
#' rownames(data1) <- seq(1,5)
#' rownames(data2) <- seq(1,5)
#' colnames(data1) <- seq(1,10)
#' colnames(data2) <- seq(1,10)
#' 
#' for(i in 1:5){
#'   data1[i,] = c(rtruncnorm(5,a=0,b=Inf,mean = 3,sd=1),rtruncnorm(5,a=0,b=Inf,mean = 4,sd=1))
#'   data2[i,] = c(rtruncnorm(5,a=0,b=Inf,mean = 3,sd=1),rtruncnorm(5,a=0,b=Inf,mean = 3.75,sd=1))
#' }
#' 
#' # perform ECI bootstrap test
#' ECIbootstrat <- ECIbootstrapTest(data1,data2,sample1,sample2)

ECIbootstrapTest <- function(data1,data2,targets1,targets2, alpha = 0.05, analysisFunc = diffExpr){
  
  # differential gene expression
  ES_list1 <- analysisFunc(data1,targets1)
  ES_list2 <- analysisFunc(data2,targets2)
  
  #get beta values and pvalues from both studies
  beta1 = ES_list1$Es
  beta2 = ES_list2$ES
  pval1 = ES_list1$pval
  pval2 = ES_list2$pval
  names(beta1) = rownames(ES_list1)
  names(beta2) = rownames(ES_list2)
  names(pval1) = rownames(ES_list1)
  names(pval2) = rownames(ES_list2)

  #get eci for all genes
  eci <- ECEA::getECI(beta1,beta2,pval1,pval2)

  n <- 1000
  len <- dim(ES_list1)[1]
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
    new2 = c(sample(group3,length(group3),replace = TRUE),sample(group4,length(group4),replace = TRUE))
    targets1new = targets1[new1,]
    targets2new = targets2[new2,]
    data1new = data1[,new1]
    data2new = data2[,new2]
    
    #diff gene expression
    ES_list_B1 = analysisFunc(data1new,targets1new)
    ES_list_B2 = analysisFunc(data2new,targets2new)
    
    #ECI
    beta_B1 = ES_list_B1$ES
    beta_B2 = ES_list_B2$ES
    pval_B1 = ES_list_B1$pval
    pval_B2 = ES_list_B2$pval
    names(beta_B1) = rownames(ES_list_B1)
    names(beta_B2) = rownames(ES_list_B2)
    names(pval_B1) = rownames(ES_list_B1)
    names(pval_B2) = rownames(ES_list_B2)

    bootstrap[,i] <- ECEA::getECI(beta_B1,beta_B2,pval_B1,pval_B2)
  }

  #confidence intervals, pvalue, and qvalue
  CI <- matrix(NA,ncol = 2, nrow = len)
  pval <- c()
  for(i in 1:len){
    CI[i,] <- coxed::bca(bootstrap[i,], conf.level = 1-alpha)
    pval[i] <- bootstrap.pval(bootstrap[i,])
  }
  names(pval) <- rownames(ES_list1)
  qvalue <- bootstrap.qval(pval)
  qvalue <- qvalue[names(pval)]

  # merging results
  result <- data.frame(ECI = eci, CI.LL = CI[,1], CI.UL = CI[,2], p_value = pval, q_value = qvalue)
  rownames(result) <- rownames(ES_list1)
  
  result_list <- list(EffectSize1 = ES_list1, EffectSize2 = ES_list2, ECIboostrap = result)
  return(result_list)
}
