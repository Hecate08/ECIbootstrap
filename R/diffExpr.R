#' Differential Expression using limma
#'
#' This function performs differential expression using Limma and
#' calculates next to the pvalue the standard deviation.
#'
#' @param data A gene expression matrix with genes as rows and subjects as columns.
#' @param targets A data.frame of subject "ID" and "Type" as columns. "Type" has to be "normal" or "tumor".
#' @return A matrix of differential gene expression values in log2FC format, standard deviation,
#' p value, and degrees of freedom.
#' @export
diffExpr <- function(data, targets){
  # getting pval
  design <- model.matrix(~targets$Type)
  colnames(design) <- levels(targets$Type)

  fit <- limma::lmFit(data, design)
  fit <- limma::eBayes(fit)
  #gene_list1 <- topTable(fit, coef=2, number=1000000, sort.by="none",confint = TRUE)

  # getting sd
  beta <- fit$coefficients[,2]
  pval <- fit$p.value[,2]
  sd <- ((sqrt(fit$s2.post)) * (fit$stdev.unscaled))[,2]
  gene_list <- data.frame(ES = beta, pval = pval, sd = sd)

  return(gene_list)
}
