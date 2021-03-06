% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ECIbootstrapTest.R
\name{ECIbootstrapTest}
\alias{ECIbootstrapTest}
\title{ECI bootstrap test}
\usage{
ECIbootstrapTest(
  data1,
  data2,
  targets1,
  targets2,
  alpha = 0.05,
  analysisFunc = diffExpr
)
}
\arguments{
\item{data1}{A matrix with attributes (e.g. genes) as rows and subjects as columns.}

\item{data2}{A matrix with attributes (e.g. genes) as rows and subjects as columns.}

\item{targets1}{A data.frame of subject "ID" and "Type" as columns. "Type" has to be "normal" or "tumor".}

\item{targets2}{A data.frame of subject "ID" and "Type" as columns. "Type" has to be "normal" or "tumor".}

\item{alpha}{Significance level for the confidence interval.}

\item{analysisFunc}{Function to analyze the data. Should return "ES" (e.g. log2FC) and "pval" as columns of a matrix}
}
\value{
The function returns a list with the returned values for the raw data analysis (e.g. differential gene expression analysis) for both studies and a matrix with the results for the ECI bootstrap test with ECI values, the Bca confidence interval, p-pvalue, and q-value of each ECI value.
}
\description{
This function performs a molecular equivalence test based on the ECI statistic using a bootstrapt test.
}
\examples{
set.seed(45)
# generate sample data
sample1 <- data.frame(ID = seq(1,10),Type = rep(c("normal","tumor"), each = 5))
sample2 <- data.frame(ID = seq(1,10),Type = rep(c("normal","tumor"), each = 5))

data1 <- matrix(NA,nrow = 5,ncol = 10)
data2 <- matrix(NA,nrow = 5,ncol = 10)
rownames(data1) <- seq(1,5)
rownames(data2) <- seq(1,5)
colnames(data1) <- seq(1,10)
colnames(data2) <- seq(1,10)

for(i in 1:5){
  data1[i,] = c(rnorm(5,mean = 3,sd=1),rnorm(5,mean = 4,sd=1))
  data2[i,] = c(rnorm(5,mean = 3,sd=1),rnorm(5,mean = 3.75,sd=1))
}

# perform ECI bootstrap test
ECIbootstrat <- ECIbootstrapTest(data1,data2,sample1,sample2)
}
