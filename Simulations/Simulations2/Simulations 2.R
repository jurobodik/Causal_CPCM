library(copula)
library(boot)
library("evd")
library("mgcv")
library("gamlss")
library(CondIndTests)
library("twosamples")
library(tidyverse)
library(kSamples)
library(independence)
library(Pareto)
library(dHSIC)
library(rje)
library(stringr)

#Function CPCM_graph_estimate can be found in file "CPCM_graph_estimate (main function)"
source('CPCM_function.R')
#Implementations of other methods in file "baselines_methods_in_R.R"  
source('baselines_methods_in_R.R')
#Data generations in file "data_generators"  
source('Data generators.R')

#The following code generates the first column in the table with X = sample_AN(n)

set.seed(0)

CPCM = function(X, family_of_distributions='Sequential choice'){
  result = CPCM_graph_estimate(X, family_of_distributions = family_of_distributions)[[1]][6]
  return(2-as.numeric(substr(result,1,1)))  #return 1 if 1-->2, return 0 if 2-->1
}

transmission_for_loci <- function(x){ #return 1 if x>0, return 0 otherwise
  return(ifelse(x>0, 1, 0))
}

transmission_for_heci <- function(x){ #return 1 if x>0, return 0 otherwise
  return(ifelse(x[[1]]==TRUE, 1, 0))
}






n=1000
r1=0;r2=0;r3=0;r4=0;r5=0;r6=0;r7=0;r8=0
for (i in 1:100) {
  
X = sample_AN(n)      #First column
#X = sample_ANs(n)     #Second column
#X = sample_MN_u(n)    #Third column
#X = sample_LS(n)      #Fourth column
#X = sample_LSs(n)     #Fifth column
#X = sample_Poisson(n)  #Sixth column
#X = sample_pareto(n)  #Seventh column
x = X[,1]; y=X[,2]

r1=r1+ CPCM(X, family_of_distributions='Sequential choice') #Possibly CPCM(X, 1) or CPCM(X, 2) if we want better performance
#r2=r2+ ResitWrap(X)$cd  #It seems that package CAM was erased or some other issue between the review rounds
r3=r3+  QCCD(X, m=3)$cd  
r4=r4+  SlopeWrap(X)$cd
r5=r5+  IGCIWrap_G(X)$cd
r6=r6+  IGCIWrap_Unif(X)$cd
r7=r7+  transmission_for_heci(my_module$HECI(as.array(x), as.array(y)))
#r8=r8+  transmission_for_loci(my_module$loci(as.array(x), as.array(y)))#We suggest running this line separately since it takes a long time to run (few hours). It seems that it is much slower in R than in Python

cat("Time remaining", 100- i, "\n")
}

r1/i 
#r2/i 
r3/i
r4/i
r5/i
r6/i
r7/i
r8/i







