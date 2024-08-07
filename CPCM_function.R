#This is the code for the main function that estimates the causal graph via CPCM, called CPCM_graph_estimate(X, family_of_distributions). 

#family_of_distributions correspond to the models we use. We implemented a joint CPCM(F1...Fk) model with the following:
#family1 = c('Gaussian with fixed sigma', 'Pareto',  'Poisson')
#family2 = c('Gaussian', 'Gamma', 'Negative_binomial')

#family_of_distributions = 1 if we use family1
#family_of_distributions = 2 if we use family2

#Rule of thumb is that we use family1 if n<=1000 and we use family2 if n>1000, but it should depend on the complexity of the dataset

#If you want to use CPCM(F) model, the choices for 'family_of_distributions' 
#are the follwing: "Gaussian", "Gaussian with fixed sigma", "Pareto", "Gamma", "Gamma with fixed scale", "Gumbel", "Gumbel with fixed scale"


#If you want different family than already coded or different method than implemented, it is easy! For that,
#simply add into function "estimate_epsilon(Y, X, family="YOUR DISTRIBUTION FUNCTION")" a code how you want to estimate your epsilons
#We always used "gam" function from "MGCV" package for the smooth parameter estimation. 


library(MASS)
library("mgcv")
library("gamlss")
library(dHSIC)
library(rje)
library(stringr)
library(EnvStats)
library(dplyr)
library(independence)

##################################   Example   ###########################################
#  n=500
#  X1 = rnorm(n)
#  X2=c();
#  for (i in 1:n) {
#    X2 = c(X2, rnorm(1, X1[i], X1[i]^2+1))
#  }
#  X= data.frame(X1, X2)
#  CPCM_graph_estimate(X, family_of_distributions = 2) 
#  
#  plot(X)
############################ CPCM estimation of the causal graph ##########################
  

CPCM_graph_estimate <- function(X, family_of_distributions = 1, force_estimate=FALSE){  
    n = length(X[,1])
    lambda = 2 #penalty for more edges
    
    
    #Estimation \hat{S}
    estimate_support_of_X_and_pair_a_family_for_X<-function(X, family_of_distributions = 1){
      
      determine_support <- function(X) {
        
        # Check if the data is discrete (integer values only)
        if (all(X == floor(X)) & length(unique(X))< length(X)/10) {
          return('discrete')  # Poisson distribution
        }
        
        # Check if the data is within the interval [0, 1]
        if (all(X >= 0 & X <= 1)) {
          return('interval')  # Beta distribution
        }
        
        # Check for Gamma distribution characteristics
        if (all(X > 0)) {
          # Check for skewness to differentiate from Gaussian
          skewness <- mean((X - mean(X))^3) / (mean((X - mean(X))^2)^(3/2))
          if (skewness > 0) {
            return('half-line')  # Gamma distribution
          }
        }
        
        # If none of the above conditions are met, assume Gaussian
        return('full support')  # Gaussian distribution
      }
      
      if(family_of_distributions == 1){
        if( determine_support(X)=='full support') return('Gaussian with fixed sigma')
        if( determine_support(X)=='discrete') return('Poisson')
        if( determine_support(X)=='interval') return('Gaussian with fixed sigma')
        if( determine_support(X)=='half-line') return('Pareto')
      }
      
      if(family_of_distributions == 2){
        if( determine_support(X)=='full support') return('Gaussian')
        if( determine_support(X)=='discrete') return('Negative_binomial')
        if( determine_support(X)=='interval') return('Gaussian')
        if( determine_support(X)=='half-line') return('Gamma')
      }
      
      
    }
    
    #Estimation \hat{epsilon} using given F
    estimate_epsilon<-function(Y, X, family="Gaussian"){ 
      
      ### ### ### ### ### ### ### ###
      rename_columns_to_1_to_q <-function(X){ #Just making sure that our variables are called X1, X2, ...
        if (is.null(nrow(X))) {X = as.data.frame(X); names(X)[1] <- "X1"; return(X) }
        return(X %>% rename_with(~ str_c("X", seq_along(.))))
      }
      X = rename_columns_to_1_to_q(X)
      names(Y)<- "Y"
      ### ### ### ### ### ### ### ###
      formula <- function(X, withY=TRUE){ #Definition of formula. We want to distinguish between discrete and continous variables for GAM. But anyway this function basically returns "Y~s(X1) + s(X2)+...+s(Xq)"
        q = ncol(X); if (is.null(nrow(X))) {q=1} 
        if(withY==TRUE){ form="Y ~ "}
        if(withY==FALSE){form=" ~ "}
        for (i in 1:q) {#Discrete variable is if it contains <=9 different values
          if (  length(unique(X[,i]))>9  ) { 
            form=paste(form, paste0("s(X", i,  ")+")) 
          }
          if (  length(unique(X[,i]))<=9  ) {  form=paste(form, paste0("as.factor(X", i,  ")+"))  }
        }
        
        form=as.formula(  substr(form,1,nchar(form)-1)  )
        return(form)
      }
      ### ### ### ### ### ### ### ### Here we start with the estimation of epsilons for each family separately.### ### ### ### ### ### ### ###
      
      
      if (family=="Gaussian") {
        
        fit=  gam(list(formula(X, TRUE),formula(X, FALSE)), data.frame(Y,X) ,family=gaulss(), link=list("identity","logb"))
        
        residuals=c()
        for (i in 1:length(Y)) {
          residuals=c(residuals,  (Y[i] -fit$fitted.values[i,1])*(fit$fitted.values[i,2]) ) #fit$fitted.values[i,2]=1/sigma[i], fit$fitted.values[i,1]=mu[i]
        }
        
        
        residuals=pnorm(residuals)
        
        return(residuals)
      }
      
      if (family=="Gaussian with fixed sigma") {
        
        fit=  gam(list(formula(X, TRUE),~1), data.frame(Y,X) ,family=gaulss(), link=list("identity","logb"))
        
        residuals=c()
        for (i in 1:length(Y)) {
          residuals=c(residuals,  (Y[i] -fit$fitted.values[i,1])*(fit$fitted.values[i,2]) ) #fit$fitted.values[i,2]=1/sigma[i], fit$fitted.values[i,1]=mu[i]
        }
        
        residuals=pnorm(residuals)
        
        
        return(residuals)
      }
      
      
      
      if (family=="Pareto") { # -log(log(Pareto) = gumbls

        Y=-log(log(Y))
        fit=  gam(list(formula(X, TRUE), ~1), data.frame(Y,X) ,family=gumbls())
        residuals=c()
        for (i in 1:n) {
          residuals=c(residuals,  pgumbel(Y[i], fit$fitted.values[i,1], 1)  )
        }
        
        
        return(residuals)
        
      }
      if (family=="Gumbel") {
        
        fit=  gam(list(formula(X, TRUE), formula(X, FALSE)), data.frame(Y,X) ,family=gumbls())
        residuals=c()
        for (i in 1:n) {
          residuals=c(residuals,  pgumbel(Y[i], fit$fitted.values[i,1], exp(fit$fitted.values[i,2]))  )
        }
        
        return(residuals)
        
      }
      
      if (family=="Gumbel with fixed scale") {
        
        fit=  gam(list(formula(X, TRUE), ~1), data.frame(Y,X) ,family=gumbls())
        residuals=c()
        for (i in 1:n) {
          residuals=c(residuals,  pgumbel(Y[i], fit$fitted.values[i,1], exp(fit$fitted.values[i,2]))  )
        }
        
        return(residuals)
        
      }
      
      
      
      if (family=="Gamma") {
        
        fit=  gam(list(formula(X, TRUE),formula(X, FALSE)), data.frame(Y,X) ,family=gammals())
        
        residuals=c()
        for (i in 1:length(Y)) {
          residuals=c(residuals,  pgamma( Y[i], shape = 1/exp(fitted(fit)[i,2]), scale = fitted(fit)[i,1]*exp(fitted(fit)[i,2])   ) ) 
        }
        
        return(residuals)
        
      }
      
      if (family=="Gamma with fixed scale") {
        
        fit=  gam(list(formula(X, TRUE),~1), data.frame(Y,X) ,family=gammals())
        
        residuals=c()
        for (i in 1:length(Y)) {
          residuals=c(residuals,  pgamma( Y[i], shape = 1/exp(fitted(fit)[i,2]), scale = fitted(fit)[i,1]*exp(fitted(fit)[i,2])   )) 
        }
        
        
        
        return(residuals)
        
      }
      
      if (family=="Negative_binomial") {
        
        X= as.numeric(X[,1])
        fit=  gam(Y~s(X), data.frame(Y, X) ,family=nb())
        
        residuals=c()
        for (i in 1:length(Y)) {
          residuals=c(residuals,  pnbinom(Y[i], mu =fitted(fit)[i],size =fit$family$getTheta(TRUE)   )  ) 
        }
        
        return(residuals)
        
      }
      
      if (family=="Poisson") {
        
        fit <- gam(formula(X, TRUE), data.frame(Y,X) ,family=poisson())
        
        residuals=c()
        for (i in 1:length(Y)) {
          residuals=c(residuals,  ppois(Y[i], lambda = fitted(fit)[i])  ) 
        }
        
        
        return(residuals)
        
      }
      
      
      
      return(  "Family not implemented, sorry.")
      
    }
    
    
    
    
  bivariate_CPCM_graph_estimate <- function(X, family_of_distributions = 1, force_estimate=FALSE){

      #Are X1 and X2 idependent? If yes, we return an empty graph as an estimate
      X1 = X[,1]; X2 = X[,2]
      
      if(force_estimate==FALSE){
      if (n>500 &   all(X1 == floor(X1))==FALSE & all(X2 == floor(X2))==FALSE )
        { indep = hoeffding.D.test(X1, X2, na.rm = TRUE, collisions = FALSE, precision = 1e-05) }else{ indep=dhsic.test(data.frame(X1, X2))}
        
      if(indep$p.value>=0.05) {res =  as.data.frame(c("-", "-", "Empty graph", "Empty graph"));     
      colnames(res)<-c("Causal results (p-values)");    
      rownames(res)<-c("p-value 1-->2", "p-value 2-->1", "Score-based graph estimate", 'Testing estimate')
      }}
      

      #X1-->X2 direction
      Y = X2; X=data.frame(X1 = X1)
      if(family_of_distributions == 1 || family_of_distributions==2)
        {family2 = estimate_support_of_X_and_pair_a_family_for_X(Y, family_of_distributions = family_of_distributions)}else{family2 = family_of_distributions}
      r2=as.numeric(estimate_epsilon(Y, X, family = family2))
      
      #X2-->X1 direction
      Y= X1; X=data.frame(X1 = X2)
      if(family_of_distributions == 1 || family_of_distributions==2)
      {family1 = estimate_support_of_X_and_pair_a_family_for_X(Y, family_of_distributions = family_of_distributions)}else{family1 = family_of_distributions}
      r1= as.numeric(estimate_epsilon(Y, X, family = family1))
      
      #Which independence test should we choose? hoeffding.D.test or dhsic.test? hoeffding.D.test does not work well for discrete variables, but is faster
      if(family1 !='Poisson' & family1 !='Negative_binomial' &family2 !='Poisson' & family2 !='Negative_binomial' ){ 
        z1=hoeffding.D.test(r2, X1, na.rm = TRUE,  precision = 1e-05)$p.value
        z2=hoeffding.D.test(r1, X2, na.rm = TRUE,  precision = 1e-05)$p.value}else{
          z1=dhsic.test(data.frame(r2, X1))$p.value
          z2=dhsic.test(data.frame(r1, X2))$p.value
        }
      
      
      #Results written in a nice way
      if (z1>=z2) res1="1 --> 2"   #Note that we do not compare it with an empty graph - we just return empty graph if and only if we do not reject H_0: X indep Y
      if (z1<z2)  res1="2 --> 1"
      
      res2=res1
      if(z1>0.05 & z2>0.05) res2="Unidentifiable  (both directions are plausible)" 
      if (z1<0.05 & z2<0.05) res2 = "Assumptions not fulfilled (both directions are not plausible)"
      
      res=as.data.frame(c( round(z1, digits=6), round(z2, digits = 6), res1, res2, paste0(family1, ";", family2)))
      colnames(res)<-c("Results")
      rownames(res)<-c("p-value 1-->2", "p-value 2-->1", "Score-based graph estimate", 'Testing estimate', 'Families used' )
      return(res)
    }
    
    
  
  trivariate_CPCM_graph_estimate <- function(X, family_of_distributions = 1){
    showing_computing_time=TRUE #Change to FALSE if annoyed by showing the remaining time when computing
    accuracy = 0.0025 #multIndepTest requires accuracy of the computed p_value
    
    make_all_dags_on_three_nodes <-function(){#I am sorry about this monstrosity. Thanks chatGPT for generating it, it was easier
      all_dags = list()
      all_dags[[length(all_dags)+1]] = matrix(c(0,0,0,   0,0,0,    0,0,0), nrow = 3, ncol = 3, byrow = TRUE)
      all_dags[[length(all_dags)+1]] = matrix(c(0,1,0,   0,0,0,    0,0,0), nrow = 3, ncol = 3, byrow = TRUE)
      all_dags[[length(all_dags)+1]] = matrix(c(0,1,0,   0,0,1,    0,0,0), nrow = 3, ncol = 3, byrow = TRUE)
      all_dags[[length(all_dags)+1]] = matrix(c(0,1,0,   0,0,0,    0,1,0), nrow = 3, ncol = 3, byrow = TRUE)
      all_dags[[length(all_dags)+1]] = matrix(c(0,1,0,   0,0,0,    1,0,0), nrow = 3, ncol = 3, byrow = TRUE)
      all_dags[[length(all_dags)+1]] = matrix(c(0,1,0,   0,0,0,    1,1,0), nrow = 3, ncol = 3, byrow = TRUE)
      all_dags[[length(all_dags)+1]] = matrix(c(0,1,1,   0,0,0,    0,0,0), nrow = 3, ncol = 3, byrow = TRUE)
      all_dags[[length(all_dags)+1]] = matrix(c(0,1,1,   0,0,1,    0,0,0), nrow = 3, ncol = 3, byrow = TRUE)
      all_dags[[length(all_dags)+1]] = matrix(c(0,1,1,   0,0,0,    0,1,0), nrow = 3, ncol = 3, byrow = TRUE)
      all_dags[[length(all_dags)+1]] = matrix(c(0,0,1,   0,0,1,    0,0,0), nrow = 3, ncol = 3, byrow = TRUE)
      
      all_dags[[length(all_dags)+1]] = matrix(c(0,0,1,   0,0,0,    0,1,0), nrow = 3, ncol = 3, byrow = TRUE)
      all_dags[[length(all_dags)+1]] = matrix(c(0,0,1,   1,0,0,    0,0,0), nrow = 3, ncol = 3, byrow = TRUE)
      all_dags[[length(all_dags)+1]] = matrix(c(0,0,1,   1,0,1,    0,0,0), nrow = 3, ncol = 3, byrow = TRUE)
      all_dags[[length(all_dags)+1]] = matrix(c(0,0,1,   0,0,0,    0,0,0), nrow = 3, ncol = 3, byrow = TRUE)
      all_dags[[length(all_dags)+1]] = matrix(c(0,0,0,   1,0,0,    0,0,0), nrow = 3, ncol = 3, byrow = TRUE)
      all_dags[[length(all_dags)+1]] = matrix(c(0,0,0,   0,0,1,    0,0,0), nrow = 3, ncol = 3, byrow = TRUE)
      all_dags[[length(all_dags)+1]] = matrix(c(0,0,0,   0,0,0,    1,0,0), nrow = 3, ncol = 3, byrow = TRUE)
      all_dags[[length(all_dags)+1]] = matrix(c(0,0,0,   0,0,0,    0,1,0), nrow = 3, ncol = 3, byrow = TRUE)
      all_dags[[length(all_dags)+1]] = matrix(c(0,0,0,   1,0,0,    0,1,0), nrow = 3, ncol = 3, byrow = TRUE)
      all_dags[[length(all_dags)+1]] = matrix(c(0,0,0,   1,0,0,    1,0,0), nrow = 3, ncol = 3, byrow = TRUE)
      
      all_dags[[length(all_dags)+1]] = matrix(c(0,0,0,   1,0,0,    1,1,0), nrow = 3, ncol = 3, byrow = TRUE)
      all_dags[[length(all_dags)+1]] = matrix(c(0,0,0,   1,0,1,    0,1,0), nrow = 3, ncol = 3, byrow = TRUE)
      all_dags[[length(all_dags)+1]] = matrix(c(0,0,0,   0,0,0,    1,1,0), nrow = 3, ncol = 3, byrow = TRUE)
      all_dags[[length(all_dags)+1]] = matrix(c(0,0,0,   1,0,1,    1,0,0), nrow = 3, ncol = 3, byrow = TRUE)
      all_dags[[length(all_dags)+1]] = matrix(c(0,0,0,   1,0,1,    0,0,0), nrow = 3, ncol = 3, byrow = TRUE)
      
      return(all_dags)}
    all_dags = make_all_dags_on_three_nodes()
    
    from_matrix_to_values<-function(variable, dag){
      if (variable ==1) {
        if (all(dag[,variable]==c(0,0,0))) {return(1)}
        else if (all(dag[,variable]==c(0,1,0))) {return(2)}
        else if (all(dag[,variable]==c(0,0,1))) {return(3)}
        else if (all(dag[,variable]==c(0,1,1))) {return(4)}
      }
      if (variable ==2) {
        if (all(dag[,variable]==c(0,0,0))) {return(1)}
        else if (all(dag[,variable]==c(1,0,0))) {return(2)}
        else if (all(dag[,variable]==c(0,0,1))) {return(3)}
        else if (all(dag[,variable]==c(1,0,1))) {return(4)}
      }  
      if (variable ==3) {
        if (all(dag[,variable]==c(0,0,0))) {return(1)}
        else if (all(dag[,variable]==c(1,0,0))) {return(2)}
        else if (all(dag[,variable]==c(0,1,0))) {return(3)}
        else if (all(dag[,variable]==c(1,1,0))) {return(4)}
      }  
      
      
    }
    res=c(); #res = final result
    ##################We coded it only for d=2 and d=3 using ugly brute-force
    
    #TRIVARIATE CASE
    if (ncol(X)==3) {
      X1 = X[,1]; X2 = X[,2]; X3 = X[,3]
      

      

      if(family_of_distributions == 1 || family_of_distributions==2)
      {family1 = estimate_support_of_X_and_pair_a_family_for_X(X1, family_of_distributions = family_of_distributions)
      family2 = estimate_support_of_X_and_pair_a_family_for_X(X2, family_of_distributions = family_of_distributions)
      family3 = estimate_support_of_X_and_pair_a_family_for_X(X3, family_of_distributions = family_of_distributions)
      }else{
        family1 = family_of_distributions
        family2 = family_of_distributions
        family3 = family_of_distributions}
      
      
      r123=as.numeric( estimate_epsilon(Y= X1, X=data.frame(X1=  X2  , X2=  X3  ), family = family1) )
      r12=as.numeric(estimate_epsilon(Y= X1, X=data.frame(X1=  X2   ), family = family1) )
      r13=as.numeric(estimate_epsilon(Y= X1, X=data.frame(X1=   X3  ), family = family1) )
      
      r213=as.numeric(estimate_epsilon(Y= X2, X=data.frame(X1=  X1  , X2=  X3  ), family = family2) )
      r21=as.numeric(estimate_epsilon(Y= X2, X=data.frame(X1=  X1   ), family = family2) )
      r23=as.numeric(estimate_epsilon(Y= X2, X=data.frame(X1=  X3  ), family = family2) )
      
      r312=as.numeric(estimate_epsilon(Y= X3, X=data.frame(X1=  X1  , X2=  X2  ), family = family3) )
      r31 =as.numeric(estimate_epsilon(Y= X3, X=data.frame(X1=  X1 ), family = family3) )
      r32 =as.numeric(estimate_epsilon(Y= X3, X=data.frame(X1=  X2 ), family = family3) )
      
      
      A = data.frame(X1, r12, r13, r123)
      B = data.frame(X2, r21, r23, r213)
      C = data.frame(X3, r31, r32, r312)
      
      p_values = c(); p_values_adjusted = c()
      k = length(all_dags)
      for (i in 1:length(all_dags)) {
        if(showing_computing_time==TRUE){cat("Time remaining:",k, "\n"); k=k-1}
        dag = all_dags[[i]]
        p_value =  multIndepTest(data.frame(A[from_matrix_to_values(1,dag)], 
                                            B[from_matrix_to_values(2, dag)], 
                                            C[from_matrix_to_values(3, dag)]) , 
                                 d=c(1,1,1), verbose = FALSE, N=as.integer(1/(2*accuracy)-1) )$fisher.pvalue
        if(p_value<=accuracy){p_value = 0}
        p_values = c(p_values,  p_value)
        p_values_adjusted = c(p_values_adjusted, -log(p_value) + lambda*sum(dag))
      }
      
      M_best = all_dags[[ which.min(p_values_adjusted)    ]]
      M_plausible = list()
      plausible = 1:length(all_dags); plausible = plausible[p_values>0.01]
      for (i in which(p_values > 0.01)[order(p_values_adjusted[plausible])]) {
        M_plausible = c(M_plausible, list(all_dags[[i]]))
      }
      
      
      return(list(plausible=M_plausible, p_values=p_values, p_values_adjusted=p_values_adjusted, result=M_best))
      
    }

  }
  

  if(ncol(X)==2) return(bivariate_CPCM_graph_estimate(X, family_of_distributions=family_of_distributions, force_estimate = force_estimate))
  if(ncol(X)==3) return(trivariate_CPCM_graph_estimate(X, family_of_distributions=family_of_distributions))
  
  if (ncol(X)!=2 && ncol(X)!=3) {
    return("Sorry, the code is avaliable only for d<=3 for now")
  }
  


































  }  
  
  
  
  
  
  
  
  
  
  
  
  
