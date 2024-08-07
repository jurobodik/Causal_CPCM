#Code taken from Natasa Tagasovska et al, https://github.com/tagas/bQCD.git
#This file contains implementations of other causal discovery methods 
#If you want to run these lines, download the zip from https://github.com/tagas/bQCD.git and set directory to it. 


###################################################################################
# causal inference methods
# 
# implementation details:
# if 'X -> Y' they output cd  = 1
# if 'Y -> X' they output cd = 0
# epsilon: confidence (or score)

#library(CAM)
library(Hmisc)
library(extraDistr)
library(quantregForest)
library(rvinecopulib)
library(statmod)
library(qrnn)

#Set directory to the downloaded zip bQCD-master.zip from https://github.com/tagas/bQCD.git
source("R/baselines/Slope/Slope.R")
source("R/baselines/Slope/utilities.R")
source("R/baselines/GRAN_EMD/kernelMethod.R")
source("R/baselines/GRAN_EMD/rkhse.R")
source("R/baselines/Slope/resit/code/startups/startupLINGAM.R", chdir = TRUE)
source("R/baselines/Slope/resit/code/startups/startupICML.R", chdir = TRUE)

#Set directory to the downloaded loci-master.zip from https://github.com/AlexImmer/loci
library(feather)
library(reticulate)
my_module <- import("__main__") #my_module$loci should work now


##### Infer the best model

## Slope: Working title for algorithm

##Params:
# a + bx + cx^2 + dx^3 + e*exp(x) + fx^-1 + gx^-2 + h*log2(x)


#!/usr/bin/Rscript

sysouts = T



### na to zero
naTo0 = function(c){
  c[is.na(c)] = 0
  return(c)
}

###### Normalization
normX = function(x, n){
  if(min(x) == max(x)){
    return( rep(n, length(x)) )
  }else{
    return( ((x-min(x)) / (max(x) - min(x))) * n )
  }
}

logg = function(x){
  if(x == 0){
    return(0)
  }else{
    return(log2(x))
  }
}
log2fac = function(n){
  sum = 0
  for(i in 2:n){
    sum = sum + logg(i)
  }
  return(sum)
}
log2nChoosek = function(n, k){
  if(k > n | k == 0){
    return(0)
  }else{
    return(log2fac(n) - log2fac(k) - log2fac(n-k))
  }
}
logN = function(z){
  z = ceiling(z)
  if(z < 1){
    return(0)
  }else{
    logstar = logg(z)
    sum = logstar
    while(logstar > 0){
      logstar = logg(logstar)
      sum = sum + logstar
    }
    return(sum + logg(2.865064))
  }
}

fitTestWMW = function(y,x,cycles=100,alpha=0.05){
  set.seed(1234)
  sigCount = 0
  for(i in 1:cycles){
    random_permutation = sample(length(x))
    cp = ceiling(length(x) / 2)
    x.train = x[random_permutation[1:cp]]
    y.train = y[random_permutation[1:cp]]
    x.test  = x[random_permutation[(cp+1):length(x)]]
    y.test  = y[random_permutation[(cp+1):length(x)]]
    res = findBestFit(y.train,x.train)
    err.train = (y.train - fofx(x.train, res$par))^2
    err.test  = (y.test  - fofx(x.test,  res$par))^2
    wmw.res = wilcox.test(err.train, err.test, alternative="two.sided")
    p.value = wmw.res$p.value
    if(!is.na(p.value) & p.value < alpha){
      sigCount = sigCount + 1
    }
  }
  return(sigCount/cycles)
}

getDuplicatePositions = function(x){
  last = x[1]
  pos = 0
  s = ""
  for(i in 2:length(x)){
    if(x[i] != last){
      last = x[i]
      if(i - (pos + 1) > 1){
        k = paste(pos, (i-1), sep=":")
        if(nchar(s) > 0){
          s = paste(s, k, sep=";")
        }else{
          s = k
        }
      }
      pos = i - 1
    }
  }
  # last elem involved
  if((pos + 1) != length(x)){
    k = paste(pos, length(x), sep=":")
    if(nchar(s) > 0){
      s = paste(s, k, sep=";")
    }else{
      s = k
    }
  }
  return(s)
}

###### Entropy and conditional entropy
entropy = function(x){
  x.unique = unique(sort(x))
  N = length(x)
  H = 0
  for(i in 1:length(x.unique)){
    frac = sum(x == x.unique[i]) / N
    H = H - (frac * logg(frac))
  }
  return(H)
}
conditionalEntropy = function(x,y){
  y.unique = unique(sort(y))
  N = length(y)
  H = 0
  for(i in 1:length(y.unique)){
    xiy = x[y == y.unique[i]]
    H = H + ( (length(xiy) / N) * entropy(xiy) )
  }
  return(H)
}

###### Discretization
equiWidthBinning = function(dataset, bins){
  newd = cut(dataset, bins)
  return(as.numeric(newd))
}

###### Helper
getMeanOccurance = function(){
  df = data.frame(I=uv, X1=rep(0,length(uv)), SD1=rep(0,length(uv)), X2=rep(0,length(uv)), SD2=rep(0,length(uv)), R=ref.uv$V6)
  for(i in 1:length(uv)){
    t = readI(uv[i])
    df$X1[i] = mean(table(t[,1]))
    df$SD1[i] = sd(table(t[,1]))
    df$X2[i] = mean(table(t[,2]))
    df$SD2[i] = sd(table(t[,2]))
  }
  return(df)
}
getMinOccurance = function(){
  df = data.frame(I=uv, X1=rep(0,length(uv)), SD1=rep(0,length(uv)), X2=rep(0,length(uv)), SD2=rep(0,length(uv)), R=ref.uv$V6)
  for(i in 1:length(uv)){
    t = readI(uv[i])
    df$X1[i] = min(table(t[,1]))
    df$SD1[i] = sd(table(t[,1]))
    df$X2[i] = min(table(t[,2]))
    df$SD2[i] = sd(table(t[,2]))
  }
  return(df)
}
getDims = function(){
  df = data.frame(I=uv, L=rep(0,length(uv)))
  for(i in 1:length(uv)){
    t = readI(uv[i])
    df$L[i] = dim(t)[1]
  }
  return(df)
}

###### Calculate decision rate (-- = 1/2)
decision_rate_w = function(res){
  return(decision_rate(res$Correct, res$Eps, res$Cds))
}
decision_rate = function(corr, eps, cds){
  df = data.frame(A=corr, B=abs(eps), C=cds)
  df = df[with(df, order(-B)),]
  sum = 0
  dr = rep(0, dim(df)[1])
  for(i in 1:dim(df)[1]){
    if(df$C[i] == "--"){
      sum = sum + 0.5
    }else{
      sum = sum + df$A[i]
    }
    dr[i] = sum / i
  }
  return(c(1,dr))
}
decision_rate_meta = function(res, uv=uv, oneHalf=T){
  corr = res$Correct
  eps = res$Eps
  cds = res$Cds
  df = data.frame(Uv=uv, A=corr, B=abs(eps), C=cds)
  df = df[with(df, order(-B)),]
  sum = 0
  total = 0
  step = rep(0, dim(df)[1])
  dr = rep(0, dim(df)[1])
  for(i in 1:dim(df)[1]){
    val = meta$V6[df$Uv[i]]
    total = total + val
    if(df$C[i] == "--"){
      sum = sum + 0.5 * val
      if(oneHalf){
        dr[i] = sum / total
      }else{
        dr[i] = dr[i-1]
      }
    }else{
      sum = sum + df$A[i] * val
      dr[i] = sum / total
    }
    step[i] = total
  }
  step = step / total
  return(list(D=c(1,dr), S=c(0,step)))
}
decision_rate_meta_pos = function(res, uv=uv){
  corr = res$Correct
  eps = res$Eps
  cds = res$Cds
  df = data.frame(Uv=uv, A=corr, B=abs(eps), C=cds)
  df = df[with(df, order(B)),]
  sum = 0
  total = 0
  step = rep(0, dim(df)[1])
  dr = rep(0, dim(df)[1])
  for(i in 1:dim(df)[1]){
    val = meta$V6[df$Uv[i]]
    total = total + val
    if(df$C[i] == "--"){
      sum = sum + 0.5 * val
    }else{
      sum = sum + df$A[i] * val
    }
    dr[i] = sum / total
    step[i] = total
  }
  step = step / total
  return(list(D=c(1,dr), S=c(0,step)))
}
decision_rate_w_pos = function(res){
  corr = res$Correct
  eps = res$Eps
  cds = res$Cds
  df = data.frame(A=corr, B=abs(eps), C=cds)
  df = df[with(df, order(B)),]
  sum = 0
  dr = rep(0, dim(df)[1])
  for(i in 1:dim(df)[1]){
    if(df$C[i] == "--"){
      sum = sum + 0.5
    }else{
      sum = sum + df$A[i]
    }
    dr[i] = sum / i
  }
  return(c(1,dr))
}

###### SSE
SSE = function(x, yhead){
  e = sum((x-yhead)^2)
  return(e)
}

###### PlotI
plotI = function(i){plot(readI(i))}

###### MDL score
resolution = 0.01
setResolution = function(val){
  resolution <<- val
}
gaussian_score_emp = function(x){
  sse = SSE(x, mean(x))
  var = sse / length(x)
  sigma = sqrt(var)
  return(gaussian_score(sigma, x))
}
gaussian_score = function(sigma, x){
  sse = SSE(x, mean(x))
  n = length(x)
  sigmasq = sigma^2
  if(sse == 0.0 | sigmasq == 0.0){
    return(0.0)
  }else{
    err = (sse / (2 * sigmasq * log(2))) + ((n/2) * logg(2 * pi * sigmasq)) - n * logg(resolution)
    return(err)
  }
}
gaussian_score_emp_sse = function(sse, n){
  var = sse / n
  sigma = sqrt(var)
  return(gaussian_score_sse(sigma, sse, n))
}

gaussian_score_sse = function(sigma, sse, n){
  sigmasq = sigma^2
  if(sse == 0.0 | sigmasq == 0.0){
    return(0.0)
  }else{
    err = (sse / (2 * sigmasq * log(2))) + ((n/2) * logg(2 * pi * sigmasq)) - n * logg(resolution)
    return(max(err,0))
  }
}
parameterScore = function(model){
  sum = 0
  coeff = naTo0(model$coefficients)
  for(c in coeff){
    if(c != 0){
      ca = abs(c)
      cdummy = ca
      prec = 1
      while(cdummy < 1000){
        cdummy = cdummy * 10
        prec = prec + 1
      }
      sum = sum + logN(cdummy) + logN(prec) + 1
    }
  }
  return(sum)
}

### igci
S = function(x.uns){
  sum = 0.0
  x = sort(x.uns)
  m = length(x)
  for(i in 1:(m-1)){
    sum = sum + logg(abs(x[i+1] - x[i]))
  }
  sum = sum * (1/(m-1))
  sum = sum + logg(m)
  return(sum)
}
normalize = function(tt){
  t = tt
  # normalize to mean 0 and sd 1
  for(j in 1:(dim(t)[2])){
    t[,j] = as.numeric(t[,j])
    t[,j] = (t[,j] - min(t[,j])) / (max(t[,j]) - min(t[,j]))
  }
  return(t)
}


normalize_pobs = function(tt){
  n <- nrow(tt)
  u1 <- rank(as.numeric(tt[,1]), ties.method = "random")/(n + 1)
  u2 <- rank(as.numeric(tt[,2]),ties.method = "random")/(n + 1)
  return(cbind(u1,u2))
}



# case 1
# % uniform reference measure
# x = (x - min(x)) / (max(x) - min(x));
# y = (y - min(y)) / (max(y) - min(y));
# case 2
# % Gaussian reference measure
# x = (x - mean(x)) ./ std(x);
# y = (y - mean(y)) ./ std(y);
# otherwise

normalize_G = function(tt){
  t = tt
  # normalize to mean 0 and sd 1
  for(j in 1:(dim(t)[2])){
    t[,j] = as.numeric(t[,j])
    #t[,j] = (t[,j] - min(t[,j])) / (max(t[,j]) - min(t[,j]))
    t[,j] = (t[,j] - mean(t[,j]) / sd(t[,j]))
  }
  return(t)
}

## mindiff
mindiff = function(x){
  xs = sort(x)
  diff = 0.01
  for(i in 1:(length(x)-1)){
    new_diff = xs[i+1] - xs[i]
    #print(new_diff)
    #print(xs[i+1])
    if(new_diff != 0 & new_diff < diff){
      diff = new_diff
    }
  }
  return(diff)
}
data_precision = function(x){
  precision = 1
  set = x != round(x)
  x = x[set]
  while(length(x) > 0){
    precision = precision / 10
    x = 10 * x
    set = x != round(x)
    x = x[set]
  }
  return(precision)
}

##### IGCI
IGCI = function(t){
  t = normalize(t)#normalize(t)
  diffs = S(t[,2]) - S(t[,1])
  causd = "--"
  if(diffs < 0){
    causd = "->"
  }else if(diffs > 0){
    causd = "<-"
  }
  return(list(cd=causd, epsilon=diffs))
}

##### IGCI
IGCI_G = function(t){
  t = normalize_G(t)
  diffs = S(t[,2]) - S(t[,1])
  causd = "--"
  if(diffs < 0){
    causd = "->"
  }else if(diffs > 0){
    causd = "<-"
  }
  return(list(cd=causd, epsilon=diffs))
}





fofx = function(x, p){
  x.neg = x
  x.neg[x == 0] = resolution
  y.head = p$a + p$b * x + p$c * (x^2) + p$d * (x^3) + p$e * (exp(x)) + p$f * (x.neg^(-1)) + p$g * (x.neg^(-2))
  return(y.head)
}
getFunctionIndex = function(p){
  if(abs(p$b) > 0){
    return(1)
  }else if(abs(p$c) > 0){
    return(2)
  }else if(abs(p$d) > 0){
    return(3)
  }else if(abs(p$e) > 0){
    return(4)
  }else if(abs(p$f) > 0){
    return(5)
  }else if(abs(p$g) > 0){
    return(6)
  }else{
    return(0)
  }
}
###### Best fit
fitLin = function(y,x){
  m = lm(y~x)
  params = list(a=0, b=0, c=0, d=0, e=0, f=0, g=0, h=0)
  sse = tail(anova(m)[,2],1)
  coeff = naTo0(m$coefficients)
  params$a = coeff[1]
  params$b = coeff[2]
  res = list(sse=sse, model=parameterScore(m), par=params)
  return(res)
}
fitQuad = function(y,x){
  m = lm(y~I(x^2))
  params = list(a=0, b=0, c=0, d=0, e=0, f=0, g=0, h=0)
  sse = tail(anova(m)[,2],1)
  coeff = naTo0(m$coefficients)
  params$a = coeff[1]
  params$c = coeff[2]
  res = list(sse=sse, model=parameterScore(m), par=params)
  return(res)
}
fitCub = function(y,x){
  m = lm(y~I(x^3))
  params = list(a=0, b=0, c=0, d=0, e=0, f=0, g=0, h=0)
  sse = tail(anova(m)[,2],1)
  coeff = naTo0(m$coefficients)
  params$a = coeff[1]
  params$d = coeff[2]
  res = list(sse=sse, model=parameterScore(m), par=params)
  return(res)
}
fitExp = function(y,x){
  xe = exp(x)
  m = lm(y~xe)
  params = list(a=0, b=0, c=0, d=0, e=0, f=0, g=0, h=0)
  sse = tail(anova(m)[,2],1)
  coeff = naTo0(m$coefficients)
  params$a = coeff[1]
  params$e = coeff[2]
  res = list(sse=sse, model=parameterScore(m), par=params)
  return(res)
}
fitNegLin = function(y,x){
  x[x == 0] = resolution
  xe = x^(-1)
  m = lm(y~xe)
  params = list(a=0, b=0, c=0, d=0, e=0, f=0, g=0, h=0)
  sse = tail(anova(m)[,2],1)
  coeff = naTo0(m$coefficients)
  params$a = coeff[1]
  params$f = coeff[2]
  res = list(sse=sse, model=parameterScore(m), par=params)
  return(res)
}
fitNegQuad = function(y,x){
  x[x == 0] = resolution
  xe = x^(-1)
  m = lm(y~xe)
  params = list(a=0, b=0, c=0, d=0, e=0, f=0, g=0, h=0)
  sse = tail(anova(m)[,2],1)
  coeff = naTo0(m$coefficients)
  params$a = coeff[1]
  params$g = coeff[2]
  res = list(sse=sse, model=parameterScore(m), par=params)
  return(res)
}

findBestFit = function(y,x){
  r1 = fitLin(y,x)
  r2 = fitQuad(y,x)
  r3 = fitCub(y,x)
  r4 = fitExp(y,x)
  r5 = fitNegLin(y,x)
  r6 = fitNegQuad(y,x)
  s = rep(0,6)
  s[1] = gaussian_score_emp_sse(r1$sse, length(x)) + r1$model
  s[2] = gaussian_score_emp_sse(r2$sse, length(x)) + r2$model
  s[3] = gaussian_score_emp_sse(r3$sse, length(x)) + r3$model
  s[4] = gaussian_score_emp_sse(r4$sse, length(x)) + r4$model
  s[5] = gaussian_score_emp_sse(r5$sse, length(x)) + r5$model
  s[6] = gaussian_score_emp_sse(r6$sse, length(x)) + r6$model
  i = which(s == min(s))
  if(i[1] == 1){
    return(r1)
  }else if(i[1] == 2){
    return(r2)
  }else if(i[1] == 3){
    return(r3)
  }else if(i[1] == 4){
    return(r4)
  }else if(i[1] == 5){
    return(r5)
  }else{
    return(r6)
  }
}

fitWrapperStirr = function(y,x){
  maxX = max(x)
  xf = round(x*100000)
  tx = table(xf)
  if(length(tx) <= 2){ ## binary data
    return(list(s=0, b=-2, p=-2))
  }
  mx = mean(tx)
  sd = sd(tx)
  costsND = 2147483647
  potential_bins = 0
  if(mx >= 10){
    potential_bins = 1
    xi = xf %in% as.numeric(names(tx)[tx >= 5])
    xfg = xf[xi]
    xg = x[xi]
    yg = y[xi]
    score1 = 0
    score2 = 0
    score3 = 0
    score4 = 0
    score5 = 0
    score6 = 0
    sse = 0.0
    if(sum(!xi) > 1){
      xl = x[!xi]
      yl = y[!xi]
      resSingletons = findBestFit(yl,xl)
      sse = gaussian_score_emp_sse(resSingletons$sse, length(xl)) + resSingletons$model
    }
    for(e in unique(sort(xfg))){
      ones = xfg == e
      yt = sort(yg[ones])
      xt = normX(1:length(yt), 10) - 5.0
      f1 = fitLin(yt, xt)
      f2 = fitQuad(yt, xt)
      f3 = fitCub(yt, xt)
      f4 = fitExp(yt, xt)
      f5 = fitNegLin(yt, xt)
      f6 = fitNegQuad(yt, xt)
      score1 = score1 + gaussian_score_emp_sse(f1$sse, length(xt)) + f1$model
      score2 = score2 + gaussian_score_emp_sse(f2$sse, length(xt)) + f2$model
      score3 = score3 + gaussian_score_emp_sse(f3$sse, length(xt)) + f3$model
      score4 = score4 + gaussian_score_emp_sse(f4$sse, length(xt)) + f4$model
      score5 = score5 + gaussian_score_emp_sse(f5$sse, length(xt)) + f5$model
      score6 = score6 + gaussian_score_emp_sse(f6$sse, length(xt)) + f6$model
    }
    score = sse + min(score1,score2,score3,score4,score5,score6)
    modelCosts = 1 + log2nChoosek(length(x)-1, length(unique(sort(xfg)))-1) + logg(length(unique(sort(xfg)))) ## cut points
    costs = score + modelCosts
    return(list(s=costs, b=c(1)))
    if(costs < costsND){
      costsND = costs
    }
  }
  res = findBestFit(y,x)
  score = gaussian_score_emp_sse(res$sse, length(x)) + res$model
  modelCosts = 1 ## additional model cost of 1 cause no bins have to be identified
  costsD = score + modelCosts
  if(costsND < costsD){
    return(list(s=costsND, b=1, p=potential_bins))
  }else{
    return(list(s=costsD, b=0, p=potential_bins))
  }
}

fitComparison = function(fit, l, lx, score, bins_old){
  newScore = gaussian_score_emp_sse(fit$sse, l) + fit$model
  delta = newScore - score + log2nChoosek(lx-1, bins_old+1) + logg(bins_old+1) - logg(bins_old) - log2nChoosek(lx-1, bins_old)
  return(delta)
}

fitWrapper = function(y,x){
  minNum = 5
  maxX = max(x)
  xf = round(x*100000)
  tx = table(xf)
  if(length(tx) <= 2){  ## binary data
    return(list(s=0, b=-2, p=-2))
  }
  mx = mean(tx)
  sd = sd(tx)
  lx = length(x)
  res = findBestFit(y,x)
  fun = getFunctionIndex(res$par)
  score = gaussian_score_emp_sse(res$sse, length(x))
  modelCosts = 1 + res$model ## additional model cost of 1 cause no bins have to be identified
  costs = score + modelCosts
  xi = xf %in% as.numeric(names(tx)[tx >= minNum])
  xfg = xf[xi]
  xg = x[xi]
  yg = y[xi]
  score1 = score
  score2 = score
  score3 = score
  score4 = score
  score5 = score
  score6 = score
  bins = rep(1,6)
  for(e in unique(sort(xfg))){
    ones = xfg == e
    yt = sort(yg[ones])
    xt = normX(1:length(yt), 10) - 5.0
    curr_l = length(xt)
    old_sse = sum((yg[ones] - fofx(xg[ones], res$par))^2)
    old_score = gaussian_score_emp_sse(old_sse, curr_l)
    f1 = fitLin(yt, xt)
    f2 = fitQuad(yt, xt)
    f3 = fitCub(yt, xt)
    f4 = fitExp(yt, xt)
    f5 = fitNegLin(yt, xt)
    f6 = fitNegQuad(yt, xt)
    s1 = fitComparison(f1, curr_l, lx, old_score, bins[1])
    if(s1 < 0){
      bins[1] = bins[1] + 1
      score1 = score1 + s1
    }
    s2 = fitComparison(f2, curr_l, lx, old_score, bins[2])
    if(s2 < 0){
      bins[2] = bins[2] + 1
      score2 = score2 + s2
    }
    s3 = fitComparison(f3, curr_l, lx, old_score, bins[3])
    if(s3 < 0){
      bins[3] = bins[3] + 1
      score3 = score3 + s3
    }
    s4 = fitComparison(f4, curr_l, lx, old_score, bins[4])
    if(s4 < 0){
      bins[4] = bins[4] + 1
      score4 = score4 + s4
    }
    s5 = fitComparison(f5, curr_l, lx, old_score, bins[5])
    if(s5 < 0){
      bins[5] = bins[5] + 1
      score5 = score5 + s5
    }
    s6 = fitComparison(f6, curr_l, lx, old_score, bins[6])
    if(s6 < 0){
      bins[6] = bins[6] + 1
      score6 = score6 + s6
    }
  }
  # correct for all divided
  for(i in 1:length(bins)){
    if(bins[i] > length(tx)){
      bins[i] = bins[i] - 1
    }
  }
  # add model costs
  score1 = score1 + 1 + log2nChoosek(length(x)-1, bins[1]-1)
  score2 = score2 + 1 + log2nChoosek(length(x)-1, bins[2]-1)
  score3 = score3 + 1 + log2nChoosek(length(x)-1, bins[3]-1)
  score4 = score4 + 1 + log2nChoosek(length(x)-1, bins[4]-1)
  score5 = score5 + 1 + log2nChoosek(length(x)-1, bins[5]-1)
  score6 = score6 + 1 + log2nChoosek(length(x)-1, bins[6]-1)
  scores = c(score1,score2,score3,score4,score5,score6)
  i = which(scores == min(scores))[1]
  score = scores[i]
  costs = score
  potential_bins = sum(tx >= minNum)
  if(potential_bins != length(tx)){
    potential_bins = potential_bins + 1
  }
  if(bins[i] > 1){
    fun = i
  }
  return(list(s=score, b=bins[i], p=potential_bins, f=fun))
}

fitFG = function(y,x){
  xf = round(x*100000)
  tx = table(xf)
  mx = mean(tx)
  ## if non-deterministic
  if(mx >= 10){
    res = findBestFit(y,x)
    model = res$model
    y.rem = y - fofx(x, res$par)
    xt = normX(1:length(y), 10) - 5.0
    res2 = findBestFit(sort(y.rem),xt)
    model = model + res$model
    data = gaussian_score_emp_sse(res2$sse, length(x))
    return(data + model)
  }else{
    return(fitWrapper(y,x))
  }
}

defaultScore = function(x){
  score = -logg(resolution) * length(x)
  return(score)
}

### Slope algorithm
Slope = function(t, prune=1.0, alpha=0.001){
  ## remove percentage of outliers
  if(prune < 1.0){
    require("ldbod")
    count = floor(dim(t)[1] * prune)
    scores = ldbod(t, k=10)
    t <- t[order(scores$lof,decreasing=F)[1:count],]
  }
  x = normX(t[,1],1)
  y = normX(t[,2],1)
  print("Calculate X->Y...")
  setResolution(mindiff(y))
  print(resolution)
  dy = defaultScore(y)
  resXtoY = fitWrapper(y,x)
  sseXtoY = resXtoY$s
  
  print("Calculate Y->X...")
  setResolution(mindiff(x))
  print(resolution)
  dx = defaultScore(x)
  resYtoX = fitWrapper(x,y)
  sseYtoX = resYtoX$s
  
  dXY = sseXtoY + dx
  dYX = sseYtoX + dy
  dXtoY = dXY / (dx + dy)
  dYtoX = dYX / (dx + dy)
  
  # Get delta
  eps = dXtoY - dYtoX
  pv = 2^(-(abs(dXY - dYX)/2))
  ## if data is binary
  if(resXtoY$p == -2 | resYtoX$p == -2){
    eps = 0.0
    pv = 1.0
  }
  
  # Determine causal direction
  causd = "--"
  if(abs(eps) > 0.0 & pv < alpha){
    if(eps < 0){
      causd = "->"
    }else{
      causd = "<-"
    }
  }
  
  r = list(epsilon = eps, cd = causd, p.value=pv)
  return(r)
}

SlopeInfo = function(t, prune=1.0, alpha=0.001){
  ## remove percentage of outliers
  if(prune < 1.0){
    require("ldbod")
    count = floor(dim(t)[1] * prune)
    scores = ldbod(t, k=10)
    t <- t[order(scores$lof,decreasing=F)[1:count],]
  }
  x = normX(t[,1],1)
  y = normX(t[,2],1)
  #print("Calculate X->Y...")
  setResolution(mindiff(y))
  #print(resolution)
  dy = defaultScore(y)
  resXtoY = fitWrapper(y,x)
  sseXtoY = resXtoY$s
  
  #print("Calculate Y->X...")
  setResolution(mindiff(x))
  #print(resolution)
  dx = defaultScore(x)
  resYtoX = fitWrapper(x,y)
  sseYtoX = resYtoX$s
  
  dXY = sseXtoY + dx
  dYX = sseYtoX + dy
  dXtoY = dXY / (dx + dy)
  dYtoX = dYX / (dx + dy)
  
  # Get delta
  eps = dXtoY - dYtoX
  pv = 2^(-(abs(dXY - dYX)/2))
  ## if data is binary
  if(resXtoY$p == -2 | resYtoX$p == -2){
    eps = 0.0
    pv = 1.0
  }
  
  # Determine causal direction
  causd = "--"
  if(abs(eps) > 0.0 & pv < alpha){
    if(eps < 0){
      causd = "->"
    }else{
      causd = "<-"
    }
  }
  
  r = list(epsilon = eps, cd = causd, p.value=pv, sc=c(dXY, dYX))
  return(r)
}




SlopeWrap <- function(pair){
  res = SlopeInfo(data.frame(pair), alpha = 1.01)
  cd = NA
  epsilon = NA
  if(res$cd == "->") cd = 1
  else if(res$cd == "<-") cd = 0
  else if(res$cd == "--") cd = NA
  list(epsilon = res$eps, cd = cd)
}

IGCIWrap_G <- function(pair){
  res = IGCI_G(data.frame(pair))
  cd = NA
  epsilon = NA
  if(res$cd == "->") cd = 1
  else if(res$cd == "<-") cd = 0
  list(epsilon = res$eps, cd = cd)
  
}

IGCIWrap_Unif <- function(pair){
  res = IGCI(data.frame(pair))
  cd = NA
  epsilon = NA
  if(res$cd == "->") cd = 1
  else if(res$cd == "<-") cd = 0
  list(epsilon = res$eps, cd = cd)
  
}

lingamWrapper = function(t) {
  res = tryCatch({
    r = lingamWrap(t)
  }, error = function(e) {
    print(e)
    return(list(
      B = matrix(c(0, 0, 0, 0), nrow = 2, ncol = 2),
      Adj = matrix(c(F, F, F, F), nrow = 2, ncol = 2)
    ))
  })
  C = 1 * res$Adj
  if (sum(C) != 1) {
    causd = "--"
    p_val = 0
  } else{
    if (C[2, 1] == 1) {
      causd = 0#"<-"
      p_val = res$B[1, 2]
    } else if (C[1, 2] == 1) {
      causd = 1#"->"
      p_val = res$B[2, 1]
    } else{
      causd = "--"
      p_val = 0
    }
  }
  print(res)
  return(list(cd = causd, epsilon = p_val))
  #causd
}

CAMWrapper <- function(X){
  
  res = CAM(X, scoreName = "SEMGAM")
  cd = NA
  if(res$Adj[3] == 1) cd = 1
  else if(res$Adj[2] == 1) cd = 0
  return(list(cd = cd, epsilon = res$Score))
}


ICMLWrapper <- function(X){
  res=ICML(X,
           model = train_gp,
           indtest = indtestHsic,
           output = FALSE)
  cd = NA
  epsilon = NA
  if(res$Cd == "->") cd = 1
  else if(res$Cd == "<-") cd = 0
  return(list(cd = cd, epsilon = res$Eps))
}


# ****synthetic data only****
# excript from original paper: Chen 2014 ; same parameters used in Gaussianity measures 2016
# For EMD we use gaussian kernel with width 1/5*dm, 
# where d_m is the median of distances among all input patterns.

RKHSEWrap <- function(t) {
  res = determine_causal_direction_rkhse(t[,1], t[,2])
  cd = NA
  epsilon = NA
  if(res$relation == "x->y") cd = 1
  else if(res$relation == "y->x") cd = 0
  list(cd = cd, epsilon = res$conf)
}


### test gaussianity measures
###
###

GaussianityWrap <- function(t) {
  res = determine_causal_direction(t[, 1], t[, 2])
  cd = NA
  epsilon = NA
  if(res$relation == "x->y") cd = 1
  else if(res$relation == "y->x") cd = 0
  list(cd = cd, epsilon = res$conf)
}


## Proper scoring rule for predicted quantiles
quantileScoring <- function(actual, pred, prob = 0.95) {
  mean((as.numeric(actual <= pred) - prob) * (pred - actual))
}

## Quantile Causal discovery method nonparametric copula
QCCD <- function(pair, m=1) {
  
  # to get reproducible jittering results for discreet data
  #set.seed(0)
  
  n <- nrow(pair)
  # Recover pseudo-observations and estimate the copula non-parametrically
  u <- apply(pair, 2, function(x) rank(x, ties.method = "random")/(n + 1))
  cop <- bicop(data = u,
               family_set = "tll",
               nonpar_method = "constant")
  
  # deal with discrete data
  pair_scaled <- qnorm(u)
  
  # integrate over quantiles
  if(n < 200){
    uw <- gauss.quad.prob(1)
  } else {
    uw <- gauss.quad.prob(m)
  }
  
  cls <- sapply(uw$nodes, function(uu) {
    u_pred <- cbind(predict(object = cop,
                            newdata = cbind(uu, u[, 2]),
                            what = "hinv2"),
                    predict(object = cop,
                            newdata = cbind(u[, 1], uu),
                            what = "hinv1"))
    
    # marginal and conditional quantiles
    marg_q <- sapply(1:2, function(i) quantile(pair_scaled[,i], uu))
    cond_q <- sapply(1:2, function(i) quantile(pair_scaled[,i], u_pred[, i]))
    
    # code lengths
    cl_marginal <- sapply(1:2, function(i)
      quantileScoring(pair_scaled[, i], marg_q[i], uu))
    cl_conditional <- sapply(1:2, function(i)
      quantileScoring(pair_scaled[, i], cond_q[, i], uu))
    
    c(cl_marginal, cl_conditional)
  })
  
  sel <- !apply(is.na(cls), 2, any)
  uw$weights <- uw$weights[sel] / sum(uw$weights[sel])
  cls <- apply(na.omit(t(cls)) * uw$weights, 2, sum)
  
  dx_to_y <- (cls[1] + cls[4])/sum(cls[1:2])
  dy_to_x <- (cls[2] + cls[3])/sum(cls[1:2])
  
  cd <- ifelse(dy_to_x > dx_to_y, 1, 0)
  
  epsilon <-  (-dx_to_y + dy_to_x )
  
  return(list(cd = cd, epsilon = epsilon))
}


## Quantile Causal discovery method - quantile forest
QCD_qrf <- function(pair, m=3) {
  
  # to get reproducible jittering results for discreet data
  #set.seed(0)
  
  n <- nrow(pair)
  # Recover pseudo-observations and estimate the copula non-parametrically
  u <- apply(pair, 2, function(x) rank(x, ties.method = "random")/(n + 1))
  # deal with discrete data
  pair_scaled <- qnorm(u)
  
  # integrate over quantiles
  if(n < 200){
    uw <- gauss.quad.prob(1)
  } else {
    uw <- gauss.quad.prob(m)
  }
  
  colnames(pair_scaled) <- c("x", "y")
  pair_scaled <- as.data.frame(pair_scaled)
  qrf_x_to_y <-  quantregForest(x=as.matrix(pair_scaled[,1]), y=as.matrix(pair_scaled[,2]), nodesize=10,sampsize=50)
  qrf_y_to_x <- quantregForest(x=as.matrix(pair_scaled[,2]), y=as.matrix(pair_scaled[,1]), nodesize=10,sampsize=50)
  
  cls <- sapply(uw$nodes, function(uu) {
    
    pred <- cbind(predict(qrf_y_to_x, as.matrix(pair_scaled[,2]), what=uu),
                  predict(qrf_x_to_y,  as.matrix(pair_scaled[,1]), what=uu))
    
    marg_q <- sapply(1:2, function(i) quantile(pair_scaled[,i], uu))
    # code lengths
    cl_marginal <- sapply(1:2, function(i)
      quantileScoring(pair_scaled[, i], marg_q[i], uu))
    cl_conditional <- sapply(1:2, function(i)
      quantileScoring(pair_scaled[, i], pred[, i], uu))
    
    c(cl_marginal, cl_conditional)
  })
  
  sel <- !apply(is.na(cls), 2, any)
  uw$weights <- uw$weights[sel] / sum(uw$weights[sel])
  cls <- apply(na.omit(t(cls)) * uw$weights, 2, sum)
  
  dx_to_y <- (cls[1] + cls[4])/sum(cls[1:2])
  dy_to_x <- (cls[2] + cls[3])/sum(cls[1:2])
  
  cd <- ifelse(dy_to_x > dx_to_y, 1, 0)
  epsilon <-  (-dx_to_y + dy_to_x )
  
  return(list(cd = cd, epsilon = epsilon))
}



## Quantile Causal discovery method - Neural net with quantile loss (simoultaneous estimation)
QCD_qnn <- function(pair, m=3) {
  # to get reproducible jittering results for discreet data
  # set.seed(0)
  
  n <- nrow(pair)
  # Recover pseudo-observations and estimate the copula non-parametrically
  u <- apply(pair, 2, function(x) rank(x, ties.method = "random")/(n + 1))
  # deal with discrete data
  pair_scaled <- qnorm(u)
  
  # integrate over quantiles
  if(n < 200){
    uw <- gauss.quad.prob(1)
  } else {
    uw <- gauss.quad.prob(m)
  }
  
  if(length(uw$nodes) == 1){
    qnn_x_to_y <- qrnn2.fit(x=as.matrix(pair_scaled[,1]), y=as.matrix(pair_scaled[,2]), tau=uw$nodes,
                            n.hidden=10, n.hidden2=5, n.trials=1,
                            iter.max=1000)
    qnn_y_to_x <- qrnn2.fit(x=as.matrix(pair_scaled[,2]), y=as.matrix(pair_scaled[,1]), tau=uw$nodes,
                            n.hidden=5, n.hidden2=10, n.trials=1,
                            iter.max=1000)
    pred.x_to_y <- qrnn2.predict(as.matrix(pair_scaled[,1]), qnn_x_to_y)
    pred.y_to_x <- qrnn2.predict(as.matrix(pair_scaled[,2]), qnn_y_to_x)
  }
  else{
    qnn_x_to_y <- mcqrnn.fit(x=as.matrix(pair_scaled[,1]), y=as.matrix(pair_scaled[,2]), tau=uw$nodes,
                             n.hidden=10, n.hidden2=5, n.trials=1,
                             iter.max=1000)
    qnn_y_to_x <- mcqrnn.fit(x=as.matrix(pair_scaled[,2]), y=as.matrix(pair_scaled[,1]), tau=uw$nodes,
                             n.hidden=5, n.hidden2=10, n.trials=1,
                             iter.max=1000)
    pred.x_to_y <- mcqrnn.predict(as.matrix(pair_scaled[,1]), qnn_x_to_y)
    pred.y_to_x <- mcqrnn.predict(as.matrix(pair_scaled[,2]), qnn_y_to_x)
  }
  
  cls <- sapply(uw$nodes, function(uu) {
    
    pred <- cbind(pred.y_to_x[ ,which(uw$nodes==uu)],
                  pred.x_to_y[ ,which(uw$nodes==uu)])
    
    marg_q <- sapply(1:2, function(i) quantile(pair_scaled[,i], uu))
    # code lengths
    cl_marginal <- sapply(1:2, function(i)
      quantileScoring(pair_scaled[, i], marg_q[i], uu))
    cl_conditional <- sapply(1:2, function(i)
      quantileScoring(pair_scaled[, i], pred[, i], uu))
    
    c(cl_marginal, cl_conditional)
  })
  
  sel <- !apply(is.na(cls), 2, any)
  uw$weights <- uw$weights[sel] / sum(uw$weights[sel])
  cls <- apply(na.omit(t(cls)) * uw$weights, 2, sum)
  
  dx_to_y <- (cls[1] + cls[4])/sum(cls[1:2])
  dy_to_x <- (cls[2] + cls[3])/sum(cls[1:2])
  
  cd <- ifelse(dy_to_x > dx_to_y, 1, 0)
  epsilon <-  (-dx_to_y + dy_to_x )
  
  return(list(cd = cd, epsilon = epsilon))
}


# wrapper used for the real data pairs
QCD_wrap <- function(X, Y, type="QCCD", m=7){
  method_to_run <- QCCD
  
  res = method_to_run(cbind(X,Y), m)
  if(!is.na(res$cd)) {
    cd = ifelse(res$cd == 1, "->", "<-")
  } else{
    cd = "--"
  }
  list(cd = cd, eps = res$eps)
}


