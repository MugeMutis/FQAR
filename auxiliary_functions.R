# Packagess
library(fda) 
library(quantreg)
library(ftsa)
library(MASS)

# Split the data into a functional response and 
# one-lagged predictor variable

var_fun = function(data,order){
  n = dim(data)[1]
  y = data[((order+1):n),]
  x=list()
  a=1
  b=order
  for(i in 1:order){
    x[[i]] = data[(a:(n-b)),]
    a = a+1
    b = b-1
  }
  return(list(x=x,y=y))
}

# Obtain the functional principal components and the 
# corresponding principal components scores
getPCA = function(data, nbasis, ncomp, rangeval, cen){
  n = dim(data)[1]
  p = dim(data)[2]
  dimnames(data)=list(as.character(1:n), as.character(1:p))
  grid_points = seq(rangeval[1], rangeval[2], length.out = p)
  bs_basis = create.bspline.basis(rangeval, nbasis = nbasis)
  evalbase = eval.basis(grid_points, bs_basis)
  fdobj = fdPar(bs_basis, int2Lfd(2), lambda=0)
  pcaobj = smooth.basisPar(grid_points, t(data), bs_basis, Lfdobj=NULL, lambda=0)$fd
  dpca = pca.fd(pcaobj, nharm = ncomp, fdobj, centerfns = cen)
  PCAscore = dpca$scores
  PCAcoef = dpca$harmonics
  mean_coef = dpca$meanfd
  return(list(PCAcoef = PCAcoef, PCAscore = PCAscore,
              meanScore = mean_coef, evalbase = evalbase))
}

# Estimate the regression parameter matrix (Quantile regression)
est_fun_qr = function(sco_Y, sco_X, tau){
  Bhat = t(rqs.fit(x = cbind(1,sco_X), y = sco_Y,
                   tau = tau))
  B0 = Bhat[1,]
  fc = dim(Bhat)[2]
  fBhat = Bhat[-1,]
  if(is.null(dim(fBhat)) == TRUE & fc >1){
    fBhat = t(as.matrix(fBhat))
  }else if(is.null(dim(fBhat)) == TRUE & fc == 1){
    fBhat = as.matrix(fBhat)
  }
  return(list(Bhat=Bhat,B0=B0,fBhat=fBhat))
}

# Estimate the regression parameter matrix (MLE)
est_fun_ml = function(sco_Y, sco_X){
  sco_X=cbind(1,sco_X)
  Bhat = ginv(t(sco_X) %*% sco_X) %*% t(sco_X) %*% sco_Y
  B0 = Bhat[1,]
  fBhat = Bhat[-1,]
  return(list(Bhat=Bhat,B0=B0,fBhat=fBhat))
}

# Obtain the forecast conditionally on the past observations
pred_fun = function(comp_Y, sco_X, Bhat){
  ncomp = dim(comp_Y$coefs)[2]
  nest = t(sco_X %*% Bhat)
  if(ncomp == 1){
    nh = nest[1] * comp_Y[1,]
  }else{
    nh = nest[1] * comp_Y[1,]
    for (j in 2:ncomp){
      nh = nh + nest[j] * comp_Y[j,]  
    }
  }
  return(nh)
}

# Compute the mia metric for the observed and forecasted values
mia = function(Y_true, Y_pred){
  a=abs(Y_true-mean(Y_true))
  b=abs(Y_pred-mean(Y_true))
  mia = 1-(sum(abs(Y_true-Y_pred))/sum(a+b))
  return(mia)
}

# Compute the nse metric for the observed and forecasted values
nse = function(Y_true, Y_pred){
  a=(Y_true-mean(Y_true))
  nse = 1-(sum((Y_true-Y_pred)^2)/sum(a^2))
  return(nse)
}

# Compute the interval score for the prediction interval
interval_score <- function(holdout, lb, ub, alpha){
  lb_ind = ifelse(holdout < lb, 1, 0)
  ub_ind = ifelse(holdout > ub, 1, 0)
  score = (ub - lb) + 2/alpha * ((lb - holdout) * lb_ind + (holdout - ub) * ub_ind)
  return(mean(score))
}


