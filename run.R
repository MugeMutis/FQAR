# load the auxilary functions
    source("auxiliary_functions.R")

    river = read.delim("mary_river.txt",header = FALSE)
    river = as.matrix(river)

    Y_true = river[dim(river)[1],]
    data = river[-dim(river)[1],]
   
    # obtain the functional response and one-lagged functional predictor variable
    order = 1
    xy = var_fun(data=data[(1:dim(data)[1]),], order=order)
    x = xy$x
    y = xy$y
    
    # run the FPCA for functional variables
    rangeval = c(1,24)
    refinement = 24
    nbasis = 10
    ncomp = 5
    
    fPCA_Y = getPCA(data=y, nbasis, ncomp, rangeval, cen=FALSE)
    fsco_Y = fPCA_Y$PCAscore
    fcomp_Y = fPCA_Y$PCAcoef
    fmean_Y = fPCA_Y$meanScore
    evaly = fPCA_Y$evalbase
    fPCA_Y2 = getPCA(data=y, nbasis, ncomp, rangeval, cen=TRUE)
    fsco_Y2 = fPCA_Y2$PCAscore
    
    fsco_X = list()
    evalx = list()
    fcomp_X = list()
    fnp = length(x)
    for(fij in 1:fnp){
      fPCA_X = getPCA(data = x[[fij]], nbasis, ncomp, rangeval, cen = TRUE)
      fsco_X[[fij]] = fPCA_X$PCAscore
      evalx[[fij]] = fPCA_X$evalbase
      fcomp_X[[fij]] = fPCA_X$PCAcoef
    }  
    
    
    # Obtain forecasts with FQAR(0.5)
    Bhat_qr_50 = est_fun_qr(sco_Y = fsco_Y, sco_X = do.call(cbind, fsco_X), tau = 0.50)
    fB0_qr_50 = Bhat_qr_50$B0 %*% t(fcomp_Y$coefs) %*% t(evaly)
    fmodel_qr_50 = pred_fun(comp_Y = fcomp_Y, sco_X = c(t(fsco_Y2[((dim(fsco_Y2)[1]-order+1):(dim(fsco_Y2)[1])),])),
                            Bhat = Bhat_qr_50$fBhat)
    fYpred_qr_50 = eval.fd(fmodel_qr_50, seq(rangeval[1], rangeval[2], length.out = dim(data)[2])) + t(fB0_qr_50)
    
    
    # Obtain forecasts with FPCR
    Bhat_ml = est_fun_ml(sco_Y = fsco_Y, sco_X = do.call(cbind, fsco_X))
    fB0_ml = Bhat_ml$B0 %*% t(fcomp_Y$coefs) %*% t(evaly)
    fmodel_ml = pred_fun(comp_Y = fcomp_Y, sco_X = c(t(fsco_Y2[((dim(fsco_Y2)[1]-order+1):(dim(fsco_Y2)[1])),])),
                         Bhat = Bhat_ml$fBhat)
    fYpred_ml = eval.fd(fmodel_ml, seq(rangeval[1], rangeval[2], length.out = dim(data)[2])) + t(fB0_ml)
    
    
    colnames(data) = 1:dim(data)[2]
    rownames(data) = 1:dim(data)[1]
    
    # Obtain forecasts with FPLS
    fYpred_pls = fplsr(data=fts(x=1:dim(data)[2],y=t(data)), order=ncomp, type="simpls")$Ypred$y
    
    # Obtain forecasts with ARIMA
    fYpred_arima = t(forecast(object=ftsm(fts(x=1:dim(data)[2],y=t(data)), order=ncomp, method="classical", mean=TRUE), 
                              h=order, method="arima")$mean$y)
    

    # mape values
    mape_qr_50 = mean(abs((fYpred_qr_50-Y_true)/Y_true))  # 0.004016443
    mape_ml = mean(abs((fYpred_ml-Y_true)/Y_true))        # 0.01792263
    mape_pls = mean(abs((fYpred_pls-Y_true)/Y_true))      # 0.01681547
    mape_arima = mean(abs((fYpred_arima-Y_true)/Y_true))  # 0.01428254
    
    # rmspe values
    rmspe_qr_50 = sqrt(mean(((fYpred_qr_50-Y_true) / Y_true)^2))  # 0.004581631
    rmspe_ml = sqrt(mean(((fYpred_ml-Y_true) / Y_true)^2))        # 0.02135923
    rmspe_pls = sqrt(mean(((fYpred_pls-Y_true) / Y_true)^2))      # 0.02008463
    rmspe_arima = sqrt(mean(((fYpred_arima-Y_true) / Y_true)^2))  # 0.01577279
    
    # rmespe values
    rmespe_qr_50 = sqrt(median(((fYpred_qr_50-Y_true) / Y_true)^2))  # 0.003796296
    rmespe_ml = sqrt(median(((fYpred_ml-Y_true) / Y_true)^2))        # 0.01874972
    rmespe_pls = sqrt(median(((fYpred_pls-Y_true) / Y_true)^2))      # 0.0145404
    rmespe_arima = sqrt(median(((fYpred_arima-Y_true) / Y_true)^2))  # 0.01344145
    
    # re values
    re_qr_50 = mean((fYpred_qr_50-Y_true)/Y_true)*100  # 0.3651902
    re_ml = mean((fYpred_ml-Y_true)/Y_true)*100        # 1.792263
    re_pls = mean((fYpred_pls-Y_true)/Y_true)*100      # 1.681547
    re_arima = mean((fYpred_arima-Y_true)/Y_true)*100  # 1.415055
    
    # nse values
    nse_qr_50 = nse(Y_true,fYpred_qr_50)  # -16.56743
    nse_ml = nse(Y_true,fYpred_ml)        # -380.8936
    nse_pls = nse(Y_true,fYpred_pls)      # -336.6714
    nse_arima = nse(Y_true,fYpred_arima)  # -207.2162
    
    # mia values
    mia_qr_50 = mia(Y_true,fYpred_qr_50)  # 0.1636696
    mia_ml = mia(Y_true,fYpred_ml)        # 0.04833109
    mia_pls = mia(Y_true,fYpred_pls)      # 0.05135226
    mia_arima = mia(Y_true,fYpred_arima)  # 0.05166005
    
    
    # interval score value for 95% nominal confidence level
    Bhat_qr_025 = est_fun_qr(sco_Y = fsco_Y, sco_X = do.call(cbind, fsco_X), tau = 0.025)
    fB0_qr_025 = Bhat_qr_025$B0 %*% t(fcomp_Y$coefs) %*% t(evaly)
    fmodel_qr_025 = pred_fun(comp_Y = fcomp_Y, sco_X = c(t(fsco_Y2[((dim(fsco_Y2)[1]-order+1):(dim(fsco_Y2)[1])),])),
                             Bhat = Bhat_qr_025$fBhat)
    fYpred_qr_025 = eval.fd(fmodel_qr_025, seq(rangeval[1], rangeval[2], length.out = dim(data)[2])) + t(fB0_qr_025)
    
    
    Bhat_qr_975 = est_fun_qr(sco_Y = fsco_Y, sco_X = do.call(cbind, fsco_X), tau = 0.975)
    fB0_qr_975 = Bhat_qr_975$B0 %*% t(fcomp_Y$coefs) %*% t(evaly)
    fmodel_qr_975 = pred_fun(comp_Y = fcomp_Y, sco_X = c(t(fsco_Y2[((dim(fsco_Y2)[1]-order+1):(dim(fsco_Y2)[1])),])),
                             Bhat = Bhat_qr_975$fBhat)
    fYpred_qr_975 = eval.fd(fmodel_qr_975, seq(rangeval[1], rangeval[2], length.out = dim(data)[2])) + t(fB0_qr_975)
    
    
    int_score = interval_score(Y_true, fYpred_qr_025, fYpred_qr_975, alpha=0.05)  # 0.2029567
