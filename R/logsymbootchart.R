
#rm(list = ls())

#====================================================
## Packages

#if("ssym" %in% rownames(installed.packages())==FALSE){install.packages("ssym");library(ssym)} else library(ssym)
#if("gbs" %in% rownames(installed.packages())==FALSE){install.packages("gbs");library(gbs)} else library(gbs)
#if("GoFKernel" %in% rownames(installed.packages())==FALSE){install.packages("GoFKernel");library(GoFKernel)} else library(GoFKernel)
#if("cubature" %in% rownames(installed.packages())==FALSE){install.packages("scubature");library(cubature)} else library(cubature)

#====================================================
## Maximum likelihood estimator - log-symmetric family

mle.logsym1 <- function(y, family, linf, lsup, by){

  if(family == "BS"){
    par <- as.numeric(mlebs(y)[c(2,3)])
    phi <- 4 # phi - fixed
    xi <- par[1] # xi - alpha
    eta <- par[2] # eta - beta
    logLik <- sum(dgbs(y, xi, eta, log = T))
    AIC <- -2*logLik+2*2
    BIC <- -2*logLik+log(length(y))*2
  }

  else{ if(family == "Normal"){xi <- 0} else {xi <- seq(linf, lsup, by)}

    n <- length(xi)
    eta <- numeric(n)
    phi <- numeric(n)
    logLik <- numeric(n)
    AIC <- BIC <- numeric(n)

    for(i in 1:n){
      fit <- ssym.l(log(y) ~ rep(1, length(y)) ,family = family, xi = xi[i])
      eta[i] <- exp(fit$theta.mu)
      phi[i] <- exp(fit$theta.phi)
      logLik[i] <- attr(logLik(fit),"log")
      AIC[i] <- -2*logLik[i]+2*2
      BIC[i] <- -2*logLik[i]+log(length(y))*2
    }

    imax <- order(logLik)[n]
    eta <- eta[imax]
    phi <- phi[imax]
    xi <- xi[imax]
    AIC <- AIC[imax]
    BIC <- BIC[imax]
  }

  return(list(eta = eta, phi = phi, xi = xi, AIC = AIC, BIC = BIC, family = family))
}

#====================================================
## Choosing the best log-symmetric model for the data (AIC)

best.logsym <- function(y){

BS <- mle.logsym1(y, family="BS")
NO <- mle.logsym1(y, family="Normal", linf=0, lsup=0)
ST <- mle.logsym1(y, family="Student", linf=2, lsup=20,by=1)
PO <- mle.logsym1(y, family="Powerexp", linf=-.5, lsup=.5, by=0.1)
SL <- mle.logsym1(y, family="Slash", linf=2, lsup=20 , by=0.1)
HY <- mle.logsym1(y, family="Hyperbolic", linf=1, lsup=20,by=.1)

family <- c("BS","Normal","Student","Powerexp","Slash","Hyperbolic")
eta      <- c(BS$eta,NO$eta,ST$eta,PO$eta,SL$eta,HY$eta)
phi      <- c(BS$phi,NO$phi,ST$phi,PO$phi,SL$phi,HY$phi)
xi       <- c(BS$xi,NO$xi,ST$xi,PO$xi,SL$xi,HY$xi)
AIC      <- c(BS$AIC,NO$AIC,ST$AIC,PO$AIC,SL$AIC,HY$AIC)
BIC      <- c(BS$BIC,NO$BIC,ST$BIC,PO$BIC,SL$BIC,HY$BIC)

ind <- which(min(AIC)==AIC)

return(list(eta = eta[ind], phi = phi[ind], xi = xi[ind], AIC = AIC[ind], BIC = BIC[ind], family = family[ind]))
}




#====================================================
## functions for each family member


quantile.function <- function(q, xi, family){

if(family=="Normal"){
  xi <- 0
  v <- function(z) rep(1,length(z))
  vp <- function(z) rep(0,length(z))
  dg <- 1
  fg <- 3
  deviance.mu <- function(z) z^2
  deviance.phi <- function(z) abs(z^2-1-log(z^2))
  cdf <- function(z) pnorm(z)
  pdf <- function(z) dnorm(z)
  qf  <- function(z) qnorm(z)
  xix <- 1
}
if(family=="Student"){
  if(xi[1]<=0) stop("the extra parameter must be positive!!",call.=FALSE)
  nu <- xi[1]
  v <- function(z) (nu + 1)/(nu + z^2)
  vp <- function(z) -2*z*(nu + 1)/((nu + z^2)^2)
  dg <- (nu + 1)/(nu + 3)
  fg <- 3*(nu + 1)/(nu + 3)
  deviance.mu <- function(z) abs((nu+1)*log(1 + z^2/nu))
  deviance.phi <- function(z) abs((nu+1)*log((nu + z^2)/(nu + 1)) -log(z^2))
  cdf <- function(z) pt(z,nu)
  pdf <- function(z) dt(z,nu)
  qf  <- function(z) qt(z,nu)
  xix <- nu/(nu-2)
  if(nu<=2) xix <- as.null(xix)
}


  if (family == "BS") {
    if (xi[1] <= 0)
      stop("the extra parameter must be positive!!", call. = FALSE)
    alpha <- xi[1]
    v <- function(z) 4 * sinh(z) * cosh(z)/(alpha^2 * z) -
      tanh(z)/z
    vp <- function(z) ((cosh(z) * z - sinh(z))/z^2) * (4 *
                                                         cosh(z)/alpha^2 - 1/cosh(z)) + (sinh(z)/z) * (4 *
                                                                                                         sinh(z)/alpha^2 + sinh(z)/(cosh(z)^2))
    dg <- 2 + 4/(alpha^2) - (sqrt(2 * pi)/alpha) * (1 -
                                                      2 * (pnorm(sqrt(2)/alpha, mean = 0, sd = sqrt(2)/2) -
                                                             0.5)) * exp(2/(alpha^2))
    dshn <- function(z) 2 * cosh(z) * exp(-(2/alpha^2) *
                                            (sinh(z))^2)/(alpha * sqrt(2 * pi))
    fgf <- function(z) dshn(z) * (4 * sinh(z) * cosh(z)/(alpha^2) -
                                    tanh(z))^2 * z^2
    fg <- 2 * integrate(fgf, 0, 60)$value
    deviance.mu <- function(z) {
      if (alpha <= 2)
        abs(4 * (sinh(z))^2/alpha^2 - log((cosh(z))^2))
      else {
        2 * log(dshn(acosh(alpha/2))/dshn(z))
      }
    }
    tau <- uniroot(function(x) 4 * x * sinh(x) * cosh(x)/(alpha^2) -
                     tanh(x) * x - 1, lower = 0, upper = 50)$root
    deviance.phi <- function(z) {
      a <- 2 * log(cosh(tau) * exp(-(2/alpha^2) * (sinh(tau))^2)) +
        log(tau^2)
      b <- 2 * log(cosh(z) * exp(-(2/alpha^2) * (sinh(z))^2)) +
        log(z^2)
      ifelse(a < b, 0, a - b)
    }
    cdf <- function(z) pnorm(2 * sinh(z)/alpha)
    pdf <- function(z) (2 * cosh(sqrt(z^2)) * exp(-2 * sinh(sqrt(z^2)) *
                                                    sinh(sqrt(z^2))/alpha^2)/(sqrt(2 * pi) * alpha))
    qf <- function(z) asinh(alpha * qnorm(z)/2)
    dzg <- function(z) sinh(sqrt(z)) * exp(-(2/alpha^2) * sinh(sqrt(z))^2) * (alpha^2 - 4 * cosh(sqrt(z))^2)/(2 * alpha^2 * sqrt(z))
    fgf <- function(z) dshn(z) * z^2
    xix <- 2 * integrate(fgf, 0, 20)$value
  }

if (family == "Powerexp") {
  if (xi[1] <= -1 | xi[1] >= 1)
    stop("the extra parameter must be within the interval (-1, 1)!!",
         call. = FALSE)
  kk <- xi[1]
  v <- function(z) abs(z)^(-(2 * kk/(1 + kk)))/(1 + kk)
  vp <- function(z) -2 * kk * ifelse(z >= 0, 1, -1) *
    abs(z)^(-((3 * kk + 1)/(1 + kk)))/(1 + kk)^2
  dg <- 2^(1 - kk) * gamma((3 - kk)/2)/((1 + kk)^2 * gamma((1 +
                                                              kk)/2))
  fg <- (kk + 3)/(kk + 1)
  deviance.mu <- function(z) abs((abs(z))^(2/(1 + kk)))
  deviance.phi <- function(z) abs((abs(z))^(2/(1 + kk)) -
                                    (1 + kk) - log(z^2/((1 + kk)^(1 + kk))))
  pp <- 2/(kk + 1)
  sigmap <- (1 + kk)^((kk + 1)/2)
  cdf <- function(z) pnormp(z, mu = 0, sigmap = sigmap,
                            p = pp)
  pdf <- function(z) dnormp(z, mu = 0, sigmap = sigmap,
                            p = pp)
  qf <- function(z) qnormp(z , mu=0, sigmap= sigmap,
                           p= pp)
  dzg <- function(z) -(1/2) * exp(-(1/2) * z^(1/(1 + kk))) * (1/(1 + kk)) * z^(1/(1 + kk) - 1)
  xix <- 2^(1 + kk) * gamma(3 * (1 + kk)/2)/gamma((1 +
                                                     kk)/2)
}

if (family == "Hyperbolic") {
  if (xi[1] <= 0)
    stop("the extra parameter must be positive!!", call. = FALSE)
  nu <- xi[1]
  v <- function(z) nu/sqrt(1 + z^2)
  vp <- function(z) -nu * z/((1 + z^2)^(3/2))
  dh <- function(z) exp(-nu * sqrt(1 + z^2))/(2 * besselK(nu,
                                                          1))
  dgd <- function(z) dh(z) * (nu * z/sqrt(1 + z^2))^2
  dg <- 2 * cubintegrate(dgd, 0, Inf, method = "pcubature")$integral
  fgf <- function(z) dh(z) * (nu * z^2/sqrt(1 + z^2))^2
  fg <- 2 * cubintegrate(fgf, 0, Inf, method = "pcubature")$integral
  deviance.mu <- function(z) abs(2 * nu * (sqrt(1 + z^2) -
                                             1))
  tau <- sqrt((1 + sqrt(1 + 4 * nu^2))/(2 * nu^2))
  deviance.phi <- function(z) abs(2 * nu * (sqrt(1 + z^2) -
                                              sqrt(1 + tau^2)) - log(z^2/tau^2))
  cdf <- function(z){
    temporal <- matrix(z, length(z), 1)
    phyperbolic <- function(z) cubintegrate(dh, method = "pcubature", lower = -Inf, upper = z)$integral
    result <- apply(X = temporal, MARGIN = 1, FUN = phyperbolic)
    return(result)
  }
  pdf <- function(z) dh(z)
  #qf <- function(z) qhyperb(z, Theta = c(pi = 0, zeta =  nu, delta = 1, mu = 0))
  qf <- function(z) {
    temporal <- matrix(z, length(z), 1)
    qhyperbolic <- function(z) uniroot(function(x) cdf(x) - z, interval = c(-100000,100000))$root
    result <- apply(X = temporal, MARGIN = 1, FUN = qhyperbolic)
    return(result)
  }
  dzg <- function(z) - nu * (2 * sqrt(1 + z)) * exp(-nu * sqrt(1 + z))
  fgf <- function(z) dh(z) * z^2
  xix <- 2 * cubintegrate(fgf, 0, Inf, method = "pcubature")$integral
}
if (family == "Slash") {
  if (xi[1] <= 0)
    stop("the extra parameter must be positive!!", call. = FALSE)
  nu <- xi[1]
  G <- function(a, x) gamma(a) * pgamma(1, shape = a,
                                        scale = 1/x)/(x^a)
  v <- function(z) G(nu + 3/2, z^2/2)/G(nu + 1/2, z^2/2)
  vp <- function(z) grad(v, z)
  ds <- function(z) nu * G(nu + 1/2, z^2/2)/sqrt(2 * pi)
  gdg <- function(z) ds(z) * (v(z))^2 * z^2
  dg <- 2 * cubintegrate(gdg, 0, Inf, method = "pcubature")$integral
  gfg <- function(z) ds(z) * (v(z))^2 * z^4
  fg <- 2 * cubintegrate(gfg, 0, Inf, method = "pcubature")$integral
  deviance.mu <- function(z) abs(2 * log(2/(2 * nu + 1)) -
                                   2 * log(G(nu + 1/2, z^2/2)))
  tau <- uniroot(function(x) v(x) * x^2 - 1, lower = 1e-04,
                 upper = 1000)$root
  deviance.phi <- function(z) {
    a <- 2 * log(G(nu + 1/2, tau^2/2)) + log(tau^2)
    b <- 2 * log(G(nu + 1/2, z^2/2)) + log(z^2)
    ifelse(a < b, 0, a - b)
  }
  cdf <- function(z) {
    temporal <- matrix(z, length(z), 1)
    pslash <-  function(z) cubintegrate(ds, method = "pcubature", lower = -Inf, upper = z)$integral
    result <- apply(X = temporal, MARGIN = 1, FUN = pslash)
    return(result)
  }
  pdf <- function(z) ds(z)
  qf <- function(z) {
    temporal <- matrix(z, length(z), 1)
    qslash <- function(z) uniroot(function(x) cdf(x) - z, interval = c(-1000000,1000000))$root
    result <- apply(X = temporal, MARGIN = 1, FUN = qslash)
    return(result)
  }
  dzg <- function(z) -2^(-nu - 1/2) * exp(-z/2) * z^(nu - 1/2)
  gfg <- function(z) {
    ds(z) * z^2
  }
  xix <- 2 * cubintegrate(gfg, 0, 60, method = "pcubature")$integral
}

return(qf(q))

}

#====================================================
## Random number generation of the log-symmetric model

sample.logsym <- function(n, par, family, xi){
  if(!is.numeric(par)){stop("parameters must be numeric")}
npas <- numeric()
  switch(family,
         "BS"         = {eta <- par[1]; phi <- par[2];    npas  <- eta * exp(sqrt(phi) * rvgs(n, family = "Sinh-normal", xi = xi))},
         "Normal"     = {eta <- par[1]; phi <- par[2];    npas  <- eta * exp(sqrt(phi) * rvgs(n, family = "Normal", xi = xi))},
         "Student"    = {eta <- par[1]; phi <- par[2];    npas  <- eta * exp(sqrt(phi) * rvgs(n, family = "Student", xi = xi))},
         "Powerexp"   = {eta <- par[1]; phi <- par[2];    npas  <- eta * exp(sqrt(phi) * rvgs(n, family = "Powerexp", xi = xi)) },
         "Hyperbolic" = {eta <- par[1]; phi <- par[2];    npas  <- eta * exp(sqrt(phi) * rvgs(n, family = "Hyperbolic", xi = xi))},
         "Slash"      = {eta <- par[1]; phi <- par[2];    npas  <- eta * exp(sqrt(phi) * rvgs(n, family = "Slash", xi = xi))}
          )
return(npas)
}

#====================================================
## Estimate the parameters of the log-symmetric model

mle.logsym2 <- function(x, family, xi){
  if(!is.numeric(x)){stop("data must be numeric")}
  switch(family,
         "BS"  = {ml.est <- ssym.l(log(x) ~ rep(1, length(x)) ,family = "Sinh-normal", xi = xi); a <- exp(ml.est$theta.mu); b <- exp(ml.est$theta.phi)},
         "Normal"  = {ml.est <- ssym.l(log(x) ~ rep(1, length(x)) ,family = "Normal", xi = xi); a <- exp(ml.est$theta.mu); b <- exp(ml.est$theta.phi)},
         "Student"   = {ml.est <- ssym.l(log(x) ~ rep(1, length(x)) ,family = "Student", xi = xi); a <- exp(ml.est$theta.mu); b <- exp(ml.est$theta.phi)},
         "Powerexp"  = {ml.est <- ssym.l(log(x) ~ rep(1, length(x)) ,family = "Powerexp", xi = xi); a <- exp(ml.est$theta.mu); b <- exp(ml.est$theta.phi)},
         "Hyperbolic"   = {ml.est <- ssym.l(log(x) ~ rep(1, length(x)) ,family = "Hyperbolic", xi = xi); a <- exp(ml.est$theta.mu); b <- exp(ml.est$theta.phi)},
         "Slash" = {ml.est <- ssym.l(log(x) ~ rep(1, length(x)) ,family = "Slash", xi = xi); a <- exp(ml.est$theta.mu); b <- exp(ml.est$theta.phi)}
         )
  par <- c(a, b)
  return(par)
}

#====================================================
## quantile functions

quantile.logsym <- function(p, par, family, xi){
  if(!is.numeric(p) | p>=1){stop("p must be 0 < p < 1")}
  switch (family,
          "BS"     = {eta <- par[1]; phi <- par[2]; w <- eta * exp(sqrt(phi) * quantile.function(q = p, xi = xi, family = family))},
          "Normal"  = {eta <- par[1]; phi <- par[2]; w <- eta * exp(sqrt(phi) * quantile.function(q = p, xi = xi, family = family))},
          "Student"   = {eta <- par[1]; phi <- par[2]; w <- eta * exp(sqrt(phi) * quantile.function(q = p, xi = xi, family = family))},
          "Powerexp"  = {eta <- par[1]; phi <- par[2]; w <- eta * exp(sqrt(phi) * quantile.function(q = p, xi = xi, family = family))},
          "Hyperbolic"   = {eta <- par[1]; phi <- par[2]; w <- eta * exp(sqrt(phi) * quantile.function(q = p, xi = xi, family = family))},
          "Slash" = {eta <- par[1]; phi <- par[2]; w <- eta * exp(sqrt(phi) * quantile.function(q = p, xi = xi, family = family)   )}
  )
  return(w)
}


#====================================================
## plot the control chart

## Fase I
bchart1 <- function(x, lcl, ucl, cl, n){
  par(mar = c(2.2, 2.2, .5, .5), bty = "n", mgp = c(1.2, .1, 0), tck = .01, family = "sans", cex = .95, col.axis = "grey30")

  plot(1:length(x), x, main="", xlab="Subgroup", ylab="Percentiles",
       pch=ifelse(x >= lcl & x <=ucl, 16, 8), type="b", lty=2, axes=F, xpd=T, ylim=c(lcl-3*sd(x),ucl+3*sd(x)))
  #axis(1, 1:length(x), labels=1:length(x))
  axis(1, 1:n, labels=1:n)
  axis(2, round(seq(lcl-3*sd(x),ucl+3*sd(x),.5),3))
  abline(h=cl, lty=2)
  abline(h=c(lcl, ucl), lty=1, lwd=2)
  par(adj = 0)
  text(1, c(lcl, cl, ucl)+.1, paste(c("LCL=","CL=", "UCL="),
                                    round(c(lcl, cl, ucl),3)), cex=.8, font=2)
}

## Fase II
bchart2 <- function(x, lcl, ucl, cl, n){
  par(mar = c(2.2, 2.2, .5, .5), bty = "n", mgp = c(1.2, .1, 0), tck = .01, family = "sans", cex = .95, col.axis = "grey30")

  plot(1:length(x), x, main="", xlab="Subgroup", ylab="Percentiles",
       pch=ifelse(x >= lcl & x <=ucl, 16, 8), type="b", lty=2, axes=F, xpd=T, ylim=c(lcl-3*sd(x),ucl+3*sd(x)))
  #axis(1, 1:length(x), labels=1:length(x))
  axis(1, 1:n, labels=(n+1):(n+n))
  axis(2, round(seq(lcl-3*sd(x),ucl+3*sd(x),.5),3))
  abline(h=cl, lty=2)
  abline(h=c(lcl, ucl), lty=1, lwd=2)
  par(adj = 0)
  text(1, c(lcl, cl, ucl)+.1, paste(c("LCL=","CL=", "UCL="),
                                    round(c(lcl, cl, ucl),3)), cex=.8, font=2)
}

#====================================================
## function to construct the control limits of the
## log-symmetric bootstrap control chart (paper)

chart.logsym.paper <- function(dataic, dataofc, paric, family, p, n, m, gamma, boot)
{

  dataicm    <- matrix(dataic,n,m,byrow = TRUE )
  dataofcm   <- matrix(dataofc,n,m,byrow = TRUE )

  etahat    <- paric[1]
  phihat    <- paric[2]
  xihat     <- paric[3]
  familyhat <- family

  W <- numeric()

  for(i in 1:boot){
    npas <- sample.logsym(n = m, par = c(etahat,phihat), family = familyhat, xi = xihat)
    par <- mle.logsym2(x = npas, family = familyhat, xi = xihat )
    W[i] <- quantile.logsym(p, par, family = familyhat, xi =xihat )
  }

  #Control limits
  lcl <- quantile(W, gamma/2)
  cl  <- quantile.logsym(p, c(etahat,phihat), family = familyhat, xi = xihat )
  ucl <- quantile(W, 1-gamma/2)

  Wt <- double()
  for(i in 1:nrow(dataofcm)){
    parp <- mle.logsym2(dataofcm[i,], family = familyhat, xi = xihat)
    Wt[i] <- quantile.logsym(p = p, parp, family = familyhat, xi = xihat)
  }

  bchart2(x=Wt, lcl=lcl, ucl=ucl, cl=cl, n=n)

return(list(family = familyhat, eta = etahat, phi = phihat, xi = xihat, lcl = lcl, cl = cl, ucl = ucl))
}



#====================================================
## function to construct the log-symmetric bootstrap
## control chart - Fase I


chart.logsym1 <- function(data, p, n, m, gamma, boot)
{

  datam <- matrix (data,n,m,byrow = TRUE )

  resultsmle <- best.logsym(y = data)
  etahat     <- resultsmle$eta
  phihat     <- resultsmle$phi
  xihat      <- resultsmle$xi
  familyhat  <- resultsmle$family

  W <- numeric()

  for(i in 1:boot){
    npas <- sample.logsym(n = m, par = c(etahat,phihat), family = familyhat, xi = xihat)
    par <- mle.logsym2(x = npas, family = familyhat, xi = xihat )
    W[i] <- quantile.logsym(p, par, family = familyhat, xi =xihat )
  }

  #Control limits
  lcl <- quantile(W, gamma/2)
  cl  <- quantile.logsym(p, c(etahat,phihat), family = familyhat, xi = xihat )
  ucl <- quantile(W, 1-gamma/2)


  Wt <- double()
  for(i in 1:nrow(datam)){
    parp <- mle.logsym2(datam[i,], family = familyhat, xi = xihat)
    Wt[i] <- quantile.logsym(p = p, parp, family = familyhat, xi = xihat)
  }

  bchart1(x=Wt, lcl=lcl, ucl=ucl, cl=cl, n=n)

  return(list(family = familyhat, eta = etahat, phi = phihat, xi = xihat, lcl = lcl, cl = cl, ucl = ucl))
}


#====================================================
## function to construct the log-symmetric bootstrap
## control chart - Fase II

chart.logsym2 <- function(data, p, n, m, family, lcl, ucl, cl,  linf, lsup, by)
{

datam <- matrix (data,n,m,byrow = TRUE )

Wt <- numeric()
for(i in 1:nrow(datam)){
  #par <- mle.logsym(data[i,], family = family, xi = xi)
  res = mle.logsym1(y = datam[i,], family =family , linf=linf, lsup=lsup, by=by)
  par = c(res$eta,res$phi,res$xi)
  Wt[i] <- quantile.logsym(p = p, par[1:2], family = family, xi = par[3])
}

bchart2(x=Wt, lcl=lcl, ucl=ucl, cl=cl, n=n)

}

#====================================================









