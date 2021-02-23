########################################################################
## REPLICATION FILES: Calonico, Cattaneo and Farrell (2018)
## nprobust: Nonparametric Kernel-Based Estimation and Robust Bias-Corrected Inference
## Last update: 17-April-2018
########################################################################
## This script contains only auxiliary functions. Please run CCF_nprobust_gentables.R

scale <- 0
n     <- 500
sim   <- 5000
kernel<-"epa"
vce   <- "nn"
level <- 5
qz    <- qnorm((100-level/2)/100)
neval <- 5
x.lb  <- 0
x.ub  <- 1
evalx <- seq(x.lb,x.ub,length.out=neval)
check <- 21


g.0.fun = function(x) {sin(2*x-1) + 2*exp(-16*(x-0.5)^2)}
g.1.fun = function(x) {cos(2 * x - 1) * 2 - 2 * (exp(-16 * (x - 0.5)^2) * (16 * (2 * (x - 0.5))))}
g.2.fun = function(x) {-(sin(2 * x - 1) * 2 * 2 + 2 * (exp(-16 * (x - 0.5)^2) * (16 * 2) - exp(-16 * (x - 0.5)^2) * (16 * (2 * (x - 0.5))) * (16 *(2 * (x - 0.5)))))}
g.3.fun = function(x) {-(cos(2 * x - 1) * 2 * 2 * 2 - 2 * (exp(-16 * (x - 0.5)^2) * (16 * (2 * (x - 0.5))) * (16 * 2) + ((exp(-16 * (x - 0.5)^2) * (16 * 2) - exp(-16 * (x - 0.5)^2) * (16 * (2 * (x - 0.5))) * (16 * (2 * (x - 0.5)))) * (16 * (2 * (x - 0.5))) + exp(-16 * (x - 0.5)^2) * (16 * (2 * (x - 0.5))) * (16 * 2))))}
g.4.fun = function(x) { sin(2 * x - 1) * 2 * 2 * 2 * 2 + 2 * ((exp(-16 * (x - 0.5)^2) * (16 * 2) - exp(-16 * (x - 0.5)^2) * (16 * (2 * (x - 0.5))) *  (16 * (2 * (x - 0.5)))) * (16 * 2) + ((exp(-16 * (x - 0.5)^2) * (16 * 2) - exp(-16 * (x - 0.5)^2) * (16 * (2 * (x - 0.5))) * (16 * (2 * (x - 0.5)))) * (16 * 2) - (exp(-16 * (x - 0.5)^2) * (16 * (2 * (x - 0.5))) * (16 * 2) + ((exp(-16 * (x - 0.5)^2) * (16 * 2) - exp(-16 * (x - 0.5)^2) * (16 * (2 * (x - 0.5))) *  (16 * (2 * (x - 0.5)))) * (16 * (2 * (x - 0.5))) + exp(-16 *(x - 0.5)^2) * (16 * (2 * (x - 0.5))) * (16 * 2))) * (16 *  (2 * (x - 0.5))) + (exp(-16 * (x - 0.5)^2) * (16 * 2) - exp(-16 * (x - 0.5)^2) * (16 * (2 * (x - 0.5))) * (16 * (2 * (x - 0.5)))) * (16 * 2)))}

sigma   = 1
x.gen   = function(n) {runif(n, min=x.lb, max=x.ub)}  
fx      = function(x) {dunif(x, min=x.lb, max=x.ub)}
u.rnd   = function(n) {rnorm(n,0,sigma)}


Cbw.lp.fun <- function(p,v,kernel) {
  m1 = function(i,j,k)   integrate(function(x) x^i*x^j*k(x),0,Inf)$value
  m2 = function(i,j,k)   integrate(function(x) x^i*x^j*(k(x))^2,0,Inf)$value
  GAMMA = function(p,k) {out=matrix(NA,p+1,p+1); for (i in 0:p) for (j in 0:p) out[i+1,j+1]=m1(i,j,k); return(out)}
  NU    = function(p,k) {out=matrix(NA,p+1,1); for (i in 0:p) out[i+1,1]=m1(i,p+1,k); return(out)}
  PSI   = function(p,k) {out=matrix(NA,p+1,p+1); for (i in 0:p) for (j in 0:p) out[i+1,j+1]=m2(i,j,k); return(out)}
  B.lp  = function(p,k) {out=solve(GAMMA(p,k))%*%NU(p,k); out[1]}
  
  C1.fun = function(p0,v,K) {
    S.inv = solve(GAMMA(p0,K))
    C1 = (S.inv%*%NU(p0,K))[v+1]
    return(C1)
  }
  
  C2.fun = function(p0,v,K) {
    S.inv = solve(GAMMA(p0,K))
    C2 = (S.inv%*%PSI(p0,K)%*%S.inv)[v+1,v+1]
    return(C2)
  }
  
  k.fun = function(u){
    if (kernel=="epa") w = 0.75*(1-u^2)*(abs(u)<=1)
    if (kernel=="uni") w =          0.5*(abs(u)<=1)
    if (kernel=="tri") w =   (1-abs(u))*(abs(u)<=1)
    return(w)
  }
  
  C1.h = C1.fun(p0=p, v=v, K=k.fun)
  C2.h = C2.fun(p0=p, v=v, K=k.fun)
  
  out = list(C1.h=C1.h, C2.h = C2.h)
  return(out)
}

C.pob <- Cbw.lp.fun(p, deriv, kernel) 
C.B=as.numeric(C.pob[1])
C.V=as.numeric(C.pob[2])

h.mse.pob = m.v.pob = 0

for (i in 1:neval) {
  m.v.pob[i] = g.0.fun(evalx[i])
  if (deriv==1) m.v.pob[i] = g.1.fun(evalx[i])
  mp1.pob = g.2.fun(evalx[i])
  mp2.pob = g.3.fun(evalx[i])
  if (p==0) {
    mp1.pob = g.1.fun(evalx[i])
    mp2.pob = g.2.fun(evalx[i])
  }
  if (p==2) {
    mp1.pob = g.3.fun(evalx[i])
    mp2.pob = g.4.fun(evalx[i])
  }
  f0.pob    = fx(evalx[i])
  
  even <- (p-deriv)%%2==0
  
  V.pob = sigma^2/f0.pob
  B1.pob = mp1.pob/factorial(p+1)
  B2.pob = mp2.pob/factorial(p+2)
  
  if (even==FALSE) {
    h.mse.pob[i] = ( ((1+2*deriv)*C.V*V.pob) / (2*(p+1-deriv)*(C.B*B1.pob)^2*n) )^(1/(2*p+3))
  }
  else {
    bw.fun  <-function(H) {abs(H^(2*(p+1-deriv))*(C.B*B1.pob + H*C.B*B2.pob)^2 + C.V*V.pob/(n*H^(1+2*deriv)))}
    h.mse.pob[i] <- optimize(bw.fun, interval=c(.Machine$double.eps, 1))$minimum
  }
}

if  (even==TRUE ) { 
  h.ce.pob <- h.mse.pob*n^(-((p+2)/((2*p+5)*(p+3))))
} else{
  h.ce.pob <- h.mse.pob*n^(-((p)/((2*p+3)*(p+3))))
}

I.V = sigma^2*integrate(function(x) 1/fx(x),      lower=x.lb, upper=x.ub)$value 
I.B1 = integrate(function(x) g.2.fun(x)^2, lower=x.lb, upper=x.ub)$value/factorial(p+1)^2
I.B2 = integrate(function(x) g.3.fun(x)^2, lower=x.lb, upper=x.ub)$value/factorial(p+1)^2
I.B1.B2 = integrate(function(x) g.2.fun(x)*g.3.fun(x), lower=x.lb, upper=x.ub)$value/factorial(p+1)^2
if (p==0) {
  I.B1 = integrate(function(x) g.1.fun(x)^2, lower=x.lb, upper=x.ub)$value/factorial(p+1)^2
  I.B2 = integrate(function(x) g.2.fun(x)^2, lower=x.lb, upper=x.ub)$value/factorial(p+1)^2
  I.B1.B2 = integrate(function(x) g.1.fun(x)*g.2.fun(x), lower=x.lb, upper=x.ub)$value/factorial(p+1)^2
}
if (p==2) {
  I.B1 = integrate(function(x) g.3.fun(x)^2, lower=x.lb, upper=x.ub)$value/factorial(p+1)^2
  I.B2 = integrate(function(x) g.4.fun(x)^2, lower=x.lb, upper=x.ub)$value/factorial(p+1)^2
  I.B1.B2 = integrate(function(x) g.3.fun(x)*g.4.fun(x), lower=x.lb, upper=x.ub)$value/factorial(p+1)^2
}          

if (even==FALSE) {
  h.imse.pob = ( ((2*deriv+1)*C.V*I.V) / (2*(p+1-deriv)*(C.B^2*I.B1)*n) )^(1/(2*p+3))
} else {
  bw.fun  <-function(H) {abs(H^(2*p+2-2*deriv)*( C.B^2*I.B1 + H^2*C.B^2*I.B2 + 2*C.B*H*C.B*I.B1.B2 ) + C.V*I.V/(n*H^(1+2*deriv)))}
  h.imse.pob <- optimize(bw.fun, interval=c(.Machine$double.eps, 1))$minimum
  #h.imse.pob <- mean(h.mse.pob)
}



gx.lf.pob = se.lf.pob = matrix(NA,sim,neval)
gx.lp.pob = se.lp.pob = matrix(NA,sim,neval)
gx.rg.pob = se.rg.pob = matrix(NA,sim,neval)
gx.np.pob = se.np.pob = matrix(NA,sim,neval)
gx.rb.pob = se.rb.pob = matrix(NA,sim,neval)
gx.us.pob = se.us.pob = matrix(NA,sim,neval)

gx.lf.ipob = se.lf.ipob = matrix(NA,sim,neval)
gx.lp.ipob = se.lp.ipob = matrix(NA,sim,neval)
gx.rg.ipob = se.rg.ipob = matrix(NA,sim,neval)
gx.np.ipob = se.np.ipob = matrix(NA,sim,neval)
gx.rb.ipob = se.rb.ipob = matrix(NA,sim,neval)
gx.us.ipob = se.us.ipob = matrix(NA,sim,neval)

gx.rb.cepob = se.rb.cepob = matrix(NA,sim,neval)
gx.us.cepob = se.us.cepob = matrix(NA,sim,neval)

gx.lf.hat   = se.lf.hat   = matrix(NA,sim,neval)
gx.lp.hat   = se.lp.hat   = matrix(NA,sim,neval)
gx.rg.hat   = se.rg.hat   = matrix(0, sim,neval)
gx.np.hat   = se.np.hat   = matrix(NA,sim,neval)
gx.rb.hat   = se.rb.hat   = matrix(NA,sim,neval)
gx.rb.ihat  = se.rb.ihat  = matrix(NA,sim,neval)
gx.us.hat   = se.us.hat   = matrix(NA,sim,neval)
gx.us.ihat  = se.us.ihat  = matrix(NA,sim,neval)
gx.us.cehat = se.us.cehat = matrix(NA,sim,neval)
gx.rb.cehat = se.rb.cehat = matrix(NA,sim,neval)

h.mse.dpi  = matrix(NA,sim,neval)
h.mse.rot  = matrix(NA,sim,neval)
h.ce.rot   = matrix(NA,sim,neval)
h.ce.dpi   = matrix(NA,sim,neval)
h.imse.dpi = h.imse.rot = 0

lpbw.mse = lpbw.ce = matrix(NA,sim,neval)
h.hat    = matrix(NA,sim,8)
# Loop
set.seed(2016)

showwhen = 1; showevery=100
for (i in 1:sim) {
  if (i==showwhen) {cat(paste("\nSimulations Completed:",i-1,"of",sim,"-", Sys.time())); showwhen=showwhen+showevery}
  
  # Generate random data
  x = x.gen(n) 
  u = u.rnd(n)
  #for (j in 1:n) u[j] = u.rnd(x[j])
  y = g.0.fun(x) + u

  if (p==0) reg.type = "lc"
  if (p==1) reg.type = "ll"
  
  ### nprobust
  lp.bw = lpbwselect(y=y, x=x, eval=evalx, p=p, deriv=deriv, kernel=kernel, bwselect="all", bwcheck=check, bwregul=scale)
  
  h.mse.dpi[i,] = lp.bw$bws[,"h.mse.dpi"]
  h.mse.rot[i,] = lp.bw$bws[,"h.mse.rot"]
  h.ce.dpi[i,]  = lp.bw$bws[,"h.ce.dpi"]
  h.ce.rot[i,]  = lp.bw$bws[,"h.ce.rot"]
  h.imse.dpi[i] = lp.bw$bws.imse[1,1]
  h.imse.rot[i] = lp.bw$bws.imse[1,2]
  
  lpbw.imse = h.imse.dpi[i]
  
  lpol.bw = np.bw = NA
  grad = FALSE
  if (deriv==1) grad = TRUE
  
  if (deriv==0 & p>0) lpol.bw   = locpol(y~x, data=data.frame(y,x), kernel=EpaK, deg=p)$bw
  if (p<2) np.bw     = npregbw(y ~ x, regtype = reg.type, gradient=grad, ckertype="epanechnikov")$bw

  h.hat[i,] = c(lpol.bw, NA, NA, np.bw, NA,   NA, lpbw.imse, NA)
  
  for (k in 1:neval) {
    c     = evalx[k] 
    lpbw.mse[i,k] = h.mse.dpi[i,k]
    lpbw.ce[i,k]  = h.ce.dpi[i,k]
    
  ### locfit
  if (deriv==0 & p>0) {
  locfit.est  = locfit(y~lp(x, deg = p, h = h.mse.pob[k]))
  locfit.pred = predict(locfit.est, c(c), se.fit=TRUE)
  gx.lf.pob[i,k] = locfit.pred$fit
  se.lf.pob[i,k] = locfit.pred$se.fit
  
  locfit.est = locfit(y~lp(x, deg = p))
  locfit.pred = predict(locfit.est,c(c),se.fit=TRUE)
  gx.lf.hat[i,k] = locfit.pred$fit
  se.lf.hat[i,k] = locfit.pred$se.fit
  }
    
  ### locpol
  if (p>0 & deriv==0) {
  lpol.est    = locpol(y~x, data=data.frame(y,x), bw=h.mse.pob[k],   kernel=EpaK, deg=p, xeval=c)
  gx.lp.pob[i,k] = as.numeric(lpol.est$lpFit[2])
  se.lp.pob[i,k] = sqrt(as.numeric(lpol.est$lpFit[5]))

  lpol.est    = locpol(y~x, data=data.frame(y,x), bw=lpol.bw, kernel=EpaK, deg=p, xeval=c)
  gx.lp.hat[i,k] = as.numeric(lpol.est$lpFit[2])
  se.lp.hat[i,k] = sqrt(as.numeric(lpol.est$lpFit[5]))
  }
  
  ### lpridge
  lpridge.est = lpridge(x, y, bandwidth=h.mse.pob[k], deriv=deriv, x.out=c, order = p, ridge = NULL, weight = kernel, mnew = 100,var = TRUE)
  gx.rg.pob[i,k] = lpridge.est$est
  se.rg.pob[i,k] = sqrt(lpridge.est$est.var)
  
  ### np
  if (p<2) {
  np.est <- npreg(y ~ x, regtype = reg.type, bws = h.mse.pob[k], exdat=c, gradients=grad, ckertype="epanechnikov")
  gx.np.pob[i,k] = np.est$mean
  se.np.pob[i,k] = np.est$merr
  if (deriv==1) {
    gx.np.pob[i,k] = np.est$grad
    se.np.pob[i,k] = np.est$gerr
  }
  
  np.est <- npreg(y ~ x, regtype = reg.type, bws = np.bw, exdat=c, gradients=grad)
  gx.np.hat[i,k] = np.est$mean
  if (deriv==1) gx.np.hat[i,k] = np.est$grad
  se.np.hat[i,k] = np.est$merr
  }
  
  ### nprobust
  # h.mse
  lp.out = lprobust(y=y, x=x, p=p, deriv=deriv, eval=c, h=h.mse.pob[k], kernel=kernel, vce=vce, bwcheck=check, bwregul=scale)
  gx.rb.pob[i,k] = lp.out$Estimate[,"tau.bc"]
  se.rb.pob[i,k] = lp.out$Estimate[,"se.rb"]
  gx.us.pob[i,k] = lp.out$Estimate[,"tau.us"]
  se.us.pob[i,k] = lp.out$Estimate[,"se.us"]
  
  lp.out = lprobust(y=y, x=x, p=p, deriv=deriv, eval=c, h=lpbw.mse[i,k], kernel=kernel, vce=vce, bwcheck=check, bwregul=scale)
  gx.rb.hat[i,k] = lp.out$Estimate[,"tau.bc"]
  se.rb.hat[i,k] = lp.out$Estimate[,"se.rb"]
  gx.us.hat[i,k] = lp.out$Estimate[,"tau.us"]
  se.us.hat[i,k] = lp.out$Estimate[,"se.us"]
  
  # h.imse
  lp.out = lprobust(y=y, x=x, p=p, deriv=deriv, eval=c, h=h.imse.pob, kernel=kernel, vce=vce, bwcheck=check, bwregul=scale)
  gx.rb.ipob[i,k] = lp.out$Estimate[,"tau.bc"]
  se.rb.ipob[i,k] = lp.out$Estimate[,"se.rb"]
  gx.us.ipob[i,k] = lp.out$Estimate[,"tau.us"]
  se.us.ipob[i,k] = lp.out$Estimate[,"se.us"]
  
  lp.out = lprobust(y=y, x=x, p=p, deriv=deriv, eval=c, h=lpbw.imse, kernel=kernel, vce=vce, bwcheck=check, bwregul=scale)
  gx.rb.ihat[i,k] = lp.out$Estimate[,"tau.bc"]
  se.rb.ihat[i,k] = lp.out$Estimate[,"se.rb"]
  gx.us.ihat[i,k] = lp.out$Estimate[,"tau.us"]
  se.us.ihat[i,k] = lp.out$Estimate[,"se.us"]
  
  # h.ce
  #if (deriv==0) {
  lp.out = lprobust(y=y, x=x, p=p, deriv=deriv, eval=c, h=h.ce.pob[k], kernel=kernel, vce=vce, bwcheck=check, bwregul=scale)
  gx.rb.cepob[i,k] = lp.out$Estimate[,"tau.bc"]
  se.rb.cepob[i,k] = lp.out$Estimate[,"se.rb"]
  gx.us.cepob[i,k] = lp.out$Estimate[,"tau.us"]
  se.us.cepob[i,k] = lp.out$Estimate[,"se.us"]
  
  lp.out = lprobust(y=y, x=x, p=p, deriv=deriv, eval=c, h=lpbw.ce[i,k], kernel=kernel, vce=vce, bwcheck=check, bwregul=scale)
  gx.rb.cehat[i,k] = lp.out$Estimate[,"tau.bc"]
  se.rb.cehat[i,k] = lp.out$Estimate[,"se.rb"]
  gx.us.cehat[i,k] = lp.out$Estimate[,"tau.us"]
  se.us.cehat[i,k] = lp.out$Estimate[,"se.us"]
  #}
}

}

m.v.pob.m = matrix(NA,sim,neval)
for (i in 1:sim) m.v.pob.m[i,] = m.v.pob

output.fun = function(gx,se) {
  var = matrix(NA,1,neval)
  for (i in 1:neval) var[,i] = var(gx[,i])
  T = (gx - m.v.pob.m) / se
  bias = abs(colMeans(gx- m.v.pob.m, na.rm=TRUE))
  mse  = var + bias^2
  ec   = colMeans(1*(abs(T)<=qz), na.rm=TRUE)
  il   = colMeans(2*qz*se,        na.rm=TRUE)
  out=list(bias=bias, var=var, mse=mse, ec=ec, il=il)
}

locpol.pob   = output.fun(gx.lp.pob, se.lp.pob)
locfit.pob   = output.fun(gx.lf.pob, se.lf.pob)
lpridge.pob  = output.fun(gx.rg.pob, se.rg.pob)
np.pob       = output.fun(gx.np.pob, se.np.pob)
nprob.pob    = output.fun(gx.rb.pob, se.rb.pob)
npus.pob     = output.fun(gx.us.pob, se.us.pob)

nprob.ipob   = output.fun(gx.rb.ipob,se.rb.ipob)
npus.ipob    = output.fun(gx.us.ipob,se.us.ipob)
nprob.cepob   = output.fun(gx.rb.cepob,se.rb.cepob)
npus.cepob    = output.fun(gx.us.cepob,se.us.cepob)

locpol.hat   = output.fun(gx.lp.hat, se.lp.hat)
locfit.hat   = output.fun(gx.lf.hat, se.lf.hat)
lpridge.hat  = output.fun(gx.rg.hat, se.rg.hat)
np.hat       = output.fun(gx.np.hat, se.np.hat)
nprob.hat    = output.fun(gx.rb.hat, se.rb.hat)
nprob.ihat   = output.fun(gx.rb.ihat,se.rb.ihat)
npus.hat     = output.fun(gx.us.hat, se.us.hat)
npus.ihat    = output.fun(gx.us.ihat,se.us.ihat)
nprob.cehat  = output.fun(gx.rb.cehat, se.rb.cehat)
npus.cehat   = output.fun(gx.us.cehat, se.us.cehat)

#h = rep(h.imse.pob, 5*neval)
h = rep(h.mse.pob,each=8)
h[seq(7,length(h),8)]=h.imse.pob
h[seq(8,length(h),8)]=h.ce.pob
bias = c(rbind(locpol.pob$bias, lpridge.pob$bias, locfit.pob$bias,  np.pob$bias,NA, npus.pob$bias,  npus.ipob$bias,  npus.cepob$bias))
var  = c(rbind(locpol.pob$var,  lpridge.pob$var,  locfit.pob$var,   np.pob$var, NA, npus.pob$var,   npus.ipob$var,   npus.cepob$var))
mse  = c(rbind(locpol.pob$mse,  lpridge.pob$mse,  locfit.pob$mse,   np.pob$mse, NA, npus.pob$mse,   npus.ipob$mse,   npus.cepob$mse))
ec   = c(rbind(locpol.pob$ec,   lpridge.pob$ec,   locfit.pob$ec,    np.pob$ec,  NA, nprob.pob$ec,   nprob.ipob$ec,   nprob.cepob$ec))
il   = c(rbind(locpol.pob$il,   lpridge.pob$il,   locfit.pob$il,    np.pob$il,  NA, nprob.pob$il,   nprob.ipob$il,   nprob.cepob$il))
out1 = cbind(h, bias, var, mse, ec, il)

h.mse = colMeans(lpbw.mse, na.rm=TRUE)
h.ce  = colMeans(lpbw.ce,  na.rm=TRUE)
h = colMeans(h.hat, na.rm=TRUE)
h = rep(h,neval)
h[seq(6,length(h),8)]=h.mse
h[seq(8,length(h),8)]=h.ce
bias = c(rbind(locpol.hat$bias, lpridge.hat$bias, locfit.hat$bias, np.hat$bias, NA, npus.hat$bias,  npus.ihat$bias,  npus.cehat$bias))
var  = c(rbind(locpol.hat$var,  lpridge.hat$var,  locfit.hat$var,  np.hat$var,  NA, npus.hat$var,   npus.ihat$var,   npus.cehat$var))
mse  = c(rbind(locpol.hat$mse,  lpridge.hat$mse,  locfit.hat$mse,  np.hat$mse,  NA, npus.hat$mse,   npus.ihat$mse,   npus.cehat$mse))
ec   = c(rbind(locpol.hat$ec,   lpridge.hat$ec,   locfit.hat$ec,   np.hat$ec,   NA, nprob.hat$ec,   nprob.ihat$ec,   nprob.cehat$ec))
il   = c(rbind(locpol.hat$il,   lpridge.hat$il,   locfit.hat$il,   np.hat$il,   NA, nprob.hat$il,   nprob.ihat$il,   nprob.cehat$il))
out2 = cbind(h, bias, var, mse, ec, il)

ev  = rep(evalx,each=8)
out = cbind(ev,out1,out2)

write.csv(out, file=paste("output/lp","_p",p,"_d",deriv,"_r",scale,"_n",n,".csv",sep=""))

table1 = as.matrix(out[,2:ncol(out)])
table=formatC(table1,     format = "f", digits = 3)
for (i in seq(2,nrow(table1),8)) table[i,7:12] = paste(rep("\\texttt{na}",6))
for (i in seq(3,nrow(table1),8)) table[i,7]    = paste("\\texttt{na}")
for (i in seq(5,nrow(table1),8)) table[i,] = paste(rep("",6))

if (deriv==1 | p==0) {
  for (i in seq(1,nrow(table1),8)) table[i,1:12]=paste(rep("\\texttt{na}",6))
  for (i in seq(3,nrow(table1),8)) table[i,1:12]=paste(rep("\\texttt{na}",6))
}

if (p>1) {
  for (i in seq(4,nrow(table1),8)) table[i,1:12]=paste(rep("\\texttt{na}",6))
}

#if (deriv==1) {
#  for (i in seq(8,nrow(table1),8)) table[i,1:12]=paste(rep("\\texttt{na}",6))
#}


colnames(table) = rep(c("$h$","Bias","Var","MSE","EC","IL"),2)
rownames(table) = rep(c("\\texttt{locpol}","\\texttt{lpridge}","\\texttt{locfit}","\\texttt{np}","\\texttt{nprobust}", "      $h_{\\texttt{MSE}}$", "     $h_{\\texttt{IMSE}}$", "     $h_{\\texttt{CE}}$"),neval)

evals=out[seq(1,nrow(table1),8),1]
eval.group = paste("$x=$",round(evals,2))
ncol = ncol(table)
table1_tex = latex(table, file = paste("output/lp","_p",p,"_d",deriv,"_r",scale,"_n",n,".txt",sep=""), landscape=FALSE,
                   outer.size='scriptsize', col.just=rep('c',ncol), center='none', title='', table.env=FALSE,
                   n.rgroup=rep(8,neval), rgroup = eval.group,
                   n.cgroup=c(6,6), cgroup = c("Population bandwidth","Estimated bandwidth")) 


### lpbwselect output

out=matrix(0,neval+1,6)

out[1:neval,1] = h.mse.pob
out[1:neval,2] = colMeans(h.mse.dpi, na.rm=TRUE)
out[1:neval,3] = colMeans(h.mse.rot, na.rm=TRUE)
out[1:neval,4] = h.ce.pob
out[1:neval,5] = colMeans(h.ce.dpi, na.rm=TRUE)
out[1:neval,6] = colMeans(h.ce.rot, na.rm=TRUE)

out[neval+1,1] = h.imse.pob
out[neval+1,2] = mean(h.imse.dpi, na.rm=TRUE)
out[neval+1,3] = mean(h.imse.rot, na.rm=TRUE)

write.csv(out, file=paste("output/lpbw","_p",p,"_d",deriv,"_r",scale,"_n",n,".csv",sep=""))

outex = formatC(out,     format = "f", digits = 3)
#outex[neval+1,4:6] = rep(paste("-"),3)

colnames(outex) = c("POP", "DPI", "ROT", "POP", "DPI", "ROT")
rownames(outex) = c(paste("$x=",evalx,"$"),"\\textbf{IMSE}")

table1_tex = latex(outex, file = paste("output/lpbw","_p",p,"_d",deriv,"_r",scale,"_n",n,".txt",sep=""),
                   landscape=FALSE, outer.size='scriptsize', col.just=rep('c',6), center='none', title='', 
                   table.env=FALSE, n.rgroup=c(neval,1), n.cgroup=c(3,3), cgroup = c("MSE","CE")) 



