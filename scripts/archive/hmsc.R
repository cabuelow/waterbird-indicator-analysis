# joint species distribution model (JSDM)

# libaries

library(Hmsc)
library(MASS)
library(tidyverse)
library(tmap)
library(sf)
library(raster)

vignette(package = "Hmsc")

set.seed(6)

# generate simulated data

n = 100
ns = 5
beta1 = c(-1,-1,0,1,2)
alpha = rep(0,ns)
beta = cbind(alpha, beta1)
x = cbind(rep(1,n), rnorm(n))
Lf = x%*%t(beta)

xycoords2 = matrix(runif(2*n), ncol=2)

colnames(xycoords2) = c('x-coordinate', 'y-coordinate')
rownames(xycoords2) = 1:n

sigma.spatial = c(2)
alpha.spatial = c(0.35)
Sigma = sigma.spatial^2*exp(-as.matrix(dist(xycoords))/alpha.spatial)
eta1 = mvrnorm(mu=rep(0,n), Sigma=Sigma)
lambda1 = c(1,2,-2,-1,0)
Lr = eta1%*%t(lambda1)
L = Lf + Lr
y = as.matrix(L + matrix(rnorm(n*ns),ncol=ns))
yprob = 1*((L +matrix(rnorm(n*ns),ncol=ns))>0)
XData = data.frame(x1=x[,2])

# plot the spatially structured data

rbPal = colorRampPalette(c('cyan','red'))
par(mfrow=c(2,3))
Col = rbPal(10)[as.numeric(cut(x[,2],breaks = 10))]
plot(xycoords[,2],xycoords[,1],pch = 20,col = Col,main=paste('x'))
for(s in 1:ns){
  Col = rbPal(10)[as.numeric(cut(y[,s],breaks = 10))]
  plot(xycoords[,2],xycoords[,1],pch = 20,col = Col,main=paste('Species',s))
}

# spatially explicit model in Hmsc
# use SData to give coordinates

studyDesign = data.frame(sample = as.factor(1:n))
rL.spatial2 = HmscRandomLevel(sData = xycoords2)
rL.spatial2 = setPriors(rL.spatial2,nfMin=1,nfMax=1) # limiting to one LV
m.spatial = Hmsc(Y=yprob, XData=XData, XFormula=~x1,
                 studyDesign=studyDesign, ranLevels = list('sample'=rL.spatial), distr = 'probit')

# set mcmc sampling parameters

nChains = 2
test.run = TRUE
if (test.run){
  # with this option, the vignette runs fast but results are not reliable
  thin = 1
  samples = 10
  transient = 5
  verbose = 0
} else {
  # with this option, the vignette evaluates slow but it reproduces the results of
  # the .pdf version
  thin = 10
  samples = 1000
  transient = 1000
  verbose = 0 # set to 0 to suppress output reporting how MCMC samplingis proceeding
}

m.spatial = sampleMcmc(m.spatial, thin = thin, samples = samples, transient = transient,
                       nChains = nChains, verbose = verbose,updater=list(GammaEta=FALSE))
