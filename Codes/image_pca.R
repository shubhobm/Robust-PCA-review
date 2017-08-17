#' Reconstructing Images using PCA
#' Examples are shown for PCA 
#' TODO: example of reconstruction using KPCA, depth-pca and other pca
rm(list=ls())
library(png)
library(EBImage)
library(ripa)
library(rpca)
library(robustbase)
library(fda.usc)
library(rrcov)
library(pixmap)

setwd("D:/Study/My projects/Robust-PCA-review/Codes")
source('misc_functions.R')

## function to do proces images
process.img = function(img){
  ## add noise
  n = nrow(img)
  p = ncol(img)
  img.noise = img
  img.noise[sample(1:n,20),sample(1:p,20)] = sample(c(0,1), 20^2, replace=T)
  
  mus = apply(img.noise, 2, mean)
  sigmas = apply(img.noise, 2, sd)
  img.noise.scaled = scale(img.noise, mus, scale=F)
  
  res <- prcomp(img.noise.scaled)
  pcnum <- 10 #Use first 20 principal components
  img.trunc <- res$x[,1:pcnum] %*% t(res$rotation[,1:pcnum])
  img.trunc = scale(img.trunc, -mus, scale=F)
  
  ## PCP
  pcpmod = rpca(img.noise.scaled)
  img.pcp = scale(pcpmod$L, -mus, scale=F)
  
  ## DPCA
  pcarank <- PcaRank(img.noise.scaled, k=pcnum)
  img.dpca  = pcarank@scores %*% t(pcarank@loadings)
  img.dpca = scale(img.dpca, -mus, scale=F)
  # img.trunc = pcarank@scores
  
  ## ROBPCA
  pcarob <- PcaHubert(img.noise.scaled, k=pcnum)
  img.rpca  = pcarob@scores %*% t(pcarob@loadings)
  img.rpca = scale(img.rpca, -mus, scale=F)
  
  plot(imagematrix(img))
  plot(imagematrix(img.noise))
  plot(imagematrix(img.dpca))
  plot(imagematrix(img.rpca))
  plot(imagematrix(img.pcp))
}

lenna = readPNG("lenna.png")
yale01 = imagematrix(getChannels(read.pnm(("yaleB01_P00A+005E+10.pgm"))))
yale02 = imagematrix(getChannels(read.pnm(("yaleB02_P00A+005E+10.pgm"))))
yale28 = imagematrix(getChannels(read.pnm(("yaleB28_P00A+005E+10.pgm"))))

lenna.small = imagematrix(resize(lenna, 96, 96))
yale01.small = imagematrix(resize(yale01, 96, 84))
yale02.small = imagematrix(resize(yale02, 96, 84))
yale28.small = imagematrix(resize(yale28, 96, 84))

set.seed(08162017)
defaultPar = par()
par(mfrow=c(4,5), mai=rep(0.1,4))
process.img(lenna.small)
process.img(yale01.small)
process.img(yale02.small)
process.img(yale28.small)
par(defaultPar)

r1 <- writePNG(img.noise,target="lenna_orig.png")
r <- writePNG(img.trunc,target="lenna_recon.png")
writePNG(pcpmod$L,target="lenna_recon_pcp.png")
#r3 <- writePNG(img.trunc.kpc,target="lenna_recon_kpc.png")
