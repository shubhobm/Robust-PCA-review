#' Reconstructing Images using PCA
#' Examples are shown for PCA 
#' TODO: example of reconstruction using KPCA, depth-pca and other pca

library(png)
library(EBImage)
library(ripa)
setwd("D:/Study/My projects/Robust-PCA-review")
lenna.orig <- readPNG("lenna.png")
lenna64 = imagematrix(resize(lenna.orig, 64, 64))

## add noise
lenna.noise = lenna64
for(i in 31:50)
  for(j in 31:50)
    lenna.noise[i,j] = runif(1)

res <- prcomp(lenna.orig, center = TRUE, scale = FALSE)
pcnum <- 20 #Use first 20 principal components
lenna.trunc <- res$x[,1:pcnum] %*% t(res$rotation[,1:pcnum])
#how to recover image in kpca?
#lenna.trunc.kpc <- res.kpca$proj%*%t(res.kpca$pc)%*%res.kpca$x
if(res$scale != FALSE){
  lenna.trunc <- scale(lenna.trunc, center = FALSE , scale=1/res$scale)
}
if(res$center != FALSE){
  lenna.trunc <- scale(lenna.trunc, center = -1 * res$center, scale=FALSE)
}


# par(mfrow=c(1,2))
# image(t(t(lenna.noise)))
# image(lenna.trunc)
# par(mfrow=c(1,1))

## PCP
library(rpca)
pcpmod = rpca(lenna.noise)

r1 <- writePNG(lenna.noise,target="lenna_orig.png")
r <- writePNG(lenna.trunc,target="lenna_recon.png")
writePNG(pcpmod$L,target="lenna_recon_pcp.png")
#r3 <- writePNG(lenna.trunc.kpc,target="lenna_recon_kpc.png")
