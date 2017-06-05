setwd("C:/Users/livmb/Dropbox/UNI/Speciale/Brød billeder/brød")

library(imager)
library(spatstat)

##Loads Original image, which is a Window of the bread image
#Note: BlåBrødPågen15,5cmTop.jpg is the picture in the largest observation window, 
#while BlåBrødPågen7cmTop.jpg is the picture, where the holes are more homogeneous
im <- load.image("BlåBrødPågen7cmTop.JPG")


##-------------------------------------------------------------
#### Picture Manipulation ####
##-------------------------------------------------------------
# Plots Original image, extracted Blue channel image and the 
# thresholded version of this. The threshold is here determined 
# by Kmeans. 
##-------------------------------------------------------------
layout(matrix(c(1:3), ncol=3))
plot(im, xlab = "", ylab = "")

BRGim <- B(im) - (R(im)+G(im)) #Extracting the blue color chanel in the Original Image
BRGim.normalize <- (BRGim - min(BRGim))/(max(BRGim)-min(BRGim)) #Normalizes to image values 
plot(BRGim.normalize, xlab = "", ylab = "")

thres.im <- threshold(BRGim.normalize, thr = "auto") #Thresholded Image
plot(thres.im, xlab = "", ylab = "")


##-------------------------------------------------------------
# Plotting the thresholded image in an observation 
# window W and then performing Lasletts test.
##-------------------------------------------------------------
layout(1)

thres.im <- cimg2im(thres.im) #Change the image type from cimg object to im object

#Defining the observation window W. Here we just look at the whole image, since we have
#defined the observation window before loading. 
W <- owin(c(thres.im$xrange[1], thres.im$xrange[2]), c(thres.im$yrange[1], thres.im$yrange[2])) 
Y <- thres.im[W]  

top.length = 7 #Set the top length of the image 
top.pix <- W$xrange[2]-W$xrange[1] 
pixel.width = top.length/top.pix #How many cm the width of one pixel is
pixel.height = pixel.width #pixel hieght = pixel length, since pixels are quadratic 
pixel.area = pixel.width * pixel.height #area of one pixel in cm^2
area.W <- area(W)*pixel.area #area of W in cm^2

laslett = laslett(Y, plotit = F) #Lasletts transformation

W.prime.prime = laslett$Rect #the observation window W''
x.point <- laslett$df$newx
y.point <- laslett$df$y
X <- ppp(x.point, y.point, W.prime.prime) #the Laslett transformed process
plot(X, main = "")

no.points.W.prime.prime <- X$n #no. of points in the Laslett transformed process
area.W.prime.prime = area(W.prime.prime)*pixel.area #area of W'' in cm^2

##------------------------------------------------##
#### Test for Complete Spacial Randomness (CSR) ####
##------------------------------------------------##
#------------------------------------#
####    Clark And Evans Test      ####
#------------------------------------#
clarkevans.test(X) #alternative="two.sided" 
clarkevans.test(X, correction="D", nsim = 10^4) #alternative="two.sided" + Donnelly correction 

##-----------------------------------------------------##
#### Estimation of the Intensity and volume fraction ####
##-----------------------------------------------------##
#Estimate for the intensity
gamma.hat = no.points.W.prime.prime/area.W.prime.prime 

#Estimate for the volume fraction p
no.black.pixels <- sum(thres.im < 1) #number of black pixels in the thresholded image
p.hat <- no.black.pixels/area(W) #Because area.W.pix = no.black.pixels + no.white.prixels = total number of pixels in W

##-----------------------##
#### Varians Estimate  ####
##-----------------------##
#Estimate of R_0
R.zero.hat <- sqrt(-(log(1-p.hat))/(gamma.hat*pi))

#Covariogram fct for B^2
F2 <- function(x){
  if(0 < x && x < 2*R.zero.hat){
    F2 <- ((4*x)/(pi*R.zero.hat^2))*(acos(x/(2*R.zero.hat)) - (x/(2*R.zero.hat))*sqrt(1-(x^2/(4*R.zero.hat^2))))
  } else{
    F2 <- 0
  }
  return(F2)
}
F2 <- Vectorize(F2);

int <- function(x){
  return((exp(gamma.hat*R.zero.hat^2*(F2(x/R.zero.hat)))-1)*x)
}
int <- integrate(int, 0, Inf)$value 

#Estimate of sigma_{2,2}^2
sigma.hat <- 2*pi*(1-p.hat)^2*int 

#95% confidence interval for the volume fraction p 
KI = c(p.hat - qnorm(0.975)*1/sqrt(area.W)*sqrt(sigma.hat), p.hat + qnorm(0.975)*1/sqrt(area.W)*sqrt(sigma.hat))












