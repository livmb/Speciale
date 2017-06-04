setwd("C:/Users/livmb/Dropbox/UNI/Speciale/Brød billeder/brød")
library(imager)
library(spatstat)

#Sets layout in plots (3x3 matrix)
layout(matrix(c(1:3), ncol=3))

#Loads Original image, which is a Window of a bread image
im <- load.image("BlåBrødPågen7cmTop.JPG")

##-------------------------------------------------------------
#### Picture Manipulation ####
##-------------------------------------------------------------
# Plots Original image, extracted Blue channel image and the 
# thresholded version of this. The threshold is here determined 
# by Kmeans. 
##-------------------------------------------------------------
layout(matrix(c(1:3), ncol=3))
plot(im, xlab = "", ylab = "")#, main = "Original")

BRGim <- B(im) - (R(im)+G(im)) #Extracting the blue in the Original Image
BRGim.normalize <- (BRGim - min(BRGim))/(max(BRGim)-min(BRGim)) #Normalizes to image values 
plot(BRGim.normalize, xlab = "", ylab = "")#,main = "Blue-(Red+Green)")

thres.im <- threshold(BRGim.normalize, thr = "auto") #Thresholded Image
plot(thres.im, xlab = "", ylab = "")#, main = "Threshold")


##-------------------------------------------------------------
# Plotting the thresholded image in a observation 
# window, W, and then performing Lasletts test.
##-------------------------------------------------------------
layout(1)

thres.im <- cimg2im(thres.im) #Change the image type
W <- owin(c(thres.im$xrange[1], thres.im$xrange[2]), c(thres.im$yrange[1], thres.im$yrange[2])) #defining the obs window W, here we just look at the whole image
Y <- thres.im[W]  

laslett = laslett(Y)

W.prime.prime = laslett$Rect
x.point <- laslett$df$newx
y.point <- laslett$df$y
X <- ppp(x.point, y.point, W.prime.prime)
plot(X, main = "")

##------------------------------------------------##
#### Test for Complete Spacial Randomness (CSR) ####
##------------------------------------------------##
#------------------------------------#
####         Method 3:            ####
#
# Clark And Evans Test
#------------------------------------#
clarkevans.test(X) #alternative="two.sided" 
clarkevans.test(X, correction="D", nsim = 10^4) #alternative="two.sided" + Donnelly correction 
#clarkevans.test(X, correction="cdf", nsim = 10^4) 


##------------------------------------------------------------##
#### Estimation of the Intensity, volume fraction and alpha ####
##------------------------------------------------------------##
#Area of W
top.length = 7.0 #cm
top.pix <- W$xrange[2]-W$xrange[1] #1381 pix
pixel.width = top.length/top.pix #How many cm the width of one pixel is
pixel.height = pixel.width #pixels are quadratic 
pixel.area = pixel.width * pixel.height #area of one pixel in cm^2

area.W.prime.prime = area(W.prime.prime)*pixel.area #cm^2

#Estimate for the intensity gamma
no.points.W.prime.prime <- X$n #135
gamma.hat = no.points.W.prime.prime/area.W.prime.prime #57.5285


#Estimate for the volume fraction p
no.black.particles <- sum(thres.im < 1) #number of black particles in the thresholded image
#Note: sum(thres.im) gives the number of white pixels
p.hat <- no.black.particles/area(W) #Because area.W.pix = no.black.particles + no.white.particles = total number of particles in W

#Estimate of \hat{\alpha}
alpha.hat <- sqrt(-(2*gamma.hat*pi)/(log(1-p.hat)))

##-----------------------##
#### Varians Estimates ####
##-----------------------##
#----------------#
#### Method 1 ####
#----------------#

#Estimate of \sigma_{2,2}^2
sigma.hat.one <- ((6*gamma.hat*pi^2*(1-p.hat)^2)/(alpha.hat)^3) #0.06277515

#Confidence interval
KI.one = c(p.hat - 1.96*sqrt(sigma.hat.one), p.hat + 1.96*sqrt(sigma.hat.one)) 
#### ???? negativt! What! [-0.1128543  0.8693006]

#----------------#
#### Method 2 ####
#----------------#
#Estimate of r_0
r.zero.hat <- sqrt(-(log(1-p.hat))/(gamma.hat*pi))

#Covariogram fct
F2 <- function(x){
  if(0 < x && x < 2*r.zero.hat){
    F2 <- ((4*x)/(pi*r.zero.hat^2))*(acos(x/(2*r.zero.hat)) - (x/(2*r.zero.hat))*sqrt(1-(x^2/(4*r.zero.hat^2))))
  } else{
    F2 <- 0
  }
  return(F2)
}
F2 <- Vectorize(F2);
curve(F2, 0-0.1, 2*r.zero.hat+0.1)

int <- function(x){
  return((exp(gamma.hat*r.zero.hat^2*(F2(x/r.zero.hat)))-1)*x)
}
int <- integrate(int, 0, Inf)$value #1.51875e-05

#Estimate of sigma
sigma.hat.two <- 2*pi*(1-p.hat)^2*int #3.689225e-05

#Confidence interval
KI.two = c(p.hat - 1.96*sqrt(sigma.hat.two), p.hat + 1.96*sqrt(sigma.hat.two)) #[0.3663183 0.3901280]




