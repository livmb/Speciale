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
plot(im, xlab = "", ylab = "")

BRGim <- B(im) - (R(im)+G(im)) #Extracting the blue colorchanel in the Original Image
BRGim.normalize <- (BRGim - min(BRGim))/(max(BRGim)-min(BRGim)) #Normalizes to image values 
plot(BRGim.normalize, xlab = "", ylab = "")

thres.im <- threshold(BRGim.normalize, thr = "auto") #Thresholded Image
plot(thres.im, xlab = "", ylab = "")


##-------------------------------------------------------------
# Plotting the thresholded image in a observation 
# window, W, and then performing Lasletts test.
##-------------------------------------------------------------
layout(1)

thres.im <- cimg2im(thres.im) #Change the image type
W <- owin(c(thres.im$xrange[1], thres.im$xrange[2]), c(thres.im$yrange[1], thres.im$yrange[2])) #defining the obs window W, here we just look at the whole image
Y <- thres.im[W]  

top.length = 7 #cm
top.pix <- W$xrange[2]-W$xrange[1] #1381 pix
pixel.width = top.length/top.pix #How many cm the width of one pixel is
pixel.height = pixel.width #pixels are quadratic 
pixel.area = pixel.width * pixel.height #area of one pixel in cm^2
area.W <- area(W)*pixel.area #area of W in cm^2

laslett = laslett(Y, plotit = F) #Lasletts transformation

W.prime.prime = laslett$Rect
x.point <- laslett$df$newx
y.point <- laslett$df$y
X <- ppp(x.point, y.point, W.prime.prime)
plot(X, main = "")

no.points.W.prime.prime <- X$n #135
area.W.prime.prime = area(W.prime.prime)*pixel.area #cm^2

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

clarkevans.test(X, alternative="clustered")
clarkevans.test(X, alternative="clustered", correction="D", nsim = 10^4) 


##-----------------------------------------------------##
#### Estimation of the Intensity and volume fraction ####
##-----------------------------------------------------##
#Estimate for the intensity
gamma.hat = no.points.W.prime.prime/area.W.prime.prime #57.5285

#Estimate for the volume fraction p
no.black.particles <- sum(thres.im < 1) #number of black particles in the thresholded image
#Note: sum(thres.im) gives the number of white pixels
p.hat <- no.black.particles/area(W) #Because area.W.pix = no.black.particles + no.white.particles = total number of particles in W

##-----------------------##
#### Varians Estimates ####
##-----------------------##

#----------------#
#### Method 1 ####
#----------------#
#Estimate of r_0
r.zero.hat <- sqrt(-(log(1-p.hat))/(gamma.hat*pi))

#Covariogram fct for B^2
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
sigma.hat <- 2*pi*(1-p.hat)^2*int #3.689225e-05

#Confidence interval
KI = c(p.hat - qnorm(0.975)*sqrt(1/area.W)*sqrt(sigma.hat), p.hat + qnorm(0.975)*sqrt(1/area.W)*sqrt(sigma.hat)) #[0.3752955 , 0.3811508]

#----------------#
#### Method 2 ####
#----------------#
t.zero <- 2*r.zero.hat

x.less = floor(690.5-r.zero.hat)
x.greater = ceiling(690.5 + r.zero.hat)
y.less = floor(233 - r.zero.hat)
y.greater = ceiling(233 + r.zero.hat)
W.one <- owin(c(W$xrange[1], x.less), c(y.greater, W$yrange[2]))
W.three <- owin(c(W$xrange[1], x.less), c(W$yrange[1], y.less))
W.two <- owin(c(x.greater, W$xrange[2]), c(y.greater, W$xrange[2]))
W.four <- owin(c(x.greater, W$xrange[2]), c(W$yrange[1], y.less))

layout(matrix(c(1:4), ncol=2))
W.one <- Y[W.one]
W.two <- Y[W.two]
W.three <- Y[W.three]
W.four <- Y[W.four]
plot(W.one)
plot(W.two)
plot(W.three)
plot(W.four)

no.black.particles.W.one <- sum(W.one < 1) 
p.hat.one <- no.black.particles.W.one/area(W.one) 
p

