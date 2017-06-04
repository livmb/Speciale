setwd("C:/Users/livmb/Dropbox/UNI/Speciale/Brød billeder/brød")
require(imager)
require(spatstat)

layout(matrix(c(1:6), ncol=2))
#Loads Original image, which is a Window of a bread image
im <- load.image("ScanBlåMindre.JPG")
plot(im, main = "Original") 

im.blur <- isoblur(im, 1, neumann = FALSE, gaussian = TRUE)
plot(im.blur, main = "Isotropic blur")

im.gray.blur <- grayscale(im.blur) # %>% grayscale 
plot(im.gray.blur, main = "Isotropic blur")

im.thres.blur <- threshold(im.gray.blur)
plot(im.thres.blur, main = "Isotropic blur")
#display(im.thres.blur, TRUE)

#mask <- imfill(5,10,val=1) #Rectangular mask that can be used in mopening
im.thres.blur = mopening_square(im.thres.blur, 5)
plot(im.thres.blur, main = "Morphological opening")

layout(1)
im.thres.blur <- cimg2im(im.thres.blur)
w <- owin()
w <- owin(c(10, 200), c(10, 200))
Y <- im.thres.blur[w]  
plot(Y, main = "Morphological opening i et observationsvindue")

laslett.blur = laslett(Y)
#plot(laslett.blur$oldX) #Original Y
#plot(laslett.blur$TanOld) #Point pattern for tangent point 
#plot(laslett.blur$TanNew) #Point pattern for new tangent points, with changed obs-window
#plot(laslett.blur$Rect) #Plots the new obs-window 
#laslett.blur$df #Dataframe with all the new values 
#laslett.blur$type #Specify that it's the lower tagent point 

layout(matrix(c(1:2), ncol=2))
W = laslett.blur$Rect
x.point <- laslett.blur$df$newx
y.point <- laslett.blur$df$y
X <- ppp(x.point, y.point, W)
plot(X, main = "Den nye proces")


##------------------------------------------------##
#### Test for Complete Spacial Randomness (CSR) ####
##------------------------------------------------##

#------------------------------------#
#             Method 1:
#
# Pearsons chi^2 goodness-of-fit test 
# of CSR using Quadrat Counts
#------------------------------------#

layout(1)
test1 <- quadrat.test(X, nx = 2, ny =2) 
test1$p.value # > 0.05, which means that we can't reject the CSR hypothesis. 
plot.ppp(X, cols = "red", main = "GOF using Quadrat Counts")
plot(test1, add = TRUE) #(observed counts at top left; expected count at top right; Pearson residuals at bottom).

#--------------------------------#
#         Method 2:
#
# Kolmogorov-Smirnov, Cramer-von 
# Mises or Anderson-Darling test 
# of CSR
#--------------------------------#

layout(1)

# We are simply comparing the observed and expected distributions of the x coordinate.
fun <- function(x,y){
  x
}

test2 <- cdf.test(X, fun, test = "ks") #Kolmogorov-Smirnov test
test3 <- cdf.test(X, fun, test = "cvm") #Cramer-von Mises test
test4 <- cdf.test(X, fun, test = "ad") #Anderson-Darling test
test2$p.value # > 0.05
test3$p.value # > 0.05
test4$p.value # > 0.05

plot(test2)
plot(test3)
plot(test4)





####bermantest 
#ppm(X, ~1)

# fit inhomogeneous Poisson model and test
model <- ppm(X, ~x)
cdf.test(model, "x")

if(interactive()) {
  # synthetic data: nonuniform Poisson process
  X <- rpoispp(function(x,y) { 100 * exp(x) }, win=square(1))
  
  # fit uniform Poisson process
  fit0 <- ppm(X ~1)
  # fit correct nonuniform Poisson process
  fit1 <- ppm(X ~x)
  
  # test wrong model
  cdf.test(fit0, "x")
  # test right model
  cdf.test(fit1, "x")
}



