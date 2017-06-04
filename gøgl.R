library(imager)
library(spatstat)

#Sets layout in plots (3x3 matrix)
layout(matrix(c(1:9), ncol=3))

#Loads Original image, which is a Window of a bread image
im <- load.image("brødwind.JPG")

##-------------------------------------------------------------
# Plots Original image, Original grayscaled image and Original 
# image with threshold determined from Kmeans
##-------------------------------------------------------------
plot(im, main = "Original") 

im.gray <- im %>% grayscale 
plot(im.gray, main = "Original")

im.thres <- threshold(im.gray, "auto", TRUE)
plot(im.thres, main = "Original")

##-------------------------------------------------------------
# Plots anisotropic blurred image with a certain amplitude and 
# sharpness, and thereafter the grayscaled version and the 
# image with threshold determined from Kmeans
##-------------------------------------------------------------
im.blur.anis <- blur_anisotropic(im,ampl=1e4,sharp=1) 
plot(im.blur.anis, main = "Anisotropic blur")

im.gray.anis <- im.blur.anis %>% grayscale 
plot(im.gray.anis, main = "Anisotropic blur")

im.thres.anis <- threshold(im.gray.anis, "auto", TRUE)
plot(im.thres.anis, main = "Anisotropic blur")

##-------------------------------------------------------------
# Plots isotropic (Gaussian) blurred image with a certain std. 
# div., and thereafter the grayscaled version and the image 
# with threshold determined from Kmeans
##-------------------------------------------------------------
im.blur <- isoblur(im, 3, neumann = FALSE, gaussian = TRUE)
plot(im.blur, main = "Isotropic blur")

im.gray.blur <- im.blur %>% grayscale 
plot(im.gray.blur, main = "Isotropic blur")

im.thres.blur <- threshold(im.gray.blur, "auto", TRUE)
plot(im.thres.blur, main = "Isotropic blur")





####Laslett####
layout(1)
im.thres = cimg2im(im.thres)
laslett = laslett(im.thres)

im.thres.anis = cimg2im(im.thres.anis)
laslett = laslett(im.thres.anis)

im.thres.blur = cimg2im(im.thres.blur)
laslett = laslett(im.thres.blur)




