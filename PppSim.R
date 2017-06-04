library(spatstat)

#### Simulating Poisson point process #####
int = 1/3 #intensity
W = owin(c(0,5),c(0,5)) #observationwindow, W=[0,5]x[0,5] 

#Simulates a homogeneous Poisson point process with intensity 1/3 in W
X = rpoispp(int, win=W) 
x = X$x #x-values in plot
y = X$y #y-values in plot
points = cbind(x,y) #matrix of coordinates of points 

plot(X, pch = 4, xlim = c(0,5), ylim = c(0,5), xlab = "", ylab = "", main = "", lwd = 2)

#Calculates the Euclidean distance between the points 
dist = as.matrix(dist(points))

#Checks if the points are different and has distance leq to 1. Returns matrix consisting of TRUE and FALSE. 
bdist = ((dist <= 1 & dist > 0) & upper.tri(dist)) 

#Matrix of indices in bdist, where bdist == TRUE
indices = which(bdist, arr.ind = TRUE)

##-------------------------------------------------------
## Functionname: draw.line
## Purpose:      Draw lines between the points in the
##               matrix points 
## Input:        v: vector of indices 
## Output:       B: lines btween points in plot
##-------------------------------------------------------
draw.line <- function(v){
  segments(points[v[1], 1], points[v[1], 2], points[v[2], 1], points[v[2], 2], col = "red", lwd = 2)
}
apply(indices, 1, draw.line)

#Number of points in Ppp, which is in a distance at most one from each other 
sum.of.points = sum(bdist) 
