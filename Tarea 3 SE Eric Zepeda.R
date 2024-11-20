

f = function(x){ #pdf de la densidad objetivo de juguete
  dnorm(x, mean=0, sd = 1)
}


Nsim <- 10^5


### EEJMPLO ####
#este ejemplo lo hice para ver como funciona realmente.
x<- double(Nsim)
x[1]<- rnorm(1)
#x<- matrix(0,nrow = Nsim, ncol = 2)
#x[1,]<- rnorm(2)

for(i in 2:Nsim - 1){
  z<- rnorm(2, mean=0, sd = 1)
  y<-rnorm(1,mean = x[i]+2*z, 0.5)
  if(runif(1) < min(1, f(y) / f(x[i]))) #Asumiendo simetría
  {
    x[i+1] <- y
  } else x[i+1] <- x[i]
}
  
muestra<- tail(x,n=5000) 
x_axis <- as.list(seq(1, 5000))



### Metropolis Hasting TAREA  #########

dev.new()

f = function(x){ #the pdf of f - the target distribution
  return((1/(20216.335877))*exp(-0.5*((x[1]^2)*(x[2]^2)+(x[1]^2)+(x[2]^2)-8*x[1]-8*x[2])))
}
Nsim <- 10^5

x<- matrix(0,nrow = Nsim, ncol = 2)
x[1,]<- rnorm(2)
for(i in 2:Nsim - 1){
  z<- rnorm(2, mean=0, sd = 1)
  #y<-rnorm(2,mean = x[i,]+2*z, 1)
  y<- x+2*z
  if(runif(1) < min(1, f(y) / f(x[i,]))) #same as above if q(x,y) == q(y,x), i.e. q is symmetric
  {
    x[i+1,1] <- y[1]
    x[i+1,2] <- y[2]
  } else x[i+1,] <- x[i,]

}

muestra_x<- tail(x[,1],n=5000) 
muestra_y<- tail(x[,2],n=5000)
x_axis <- as.list(seq(1, 5000))

plot(muestra_x, type = "l", main = "Muestra", xlab="Indice", ylab="X1")
plot(muestra_y, type = "l", main = "Muestra", xlab="Indice", ylab="X2")
print(sum(muestra_x)/5000)
print(sum(muestra_y)/5000)


##########

dev.new()
Nsim <- 10^5
x<- rep(0,Nsim)
y<- rep(0,Nsim)
for(i in 2:Nsim - 1){
z <-  rnorm(1, mean=0, sd = 1)
a <- 1/(1+x[i]^2)
y[i] <- 4*a + z*sqrt(a)
z <-  rnorm(1, mean=0, sd = 1)
b <- 1/(1+y[i]^2)
x[i]<-4*b+z*sqrt(b)
}


muestra_x<- tail(x, n=5000) 
muestra_y<- tail(y, n=5000)
x_axis <- as.list(seq(1, 5000))

plot(muestra_x, type = "l", main = "Muestra", xlab="Indice", ylab="X")

print(sum(muestra_x)/5000)
print(sum(muestra_y)/5000)
