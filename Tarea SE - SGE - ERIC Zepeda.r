# Librerias #


# Funciones #
# h #
h <- function(x,y){
    valor <- ((x*sin(20*y)+y*sin(20*x))^2)*cosh(sin(10*x)*x)+((x*cos(10*y)-y*sin(10*x))^2)*cosh(cos(20*y)*y)
    return(valor)
}
# Gammas #
#   1    #
#gamma <- function(k){
#  return(1/log(1+k))
#}


#   2    #
gamma <- function(k){
  return(1/(100*log(1+k)))
         }


#   3    #
#gamma <- function(k){
#  return(1/(k+1))
#}

#### LAMBDA #####
#      1     ####
#lambda<-function(k){
#  return(1/(log(1+k))^(1/10))
#}

####   2     #####
lambda<-function(k){
  return(1/(100*log(1+k)))
}

#### 3 ########
#lambda <- function(k){
#  return(1/(k+1)^0.5)
#}
#### 4 ########
#lambda <- function(k){
#  return(1/(k+1)^0.1)
#}




#### SGD TEST   ####3
SGD_2D <- function(h,theta,max_iter,tol,gamma, lambda){
    iter <- 0
    lista_iter<-list()
    lista_theta<-list()
    lista_error<-list()
    error_sgd <-sqrt(theta[1]^2+theta[2]^2)
    lista_theta<-append(lista_theta,theta)
    lista_error<-append(lista_error,error_sgd)
    k<-0
    while (iter < max_iter & tol < error_sgd ) {
      k <- k+1
      iter <- iter+1
      #z_k <- runif(2)
      z_k <-rnorm(2)
      z_k <- z_k/sqrt(z_k[1]^2+z_k[2]^2)
      x_p<-theta[1] + lambda(k)*z_k[1]
      y_p<- theta[2] + lambda(k)*z_k[2]
      x_n<-theta[1] - lambda(k)*z_k[1]
      y_n<-theta[2] - lambda(k)*z_k[2]
      delta <- h(x_p,y_p) - h(x_n,y_n)
      theta<- theta - (gamma(k)/(2*lambda(k)))*delta*z_k
      error_sgd<-sqrt(theta[1]^2+theta[2]^2)
      lista_theta<-append(lista_theta,theta)
      lista_error<-append(lista_error,error_sgd)
      if(error_sgd>100){
        break
      }
    }
lista_iter<- as.list(seq(1, iter))
return(list(theta = theta,iteraciones = lista_iter,error = lista_error,thetas =lista_theta))
}



simulation_pair <- function(theta,max_iter,max_sim,tol,gamma,lambda){
  lista_sim <- as.list(seq(1, max_sim))
  lista_tiempos<- list()
  lista_error_X <- list()
  lista_error_Y <- list()
  lista_error <- list()
  lista_iter_sim <- list()
  for(k in 1:max_sim){
    t <- proc.time()
    min_est <- SGD_2D(h,theta,max_iter,tol,gamma,lambda)
    tiempo <- proc.time() - t
    theta_test<- min_est$theta
    error_X<- abs(theta_test[1])
    error_Y<- abs(theta_test[2])
    error <- sqrt(theta_test[1]^2+theta_test[2]^2)
    lista_tiempos<- append(lista_tiempos, tiempo)
    lista_error_X <-append(lista_error_X, error_X)
    lista_error_Y <- append(lista_error_Y, error_Y)
    lista_error <- append(lista_error,error)
    lista_iter_sim <- append(lista_iter_sim, min_est$lista_iter[-1])
  }
  return(list(sim = lista_sim, tiempos = lista_tiempos, error_x = lista_error_X,error_Y = lista_error_Y, error = lista_error))
  
  
}

theta <- c(0.65,0.80)
max_iter <- 100
tol<-0.1
#ttt <- SGD_2D(h,theta,max_iter,tol,gamma,lambda)
max_sim = 1000
Resultado <- simulation_pair(theta,max_iter,max_sim,tol,gamma,lambda)

exito <- Resultado$error[sapply(Resultado$error, function(x) x < tol)]
x<- as.list(seq(1,length(exito)))

plot(x,exito, main = "Casos exitosos del algoritmo", xlab = "Casos", ylab = "Error")
print(length(exito)/max_sim)
