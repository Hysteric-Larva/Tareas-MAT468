# Homework 1, stocastic simulation

#Matrix and vector definition:
A <- matrix(c(1, 1, 1, 2, 3, 1, 1, -1, -1), nrow = 3, byrow = TRUE)

b <- c(4, 9, -1)

B <- matrix(c(30, 16, 46, 16, 10, 26, 46, 26, 72), nrow = 3, byrow = TRUE)

#Sección LU
# Función para factorizar una matriz A en LU
lu_factorization <- function(A) {
  n <- nrow(A)
  L <- diag(1, n)  # Inicializamos L como matriz identidad
  U <- A           # U inicialmente es igual a A
  
  for (j in 1:(n-1)) {
    for (i in (j+1):n) {
      L[i, j] <- U[i, j] / U[j, j]  # Calculamos el factor de L
      U[i, j:n] <- U[i, j:n] - L[i, j] * U[j, j:n]  # Se actualiza U para que sea triangular superior
    }
  }
  
  return(list(L = L, U = U))
}

forward_substitution <- function(L, b) {
  n <- nrow(L)
  y <- numeric(n)
  y[1] <- b[1]
  for (i in 2:n) {
    y[i] <- b[i] - sum(L[i, 1:(i-1)] * y[1:(i-1)])
  }
  return(y)
}

backward_substitution <- function(U, y) {
  n <- nrow(U)
  x <- numeric(n)
  x[n] <- y[n] / U[n, n]
  for (i in (n-1):1) {
    x[i] <- (y[i] - sum(U[i, (i+1):n] * x[(i+1):n])) / U[i, i]
  }
  return(x)
}





lu_solver <- function(A,b){
result <- lu_factorization(A)
L <- result$L
U <- result$U
# Resolvemos Ly = b
y <- forward_substitution(L, b)

# Resolvemos Ux = y
x <- backward_substitution(U, y)
return(x)
}
## Prueba
LU <- lu_factorization(A)

# Mostrar las matrices L y U
L <- LU$L
U <- LU$U

###### FIN DE LA PARTE LU    ########

#### INICIO ITERATIVO ###

refinamiento_iterativo <- function(A, b, x0, tau, kmax) {
  k <- 0
  x <- x0  # Inicializamos x con la solución inicial
  
  # Cálculo de la norma infinita
  norm_inf <- function(v) {
    return(max(abs(v)))
  }
  
  # Bucle infinito
  while (TRUE) {
    # Paso 3: Calcular el residuo r(k)
    r_k <- b - A %*% x
    
    # Paso 4: Resolver Aδ = r(k) para obtener δ(k)
    delta_k <- solve(A, r_k)  # Usamos la función solve para resolver el sistema
    
    # Paso 5: Actualizar x
    x_new <- x + delta_k
    
    # Paso 6: Verificar la condición de parada
    if (norm_inf(delta_k) <= tau * norm_inf(x_new)) {
      return(x_new)  # Paso 7: Retornar la solución
    } else if (k < kmax) {
      # Paso 9: Actualizar k y x
      k <- k + 1
      x <- x_new
    } else {
      # Paso 12: Indicar que el algoritmo no converge
      stop("El algoritmo no converge después de kmax iteraciones.")
    }
  }
}

#################    EJERCICIO 1     ########################
x_lu <- lu_solver(A,b)
x_ref <- refinamiento_iterativo(A, b, rep(0, nrow(A)), 0.0005, 500)

###### Ejercicio 2:     OPERADOR SWEEP       ####################

sweep_operator <- function(A, k) {
  n <- nrow(A)
  # Verificamos que el elemento diagonal no sea cero
  if (A[k, k] == 0) {
    stop("El elemento diagonal a_{kk} debe ser diferente de cero.")
  }
  
  # Inicializamos la matriz B
  B <- matrix(0, n, n)
  
  # Calculamos bkk, bik, bkj y bij según la definición
  akk <- 1/A[k, k]
  
  # Elemento diagonal
  B[k, k] <- akk
  
  # Ajuste de filas
  for (i in 1:n) {
    if (i != k) {
      B[i, k] <- -A[i, k] * akk  # bik
      B[k, i] <- A[k, i] * akk  # bki
    for (j in 1:n) {
    if (j != k ) {
        B[i, j] <- A[i, j] - (A[i, k] * A[k, j]) * akk  # bij
      
    }
  }
}
  
}
return(B)
}






# Función para calcular la inversa usando el operador sweep
inverse_sweep <- function(A) {
  n <- nrow(A)
  
  for (k in 1:n) {
    # Aplicamos el operador sweep en la k-ésima posición
    A <- sweep_operator(A, k)
  }
  
  return(A)
}
#Ejercicio 3
resultado <-inverse_sweep(B)


# Ejercicio 4
res <- eigen(B)
print(res$values)

# Definir la función de factorización de Cholesky
cholesky_decomposition <- function(A) {
  if (nrow(A) != ncol(A)) {
    stop("La matriz debe ser cuadrada.")
  }
  
  if (!all(A == t(A))) {
    stop("La matriz debe ser simétrica.")
  }
  
  n <- nrow(A)
  L <- matrix(0, n, n)
  
  for (i in 1:n) {
    for (j in i:n) {
      if (i == j) {
        # Cálculo de la diagonal
        L[i, j] <- sqrt(A[i, i] - sum(L[i, 1:(i-1)]^2))
        
      } else {
        # Cálculo de los elementos no diagonales
        L[j, i] <- (A[j, i] - sum(L[j, 1:(i-1)] * L[i, 1:(i-1)])) / L[i, i]
        print(L[j,i])
      }
    }
  }
 # L[n, n] <- sqrt(A[n, n] - sum(L[n, 1:(n-1)]^2))
  
  return(L)
}

# Prueba

C = matrix(c(5,2,3,2,7,-4,3,-4,9), nrow = 3, byrow = TRUE)
L<- cholesky_decomposition(C)
# Ejercicio 4
LL<- cholesky_decomposition(B)
print(kappa(LL))