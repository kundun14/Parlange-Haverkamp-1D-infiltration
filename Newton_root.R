# implementacion del metodo de newton para hallar
# la infiltraiocn acumulada segun el modelo simplificado de haverkamp 1994
# los parametros fueron estimados con modelacion inversa con dream y 
# data de campo de un ensayo de infiltracion de doble anillo

library(numDeriv)

# implementacion newton- raphson

newton = function(fun, x_0, iter = 1000, tolerancia = 1e-5){
  
  require(numDeriv) #paquete de derivacion numerica
  
  #x_0 = x_min   #valor raiz tentativo de inicio de iteracion 
  k = iter # inicializar el numero de iteraciones, aqui 1000
  
  for (i in 1:iter) {
    
    # para cada iteracion hallar una solucion, actulizar hasta que la diferencia sea menor a la toler
    
    dx <- numDeriv::genD(func = fun, x = x_0)$D[1] # primera derivada en x_0
    x_1 <- x_0 - (fun(x_0)/dx)
    k[i] = x_1 # guardar raiz actual 
    
    # si la diferencia entre x actual y x anterior es menor a tolerancia
    # ese es el resultado
    
    if (abs(x_1 - x_0) < tolerancia ){
      
      root <- tail(k, n = 1)
      #res <- list("Raiz_aprox" = root, "Iteraciones" = k ) # lista de resultado e iteraciones 
      return(root)
      
    }
    
    x_0 <- x_1
    
  }
  
  print("No ocurrio convergencia")
  
}


# ecuacion de haverkamp de infiltracion


haver <- function(t,x){
  #parametros 
  
  S <- x[1]
  K_s <- x[2]
  beta <- x[3]
  
  if ( length(x) == 4) {
    
    K_i <- x[4]
    
  } else {K_i <- 0
  
  }
  
  dK <- K_s - K_i
  n <- length(t)
  I <- numeric(length = n)
  I_0 <- 0.2  # valor inicial para hallar la raiz
  
  # definicion de la ecuacion por partes
  
  lhs <- function(tj, dK, S, beta) (dK/S)^2 * (1 - beta) * tj
  
  rhs <- function(I, dk, S, beta) (K_s * I) / S^2 - 
                                                    (1/2) * (
                                                      
                                                      log( 
                                                        
                                                        (1/ beta) * exp (
                                                          
                                                          (2 * beta * K_s * I )/ S^2
                                                          
                                                        ) +
                                                          
                                                          (beta - 1) / beta
                                                        
                                                      )
                                                    )
                                                  
  # I implicita en dhs
  # debemos resolver para I
  
  # parte dinamica , calcular I en cada timestep
  
  for (j in 1:n){
    
    L = lhs(t[j], dK, S, beta)
    
    if ( j == 1 && t[j] == 0) { # infiltracion acumulada en el tiempo cero
      
      I[j] <- 0
      
    } else  {
      
      I[j] <- newton((function(z) rhs(z, dk, S, beta) - L), I_0) 
      
      # inicializa con I_0 y 
      # z reemplaza a I en una funcion explicita
      
    }
    
  } # cierra for
  
  return(I)
  
} # cierra la funcion


# parametros S, K_s y beta ( tiempo en horas y longitud en cm), beta es adimensional
# para clay segun paper pag 

x <- c(1.0308,0.2101,1.5219)  

# tiempos en horas
t <- seq(0, 2, 0.1) 

#maxiter <- 20 # implmentado directamente en funcion newton

I <- haver(t, x)
I

plot(x = t, y = I) 





