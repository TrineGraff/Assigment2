library(matrixcalc)
norm2 = function(x) norm(as.matrix(x), type="2")

#BACKTRACKING LINE SEARCH
#Here will a acceptable step lenth be found after a finite number of trials, becauce a_k will become small enough that the sufficient decrease hold
backtracking_line_search = function(x_k, p_k, g_k, f_k, f)
{
  #Intialiasation
  a_bar = 1  #Always one in Newton
  p = 0.5
  c = 1e-4
  a_k = a_bar 
  directional_derivative = drop(crossprod(g_k, p_k))
  k = 0L; k_max = 1000L
  done = FALSE
  
  #Iteration
  while(!done)
  {
    k = k + 1L
    a_k = p * a_k #a_{k+1}=< a_k
    #done when suffiecient decrease hold or k >=k_max; 
    done = (f(x_k + a_k * p_k) < (f_k + c * a_k * directional_derivative)) | (k >= k_max)
  }
  a_k
}

f = function(x) 100*(x[2]-x[1]^2)^2+(1-x[1])^2
g = function(x) {
  c(2 * (200 * x[1]^3 - 200 * x[1] * x[2] + x[1] - 1),
    200 * (x[2] - x[1]^2))
}
h = function(x) rbind(c(1200*x[1]^2-400*x[2]+2,-400*x[1]),c(-400*x[1],200))

Newton_Method = function(f, g, h, x_k, g_tolerance = 1e-5, k_max = 1000L)
{
  #INTIALISATION
  k = 0L;
  f_k = f(x_k)
  done = FALSE
  
  #ITERATION
  while (!done)
  {
    k = k + 1L
    g_k = g(x_k)
    p_kN = - solve(h(x_k)) %*% g_k  
    directional_derivative = crossprod(p_kN, g_k) 
    is_descent_direction = directional_derivative < 0 #gauarantess that the function f can be reduced along this direction
    p_k = if(is_descent_direction) p_kN else -g_k
    a_k = backtracking_line_search(x_k, p_k, g_k, f_k, f)
    s_k = a_k * p_k
    #Updates, new value of x_k
    x_k = x_k + s_k
    f_k = f(x_k)
    #When the norm of g is near zero we have that x_k=x* 
    magnitude_g = norm2(g_k) 
    done = (k >= k_max) || (magnitude_g < g_tolerance) 
    #Output the result
    cat("k = ", k, " ")
    cat("f = ", f_k, " ")
    cat("|g_k|", magnitude_g, " " )
    cat("x=", x_k, " ")
    cat("a = ", a_k, "\n ")
  }
}

Newton_Method(f, g, h, c(1.2,1.2))
Newton_Method(f, g, h, c(-1.2,1))

