# Step 0: restoring the workplace and importing required packages ------------------------------------
library(readxl)
library(MASS)
library(lpSolve)
library(pracma)
library(tidyverse)
library(ggplot2)

# Step 1: Data generation ------------------------------------------------------
set.seed(123) # for reproducibility 
m = 4 # number of financial assets
N = 5000 # number of simulations 
points = 10 # for grid that will be used in the EF
alpha = 0.05 # expected shortfall's coverage level 

# Generate simulated returns:
# Note: for simplicity, we simulate Gaussian returns.
mu <- c(0.053,0.141,-0.053,-0.001)
sigma <- matrix(c(1.34,1.34,1.84,0.07,
                  1.34,4.43,1.28,0.15,
                  1.84,1.28,11.3,0.20,
                  0.07,0.15,0.20,3.45),4)
ret <- mvrnorm(n = N, mu, sigma)

# Step 2: Efficient Frontier using LP ----------------------------------------------
# Linear Programming

# Set coefficients of the objective function:
f.obj <- matrix(c(1/(alpha*N)*matrix(rep(1, N)),  1,  matrix(rep(0, m))),1)

# Set matrix corresponding to coefficients of constraints by rows
# Do not consider the non-negative constraint; it is automatically assumed
ineqconst <- cbind(matrix(diag(N),nrow = N), matrix(rep(1, N)), matrix(ret,N)) 
eq1 <- cbind(matrix(rep(0, N+1),nrow = 1), matrix(mu       , 1) )
eq2 <- cbind(matrix(rep(0, N+1),nrow = 1), matrix(rep(1, m), 1) )
f.con <- rbind(ineqconst, eq1, eq2)

# Constraints
f.dir <-  cbind(matrix(rep(">=", N+1),1), matrix(rep("=", 1),1))

# Set right hand side coefficients
C = min(mu)
f.rhs <-  cbind(matrix(rep(0, N),1), matrix(c(C,1),1) )

# Global minimum risk portfolio:
start_time <- Sys.time()
ProblemSolution <- lp("min", f.obj, f.con, f.dir, f.rhs) 
end_time <- Sys.time()
timeEF_LP <- end_time - start_time 
print(paste0("LP lasts ",round(timeEF_LP,2), " seconds finding one solution"))

ES <- - ProblemSolution$objval 

# Variables final values
x1 <- ProblemSolution$solution
w <- tail(x1, n=m) # corresponding weights for the minimum E(r) portfolio
ExpRet <- as.numeric(w%*%mu)

gridER <- linspace(ExpRet, max(mu)-0.0001, points)

for (i in 2:points) {
  # Set right hand side coefficients
  C <- gridER[i]
  f.rhs <-  cbind(matrix(rep(0, N),1), matrix(c(C,1),1) )
  ProblemSolution <- lp("min", f.obj, f.con, f.dir, f.rhs) 
  ES <- rbind(ES,-ProblemSolution$objval)
  x1 <- ProblemSolution$solution
  w <- rbind(w,tail(x1, n=m))
  ExpRet <- rbind(ExpRet,w[i,]%*%mu) 
  print(i)
}

# Efficient Frontier's plot:
EFdf_LP <- data.frame(matrix(c(-ES,ExpRet),ncol = 2))
colnames(EFdf_LP) <-  c("ES", "Er")

# Plot:
win.graph(width=8,height=6,pointsize=8)
EFplot <- ggplot(data=EFdf_LP, aes(x = ES, y = Er)) +
  geom_line(linetype = "dashed", color="red")+
  geom_point(size=2)+
  labs(title = "Efficient Frontier",
       subtitle = "Linnear Programming",
       x = "-ES",
       y = "E(r)")
EFplot + theme_bw(base_line_size = 0.3)

# Step 3: Efficient Frontier using NLP -----------------------------------------------------------------
# Non-Linear Formulation (incorporating auxiliar v scalar). 
rm(w,ExpRet,ES)

# Starting values: 
v0 <- -0.1
w0 <- rep(1,m)*(1/m)  
x0 <- matrix(c(v0,w0),1)

ES_RockafellarUryasev2000 <- function(x0, r_hat, alpha) {
  m <- ncol(r_hat)
  v <- x0[1]
  w <- matrix(x0[1:m+1]) 
  rp <- t(w)%*%t(r_hat)
  di <- pmax(0, v-rp)
  ES <- - (1/alpha)*mean(di) + v
  ES <- - ES
  return(ES)
}

# Optimization restrictions:
lb <- c(-100, rep(0,m)) # lower bounds for x = (v,w1, ... ,wm)
ub <- c( 100, rep(1,m))  # upper bounds for x = (v,w1, ... ,wm)
A  <- - matrix(c(0, mu),1) # portfolio expected return must be at least as greater as the  minimum marginal expected return.
b  <- - min(mu)  
Aeq  <- matrix(c(0, rep(1,m)),1) # equality restriction: weights must add to 1.
beq   <- 1 

# Global minimum risk portfolio:
start_time <- Sys.time()
x1 <- fmincon(x0=x0,r_hat=ret,alpha=alpha,ES_RockafellarUryasev2000,ub=ub,lb=lb,A=A,b=b,Aeq=Aeq,beq=beq,tol = 1e-05)
end_time <- Sys.time()
time_NLP = end_time - start_time 
print(paste0("NLP lasts ",round(time_NLP,2), " seconds finding one solution"))

# GMR expected shortfall:  
ES <- - x1$value 

# GMR weigths:
xw <- x1$par
w <- tail(xw, m) # corresponding weights for the minimum E(r) portfolio

# Expected return associated to the GMR portfolio:
ExpRet <- as.numeric(w%*%mu)

# Efficient frontier:
gridER <- - linspace(ExpRet, max(mu)-0.0001, points)

for (i in 2 : points) {
  # Set right hand side coefficients
  b  <- gridER[i]
  x1 <- fmincon(x0=x0,r_hat=ret,alpha=alpha,ES_RockafellarUryasev2000,ub=ub,lb=lb,A=A,b=b,Aeq=Aeq,beq=beq,tol=1e-05)
  ES <- rbind(ES,-x1$value)
  xw <- x1$par
  w <- rbind(w,tail(xw, m))
  ExpRet <- rbind(ExpRet,w[i, ]%*%mu) 
  print(i)
}

wdf_NLP <- data.frame(w)
EFdf_NLP <- data.frame(matrix(c(-ES,ExpRet),ncol = 2))
colnames(EFdf_NLP) <-  c("ES", "Er")

win.graph(width=8,height=6,pointsize=8)
EFplot <- ggplot(data=EFdf_NLP, aes(x = ES, y = Er)) + 
  geom_line(linetype = "dashed", color="blue")+
  geom_point(size=2)+
  labs(title = "Efficient Frontier",
       subtitle = "Non-Linnear Programming",
       #caption = " ",
       #tag = "Figure 1",
       x = "-ES",
       y = "E(r)")
EFplot + theme_bw(base_line_size = 0.3)

# Step 4: Comparison between LP and NLP  ----------------------------------
win.graph(width=8,height=6,pointsize=8)
colors <- c("LP" = "red", "NLP" = "blue")
both <- ggplot(NULL, aes(x = ES)) + 
  geom_point(data = EFdf_LP, aes(y = Er, color="LP"), size = 3) +
  geom_point(data = EFdf_NLP, aes(y = Er, color="NLP"), size = 1.5) +
  labs(title = "Efficient Frontier", color = "Method", 
       x = "-ES",
       y = "E(r)") +
  scale_color_manual(values = colors)
both +  theme_bw(base_line_size = 0.1)
