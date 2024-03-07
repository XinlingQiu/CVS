#initialization
library('abind')
seed=sample(100000,1)
set.seed(seed)
n <- 1000#the number of sample

#parameter of the Gaussian distribution 1(G1)
alpha1 <- 0.4
miu1   <- 3
sigma1 <- 3

#parameter of the Gaussian distribution 2(G2)
alpha2 <- 0.6
miu2   <- -4
sigma2 <- 2


kk <- 6#the number of parameters
n1 <- floor(n*alpha1)#the number of sample G1
n2 <- n-n1##the number of sample G2

#Randomly generate samples
samp <-numeric(n)
samp[1:n1] <- rnorm(n1, miu1, sigma1)
samp[(n1+1):n] <- rnorm(n2, miu2, sigma2)


#The inverse of a vector x
INV<-function(x){
  return(x/sum(x^2))
}
#m(x,y),measure the distance between x and y
metrics<-function(x,y){
  return(log((sum((x-y)^2)),10))
}

#The EM algorithm
em<-function(theta){
  alpha=theta[1:(kk/3)]
  mu=theta[(kk/3+1):(2*kk/3)]
  sigma=theta[(2*kk/3+1):kk]
  
  prob <- matrix(rep(0, kk/3*n), nrow = n)
  weight <- matrix(rep(0, kk/3*n), nrow = n)
  
  # E-step
  for (i in 1:(kk/3)) {
    prob[, i]   <- sapply(samp, dnorm, mu[i], sigma[i])
    weight[, i] <- alpha[i] * prob[, i]
  }
  row_sum <- rowSums(weight)
  prob    <- weight/row_sum
  
  
  # M-step
  for (j in 1:(kk/3)) {
    sum1     <- sum(prob[, j])
    sum2     <- sum(samp*prob[, j])
    alpha[j] <- sum1/n
    mu[j]   <- sum2/sum1
    sum3     <- sum(prob[, j]*(samp-mu[j])^2)
    sigma[j] <- sqrt(sum3/sum1)
  }
  theta<-c(alpha,mu,sigma)
  return(theta)
}


#The epsilon-step of the epsilon-accelerated EM algorithm
eps_step<-function(theta1,theta2,theta3){
  a=theta3-theta2
  b=theta2-theta1
  return(theta2+INV(INV(a)-INV(b)))
}

#The vector epsilon-algorithm
eps<-function(theta1,theta2,theta3){
  return(theta1+INV(theta3-theta2))
}


#CVS
CVS<-function(init_theta,row,tau1,tau2){
  theta<-matrix(0,row,kk)
  step=1
  theta[1,]=init_theta
  theta[2,]=em(init_theta)
  start_time=proc.time()[1]
  while(metrics(theta[step,],theta[step+1,])>tau1){
    step=step+1
    theta[step+1,]=em(theta[step,])
  }
  row=step+1
  t<-matrix(0,row,kk)
  row=row-1
  while(row>1){
    for(k in 1:row){
      t[k,]=t[k+1,]+INV(theta[k+1,]-theta[k,])
    }
    row=row-1
    if(row>0){
      for(k in 1:row){
        theta[k,]=theta[k+1,]+INV(t[k+1,]-t[k,])
        if(k>1 && metrics(theta[k,],theta[k-1,])<tau2){
          time=proc.time()[1]-start_time
          return(c(time,step))#Output parameters that first meet convergence conditions,you can also return theta[k,]
        }
      }
      row=row-1
    }
  }
  time=proc.time()[1]-start_time
  return(c(time,step))#or return the last updated parameter,you can also return theta[1,]
}
############################### table 1 #####################################
#initialization
tau2=-20
init_alpha<-runif(kk/3)
init_alpha<-init_alpha/sum(init_alpha)
init_mu<-runif(kk/3,min=-10,max=10)
init_sigma<-runif(kk/3,min=0,max=10)
init_theta<-c(init_alpha,init_mu,init_sigma)



#Using EM algorithm to obtain the limit value due to the local convergence of EM
theta1<-init_theta
theta2<-em(theta1)
while(metrics(theta1,theta2)>tau2){
  theta1=theta2
  theta2=em(theta1)
}
#The last one as a limit
limit_theta=theta2

#The EM algorithm
len=100#the length of EM sequence
theta<-matrix(0,len,kk)
theta[1,]<-init_theta
for(step in 2:len){
  theta[step,]=em(theta[step-1,])
}
View(theta)

col=ceiling(len/2)##the number of columns in CVA
theta_CVA=array(rep(0,kk*len*(len+1)),dim=c(len,kk,len+1))
theta_CVA[,,1]=0
theta_CVA[,,2]=theta
metrics_matrix<-matrix(0,len,col)
for(m in 3:(len+1)){
  for(j in 1:(len-m+2)){
    theta_CVA[j,,m]=eps(theta_CVA[j+1,,m-2],theta_CVA[j,,m-1],theta_CVA[j+1,,m-1])
  }
}
for(j in 1:col){
  for(i in 1:(len-2*(j-1))){
    metrics_matrix[i,j]=metrics(theta_CVA[i,,2*j],limit_theta)
  }
}

View(metrics_matrix)


############################### table 1 end#####################################

############################### table 2##################################
num_test=10
tau1=-2
tau2=-5
iternum<-matrix(0,num_test,3)
time<-matrix(0,num_test,3)
for(i in 1:num_test){
  init_alpha<-runif(kk/3)
  init_alpha<-init_alpha/sum(init_alpha)
  init_mu<-runif(kk/3,min=-10,max=10)
  init_sigma<-runif(kk/3,min=0,max=10)
  init_theta<-c(init_alpha,init_mu,init_sigma)
  ############################# em #################################################
  len=1000
  theta_em<-matrix(0,len,kk)
  theta_em[1,]<-init_theta
  step=1
  start_time=proc.time()[1]
  while(step<len){
    theta_em[step+1,]=em(theta_em[step,])
    if(metrics(theta_em[step,],theta_em[step+1,])<tau2){
      break
    }
    step=step+1
  }
  time[i,1]=proc.time()[1]-start_time
  iternum[i,1]=step
  len=step+1#the length of EM sequence
  ############################# epsilon-accelerated EM sequence #######################################
  theta_eps_acc<-matrix(0,len,kk)
  theta_eps_acc[1,]<-init_theta
  theta_eps_acc[2,]<-em(init_theta)
  step=2
  start_time=proc.time()[1]
  while(step<len){
    theta_eps_acc[step+1,]=em(theta_eps_acc[step,])
    theta_eps_acc[step-1,]=eps_step(theta_eps_acc[step-1,],theta_eps_acc[step,],theta_eps_acc[step+1,])
    if(step>=3&&(metrics(theta_eps_acc[step-1,],theta_eps_acc[step-2,])<tau2)){
      break
    }
    step=step+1
  }
  time[i,2]=proc.time()[1]-start_time
  iternum[i,2]=step
  ###############################################CVS#################################################
  return_value=CVS(init_theta,len,tau1,tau2)
  time[i,3]=return_value[1]
  iternum[i,3]=return_value[2]
}
summary(iternum)
summary(time)


####################################table 2 end #########################################