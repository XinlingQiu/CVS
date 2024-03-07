#initialization
library("MASS")
library("abind")
seed=sample(100000,1)
set.seed(seed)
N=100
p=10
kk=p+1
X<-matrix(0,N,kk)
X[,1]=1
S <- toeplitz((p:1)/p)
R <- rWishart(1, p, S)
R=R[,,1]
miu=matrix(0,p,1)
X[,2:(p+1)]=mvrnorm(N,miu,R)
inv=solve(t(X)%*%X)
Y=X%*%rnorm(kk,0,1)
def=sample(1:N,N/2)
Y[def]=NA



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
  theta=matrix(theta,kk,1)
  #E step
  for(i in def){
    Y[i]=X[i,]%*%theta
  }
  ##M step
  theta=inv%*%t(X)%*%Y
  theta=matrix(theta,1,kk)
  theta
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
init_theta<-rnorm(kk,0,1)

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
len=50#the length of EM sequence
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