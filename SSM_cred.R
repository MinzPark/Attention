#[1] Data-----------------------------------------------------------------------------------------------------------
#### Data 

train.x <- read.csv("https://raw.githubusercontent.com/MinzPark/Attention/main/test_train_data/train_x.csv")
train.y <- read.csv("https://raw.githubusercontent.com/MinzPark/Attention/main/test_train_data/train_y.csv")
train.x <- train.x[2:9]
train.y <- train.y[2:6]

test.x <- read.csv("https://raw.githubusercontent.com/MinzPark/Attention/main/test_train_data/test_x.csv")
test.y <- read.csv("https://raw.githubusercontent.com/MinzPark/Attention/main/test_train_data/test_y.csv")
test.x <- test.x[2:9]
test.y <- test.y[3:8]


head(train.x)
head(test.x)
head(train.y)
head(test.y)

X <- model.matrix(~., data = train.x[,2:8])
dim(X)
head(X)


#[2]모수추정-----------------------------------------------------------------------------------------------------------
#### Bayesian_state space model :모수 추정 
library(nimble)
library(MCMCvis)
library(MCMCvis)


model0 <- nimbleCode({
  # likelihood
  for(i in 1:n){
    lambda[i] <-exp(inprod(X[i,1:k], beta[1:k]))
    for(j in 1:tt){
      y[i,j] ~ dpois(exp(R[i,j]) * lambda[i])
    }
  }
  # prior
  for (i in 1:k){
    beta[i] ~ dnorm(mean=0, sd=1)
  }
  
  # random effect
  for(i in 1:n){
    R[i,1] ~ dnorm(mean=0, sd= sig_e/sqrt((1-phi^2)))
    for(t in 2:tt){
      R[i,t]~ dnorm(mean=phi*R[i,t-1], sd=sig_e)
    }
  }
  
  
  sig_e ~ dgamma(1,1)
  #tau_r ~ dgamma(shape=0.01,scale=1/0.01)
  phi ~ dbeta(1,1)
  
})



my.data <- list(X = X, y = as.matrix(train.y))
my.constants <- list(n=497, k=8, tt=5)

initial.values <- function() list(beta = rnorm(n= 8, mean = 0, sd=1),
                                  sig_e = rgamma(1,1,1),
                                  phi = rbeta(1, 1,1))

# model0$initializeInfo()

parameters.to.save <- c("beta","phi","sig_e")

### Run MCMC
n.iter <- 30000
n.burnin <-10000
n.chains <- 2
n.thin<-5

mcmc.output <- nimbleMCMC(code = model0,
                          constants = my.constants,
                          data = my.data,
                          inits = initial.values,
                          monitors = parameters.to.save,
                          niter = n.iter,
                          nburnin = n.burnin,
                          nchains = n.chains, 
                          thin = n.thin)



MCMCplot(object = mcmc.output, 
         params = 'beta')

MCMCtrace(object = mcmc.output,
          pdf = TRUE, # no export to PDF
          ind = TRUE, # separate density lines per chain
          params = "beta")

MCMCtrace(object = mcmc.output,
          pdf = TRUE, # no export to PDF
          ind = TRUE, # separate density lines per chain
          params = "phi")


MCMCtrace(object = mcmc.output,
          pdf = TRUE, # no export to PDF
          ind = TRUE, # separate density lines per chain
          params = "sig_e")

result <- MCMCsummary(object = mcmc.output, round = 2)
result

#Credibility-----------------------------------------------------------------------------------------------------------
#### Credibility
#data
XX<-model.matrix(~., data = test.x[2:8])
NN<-test.y

Beta.hat=c(-3.11,  0.30,  1.62, 0.28, -1.61,-0.22,1.30,2.25)
Phi.hat=0.96
sig_e.hat = 0.25
Sigsq_r.hat=sig_e.hat^2
Sigsq.hat=Sigsq_r.hat/(1-Phi.hat^2)

start_time<-Sys.time()    #시작시간

#variance, covariance function

cov_f <- function(t1,t2,a,b,lambda,sigsq,phi){
  z=abs(t1-t2)
  t=ifelse(t1<=t2,t1,t2)
  if(t1!=t2){
    temp1=exp((1+phi^z)*phi^(t-1)*a+(1+phi^z)^2*phi^((t-1)*2)*b/2+(1+phi^z)^2*(1-phi^((t-1)*2))/(1-phi^2)*sigsq/2)
    temp2=exp((1-phi^(z*2))/(1-phi^2)*sigsq/2)
    temp3=exp(phi^(t2-1)*a+phi^(2*(t2-1))*b/2+(1-phi^((t2-1)*2))/(1-phi^2)*sigsq/2)
    temp4=exp(phi^(t1-1)*a+phi^(2*(t1-1))*b/2+(1-phi^((t1-1)*2))/(1-phi^2)*sigsq/2)
    ret=lambda^2*(temp1*temp2-temp3*temp4)
  }else{
    temp5=exp(phi^(t1-1)*a+phi^((t1-1)*2)*b/2+(1-phi^((t1-1)*2))/(1-phi^2)*sigsq/2)
    temp6=exp(2*phi^(t1-1)*a+2*phi^((t1-1)*2)*b+2*(1-phi^((t1-1)*2))/(1-phi^2)*sigsq)
    temp7=exp(2*phi^(t1-1)*a+phi^((t1-1)*2)*b+(1-phi^((t1-1)*2))/(1-phi^2)*sigsq)
    ret=lambda*temp5+lambda^2*(temp6-temp7)
  }
  ret
}

#variance, covariance matrix
cov_M <- function(t,t_pred,a,b,lambda,sigsq,phi){
  #length(t) X length(t)
  n=length(t)
  M=matrix(NA,n,n)
  for(i in 1:n){
    for(j in 1:n){
      M[i,j]=cov_f(t[i],t[j],a,b,lambda,sigsq,phi)
    }
  }
  M
}

#covariance vector
cov_CV <- function(t,t_pred,a,b,lambda,sigsq,phi){
  #length(t) X 1
  n=length(t)
  Cov=matrix(NA,n,1)
  for(i in 1:n){
    Cov[i,1]=cov_f(t[i],t_pred,a,b,lambda,sigsq,phi)
  }
  Cov
}

#alpha_hat 함수 
alpha_f <- function(t,t_pred,a,b,lambda,sigsq,phi){
  En=matrix(NA,1,t_pred)
  for(i in 1:t_pred){
    En[i]<-lambda*exp(phi^(i-1)*a+phi^((i-1)*2)*b/2+(1-phi^((i-1)*2))/(1-phi^2)*sigsq/2)
  }
  CovM<-cov_M(t,t_pred,a,b,lambda,sigsq,phi)
  CovC<-cov_CV(t,t_pred,a,b,lambda,sigsq,phi)
  alphas=solve(CovM)%*%CovC[,1]
  alpha0=En[t_pred]-sum(En[t]%*%alphas)
  (ret=c(alpha0,alphas))
}

#credibility 함수
predict_ssm<-function(ns,t,t_pred,a,b,lambda,sigsq,phi){
  alpha_hat<-alpha_f(t,t_pred,a,b,lambda,sigsq,phi)
  (ret=alpha_hat[1]+sum(ns*alpha_hat[-1]))
}

#credibility 계산
nr=dim(test.y)[1]  #1000
nc=dim(test.y)[2]-1  #5

bul=matrix(NA,nrow=nr,ncol=(nc+1))
for (i in 1:nr){
  lamb_hat=exp(XX[i,1:8] %*% Beta.hat)
  a=0
  b=Sigsq_r.hat/(1-Phi.hat^2)
  for (j in 1:nc){
    ns=NN[i,1:j]
    bul[i,(j+1)]=predict_ssm(ns=ns, t=c(1:j), t_pred=c(j+1), a=a, b=b, lambda=lamb_hat, sigsq=Sigsq_r.hat, phi=Phi.hat)
  }
}
   
cred_mse=c()
for (t in 1:(nc+1)){
  cred_mse[t]=mean((NN[,t] - bul[,t])^2)
}

end_time <-Sys.time()   #종료시간
end_time - start_time   #Time difference of 9.117023 secs

head(bul,5)    
tail(bul,5) 

print(cred_mse,7)

lamb_hat=exp(XX[i,1:8] %*% Beta.hat)
cov_M(1,2,a,b,lambda = lamb_hat,sigsq = Sigsq_r.hat,phi = Phi.hat)

