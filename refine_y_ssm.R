
#####################################################
# Load data
#####################################################

#setwd("C:/Users/none/Desktop/대학원/논문/data")

# setwd("C:/Users/parkminji/Downloads/민지/data")
# source("repeated_data_management.R")

#setwd("C:/Users/none/Desktop/대학원/논문/data")

path <- "C:/Users/none/Desktop/대학원/논문" # Download path 

setwd(path)

path_ftn <- paste0(path, "/Unbiased insurance preminum")

#setwd("C:/Users/parkminji/Downloads/민지/data")
source(paste0(path,"/data/repeated_data_management.R"))
source(paste0(path_ftn,"/lamb_kappa.R"))
summary(data.train)

dim(data.train)
dim(data.valid)

#####################################################
# Import package
#####################################################
library(tidyverse)
library(nimble)
library(MCMCvis)
library(plotrix)
library(dplyr)


#---------------------------------------------------
# find parameters :  betas, phi, sigma
#---------------------------------------------------

#####################################################
# data preprocessing
#####################################################
# (1) uniq_id for PolicyNum

id_uniq = unique(data.train$PolicyNum)

length(id_uniq) # 497

# (2) idx per id per PolicyNum's Year since each ID has different year of frequency (i.e not the same observations, 5,5,4,...)
id_idx=c(NULL)

for(i in 1:length(id_uniq)){
  ind = id_uniq[i] == data.train$PolicyNum 
  id_idx=c(id_idx, rep(i, sum(ind)))
}
Year = data.train$Year-2005 # Year[j]

head(id_idx, 20)    
tail(id_idx, 20) 
length(id_idx)


data.train <- combine_kappa_level(data.train)
data.valid <- combine_kappa_level(data.valid)

#####################################################
# estimate a.hat, lambda.hat, R, beta.hat using nimble
#####################################################


round(apply(data.train, 2, mean)*100,2) 

# Make design matrix 
X = model.matrix(n~TypeCity+TypeCounty+TypeSchool+TypeTown+TypeVillage+C2+C3, data = data.train)
#X = model.matrix(n~TypeCity+TypeCounty+TypeMisc+TypeTown+TypeSchool+C2+C3, data = data.train)


model.param <- nimbleCode({
  # likelihood
  for(j in 1:JJ){
    lambda[j] <- exp(inprod(X[j,1:KK], beta.hat[1:KK]))
    Y[j]~dpois(lambda[j]*exp(R[id_idx[j], Year[j]]))
    
  }
  
  # prior
  for(i in 1:KK){
    beta.hat[i] ~ dnorm(0,1)
  }
  
  # random effect
  for(n in 1:N.id_idx){
    R[n, 1]~dnorm(0, sigma_e/sqrt(1-phi^2))
    for(t in 2:TT) {
      R[n, t]~dnorm(phi * R[n, t-1], sigma_e)
    }
    R6[n]~dnorm(phi * R[n, TT], sigma_e)
  }
  
  # prior of param
  phi ~ dbeta(1,1)
  sigma_e ~ dgamma(1,1)
 })
 


my.data <- list(X = X, Y = data.train$n)
my.constants <- list(JJ=dim(X)[1], KK=dim(X)[2],id_idx=id_idx, N.id_idx=max(id_idx), Year = Year, TT = max(Year))

initial.values <- function() list(beta.hat = rnorm(n= 8, mean = 0, sd=1),
                                  sigma_e = rgamma(1,1,1),
                                  phi = rbeta(1,1,1))

parameters.to.save <- c("beta.hat","phi","sigma_e")

### Run MCMC
n.iter <- 20000
n.burnin <-5000
n.chains <- 2
n.thin<-5

mcmc.output <- nimbleMCMC(code = model.param,
                          constants = my.constants,
                          data = my.data,
                          inits = initial.values,
                          monitors = parameters.to.save,
                          niter = n.iter,
                          nburnin = n.burnin,
                          nchains = n.chains, 
                          thin = n.thin)


beta.plot <- MCMCtrace(object = mcmc.output,
          pdf = TRUE, # no export to PDF
          ind = TRUE, # separate density lines per chain
          params = "beta.hat")

phi.plot <- MCMCtrace(object = mcmc.output,
          pdf = TRUE, # no export to PDF
          ind = TRUE, # separate density lines per chain
          params = "phi")


sig_e.plot <- MCMCtrace(object = mcmc.output,
          pdf = TRUE, # no export to PDF
          ind = TRUE, # separate density lines per chain
          params = "sigma_e")


# save each estimation ( ex, R.hat, a.hat, beta.hat, k.hat, lambda.hat)
# total num of estimation = 497 + 1 + 8 + 1 + 2057
# round 2 point for each parameters of 10000 samples 

result <- MCMCsummary(object = mcmc.output, round = 2)
resresult

head(result); dim(result)

# save estimated params
# beta.hat
for(k in 1:KK) {betas.hat[k] = result[paste0('beta.hat[',k,']'),1]}
phi.hat <- result['phi',1]
sigma_e.hat <- result['sigma_e',1]


#####################################################
# estimation of beta.hat[1]-beta.hat[8]
# combine parameter all estimation( beta.hat,  a.hat)
#####################################################
beta.est.mean <-betas.hat;
for(k in 1:KK) {beta.est.sd[k] = result[paste0('beta.hat[',k,']'),2]}
mu <- rep(0,8)
pval <- pnorm(-abs((beta.est.mean - mu)/beta.est.sd))
beta.pval <-round(pval*2, 4)

beta.df <- data.frame(cbind(beta.est.mean= beta.est.mean, beta.est.sd = beta.est.sd,beta.pval = beta.pval))

param_result <- rbind(beta.df, phi.hat,sigma_e.hat)
rownames(param_result) <- c('beta0','TypeCity','TypeCounty','TypeSchool','TypeTown','TypeVillage',"C2","C3",'phi.hat','sigma_e.hat')

param_result




#####################################################
# Fill missing data for Y's using Y's distn assumption(ex rpois(mean= lambda * R, phi) with sampled R.hat , lambda.hat
#####################################################

# Y 구성
tau = 5 ; Y = data.frame(matrix(nrow = length(id_uniq), ncol = tau))


for(i in 1:length(id_uniq)){
  Y[i,tau+1] = id_uniq[i]
  ys = data.train[data.train$PolicyNum == id_uniq[i], ]$n
  for(j in 1: length(ys)){
    Y[i,j] = ys[j]
  }
}
#Y$n <- data.valid$n # target은 이름이 n인 열로 구성

colnames(Y) <- c("y_1", "y_2", "y_3",  "y_4", "y_5", "ID")

# set TT : target: TT +1

#-----------------------------------------------
# fill missing R2
#-----------------------------------------------

TT <- 2

#####################################################
# data preprocessing
#####################################################
# (1) uniq_id for PolicyNum

id_uniq = unique(data.train$PolicyNum)

length(id_uniq) # 497

# (2) idx per id per PolicyNum's Year since each ID has different year of frequency (i.e not the same observations, 5,5,4,...)
id_idx=c(NULL)

for(i in 1:length(id_uniq)){
  ind = id_uniq[i] == data.train$PolicyNum 
  id_idx=c(id_idx, rep(i, sum(ind)))
}
Year = data.train$Year-2005 # Year[j]

length(id_idx)

data.train <- combine_kappa_level(data.train)
data.valid <- combine_kappa_level(data.valid)


#####################################################
# estimate a.hat, lambda.hat, R, beta.hat using nimble
#####################################################


round(apply(data.train, 2, mean)*100,2) 

# Make design matrix 
X = model.matrix(n~TypeCity+TypeCounty+TypeSchool+TypeTown+TypeVillage+C2+C3, data = data.train)
#X = model.matrix(n~TypeCity+TypeCounty+TypeMisc+TypeTown+TypeSchool+C2+C3, data = data.train)


model.param <- nimbleCode({
  # likelihood
  for(j in 1:JJ){
    lambda[j] <- exp(inprod(X[j,1:KK], beta.hat[1:KK]))
    Y[j]~dpois(lambda[j]*exp(R[id_idx[j], Year[j]]))
    
  }
  
  # prior
  for(i in 1:KK){
    beta.hat[i] ~ dnorm(0,1)
  }
  
  # random effect
  for(n in 1:N.id_idx){
    R[n, 1]~dnorm(0, sigma_e/sqrt(1-phi^2))
    for(t in 2:TT) {
      R[n, t]~dnorm(phi * R[n, t-1], sigma_e)
    }
    R6[n]~dnorm(phi * R[n, TT], sigma_e)
  }
  
  # prior of param
  phi ~ dbeta(1,1)
  sigma_e ~ dgamma(1,1)
})



my.data <- list(X = X, Y = data.train$n)
my.constants <- list(JJ=dim(X)[1], KK=dim(X)[2],id_idx=id_idx, N.id_idx=max(id_idx), Year = Year, TT = max(Year))

initial.values <- function() list(beta.hat = rnorm(n= 8, mean = 0, sd=1),
                                  sigma_e = rgamma(1,1,1),
                                  phi = rbeta(1,1,1))

parameters.to.save <- c("beta.hat","phi","sigma_e")

### Run MCMC
n.iter <- 20000
n.burnin <-5000
n.chains <- 2
n.thin<-5

mcmc.output <- nimbleMCMC(code = model.param,
                          constants = my.constants,
                          data = my.data,
                          inits = initial.values,
                          monitors = parameters.to.save,
                          niter = n.iter,
                          nburnin = n.burnin,
                          nchains = n.chains, 
                          thin = n.thin)


beta.plot <- MCMCtrace(object = mcmc.output,
                       pdf = TRUE, # no export to PDF
                       ind = TRUE, # separate density lines per chain
                       params = "beta.hat")

phi.plot <- MCMCtrace(object = mcmc.output,
                      pdf = TRUE, # no export to PDF
                      ind = TRUE, # separate density lines per chain
                      params = "phi")


sig_e.plot <- MCMCtrace(object = mcmc.output,
                        pdf = TRUE, # no export to PDF
                        ind = TRUE, # separate density lines per chain
                        params = "sigma_e")


# save each estimation ( ex, R.hat, a.hat, beta.hat, k.hat, lambda.hat)
# total num of estimation = 497 + 1 + 8 + 1 + 2057
# round 2 point for each parameters of 10000 samples 

result <- MCMCsummary(object = mcmc.output, round = 2)
resresult

head(result); dim(result)

# save estimated params
# beta.hat
for(k in 1:KK) {betas.hat[k] = result[paste0('beta.hat[',k,']'),1]}
phi.hat <- result['phi',1]
sigma_e.hat <- result['sigma_e',1]






############################################
# make y filled with absent value using distn assumption
###########################################
lambda.hat <- exp(X%*%betas.hat)

# 결측치 채우기
for (i in 1:nrow(Y)){
  for (j in 1:tau){
    if (is.na(Y[i,j])==TRUE)
      Y[i,j] = rpois(1, lambda.hat[id_idx[i]] * exp(get(paste0("R",TT+1))[i,]) )}
}
head(Y)

# total Y : columns = ID, y_1,y_2,y_3,y_4, y_5, lambda.hat_per_ID
train_y <- Y
train_x <- merge(x_id, lambda_id, by.x = 'PolicyNum', by.y = 'ID')
train_x<- train_x[,-9] # drop lambda_kappa_level

# tau+1의 Y값이 있는 data.valid를 이용해 ID로 병합함. 
test_y <- merge(Y, data.valid[,c('PolicyNum','n')], by.x = "ID", by.y = 'PolicyNum') 

# concat test x
testx_id = unique(data.valid[,c("PolicyNum",'TypeCity','TypeCounty','TypeSchool','TypeTown','TypeVillage',"C2","C3")]) # 379, 8 
colnames(testx_id) <- c("PolicyNum",'TypeCity','TypeCounty','TypeSchool','TypeTown','TypeVillage',"C2","C3")

test_x <- merge(testx_id, lambda_id, by.x = 'PolicyNum', by.y = 'ID')
test_x<- test_x[,-9] # drop lambda_kappa_level


# save the file train_y, train_x, test_y, test_x
filename = paste0(path,'/data')
ifelse(dir.exists(paste0(filename,"\\test_train_data")),F, dir.create(paste0(filename,"\\test_train_data")) ) 
data_folder <- paste0(filename,"\\test_train_data")

write.csv(train_y, file = paste0(data_folder,"/train_y.csv"))
write.csv(train_x, file = paste0(data_folder,"/train_x.csv"))
write.csv(test_y, file = paste0(data_folder,"/test_y.csv"))
write.csv(test_x, file = paste0(data_folder,"/test_x.csv"))


