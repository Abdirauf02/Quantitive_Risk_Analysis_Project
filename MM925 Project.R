
set.seed(24)
par(mfrow= c(1,1))
##Question 1 

#B which is the infection pressure (weeks) is between 0.00001 to 0.0006 and likely 0,000187,
#B is a Poisson distribution. but beta is uncertain

beta_mode <- 0.000187
x.range <- seq(0,0.0006,0.000005)

calc_b_gamma <- function(a, m){
  b<- m/(a-1)
  print(b)
}

a<- 7
b<- calc_b_gamma(a,beta_mode)

plot(x.range,dgamma(x.range, shape = a, scale = b), type = "l", lwd= 2,
     main = "Gamma distribution for uncertainty of Beta", xlab = "Infection Pressure (weeks)",
     ylab= "Density")
#points(0.000187,4000)
#points(0.00001,0)
#points(0.0006,0)


mu<- 0.078

#second order model with Beta as uncertian and p 

simulate_p <- function(n_uncert,n_var,t){
  p<- array(dim = c(n_uncert,n_var))
  for (i in 1:n_uncert) {
    beta <- rgamma(n_var,shape = a, scale = b)
    for (j in 1:n_var) {
      p[i,j] <- (beta[j]*(1-exp(-(mu+beta[j])*t)))/(beta[j]+ mu)
    }
  }
  p
}

#simulate p
sim_p<- simulate_p(10,1000,2)

#Plot of ecdf
plot(ecdf(sim_p[1,]),do.p=FALSE,verticals=TRUE,xlab="Probability of infection during 2 weeks",
     main="ECDF Probabilty of Dog Infection during 2 weeks")
for (i in 2:10) {
  lines(ecdf(sim_p[i,]),do.p=FALSE,verticals=TRUE)
}

#summary statistics for each row of simulated p
for(i in 1:10){
  new <- quantile(sim_p[i,], c(0,0.025,0.5,0.975,1))
  print(new)
}


#Question 2

#since we have no prior knowledge of compliance levels 
x <- seq(0,1,0.001)
prior <-  dunif(x,0,1)
plot(x, prior, type="l", ylab="Density", lwd=2, col="blue", lty=2, 
     main="")

#likelihood modelled as binomal
n= 200
s= 20

likelihood <- dbinom(x,n,s/n)
plot(x, likelihood, main = "Likelihood",xlab=" Probabiltiy of Dog not treated ", 
     type="l", ylab="Density", col="orange", lty=2)

#posterior modelled as beta 
x <- seq(0,1,0.001)
n= 200
s= 20
posterior <- dbeta(x, s+1,n-s+1)
plot(x, posterior, main = "Posterior",xlab=" Probabiltiy of Dog not treated ", 
     type="l", ylab="Density", col="black", lty=2)

lower_interval_old<-qbeta(0.025,s+1,n-s+1)
upper_interval_old<-qbeta(0.975, s+1,n-s+1)
median_old<-qbeta(0.5, s+1,n-s+1)


#Question 3

#use posterior from question 2 as now new prior

new_prior <- dbeta(x, s+1,n-s+1)
plot(x, new_prior, main = "New Prior",xlab=" Probabiltiy of Dog not treated ", 
     type="l", ylab="Density", col="black", lty=2)

#likelihood  
likelihood <- dbinom(s,n,x)
plot(x,likelihood/max(likelihood), main = "Likelihood",xlab=" Probabiltiy of Dog not treated ", 
     type="l", ylab="Density", col="orange", lty=4)


new_s <- 1
new_n <- 26

new_posterior <- dbeta(x, new_s+1,new_n-new_s+1)
plot(x, new_posterior, main = "Likelihood",xlab=" Probabiltiy of Dog not treated ", 
     type="l", ylab="Density", col="red", lty=1)


plot(x, new_prior/max(new_prior), xlab="Probabiltiy of Dog not treated", 
     main= "Scaled Probability of Dogs not treated ",
     type="l", ylab="Scaled density", lwd=2, col="black", lty=2)
lines(x,likelihood/max(likelihood),lwd=2, col="orange", lty=4)
lines(x,new_posterior/max(new_posterior), lwd=2, col="red", lty=1)
legend("topright", c("Prior: Beta(21,181)", 
                     "Likelihood: Binomial(200,0.1)",
                     "Posterior: Beta(1,26)"), 
       col=c("black","orange", "red"), lty=c(2,4,1), cex=0.8, lwd=2)

lower_interval_new<-qbeta(0.025,new_s+1,new_n-new_s+1)
upper_interval_new<-qbeta(0.975, new_s+1,new_n-new_s+1)
median_new<-qbeta(0.5, new_s+1,new_n-new_s+1)


##Question 4


dogs_pcb <- c(102, 76, 143, 50, 103, 36, 132, 38, 118, 58,142, 86, 48, 32, 4, 
              39, 32, 112, 73, 28, 73, 90, 61, 42, 71)

pcb <- function(x, b){
  a <- array(dim = c(b,2))
  z <- c()
  m<-mean(x)
  s <- sd(x)
  for(i in 1:b){
    z<- rnorm(length(x), m,s)
    a[i,1] <- mean(z)
    a[i,2] <- sd(z)
    
  }
  a
}
#Simulate
sim_pcb <- pcb(dogs_pcb,1000)

#Get under 24 hours and over 120 hours sum
under_24 <- c()
over_120 <- c()
outside_24_120 <- c()
for (i in 1:1000) {
  under_24[i] <- pnorm(24,sim_pcb[i,1],sim_pcb[i,2])
  over_120[i] <- 1-pnorm(120,sim_pcb[i,1],sim_pcb[i,2])
  outside_24_120[i] <- under_24[i]+over_120[i]
}

#Plots
plot(ecdf(outside_24_120), verticals = TRUE, do.points= FALSE, 
     xlab = "Probability of treatment time outside 24-120 hours ", 
     main= " ECDF of Probability of treatment time outside 24-120 hours ")
plot(density(outside_24_120, bw=0.05),  xlab = "Probability of treatment time outside 24-120 hours ", 
     main= " Probability of treatment time outside 24-120 hours ")

quantile(outside_24_120, c(0,0.025 ,0.25,0.5,0.75,0.975,1))



##Question 5

# dog infected 2 weeks 
# not treated 
# outside the treatment timing 24-120 hours

#any dog is infected  -> (Treated and outside*treated + not treated)* dog infection

treated_outside <- outside_24_120
not_treated <- new_posterior/max(posterior)
dog_infection <- sim_p

any_dog_is_infected <- (treated_outside*(1-not_treated)+ not_treated)*dog_infection
#check
mean(any_dog_is_infected)
for (i in 1:10) {
  summarystat <-quantile(any_dog_is_infected[i,], c(0,0.025 ,0.25,0.5,0.75,0.975,1))
  print(summarystat)
}


#Plot of ecdf
plot(ecdf(any_dog_is_infected[1,]),do.p=FALSE,verticals=TRUE,xlab="Probability of any dog is infected",
     main="ECDF of Probability of any dog is infected")
for (i in 2:10) {
  lines(ecdf(any_dog_is_infected[i,]),do.p=FALSE,verticals=TRUE)
}

plot(density(any_dog_is_infected[1,], bw= 0.0001))
for (i in 2:10) {
  lines(density(any_dog_is_infected[i,],bw= 0.0001))
}

#Question 6 


simulate_p_yearly <- function(n_uncert,n_var,t){
  p<- array(dim = c(n_uncert,n_var))
  for (i in 1:n_uncert) {
    beta <- rgamma(n_var,shape = a, scale = b)
    for (j in 1:n_var) {
      p[i,j] <- (beta[j]*(1-exp(-(mu+beta[j])*(t/52))))/(beta[j]+ mu)
    }
  }
  p
}
#simulate p
num_of_dogs <- 1000
sim_p_yearly <- c()
for(i in 1:num_of_dogs){
  sim_p_yearly[i]<- simulate_p_yearly(10,1000,2)
}

n_dog_is_infected <- sim_p_yearly
n_dog_is_infected_yearly<-(treated_outside*(1-not_treated)+ not_treated)*n_dog_is_infected


plot(density(n_dog_is_infected_yearly, bw= 0.000001),xlab="Annual Probabilty of 1000 Dogs infection during 2 weeks",
     main="Density of Annual Probabilty of 1000 Dogs Infection during 2 weeks")

quantile(n_dog_is_infected_yearly, c(0,0.025 ,0.25,0.5,0.75,0.975,1))

#Question 7

#change beta to 0.0003, also increase the maximum to 0.0008, increase in infection pressure

beta_mode_new <- 0.0003
x.range <- seq(0,0.0008,0.000005)

calc_b_gamma <- function(a, m){
  b<- m/(a-1)
  print(b)
}

a_new<- 7
b_new<- calc_b_gamma(a,beta_mode_new)

plot(x.range,dgamma(x.range, shape = a_new, scale = b_new), type = "l", lwd= 2,
     main = "Gamma distribution for uncertainty of Beta", xlab = "Infection Pressure (weeks)",
     ylab= "Density")

simulate_p_yearly_new <- function(n_uncert,n_var,t){
  p<- array(dim = c(n_uncert,n_var))
  for (i in 1:n_uncert) {
    beta <- rgamma(n_var,shape = a_new, scale = b_new)
    for (j in 1:n_var) {
      p[i,j] <- (beta[j]*(1-exp(-(mu+beta[j])*(t/52))))/(beta[j]+ mu)
    }
  }
  p
}
#simulate p
num_of_dogs <- 1000
sim_p_yearly_new <- c()
for(i in 1:num_of_dogs){
  sim_p_yearly_new[i]<- simulate_p_yearly_new(10,1000,2)
}

n_dog_is_infected_new <- sim_p_yearly_new
n_dog_is_infected_yearly_new<-(treated_outside*(1-not_treated)+ not_treated)*n_dog_is_infected_new

plot(ecdf(n_dog_is_infected_yearly_new),do.p=FALSE,verticals=TRUE,xlab="Annual Probabilty of 1000 Dogs infection Beta = 0.0003",
     main="ECDF of Annual Probabilty of 1000 Dogs Infection during 2 weeks")
plot(density(n_dog_is_infected_yearly_new, bw= 0.000001),xlab="Annual Probabilty of 1000 Dogs infection during 2 weeks",
     main="Density of Annual Probabilty of 1000 Dogs Infection during 2 weeks")

quantile(n_dog_is_infected_new, c(0,0.025 ,0.25,0.5,0.75,0.975,1))
#increase in B increases the probabilty

# change complaince times 
#make window smaller


under_10 <- c()
over_72 <- c()
outside_10_72 <- c()
for (i in 1:1000) {
  under_10[i] <- pnorm(10,sim_pcb[i,1],sim_pcb[i,2])
  over_72[i] <- 1-pnorm(72,sim_pcb[i,1],sim_pcb[i,2])
  outside_10_72[i] <- under_10[i]+over_72[i]
}

plot(density(outside_10_72, bw=0.05),  xlab = "Probability of 10-72 hours ", 
     main= "10-72 hours ")

quantile(outside_10_72, c(0,0.025 ,0.25,0.5,0.75,0.975,1))


treated_outside_new <- outside_10_72

n_dog_is_infected_yearly_new_outside<-(treated_outside_new*(1-not_treated)+ not_treated)*n_dog_is_infected


plot(density(n_dog_is_infected_yearly_new_outside, bw= 0.000001),xlab="Annual Probabilty of 1000 Dogs infection during 2 weeks",
     main="Density of Annual Probabilty of 1000 Dogs Infection during 2 weeks")

quantile(n_dog_is_infected_yearly_new_outside, c(0,0.025 ,0.25,0.5,0.75,0.975,1))

#proababilty increase , smaller window


#increase in dogs to 2000


num_of_dogs_new <- 2000
sim_p_yearly_new <- c()
for(i in 1:num_of_dogs_new){
  sim_p_yearly_new[i]<- simulate_p_yearly(10,1000,2)
}

n_dog_is_infected_new_num_dog <- sim_p_yearly_new
n_dog_is_infected_yearly_new_num_dog<-(treated_outside*(1-not_treated)+ not_treated)*n_dog_is_infected_new_num_dog


plot(density(n_dog_is_infected_yearly_new_num_dog,bw= 0.000001),xlab="Annual Probabilty of 1000 Dogs infection during 2 weeks",
     main="Density of Annual Probabilty of 1000 Dogs Infection during 2 weeks")

quantile(n_dog_is_infected_yearly_new_num_dog, c(0,0.025 ,0.25,0.5,0.75,0.975,1))
#increased prob


#COLLECTION OF Q7 GRAPHS AND QUANTILES

par(mfrow= c(3,2))

#plot Beta = 0.0003
plot(x.range,dgamma(x.range, shape = a_new, scale = b_new), type = "l", lwd= 2,
     main = "Beta = 0.0003", xlab = "Infection Pressure (weeks)",
     ylab= "Density")
#plot annual beta = 0.0003
plot(density(n_dog_is_infected_yearly_new, bw= 0.000001),xlab="Annual Probabilty of 1000 Dogs",
     main="Beta = 0.0003")
#plot prob of 10 - 72
plot(density(outside_10_72, bw=0.05),  xlab = "Probability of 10-72 hours ", 
     main= "10-72 hours ")
#plot annual 10 -72
plot(density(n_dog_is_infected_yearly_new_outside, bw= 0.000001),xlab="Annual Probabilty of 1000 Dogs",
     main="10-72 hours")
#plot of 2000 annual dogs
plot(density(n_dog_is_infected_yearly_new_num_dog, bw= 0.000001),xlab="Annual Probabilty of 2000 Dogs",
     main="Density of Annual Probabilty of 2000 Dogs ")
#annual 

quantile(n_dog_is_infected_yearly, c(0,0.025 ,0.25,0.5,0.75,0.975,1))
#beta = 0.0003
quantile(n_dog_is_infected_yearly_new, c(0,0.025 ,0.25,0.5,0.75,0.975,1))
#outside 10-72
quantile(n_dog_is_infected_yearly_new_outside, c(0,0.025 ,0.25,0.5,0.75,0.975,1))
# 2000 dogs
quantile(n_dog_is_infected_yearly_new_num_dog, c(0,0.025 ,0.25,0.5,0.75,0.975,1))
