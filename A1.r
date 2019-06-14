#question1
options(digits=3)
set.seed(0)
N=100
CI=NULL
count=0
CI=matrix(NA, nrow = N, ncol = 2)
beta_hat_collection=c()
stats_collection=c()
for (i in 1:N){
  x=seq(-10,10,len=40)
  off=rep(c(1,-1),c(25,length(x)-25))
  y=rpois(length(x),exp(off+0.5+0.2*x))
  coeff=summary(glm(y~x+offset(off),family='poisson'))[['coefficients']]
  CI[i,1]=coeff[2,1]-2*coeff[2,2]
  CI[i,2]=coeff[2,1]+2*coeff[2,2]
  beta_hat=coeff[2,1]
  beta_hat_collection=append(beta_hat_collection, beta_hat)
  #calculate coverage probability
  if(0.2<CI[i,1]&0.2<CI[i,2]){
    count=1+count
  }
  model=glm(y~x+offset(off),family='poisson')
  #restricted model
  x02 = 0.2*x
  rmodel=glm(y ~ offset(off) + offset(x02), family='poisson')
  stats=rmodel$deviance-model$deviance
  stats_collection=append(stats_collection, stats)
  
}
#qqplot for testing chi-square distribution
qqplot(qchisq(ppoints(100), df = 1), stats_collection,
       main = expression("Q-Q plot for" ~~ {chi^2}[nu == 1]))
qqline(stats_collection, distribution = function(p) qchisq(p, df = 1))

#qqpplot for testing normal distribution
coverage_prob=(100-2*count)/100
centerialized=(beta_hat_collection-mean(beta_hat_collection))/sqrt(var(beta_hat_collection))
qqnorm(centerialized)
qqline(centerialized)

#question2
#1.log-normal distribution
mu=log(2/sqrt(7/4))
sigma=sqrt(log(7/4))
logn=function(x, mu, sigma){
    z=1/x*dnorm(log(x),mu,sigma)
    return(z)
  }
#2. weibull distribution

library(nleqslv)
params=function(x){
  y=numeric(2)
  y[1]=x[1]*gamma(1+1/x[2])-2
  y[2]=x[1]^2*(gamma(1+2/x[2])-gamma(1+1/x[2])^2)-3
  y
}
weibull_param=nleqslv(c(0.5,0.5), params)
#3. zero inflated poisson distribution
zip=function(x){
  if(x==0){
    z=1/5+(1-1/5)*exp(-2.5)
  }else{
    z=(1-1/5)*dpois(x,2.5)
    
  }
  return(z)
}
#------------sampling------------
size=1000
cts_x=seq(0,50,length=1000)
discrete_x=c(0:999)
collection=matrix(NA, nrow = 5, ncol = 1000)
for(j in 1:size){
  #log-normal
  collection[1,j]=logn(cts_x[j], mu, sigma)
  #gamma
  collection[2,j]=dgamma(cts_x[j], shape=4/3, rate=2/3)
  #weibull
  collection[3,j]=dweibull(cts_x[j], shape=weibull_param$x[2],scale=weibull_param$x[1])
  #negative binomial
  collection[4,j]=dnbinom(discrete_x[j], 4, 2/3)
  #zero_inflated poisson
  collection[5,j]=zip(discrete_x[j])
}
#------------plot-------------------
plot(cts_x, collection[1,], xlab="x", ylab="y",type="l", col="red" )
lines(cts_x, collection[2,], type="l", col="green" )
lines(cts_x, collection[3,], type="l", col="blue" )
lines(discrete_x, collection[4,], type="p", col="orange")
lines(discrete_x, collection[5,], type="p", col="black")
#---------find 99% upper quantile for each function--------
#zero-inflated model
quant=1/5+(1-1/5)*exp(-2.5); i=1
while (quant<0.99){
  quant=quant+(1-1/5)*dpois(i,2.5)
  i=i+1
}
print (i-1)
#gamma
qgamma(0.99, shape=4/3, rate=2/3)
#weibull
qweibull(0.99, shape=weibull_param$x[2],scale=weibull_param$x[1])
#log-normal
exp(qnorm(0.99)*sigma+mu)
#negative binomial
qnbinom(0.99, 4, 2/3)
#-----simulate 20 realisations from each distrbution------------
#zero-inflated poisson
zip_sample=rbinom(20,1,4/5)*rpois(20,2.5)
mean(zip_sample)
var(zip_sample)
#gamma
gamma_sample=rgamma(20,shape=4/3, rate = 2/3)
mean(gamma_sample)
var(gamma_sample)
#weibull
weibull_sample=rweibull(20, shape=weibull_param$x[2],scale=weibull_param$x[1])
mean(weibull_sample)
var(weibull_sample)
#
#generate lognormal
lognormal_sample=exp(rnorm(20, mean=mu, sd=sigma))
mean(lognormal_sample)
var(lognormal_sample)
#negative binomial
negativebinomial_sample=rnbinom(20,4,2/3)
mean(negativebinomial_sample)
var(negativebinomial_sample)
#--------------DATA ANALYSIS--------------------
data("fruitfly", package="faraway")
summary(fruitfly)
attach(fruitfly)
contrasts(activity)
thorax_center=thorax-mean(thorax)
model=glm(longevity~activity+thorax_center, data=fruitfly, family=Gamma(link="log"))
knitr::kable(rbind(summary(model)$coef, shape=c(1/summary(model)$dispersion, NA,NA,NA)),digits=4)
summary(model)
#coefficient
exp(summary(model)$coef[,"Estimate"])
shape=1/summary(model)$dispersion
scale=exp(model$coef["(Intercept)"])/shape
#--isolated
hist(fruitfly$longevity[fruitfly$activity=="isolated"&abs(thorax_center)<0.1], prob=TRUE, main="isolated", xlab="longevity", xlim=c(0,150), ylim=c(0,0.05))
xSeq=seq(0,150, len=1000)
lines(xSeq, dgamma(xSeq, shape=shape, scale=scale), col="red")
#--one
hist(fruitfly$longevity[fruitfly$activity=="one"&abs(thorax_center)<0.1], prob=TRUE, main="one", xlab="longevity", xlim=c(0,150), ylim=c(0,0.05))
xSeq=seq(0,150, len=1000)
lines(xSeq, dgamma(xSeq, shape=shape, scale=scale), col="red")
#--low
hist(fruitfly$longevity[fruitfly$activity=="low"&abs(thorax_center)<0.1], prob=TRUE, main="low", xlab="longevity", xlim=c(0,150), ylim=c(0,0.05))
xSeq=seq(0,150, len=1000)
lines(xSeq, dgamma(xSeq, shape=shape, scale=scale), col="red")
#--many
hist(fruitfly$longevity[fruitfly$activity=="many"&abs(thorax_center)<0.1], prob=TRUE, main="many", xlab="longevity", xlim=c(0,150), ylim=c(0,0.05))
xSeq=seq(0,150, len=1000)
lines(xSeq, dgamma(xSeq, shape=shape, scale=scale), col="red")
#--high
hist(fruitfly$longevity[fruitfly$activity=="high"&abs(thorax_center)<0.1], prob=TRUE, main="high", xlab="longevity", xlim=c(0,150), ylim=c(0,0.05))
xSeq=seq(0,150, len=1000)
lines(xSeq, dgamma(xSeq, shape=shape, scale=scale), col="red")

#-------------------------REPORT------------------------------
load("/Users/sanansakura/Desktop/STA2201/smoke.rdata")
tobacco=smoke$chewing_tobacco_snuff_or
race=as.vector(smoke$Race)
rural=as.vector(smoke$Rural)
tobacco[is.na(tobacco)]=7
race[is.na(race)]=8
new_tobacco=NULL
new_race=NULL
new_rural=NULL
#made up a new data set
for (i in 1:length(race)){
  if (race[i]!=8 && race[i]!='pacific' &&race[i]!='asian'
      &&race[i]!='native' && tobacco[i]!=7){
    new_tobacco=append(new_tobacco, as.numeric(tobacco[i]))
    new_race=append(new_race,race[i])
    new_rural=append(new_rural, rural[i])
  }
}
new_race=factor(new_race)
new_rural=factor(new_rural)
model1 = glm(new_tobacco~ new_rural+new_race+new_rural:new_race, family=binomial)
summary(model1)
knitr::kable(rbind(summary(model1)$coef, shape=c(1/summary(model1)$dispersion, NA,NA,NA)),digits=4)
#construct confidence interval
stderror=summary(model1)$coef[, "Std. Error"]
parTable=outer(stderror, c(MLE = 0, lower = -1.96, upper = 1.96)) +
summary(model1)$coef[, "Estimate"]
rownames(parTable) = gsub("Intercept", "baseline", rownames(parTable))
knitr::kable(exp(parTable), digits = 3)

#------hypothese 2-----------
waterpipe=smoke$ever_tobacco_hookah_or_wa
age=smoke$Age
sex=as.vector(smoke$Sex)
age[is.na(age)]=6
sex[is.na(sex)]=5
waterpipe[is.na(waterpipe)]=4
new1_waterpipe=NULL
new1_race=NULL
new1_sex=NULL
new1_age=NULL
new1_rural=NULL
for (i in 1:length(race)){
  if (race[i]!=8 && waterpipe[i]!=4 && age[i]!=6 && sex[i]!=5 ){
    new1_waterpipe=append(new1_waterpipe, as.numeric(waterpipe[i]))
    new1_race=append(new1_race,race[i])
    new1_sex=append(new1_sex, sex[i])
    new1_age=append(new1_age, age[i])
    new1_rural=append(new1_rural,rural[i])
  } }
model2 = glm(new1_waterpipe ~ new1_age+new1_race+new1_sex+new1_rural, family=binomial) 
summary(model2)
knitr::kable(rbind(summary(model2)$coef, shape=c(1/summary(model2)$dispersion, NA,NA,NA)),digits=4)
stderror=summary(model2)$coef[, "Std. Error"]
parTable=outer(stderror, c(MLE = 0, lower = -1.96, upper = 1.96)) +
  summary(model2)$coef[, "Estimate"]
rownames(parTable) = gsub("Intercept", "baseline", rownames(parTable))
knitr::kable(exp(parTable), digits = 3)



tobacco=smoke$chewing_tobacco_snuff_or
race=as.vector(smoke$Race)
tobacco[is.na(tobacco)]=7
race[is.na(race)]=8
age[is.na(age)]=6
sex[is.na(sex)]=5
new2_tobacco=NULL
new2_race=NULL
new2_age=NULL
new2_sex=NULL
for (i in 1:length(race)){
  if (race[i]!=8 && tobacco[i]!=7 && age[i]!=6 && sex[i]!=5 ){
    new2_tobacco=append(new2_tobacco, as.numeric(tobacco[i]))
    new2_race=append(new2_race,race[i])
    new2_sex=append(new2_sex, sex[i])
    new2_age=append(new2_age, age[i])
  } }

model3 = glm(new2_tobacco ~ new2_age+new2_race+new2_sex, family=binomial) 
knitr::kable(rbind(summary(model3)$coef, shape=c(1/summary(model3)$dispersion, NA,NA,NA)),digits=4)
stderror=summary(model3)$coef[, "Std. Error"]
parTable=outer(stderror, c(MLE = 0, lower = -1.96, upper = 1.96)) +
summary(model3)$coef[, "Estimate"]
rownames(parTable) = gsub("Intercept", "baseline", rownames(parTable))
knitr::kable(exp(parTable), digits = 3)

