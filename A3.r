####################question 1.1####################
cUrl = paste("http://scrippsco2.ucsd.edu/sites/default/files", "/data/flask_co2_and_isotopic/daily_co2/fldav_spo.csv",
             sep = "")
cFile = basename(cUrl)
if (!file.exists(cFile)) download.file(cUrl, cFile)
co2s = read.table(cFile, header = FALSE, sep = ",", skip = 69,
                  stringsAsFactors = FALSE)
co2s[co2s[, 6] > 0, 7] = NA
co2s = data.frame(date = strptime(co2s[, 1], format = "%Y-%m-%d", tz = "UTC"), co2 = co2s[, 7]) 
plot(co2s$date, co2s$co2, xlab="time", ylab="CO2 concentration", main="CO2 concentrations vs year")
timeOrigin = ISOdate(1980, 1, 1, 0, 0, 0, tz = "UTC")
co2s$days = as.numeric(difftime(co2s$date, timeOrigin, units = "days")) 
co2s$cos12 = cos(2 * pi * co2s$days/365.25)
co2s$sin12 = sin(2 * pi * co2s$days/365.25)
co2s$cos6 = cos(2 * 2 * pi * co2s$days/365.25)
co2s$sin6 = sin(2 * 2 * pi * co2s$days/365.25)
library("mgcv")
co2model=gam(co2~s(days)+sin12+cos6+sin6+cos12, data=co2s)
summary(co2model)
coPred = predict(co2model, co2s, se.fit = TRUE)
res=as.matrix(co2s$co2-coPred$fit)
qqnorm(res[!is.na(res)])
qqline(res)
plot(co2model, xlab="time (year)", ylab="f(d_i)", main="f(d_i) versus year")
#new data from 1970-1-1 to 2016-12-31
newX = data.frame(date = seq(ISOdate(1970, 1, 1, 0, 0, 0, tz = "UTC"), by = "1 days", length.out = 365 * 46))
newX$days = as.numeric(difftime(newX$date, timeOrigin, units = "days")) 
newX$cos12 = cos(2 * pi * newX$days/365.25)
newX$sin12 = sin(2 * pi * newX$days/365.25)
newX$cos6 = cos(2 * 2 * pi * newX$days/365.25)
newX$sin6 = sin(2 * 2 * pi * newX$days/365.25)
#calculate the slope from 1970-1-1 to 2016-12-31
dev = predict(co2model, newX, type="lpmatrix")
for(i in 6:14){
  if(i==6){
  sumdev=co2model$coef[i]*dev[,i]}
  else{sumdev=sumdev+co2model$coef[i]*dev[,i]}
}
slope=diff(sumdev)/diff(newX$days)
index=1:length(slope)
index_date=as.Date("1970-01-01") + index
plot(y=slope,x=index_date, xlab="year", main="slope from 1970-1-1 to 2016-12-31")
plot(y=slope,x=index_date, xlim=c(as.Date("1980-01-01"), as.Date("1982-12-31")),xlab="year", main="slope from 1980-1-1 to 1982-12-31")
plot(y=slope,x=index_date, xlim=c(as.Date("1990-01-01"), as.Date("1990-12-31")),xlab="year", main="slope from 1990-1-1 to 1990-12-31")


plot(newX[2190:2555,]$date, predict(co2model, newX[2190:2555,]),xlab="", main="CO2 concentration in 1976", ylab='fitted value')
plot(newX[5840:6205,]$date, predict(co2model, newX[5840:6205,]),xlab="", main="CO2 concentration in 1986", ylab='fitted value')
plot(newX[9490:9855,]$date, predict(co2model, newX[9490:9855,]),xlab="", main='CO2 concentration in 1996', ylab='fitted value')
plot(newX[13140:13505,]$date, predict(co2model, newX[13140:13505,]),xlab="", main='CO2 concentration in 2006', ylab='fitted value')

newX = data.frame(date = seq(ISOdate(1970, 1, 1, 0, 0, 0, tz = "UTC"), by = "1 days", length.out = 365 * 50))
newX$days = as.numeric(difftime(newX$date, timeOrigin, units = "days")) 
newX$cos12 = cos(2 * pi * newX$days/365.25)
newX$sin12 = sin(2 * pi * newX$days/365.25)
newX$cos6 = cos(2 * 2 * pi * newX$days/365.25)
newX$sin6 = sin(2 * 2 * pi * newX$days/365.25)
coPred = predict(co2model, newX, se.fit = TRUE)
coPred = data.frame(est = coPred$fit, lower = coPred$fit -
                      2 * coPred$se.fit, upper = coPred$fit + 2 * coPred$se.fit) 
plot(main="predicted CO2 concentration",newX$date, coPred$est,type = "l", xlab="year", ylab="CO2 concentration") 
matlines(as.numeric(newX$date), coPred[, c("lower", "upper","est")],lty = 1, col = c("blue", "blue", "red")) 
abline(h=400)

####################question 1.2####################
library(nlme)
library("INLA")
data("MathAchieve",package="MEMSS")
head(MathAchieve)
attach(MathAchieve)
#remove the data with negative math scores 
mathnew=MathAchieve[MathAchieve$MathAch>0,]
attach(mathnew)
#plot(density(MathAch))
opt = function(parameter) {sum((qgamma(c(0.05,0.95),shape=parameter[1],rate=parameter[2])-c(1, 100))^2)}
best = optim(c(1, 1), opt)
best_shape=best$par[1]
best_rate=best$par[2]
model1_2=inla(MathAch~Minority*SES+f(School,model='iid',hyper=list(theta=list(param=c( best_shape, best_rate)))),
           control.fixed=list(mean.intercept=0,prec.intercept=1/(10^2), mean=0,prec=1/(3^2)),
           family='Gamma',
           data=mathnew,
           control.family=list(hyper=list(theta=list(prior='loggamma',param=c(1,1/100)))),
           control.inla=list(fast=FALSE,h=0.2,strategy='laplace'),
           control.mode=list(theta=c(log(100),-2),restart=TRUE),
           control.compute=list(return.marginals=TRUE))
summary(model1_2)
#within school variance is:
r=qgamma(c(0.025,0.975), shape=3.071, scale = exp(2.6010)/3.071)
r[2]-r[1]
prec_post=model1_2$marginals.hyperpar$'Precision for School'
sigma_post=inla.tmarginal(function(x) sqrt(1/x),prec_post)
sigma.e=inla.emarginal(function(x) x,sigma_post)
#between school:
exp(2.6010+2*sigma.e)-exp(2.6010-2*sigma.e)
#xtable::xtable(summary(model1_2$summary.fixed,digits=2))
####################question 1.3####################
sUrl='http://www.lancaster.ac.uk/staff/diggle/APTS-data-sets/lead2000_data.txt'
download.file(sUrl, "lead2000_data.txt")
x=read.table("lead2000_data.txt", header=TRUE, skip=3)
hist(x[,"z"], breaks=50, freq=FALSE,ylim=c(0,3), main="Comparison to the empirical distribution", xlab="")
#homework2:
xNorm=rnorm(100000, mean=0.34+log(500)*0.06, sd=0.04)
xBc=(xNorm*(-0.62)+1)^(-1/0.62)
hist(xBc[xBc<10], breaks=50, main="homework2", ylim=c(0,15000))
#homework3:
xGamma=rgamma(100000,shape=31.1251, scale=exp(-0.2975+0.0881*(log(500)))/31.1251)
hist(xGamma[xGamma<10], breaks=50, main="homework3", ylim=c(0,10000))
lines(density(xBc))
lines(density(xGamma),col='blue')

####################question 2####################
library("INLA")
data="smoke.RData"
load(data)

opt_new = function(par) {sum((qgamma(c(0.1,0.9),shape=par[1],scale=par[2])-
                                c(4*qnorm(0.9)^2/(log(3))^2, 4*qnorm(0.9)^2/(log(2))^2))^2)}
state = optim(c(1, 1), opt_new)
state_param=state$par
state_param
#first is shape, second is scale
#check
#qgamma(c(0.1,0.9),shape=state_param[1],scale=state_param[2])
#prior for precision of school
opt_new2 = function(par) {sum((qgamma(c(0.1,0.9),shape=par[1],scale=par[2])-
         c(4*qnorm(0.9)^2/(log(1.2))^2, 4*qnorm(0.9)^2/(log(1.1))^2))^2)}
school = optim(c(1, 2), opt_new2)
school_param=school$par
#check
#qgamma(c(0.1, 0.9),shape=school_param[1],scale=school_param[2])
school_param
#the prior for kappa
opt_k = function(par) {sum((qnorm(c(0.8,0.9),log(3.5),par)-c(log(3),log(4)))^2)}
kappa_param=optim(2, opt_k)
kappa_scale=kappa_param$par
forInla=smoke[,c("Age", "Age_first_tried_cigt_smkg", "Sex", "Race", "state", "school", "RuralUrban")]
forInla=na.omit(forInla)
forInla=as.list(forInla)
forSurv=data.frame(time=(pmin(forInla$Age_first_tried_cigt_smkg, forInla$Age)-4)/10, event=forInla$Age_first_tried_cigt_smkg<=forInla$Age)
forInla[forInla$Age_first_tried_cigt_smkg ==8,"event"]=2
forInla$y=inla.surv(forSurv$time, forSurv$event)
model1 = inla(y ~ Race + Sex + RuralUrban +
         f(school,model = "iid", hyper = list(prec = list(prior = "loggamma",param = c(school_param[1],1/school_param[2]))))
         + f(state,model = "iid", hyper = list(prec = list(prior = "loggamma",  param = c(state_param[1],1/state_param[2])))), 
         family = "weibullsurv", data = forInla,
         control.family = list(hyper = list(alpha = list(prior = "normal",param = c(log(3.5),1/kappa_scale)))))
summary(model1)
#kappa
lnormalsample=rnorm(10^5, mean=log(3.5), sd=kappa_scale)
plot(density(exp(lnormalsample)),xlim=c(3,4),ylim=c(0,10),col='red', main="Prior and Posterior plot for Kappa",xlab="")
kappa=as.matrix(model1$marginals.hyper$'alpha parameter for weibull')
lines(kappa[,'x'],kappa[,'y'])
#school
rgamma1=rgamma(10^5, shape=school_param[1], rate=1/school_param[2])
tgamma1=1/sqrt(rgamma1)
PrecPost=model1$marginals.hyper$'Precision for school'
sigmaPost=cbind(PrecPost,sigma=1/sqrt(PrecPost[,'x']),dSigma=PrecPost[,'y']*2*PrecPost[,'x']^(3/2),
                priorSigma=2*PrecPost[,"x"]^(3/2)*dgamma(PrecPost[,"x"],shape=school_param[1], rate=1/school_param[2]))
plot(sigmaPost[,'sigma'], sigmaPost[,'dSigma'],type='l',xlim=c(-1,1),ylim=c(0,40),
     main="Prior and Posterior for school",ylab="Density",xlab="")
lines(density(tgamma1),col='red')

#state
rgamma_state=rgamma(10^5, shape=state_param[1], rate=1/state_param[2])
rgamma_state=1/sqrt(rgamma_state)
PrecPost=model1$marginals.hyper$'Precision for state'
sigmaPost=cbind(PrecPost,sigma=1/sqrt(PrecPost[,'x']),dSigma=PrecPost[,'y']*2*PrecPost[,'x']^(3/2),
                priorSigma=2*PrecPost[,"x"]^(3/2)*dgamma(PrecPost[,"x"],shape=state_param[1], rate=1/state_param[2]))
plot(sigmaPost[,'sigma'], sigmaPost[,'dSigma'],type='l',xlim=c(0,1),ylim=c(0,15),
     main="Prior and Posterior for state",ylab="Density",xlab="")
lines(density(rgamma_state),col='red')
#with interaction
model2 = inla(y ~ Race * Sex * RuralUrban +
                f(school,model = "iid", hyper = list(prec = list(prior = "loggamma",param = c(school_param[1],1/school_param[2]))))
              + f(state,model = "iid", hyper = list(prec = list(prior = "loggamma",  param = c(state_param[1],1/state_param[2])))), 
              family = "weibullsurv", data = forInla,
              control.family = list(hyper = list(alpha = list(prior = "normal",param = c(log(3.5),1/kappa_scale)))))
summary(model2)
x=seq(0,100,len=50000)
plot(x, exp(-1.8324+0.2908)*3.4*x^(2.4),type='l', xlab='age',ylab='hazard',
     main='')
lines(x, exp(-1.8324)*3.4*x^(2.4),col='blue')
legend("topleft", lty = 1, col = c("black", "blue"), legend = c("white rural male","white urban male"))
