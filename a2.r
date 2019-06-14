library(nlme)

###########question2###########
data("MathAchieve",package="MEMSS")
head(MathAchieve)
attach(MathAchieve)

q2=lme(MathAch~Minority*SES, random = ~1|School, data=MathAchieve)
knitr::kable(summary(q2)$tTable, digits = 3)
summary(q2)
library(plyr)
freq=count(MathAchieve$School)
max(freq$freq)
min(freq$freq)
###########question3###########
load("CF.Rdata")
library("nlme")
#----------first model, full model
x$ageC = x$AGE - 18
resS = lme(FEV1 ~ GENDER * F508 * ageC + PSEUDOA, random = ~1 |ID, data = x) 
knitr::kable(summary(resS)$tTable, digits = 3)
summary(resS)
#----------random slope model, full model
resRS = lme(FEV1 ~ GENDER * F508 * ageC + PSEUDOA, random = ~1 + ageC |ID, data = x) 
knitr::kable(summary(resRS)$tTable, digits = 3)
summary(resRS)
#----------models with serial correlation, full model
resC_full = lme(FEV1 ~ GENDER * F508 * ageC + PSEUDOA, random = ~1|ID, data = x, correlation=corExp(form=~ageC|ID, nugget=T)) 
knitr::kable(summary(resC_full)$tTable, digits = 3)
summary(resC_full)
#----------models without F508-----------
resC_gene = lme(FEV1 ~ GENDER * ageC + PSEUDOA, random = ~1|ID, data = x, correlation=corExp(form=~ageC|ID, nugget=T)) 
knitr::kable(summary(resC_gene)$tTable, digits = 3)
summary(resC_gene)
library(lmtest)
lrtest(resC_full, resC_gene)
#----------models without interaction-----------
resC_inter = lme(FEV1 ~ GENDER * ageC + F508*ageC + PSEUDOA, random = ~1|ID, data = x, correlation=corExp(form=~ageC|ID, nugget=T)) 
knitr::kable(summary(resC_inter)$tTable, digits = 3)
summary(resC_inter)
lrtest(resC_full, resC_inter)
###########question4###########
