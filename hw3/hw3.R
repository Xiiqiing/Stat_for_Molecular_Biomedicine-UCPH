# load data for part1
t2d <- read.table("T2D.txt", header = T)
#str(t2d)

# 1A
g1 <- glm(T2D ~ 1, data = t2d, family = "binomial")
g2 <- glm(T2D ~ OGTT, data = t2d, family = "binomial")
test1a<-anova(g1,g2,test="Chisq")
cat("1A. The test is Wald test from summary() function, the p-value is:",summary(g2)$coefficients["OGTT","Pr(>|z|)"], ". So there is a difference in diagnosis outcome from the two different diagnosis methods. And the p-value of ANOVA test (", test1a$"Pr(>Chi)"[2],") also suggests model T2D ~ OGTT is better. \n")

# 1B
#levels(as.factor(t2d$T2D))
#0 is the "unaffected" and 1 is the "affacted"
OR_T2D<-exp(coef(g2))[2]
T2D_OGTT<-coef(g2)[1]+coef(g2)[2]
T2D_HbA1c<-coef(g2)[1]
p_T2D_OGTT<-1/(1+exp(-T2D_OGTT))
p_T2D_HbA1c<-1/(1+exp(-T2D_HbA1c))
cat("1B. The odds ratio (OR) is:",OR_T2D,", and the relative risk of getting T2D from an OGTT compared to HbA1c is:",p_T2D_OGTT/p_T2D_HbA1c,"\n")

# 1C
g3 <- glm(T2D ~ OGTT + age + sex + bmi, data = t2d, family = "binomial")
cat("1C. The odds ratio of getting a T2D from an OGTT compared to HbA1c is:", exp(coef(g3))["OGTT"],"\n")


# 1D
# no OGTT
g4 <- glm(T2D ~ age + sex + bmi, data = t2d, family = "binomial")
noOGTT<-anova(g4,g3,test="Chisq") # 0.007343 **
# no age
g5 <- glm(T2D ~ OGTT + sex + bmi, data = t2d, family = "binomial")
noage1<-anova(g5,g3,test="Chisq") # 1.289015e-26 ***
# no sex
g6 <- glm(T2D ~ OGTT + age + bmi, data = t2d, family = "binomial")
nosex<-anova(g6,g3,test="Chisq") # 1.527991e-09 ***
# no bmi
g7 <- glm(T2D ~ OGTT + age + sex, data = t2d, family = "binomial")
nobmi<-anova(g7,g3,test="Chisq") # 2.239259e-18 ***
cat("1D. I did four tests, each was deducted OGTT, age, sex and bmi. The p-value is", noOGTT$"Pr(>Chi)"[2],noage1$"Pr(>Chi)"[2],nosex$"Pr(>Chi)"[2],nobmi$"Pr(>Chi)"[2],", respectively. All p-value are below 0.05, thus all factors affect T2D.\n")

# Automated F-test-based backward selection using rms::fastbw()
# lm.full <- rms::ols(A ~ B + C, data = data)
# glm.full <- rms::lrm(A ~ B + C, data = data)
# fastbw(modle, rule = "p", sls = 0.1)

# 1E
Zalpha<- function(alpha){
  qnorm(1-alpha/2)
}
my_confi<- function(estimate, SD, alpha){
  return(c(estimate-Zalpha(alpha)*SD,estimate+Zalpha(alpha)*SD))
}
coefg3<-as.data.frame(summary(g3)$coefficients)
# conf int
confintg3<-as.data.frame(rbind(my_confi(coefg3["OGTT","Estimate"],coefg3["OGTT","Std. Error"] ,0.05),my_confi(coefg3["age","Estimate"],coefg3["age","Std. Error"] ,0.05),my_confi(coefg3["sex","Estimate"],coefg3["sex","Std. Error"] ,0.05),my_confi(coefg3["bmi","Estimate"],coefg3["bmi","Std. Error"] ,0.05)),)
# OR
df1e_1<-as.data.frame(cbind(exp(coef(g3))))
df1e_1=df1e_1[-1,]
# merge df
df1e_2<-cbind(df1e_1,confintg3)
colnames(df1e_2)<-c("OR","L 95% CI","U 95% CI")
rownames(df1e_2)<-c("OGTT","age","sex","bmi")
writeLines("1E.")
print(df1e_2)

# 1F
# the probability of a 36 old male, BMI of 24.5, diagnosed with T2D, using a OGTT
res1f<-coef(g3)["age"] * 36 + coef(g3)["sex"] * 1 + coef(g3)["bmi"] * 24.5 + coef(g3)["OGTT"] * 1 + coef(g3)[1]
cat("1F. The probability of a 36 old male, BMI of 24.5, diagnosed with T2D, using a OGTT is",1/(1+exp(-res1f)),"\n")

# 1G
# coef(g3)["sex"] * 1 + coef(g3)["bmi"] * 24.5  = coef(g3)["sex"] * 2 + coef(g3)["bmi"] * BMI
res1g<-(coef(g3)["sex"] * 1 + coef(g3)["bmi"] * 24.5 - coef(g3)["sex"] * 2)/coef(g3)["bmi"]
cat("1G. Her BMI is", res1g,"\n")

# 1H
t2d1 <- t2d
t2d1$bmi_group[t2d1$bmi<25]="<25"
t2d1$bmi_group[t2d1$bmi>=25 & t2d1$bmi<=30]="25-30"
t2d1$bmi_group[t2d1$bmi>30]=">30"
g8<-glm(T2D ~ bmi_group, data = t2d1, family = "binomial")
test1h<-anova(g1,g8,test="Chisq")
cat("1H. The p-value for each group (<25, >30, 25-30): ",summary(g8)$coefficients[1:3,4], ",so BMI groups have an effect on T2D. The p-value from ANOVA (", test1h$"Pr(>Chi)"[2], ") also suggests model T2D ~ bmi_group is better.\n")

# 1I
# code
t2d1$bmi_groupcode[t2d1$bmi_group=="<25"]=0
t2d1$bmi_groupcode[t2d1$bmi_group=="25-30"]=1
t2d1$bmi_groupcode[t2d1$bmi_group==">30"]=4
t2d1$bmi_groupcode<-as.numeric(t2d1$bmi_groupcode)
# model
g9<-glm(T2D ~ bmi_groupcode, data = t2d1, family = "binomial")
# null: bmi_groupcode, alt: original bmi group
test1i<-anova(g9,g8,test="Chisq")
cat("1I. I encode 'normal weight', 'overweight' and 'obese' as 0, 1 and 4. The p-value from ANOVA is", test1i$"Pr(>Chi)"[2],", which means we cannot reject the hypothesis that the effect of being obese is 4 times bigger than the effect of being overweight.\n")

# Part 2
# load data for part2
covid <- read.table("corona.txt", header = T)

# 2A
# censored = alive + deadUnrelated + lost
cat("2A. ",as.data.frame(table(covid$event))[1,2]+as.data.frame(table(covid$event))[3,2]+as.data.frame(table(covid$event))[4,2], "individuals in the study were censored.\n")

# 2B 
library(survival)
# add observation time
covid2<-covid
covid2$days<-covid$eventDay-covid$fistSymptomDay
# model
surv.all <- survfit(Surv(days,event=="diedCorona") ~ 1,data=covid2, conf.int=0.99)
plot(surv.all, ylim=c(0.88,1), ylab="Cumulative survival",xlab = "Time (Days)", main="2B. Overall survival")
legend("topright",c("Confidence interval","Cumulative survival"),lty=c(2,1))
writeLines("2B. Please see plot\n")

# 2C
# mean survival of the non-censored individuals
time<-summary(surv.all, censored=F)$time
time<-c(0,time)
P<-summary(surv.all,censored=F)$surv
P<-c(1,P)
non_cen_mean<-sum(diff(time)*P[-length(P)])
cat("2C. The mean restricted survival time is", as.numeric(survival:::survmean(surv.all, rmean=max(covid2$days))[[1]]["*rmean"]),", the mean survival of the non-censored individuals is",non_cen_mean,"\n" )

# 2D
res2d<-summary(surv.all, times=c(14))
cat("2D. The probability of surviving the first 14 days after first symptoms is", res2d$surv,". The 99% CI is",res2d$lower, res2d$upper,"\n")

# 2E
res2e<-summary(surv.all, times=c(28))
cat("2E. The probability of surviving the next 14 days is", res2e$surv,". The 99% CI is",res2e$lower, res2e$upper,"\n")

# 2F
cox.diabetes<-coxph(Surv(days,event=="diedCorona")~DIABETES,data=covid2)
test2f1<-summary(cox.diabetes)
test2f2<-survdiff(Surv(days,event=="diedCorona")~DIABETES,data=covid2)
cat("2F. I use Log rank test, the p-value is", pchisq(test2f2$chisq, length(test2f2$n)-1, lower.tail = FALSE),", so diabetes affects the survival. The p-value from Cox regression (1.272781e-41) also shows the same conclusion.\n")
# test2f1$coefficients[5]

# 2G
cox.age<-coxph(Surv(days,event=="diedCorona")~ageGroup,data=covid2)
test2g<-summary(cox.age)
cat("2G. The hazard rates of all other age groups are smaller than group >80 (",test2g$coefficients[7,2],"), so all other ages group have lower risk compared to the group aged >80.\n")
#survdiff(Surv(days,event=="diedCorona")~ageGroup,data=covid2)

# Part 3
# load data for part3
birthWeight <- read.table("birthWeight.txt", header = T)
birthWeight$low <- as.factor(birthWeight$low)
#levels(birthWeight$low)
birthWeight$low <- factor(birthWeight$low,levels=c("normalBirthWeight","lowBirthWeight"))

# 3A
g10 <- glm(low ~ 1, data = birthWeight, family = "binomial")
g11 <- glm(low ~ smoke, data = birthWeight, family = "binomial")
test3a<-anova(g10,g11,test="Chisq")
cat("3A. The test is Wald test from summary() function, the p-value is:",summary(g11)$coefficients[2,4], ". So smoking affects the risk of low birth weight. And the p-value of ANOVA test (", test3a$"Pr(>Chi)"[2],") also suggests model low ~ smoke is better. \n")

# 3B
# normalBirthWeight is the "unaffected" and lowBirthWeight is the "affacted"
OR_low<-exp(coef(g11))[2]
# low in smokers
low_smoke<-coef(g11)[1]+coef(g11)[2]
low_nonesmoke<-coef(g11)[1]
p_low_smoke<-1/(1+exp(-low_smoke))
p_low_nonsmoke<-1/(1+exp(-low_nonesmoke))
cat("3B. The odds ratio (OR) is:",OR_low,", and the relative risk for low birth weight in smokers is:",p_low_smoke/p_low_nonsmoke,"\n")

# 3C
g12 <- glm(low ~ race + smoke + ht + age, data = birthWeight, family = "binomial")
cat("3C. The odds ratio for low birth weight in smokers is:", exp(coef(g12))["smoke"],". The 95% CI is", my_confi(summary(g12)$coefficients["smoke","Estimate"], summary(g12)$coefficients["smoke","Std. Error"],0.05),"\n")

# 3D
# white, 36 old smoker, no history of hypertension, low weight
res3d<-coef(g12)["race"] * 1 + coef(g12)["age"] * 36 + coef(g12)["smoke"] * 1 + coef(g12)["ht"] * 0 + coef(g12)[1]
cat("3D. The probability of a 36 old male, BMI of 24.5, diagnosed with T2D, using a OGTT is",1/(1+exp(-res3d)),"\n")

# 3E
# a smoking white, has a history of hypertension,have a risk of low birth weight, 35%.
# 0.35 = 1/{1+exp[-(coef(g12)["race"] * 1 + coef(g12)["age"] * AGE + coef(g12)["smoke"] * 1 + coef(g12)["ht"] * 1 + coef(g12)[1])]}
# 13/7 = exp( -(coef(g12)["race"] * 1 + coef(g12)["age"] * AGE + coef(g12)["smoke"] * 1 + coef(g12)["ht"] * 1 + coef(g12)[1]) )
# -(coef(g12)["race"] * 1 + coef(g12)["age"] * AGE + coef(g12)["smoke"] * 1 + coef(g12)["ht"] * 1 + coef(g12)[1]) = log(13/7)
# -coef(g12)["race"] * 1 - coef(g12)["age"] * AGE - coef(g12)["smoke"] * 1 - coef(g12)["ht"] * 1 - coef(g12)[1] = log(13/7)
# -coef(g12)["age"] * AGE = log(13/7) + coef(g12)["race"] * 1 + coef(g12)["smoke"] * 1 + coef(g12)["ht"] * 1 + coef(g12)[1]
AGE <- (log(13/7) + coef(g12)["race"] * 1 + coef(g12)["smoke"] * 1 + coef(g12)["ht"] * 1 + coef(g12)[1])/(- coef(g12)["age"])
cat("3E. The age is", AGE,"\n")

# 3F
# no age
g13 <- glm(low ~ race + smoke + ht, data = birthWeight, family = "binomial")
noage<-anova(g13,g12,test="Chisq") # 0.2419 reject g8, age no effect
# no ht
g14 <- glm(low ~ race + smoke + age, data = birthWeight, family = "binomial")
noht<-anova(g14,g12,test="Chisq") # 0.04623 *
# no race
g15 <- glm(low ~ smoke + ht + age, data = birthWeight, family = "binomial")
norace<-anova(g15,g12,test="Chisq") # 0.00708 **
# no smoke
g16 <- glm(low ~ race + ht + age, data = birthWeight, family = "binomial")
nosm<-anova(g16,g12,test="Chisq") # 0.002389 **
cat("3F. I did four tests, each was deducted age, ht, race and smoke. The p-value is", noage$"Pr(>Chi)"[2],noht$"Pr(>Chi)"[2],norace$"Pr(>Chi)"[2],nosm$"Pr(>Chi)"[2],", respectively. The tests turn out that hypertension, race and smoking affect low birth weight. Age has no effect.\n")

# 3G
# "2" = African American, "3" = other.
data3G<-birthWeight[birthWeight[,"race"]==2 | birthWeight[,"race"]==3,]
g17<-glm(low ~ smoke + ht + race, data = data3G, family = "binomial")
g18<-glm(low ~ smoke + ht , data = data3G, family = "binomial")
test3g<-anova(g18,g17,test="Chisq")
cat("3G. Given that 'age' is not relevent to weight, it was not included in the model. The p-value of test is", test3g$"Pr(>Chi)"[2],", so we can not reject model low ~ smoke + ht, which means there is no difference in low birth risk between African American individuals and individuals in the category 'other'.")
