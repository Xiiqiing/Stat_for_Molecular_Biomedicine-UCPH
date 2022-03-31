par(mfrow=c(2,2))
# load data for part1
intron_gc<-read.table("intron_gc.txt",header=F)
exon_gc<-read.table("exon_gc.txt",header=F)

# 1a
# Parametric test null hypothesis: no difference, mean, same gene.
test_1a<-t.test(intron_gc[,1], exon_gc[,1], paired = T)
cat("1a. Since each line corrseponds to the same gene, I chosen paired t test. 
    The p-value is", test_1a$p.value, ", reject null hypothesis, so they are different.\n" )

# 1b
# Non-parametric test for same question 1a: no difference, median, same gene
test_1b<-wilcox.test(intron_gc[,1], exon_gc[,1], paired = T)
cat("1b. Since each line corrseponds to the same gene, I chosen Wilcoxon signed rank test. 
    The p-value is", test_1b$p.value, ", still reject null hypothesis, so they are different.\n" )

# 1c
writeLines("1c. Advantages: 
    1. Useful when lack/less assumptions. 
    2. Accept small sample sizes. 
    3. Can be used for almost all data types, even unknown distribution.
    Disadvantages: Less powerful than parametric tests when assumptions hold.")

# 1d
O = c(390,1000,1149,872,481,205,90)
totalregions = sum(O*c(0,1,2,3,4,5,6.5))
meanregions = totalregions/4187
cat("1d. The mean number of regions with GC repeats in a gene is" ,meanregions,"\n")

# 1e
# according 1d, 1 gene has 2.256508 regions with CG repeats in average
# each exonic region may has a GC repeat
# there are 15 exonic regions in a gene, so max: 15 repeats
# a single exomic region has a GC repeat: mean/15
# mean = np
prob = meanregions/15
cat("1e. The probability of a single exomic region has a GC repeat is",meanregions/15,".
    The probability of observing more than 6 exonic GC regions in a gene is", 1 - pbinom(6,15,prob),"\n")

# 1f
# expected
E <- c(dbinom(0:5,15,prob)*4187,(1-sum(dbinom(0:5,15,prob)))*4187)
res_1f<-as.data.frame(rbind(O,E))
names(res_1f)<-c(0,1,2,3,4,5,'<=6')
writeLines("1f. Table of 1f:")
print(res_1f)

# 1g
X2 = sum((O-E)^2/E)
# df = K - 1 - number of estimated parameters = 7-1-1 (one para: mean)
df = 7-1-1
cat("1g. Test statics:",X2,", p-value:",1-pchisq(X2,df=df),", the number of degrees of freedom:", df,"\n")

# 1h
writeLines("1h. P-value is smaller than 0.05, thus we reject the null hypothesis: the data appear to follow binomial distribution.\n")

# load data for part2
condi<-read.table("conditions.txt",header=T)
condi_2<-read.table("conditions.txt",header=T)
# 2a.
writeLines("2a. Test: One way analysis of variance. 
    Null hypothesis: all means are equal, exercise intensity has no effect on condition. 
    Alternative hypothesis: not all means are equal, they are dependent.\n")

# 2b.
test_2b<-anova(lm(change~1,data=condi),lm(change~level,data=condi))
cat("2b. P-value:", test_2b$"Pr(>F)"[2],". Reject null hypothesis, change of physical condition does depend on exercise intensity.\n")

# 2c.
condi$group <- factor(paste(condi$age,condi$level))
means <- tapply(condi$change,condi$group,mean)
#print(means)
plot(means[1:3],pch=16,col="black", cex=1.2, main="2c. Group mean of physical condition for nine combinations", xlab="Exercise intensity", ylab = "Mean of changes of physical conditions", ylim = c(2,7), xaxt = "n")
axis(1, c(1:3), c("high","low","moderate"))
points(means[4:6],pch=16,cex=1.2,col="red")
points(means[7:9],pch=16,cex=1.2,col="blue")
legend("topright",c("age: <25","age: 25-40","age: >40"), col=c("red","black","blue"),pch=16)
writeLines("2c. Please see plot\n")

# 2d.
a<-lm(formula = change ~ age + level, data=condi)
cat("2d. The Multiple R-squared (Adjusted):", summary(a)$adj.r.squared,"\n")

# 2e.
# writeLines("2e. Design matrix:")
# print(head(cbind(model.matrix(a),condi_2)))
# condi_change=5.9577984-0.9286240(if age>40)-0.6320964(if age25-40)-2.4252735(if l)-1.5180187(if m)
cat("2e. The predicted change of physical condition for a 36 year old after moderate exercise intensity is", a$coefficients[1]+a$coefficients[3]+a$coefficients[5],"\n")

# 2f.
b <- lm(formula = change ~ level * age, data=condi)
test_2f <- anova(a,b)
cat("2f. I use the two way ANOVA with interaction test. The p-value is", test_2f$"Pr(>F)"[2],"\n")

# 2g.
writeLines("2g. Thus, there is a significant interaction between the effect of exercise intensity and age. 
    They are not independent of each other. The effect of exercise intensity on physical condition depends on age.\n")

# 2h.
plot(b,which=c(2,1), main = "2g.")
writeLines("2h. The qq plot and residuals plot looks not bad. So the basic assumptions of this ANOVA test are hold. 
    From of plot 2c, the previous additive assumption is rejectd, and interaction assumption of this test is satisffed. \n")
