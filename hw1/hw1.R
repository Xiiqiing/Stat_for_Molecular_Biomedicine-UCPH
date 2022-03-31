#In order for us to see all the plots you produce, please put this command in the beginning of your file:
par(mfrow = c(2, 2))

# 1A
# at least one time in the 18 sequenced samples -> 1 - 0 times observation
cat("1A. The probability of observing the stop-gain variant at least one time in the 18 sequenced samples is", 1 - dbinom(0, 18 * 2, 0.23))

# # 1B 
# # the probability of observing the mutation in at least 2 individuals
# cat("1B. The probability of observing the mutation in at least 2 individuals is", 1 - pbinom(1, 18, 0.23), "\n")
# 
# # 1C
# test_prob_1C <- (1 - pbinom(1, 0:18, 0.23))
# count = 0
# cat("1C. The minimum number of individuals is",
#     for (i in test_prob_1C) {
#       count = count + 1
#       if (i > 0.9) {
#         break
#       }
#     }
#     , count, ".\n")


# 1B
prob <- 0.23
# the probability of observing the mutation in 0 individual
p_0 <- (1-prob)^36 * prob^0
# the probability of observing the mutation in 1 individual (1 or 2 risk alleles in same individual)
p_1 <- (1-prob)^35 * prob^1 * 36 + (1-prob)^34 * prob^2 * 18
writeLines("\n")
cat("1B. The probability of observing the mutation in at least 2 individuals is", 1 - p_0 - p_1)

# 1C
num_sample = 1:18
count = 0
writeLines("\n")
cat("1C. The minimum number of individuals is",
    for (N in num_sample) {
      count = count + 1
      p0 <- (1-prob)^(N*2) * prob^0
      p1 <- (1-prob)^(N*2-1) * prob^1 * N*2 + (1-prob)^(N*2-2) * prob^2 * N
      p2 = 1 - p0 - p1
      if (p2 > 0.9) {
        break
      }
    }
    , count)

#1D
all_freq <- (14 / (18 * 2))
writeLines("\n")
test_1D <- binom.test(x = 14,n = 36 , p = all_freq, conf.level = 0.99)
cat("1D. The estimated allele frequency is", all_freq, ". The 99% confidence interval is", test_1D$conf.int)

#1E
P_ill <- 0.1
P_ill_homo <- 0.5
P_homo <- 0.23 * 0.23
P_homo_ill = P_ill_homo * P_homo / P_ill
writeLines("\n")
cat("1E. In Greenlandic Inuit the probability of being homozygous if having type 2 diabetes is", P_homo_ill)

#1F
P_non_homo <- 1 - P_homo
P_non_homo_ill <- 1 - P_homo_ill
P_ill_non_homo = P_non_homo_ill * P_ill / P_non_homo
writeLines("\n")
cat("1F. The probability of having type 2 diabetes if not homozygous for the variant is", P_ill_non_homo)

#1G
seeder<-function(name){
  n<-sapply(strsplit(name,""),function(x)
    sum(match(x,c(LETTERS,letters)),na.rm=T))
  set.seed(n)
  N<-2700
  genotype<-rbinom(N,2,0.23)
  glucose <- rnorm(N,5.7,sd=1.4) +
    ifelse(genotype==2,rnorm(N,4,sd=1),0)
  data.frame(stopGain=genotype,glucose=glucose)
}

genData <- seeder("tls868")
myData <- genData
myData$type[myData$stopGain == 0 | myData$stopGain == 1] <- 'nonhomozygous'
myData$type[myData$stopGain == 2] <- 'homozygous'

boxplot(
  glucose ~ type,
  data = myData,
  col = c("mistyrose","lavender"),
  main = "1G. Glucose level of homozygous and nonhomozygous mutation",
  xlab = "Type",
  ylab = "Glucose level"
)
writeLines("\n")
writeLines('1G. Please see plot.')

#1H
type0 <- myData[myData[, 'stopGain'] == 0, ]
type1 <- myData[myData[, 'stopGain'] == 1, ]
type2 <- myData[myData[, 'stopGain'] == 2, ]

calc_CI <- function(data, alpha) {
  xbar <- mean(data)
  SD <- sd(data)
  n <- length(data)
  alpha <- alpha
  Z <- qt(1 - alpha/2, df = n - 1)
  SEM <- SD/sqrt(n) 
  CI1 <- xbar - Z * SEM
  CI2 <- xbar + Z * SEM
  return(c(mean = xbar, SEM = SEM, CI1 = CI1, CI2 = CI2))
}
#or use t.test
#t.test(typeX$glucose, conf.level = 0.99)

result <- as.matrix(cbind(
  calc_CI(type0$glucose, 0.01),
  calc_CI(type1$glucose, 0.01), 
  calc_CI(type2$glucose, 0.01)))
colnames(result) <- c('non-homo(0 mut)','non-homo(1 mut)','homozygous')
writeLines("\n")
writeLines('1H.')
print(result)


#1I
library(gplots)
ci.l <- result['CI1',]
ci.u <- result['CI2',]
col = c("lightblue", "mistyrose",
        "lightcyan")
barplot2(result['mean',], xlab = "Type", ylab = "Mean of glucose levels",
         beside = TRUE, col = col,
         main = "1I. Mean of glucose levels with 99% confidence",
         plot.ci = TRUE, ci.l = ci.l, ci.u = ci.u,
         plot.grid = TRUE)
legend("topleft", legend = colnames(result), fill = col)
writeLines("\n")
writeLines('1I. Please see plot')

#1J
nonhomo <- genData[genData[, 'stopGain'] == 0 | 1,]
homo <- genData[genData[, 'stopGain'] == 2, ]
writeLines("\n")
res_1J <- t.test(homo$glucose, nonhomo$glucose, paired = F, alt = "two.sided")
cat('1J. Null hypothesis is that the difference in means is 0 (H_0: mu_1 - mu2 = 0). Alternative hypothesis is that the difference in means in not 0 (H_1: mu_1 != mu2). P-value (', res_1J$p.value,') is below signiffcance threshold (0.05), reject the null hypothesis. Test statistic is', res_1J$statistic,'.\n')

#1K
qqnorm(homo$glucose,main="1K. Normal QQ-plot for homozygous glucose level")
qqline(homo$glucose,col="red")
qqnorm(nonhomo$glucose,main="1K. Normal QQ-plot for non-homozygous glucose level")
qqline(nonhomo$glucose,col="red")
writeLines("\n")
writeLines('1K. From two qqplots we can see that glucose levels between such two types of individuals are quite different. The results of test hold.')

#1L
t_test_known_var <- function(data1, data2, var1, var2, alpha = 0.05, alt = "2side") {
  mean1 <- mean(data1); mean2 <- mean(data2)
  N1 <- length(data1); N2 <- length(data2)
  S <- sqrt((var1 / N1) + (var2 / N2))
  Z <- (mean1 - mean2 - 0) / S
  p <- if (alt == "2side") {
    2 * pnorm(abs(Z), lower.tail = FALSE)
  } else if (alt == "less") {
    pnorm(Z, lower.tail = TRUE)
  } else { #greater
    pnorm(Z, lower.tail = FALSE)
  }
  CL1 <- (mean1 - mean2 - S * qnorm(1 - alpha / 2))
  CL2 <- (mean1 - mean2 + S * qnorm(1 - alpha / 2))
  value <- c(p.value = p, Stat = Z, CL1 = CL1, CL2 = CL2)
  return(as.data.frame(value))
}
writeLines("\n")
res_1L <- t_test_known_var(homo$glucose, nonhomo$glucose, 1.7, 1.4)
cat('1L. The test statistic is', res_1L['Stat',] ,'. New p-value (', res_1L['p.value',], ') below signiffcance threshold (0.05), and far small than previous one. The results still shows that there is a difference between them.')
