#In order for us to see all the plots you produce, please put this command in the beginning of your file:
par(mfrow=c(3,2))

# 1A
writeLines("1A. The mean and median of Sepal width of 3 species ")
aggregate(Sepal.Length ~ Species, data = iris, 
          FUN = function(x) c(mean = mean(x), med = median(x)))

#1B 
var_1b <- aggregate(Sepal.Length ~ Species, data = iris, var)
cat("1B. The highest variance of Sepal Length among 3 species is",
    levels(iris$Species)
    [var_1b[var_1b$Sepal.Length == max(var_1b$Sepal.Length), 1]], "\n")

#1C
boxplot(
  Petal.Width ~ Species,
  data = iris,
  #range = 0
  main = "1C. Petal Width of 3 species",
  xlab = "Species",
  ylab = "Petal Width"
)  

#1D
colours <- c("orange", "black", "grey")
pch_ind <- c(15, 0, 0)
plot(
  iris$Sepal.Width,
  iris$Petal.Width,
  main = "1D. Width",
  xlab = "Sepal(cm)",
  ylab = "Petal(cm)",
  col = colours[iris$Species],
  pch = pch_ind[iris$Species]
)
legend("topright", legend = levels(iris$Species), 
       pch = c(15, 0, 0), col = colours)

# add mean and median
col_ind = 0
col_bor = c("black", "white", "black")
for (i in levels(iris$Species)) {
  iris_new <- iris[iris$Species == i, ]
  col_ind = col_ind + 1
  points(mean(iris_new$Sepal.Width),mean(iris_new$Petal.Width), 
           pch = 21, cex = 2.5, lwd = 2, bg = colours[col_ind])
  points(median(iris_new$Sepal.Width),median(iris_new$Petal.Width), 
         pch = 24, cex = 1.5, lwd = 2, col = col_bor[col_ind], bg = colours[col_ind])
}
legend("topleft", legend = c("mean", "median"), 
       pch = c(21, 17), col = "black")
