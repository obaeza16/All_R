## TASK 2 SCRIPT ##

rm(list=ls())
setwd("D:/Master/All_R-master/All_R-master")
bladder <- read.csv(file = "bladder.1.txt",header = TRUE,sep = "")
head(bladder,n = 3)
options(digits = 4)

#1.	Convert variables "y" and "gender" to factor variables
#2.	Define the labels for variables "y" and "gender".

bladder$y <- factor(bladder$y, levels = c(0,1),labels = c("control", "case"))
bladder$gender <- factor(bladder$gender, levels = c(1,2), labels = c("male","female"))
head(bladder,n = 3)

#3.	Build a frequency table for "y" that contains, both, the 
# absolute and the relative frequencies. The same for "gender".
#Table for variable y:

freq.y <- table(bladder$y)  #absolute frequency
freq.y   
relfreq.y <- prop.table(table(bladder$y))   #relative frequency
relfreq.y   
totfreq.y <- cbind(freq.y,relfreq.y)   #unite both tables
totfreq.y

#Table for gender:

freq.gen <- table(bladder$gender)  #absolute frequency
freq.gen   
relfreq.gen <- prop.table(table(bladder$gender))   #relative frequency
relfreq.gen   
totfreq.gen <- cbind(freq.gen,relfreq.gen)   #unite both tables
totfreq.gen

#4.	Plot gene1 levels as a function of "gender":

plot(bladder$gender, bladder$gene1,xlab="Gender",ylab="Gene1 Expression",
     main="Gene1 Expression Levels as a function of gender")

#5.	Test the normality of the expression levels of the 20 genes (use function apply). 
# How many genes are not normally distributed and which are these genes?


pval <- do.call(rbind, lapply(bladder[4:23], function(x) shapiro.test(x)[c("p.value")]))
nonorm <- (pval<0.05)
nonorm <- subset(nonorm, nonorm == "TRUE")
nonorm   #The genes not normally distributed are 2: gene 8 and gene 18

#6.	Test whether mean expression levels of gene1 and gene2 are equal. 
# Note: this is a test for equality of means with paired samples. 
# This test is performed by computing the difference of the two variables 
# (gene1-gene2) and testing whether the mean of the difference is equal to zero.

diff <- bladder$gene1-bladder$gene2
t.test(diff,mu = 0)   #The mean expression of gene1 and gene2 are equal

#7.	Test if the mean expression levels of gene1 are equal between cases 
# and controls

var.test(bladder$gene1~bladder$y)
t.test(bladder$gene1~bladder$y, var.equal=T) 
#The mean expression of gene1 and gene2 are equal

#8.	Obtain a 95% bootstrap confidence interval for the 20th percentile of gene1 

library(boot)
bootfun<-function(data,N){
  quant_20 <- NULL
  for (i in 1:N) {
    newsample <- sample(data,length(data),replace = T)
    quant_20 <- c(quant_20, quantile(newsample,probs = 0.2))
  }
return(quant_20)  
}

boot20 <- bootfun(bladder$gene1,10000)
confidn <- c(quantile(boot20,0.025),quantile(boot20,0.975))
confidn
quantile(bladder$gene1,.2)


#9.	Contrast using a permutation test whether the 20th percentile of gene1 for 
# male and female are equal.

bootfun_null<-function(data,N){
  data_s <- sample(data)
  quant_20 <- NULL
  for (i in 1:N) {
    newsample <- sample(data_s,length(data_s),replace = T)
    quant_20 <- c(quant_20, quantile(newsample,probs = 0.2))
  }
  return(quant_20)  
}

boot20perm <- bootfun_null(bladder$gene1,10000)
hist(boot20perm)
p.value<-(sum(boot20perm>boot20)+1)/(length(boot20perm)+1)
p.value


#10.	Perform a nonparametric test for association of gender and the risk of 
# disease. Provide the OR (change the levels of "gender" if necessary in order 
# that the given OR is larger than 1)

install.packages("epitools")
library(epitools)
table_OR <- table(bladder$gender,bladder$y)
table_OR
oddsratio(table_OR)

gen.dis.ass <- glm(bladder$gender~bladder$y, family=binomial()) 
gen.dis.ass
summary(gen.dis.ass)
exp(summary(gen.dis.ass)$coefficients[2,1])
#The OR is 1.263
#No te sentit ja que no estan associades


#11.	Explore for possible relationship between methylation and gene expression.

meth.gene.ass <- glm(bladder$methyl~bladder$gene1+bladder$gene2+bladder$gene3
                     +bladder$gene4+bladder$gene5+bladder$gene6+bladder$gene7
                     +bladder$gene8+bladder$gene9+bladder$gene10+bladder$gene11
                     +bladder$gene12+bladder$gene13+bladder$gene14+bladder$gene15
                     +bladder$gene16+bladder$gene17+bladder$gene18+bladder$gene19
                     +bladder$gene20, family=binomial())
summary(meth.gene.ass)


#12.	Identify genes that are related to the risk of bladder cancer using a multivariate 
# logistic regression model with stepwise variable selection. Denote the selected model as 
# "best.model". Interpret the obtained model.

dis.gene.ass <- glm(bladder$y~bladder$gene1+bladder$gene2+bladder$gene3
                    +bladder$gene4+bladder$gene5+bladder$gene6+bladder$gene7
                    +bladder$gene8+bladder$gene9+bladder$gene10+bladder$gene11
                    +bladder$gene12+bladder$gene13+bladder$gene14+bladder$gene15
                    +bladder$gene16+bladder$gene17+bladder$gene18+bladder$gene19
                    +bladder$gene20, family=binomial())

forwardmodel<-step(dis.gene.ass,direction="forward")
summary(forwardmodel)

backwardmodel<-step(dis.gene.ass,direction="backward")
summary(backwardmodel)

bothmodel<-step(dis.gene.ass,direction="both")
summary(bothmodel)

best.model <- glm(bladder$y~bladder$gene3+bladder$gene5+bladder$gene8+
                    bladder$gene17+bladder$gene19,family = binomial())

summary(best.model)

# 13.	Analyze the classification ability of "best.model" (ROC curve and AUC) 
# according to the following schemes:
# c.	Explain why both the "a" and "b" overestimate the classification accuracy of "best.model". 
# Suggest a scheme for determining an unbiased estimate of the classification accuracy of "best.model".

#Se sobreesima pq utilitzes les mateixes dades per constriur el model i per testejarlo
#Podriem reservarnos unes dades per constriurlo i unes altres per testejar-lo.

# a.	Apparent validation of "best.model" using the same data that was used for model building.

install.packages("ROCR")
library(ROCR)
lp<-best.model$linear.predictors
pred <- prediction(-lp, bladder$y)
perf <- performance(pred, "tpr", "fpr" )
plot(perf)
abline(a=0, b= 1)
title("ROC curve for best.model; bladder.1 data")

auc<-slot(performance(pred,"auc"), "y.values")[[1]]
auc

# b.	Crossvalidation with k = 5 for "best.model".

K <- 5
n <- nrow(bladder)  #number of individuals
fold <- sample(as.numeric(cut((1:n),breaks = K)))  # random assignment of each individual into one fold 

pred <- NULL # vector of predictions

for(i in 1:K){
  indTest <- which(fold==i)   # Test indices 
  indTrain <- which(fold!=i)  # Train indices
  model.i <- glm(bladder$y~bladder$gene3+bladder$gene5+bladder$gene8+
                   bladder$gene17+bladder$gene19, data=bladder[indTrain,], family=binomial())  # Adjust the model with training data
  pred.i <- predict(model.i, newdata=bladder[indTest, ])   # Predicts test data at step i
  pred[indTest] <- pred.i   # Store predicted values for test data at step i 
}  

prediction <- prediction(pred, bladder$y)
perf <- performance(prediction, "tpr", "fpr" )
plot(perf,main="Crossvalidation for best.model with k=5")
abline(a=0, b= 1)

auc<-slot(performance(prediction,"auc"), "y.values")[[1]]
auc

# 14.	For each variable "genei", i=1, ., 20, define a new factor 
# variable "levels.genei" with values "high" if gene levels are positive 
# or zero and "low" if gene levels are negative. 

bladder2 <- bladder[,4:23]
genes <- ifelse(bladder2>=0,"high","low")
colnames(genes) <- paste(rep("levels.gene",20),seq(1,20,1),sep="")
bladder2 <- data.frame(genes)
bladder2$y <- bladder[,2]
head(bladder2,n = 3)
str(bladder2$levels.gene1)

# 15.	Perform a nonparametric test of association between each variable "levels.genei" 
# and the risk of disease. Adjust the p-values for multiple testing according to an fdr 
# threshold equal to 0.1. Interpret the results.

xy <- data.frame(x1=c(3,2,3,4,5,6),x2=c(1,2,3,4,5,6))
result <- apply(X = xy,MARGIN = 1,FUN = chisq.test)

chisq.test(bladder2$levels.gene1,bladder2$y,correct = TRUE)
result <- apply(bladder2, 1, FUN = chisq.test(,bladder2$y))
bladder2 <- as.matrix(bladder2)
chisq.test(table(bladder2))

str(bladder2)
chisq.test(bladder2)
result <- apply(X = table(bladder2),MARGIN = 1,FUN = chisq.test)

# 16.	Using the last 30 variables corresponding to gene expression levels
# in three different pathways, perform a clustering analysis (hierarchical
# and k-means) and explore groups of genes that have similar expression. 
# Explain the results. 

dist1 <- dist(bladder[, c(25:50)],method="euclidian")
diam <- as.matrix(dist1)
dim(diam)

plot(hclust(dist1,method="average"))
plot(hclust(dist1,method="ward.D2"))

kmeans2<-kmeans(dist1,2)
plot(dist1, col = kmeans2$cluster)
points(kmeans2$centers,col = 1:2, pch = 8, cex=50)

data16<-t(bladder[, c(25:50)])
distind <- dist(data16,method="euclidian")
distind <- as.matrix(distind)
dim(distind)

plot(hclust(dist(data16,method="euclidian"),method="ward.D2"))
plot(hclust(dist(data16,method="euclidian"),method="average"))

kmeans3<-kmeans(data16,3)
plot(data16, col = kmeans3$cluster)
points(kmeans3$centers,col = 1:3, pch = 8, cex=2)
