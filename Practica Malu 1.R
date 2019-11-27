## TASK 2 SCRIPT ##

rm(list=ls())
setwd("~/UNI/Master/Estadistica/Task 2 practice statistics-20191117")
bladder <- read.csv(file = "bladder.1.txt",header = TRUE,sep = "")
head(bladder)

#1.	Convert variables "y" and "gender" to factor variables
#2.	Define the labels for variables "y" and "gender".

bladder$y <- factor(bladder$y, levels = c(0,1),labels = c("control", "case"))
bladder$gender <- factor(bladder$gender, levels = c(1,2), labels = c("male","female"))
head(bladder)

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

plot(bladder$gender, bladder$gene1)

#5.	Test the normality of the expression levels of the 20 genes (use function apply). 
# How many genes are not normally distributed and which are these genes?

# normality <- sapply(bladder[4:23],2, shapiro.test)

# normality <- (sapply(bladder[4:23], shapiro.test))
# nonorm2 <- normality<0.05
# nonorm2 <- subset(nonorm2, nonorm2["p.value"=="TRUE",])
# nonorm2

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
M<-function(data,i){
  m<-mean(data[i])
  m
}
M(bladder$gene1)
boot.out<-boot(bladder$gene1, M, R=1000) 
boot.out
boot.ci(boot.out, conf = 0.95, type = "perc")
quantile(M(bladder$gene1),c(0.20,0.80))

#9.	Contrast using a permutation test whether the 20th percentile of gene1 for 
# male and female are equal.


#10.	Perform a nonparametric test for association of gender and the risk of 
# disease. Provide the OR (change the levels of "gender" if necessary in order 
# that the given OR is larger than 1)

gen.dis.ass<-glm(bladder$gender~bladder$y, family=binomial()) 
gen.dis.ass
summary(gen.dis.ass)
exp(summary(gen.dis.ass)$coefficients[2,1])
#The OR is 1.263

tableOR<-table(bladder$gender,bladder$y)
tableOR
oddsratio<-231*84/(226*68)
oddsratio
#We confirm the OR found before

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



library(ROCR)
