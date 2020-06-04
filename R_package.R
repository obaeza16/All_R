GSAfisher <- function(obj,...) {
  class(obj) <- "GSAfisher"
  UseMethod("GSAfisher")
}

GSAfisher.default <- function(obj,...) {
  class(obj) <- "GSAfisher"
  fish <- sum(-2*log(obj))
  return(fish)
}

print.GSAfisher <- function(obj,...) {
  value <- sum(-2*log(obj))
  cat("The agrupated p-value from your object is",value)
}

summary.GSAfisher <- function(obj,...) {
  media <- mean(obj)
  sdev <- sd(obj)
  len <- length(obj)
  in_max <- which.max(obj)
  va_max <- max(obj)
  in_min <- which.min(obj)
  va_min <- min(obj)
  cat("The mean of your p values is",media,"\n")
  cat("The standard deviation of your p values is",sdev,"\n")
  cat("The length of your object is",len,"\n")
  cat("The maximum value of your p values is ",va_max,"in position",in_max,"\n")
  cat("The minimum value of your p values is ",va_min,"in position",in_min,"\n")
  boxplot(obj,main="Boxplot of p-values in your object",ylab="p values")
}

p.val <- runif(1000,0,0.25)


class(p.val) <- "GSAfisher"
summary(p.val<0.05)
GSAfisher(p.val)
print(p.val)
summary(p.val)



output <- data.frame()
for (i in 1:length(p.val)) {
  newpvalue <- 1 - pchisq(-2 * sum(log(p.val[i])), df = 2 * length(p.val[i]))
  output <- rbind(output,newpvalue)
}


nrSamples=5000
yourlist <- sapply(1:nrSamples, function(x) runif(1000,0,0.25))
class(yourlist) <- "GSAfisher"
yourmeans <- sapply(yourlist, GSAfisher)

summary(yourlist)
GSAfisher(yourlist)
