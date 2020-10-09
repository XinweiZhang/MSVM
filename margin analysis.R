generate_noise_circle <- function(n){
  
  t <- sample(c(1,2,3), size = n, replace = T)
  dat <- t(sapply(t, function(z){
    if(z==1){
      xz <- runif(1, min = 0, max = 3)
    }else if(z==2){
      xz <- runif(1, min = 3 , max = 11)
    }else if(z==3){
      xz <- runif(1, min = 11, max = 20)
    }
    x = c(cos(xz*pi/10), sin(xz*pi/10))

    if(z==1){
      y = sample(c(1,2,3), size=1,prob= c(.5, 1/3, 1/6))
      # y = sample(c(1,2,3), size=1,prob= c(.9, 0, .1))
    }else if(z==2){
      y = sample(c(1,2,3), size=1,prob= c(1/3, .5, 1/6))
      # y = sample(c(1,2,3), size=1,prob= c(.2, .7, .1))
    }else if(z==3){
      y = sample(c(1,2,3), size=1,prob= c(.05, .05, .9))
      # y = sample(c(1,2,3), size=1,prob= c(1/6, 1/3, .5))
    }
   
    return(c(x,y))
  }))
  X <- dat[,1:2]
  y <- dat[,3]
  list(X=X,y=y, oracle_acc = mean(t==y), t=t)
}


##############Three Class################

set.seed(123)
n <- 1000
oracle_data <- generate_noise_circle(10^5)
train_data <- generate_noise_circle(n)
train_data$oracle_acc
WW_fit <- WW_pri_opt(train_data$X,train_data$y, C=1, lambda=0, intercept = F)
CS_fit <- CS_pri_opt(train_data$X,train_data$y,  C=1, lambda=0, intercept = F)
Duchi_fit <- Duchi_pri_opt(train_data$X,train_data$y,  C=1, lambda=0, intercept = F)
MDuchi_fit <- MDuchi_pri_opt(train_data$X,train_data$y, C=1, lambda=0, intercept = F)

y <- train_data$y
X <- train_data$X
class_idx <- sort(unique(y))
Y <- sapply(class_idx, function(id){as.numeric(y==id)})
m <- length(class_idx)
t <- train_data$t

c(mean(predict(WW_fit, oracle_data$X) == oracle_data$y),
mean(predict(CS_fit, oracle_data$X) == oracle_data$y),
mean(predict(Duchi_fit, oracle_data$X) == oracle_data$y),
mean(predict(MDuchi_fit, oracle_data$X) == oracle_data$y))



WW_slack <- pmax(1-((cbind(train_data$X,1)%*%t(WW_fit$beta)*Y)%*%matrix(1, nrow = m, ncol = m) - cbind(train_data$X,1)%*%t(WW_fit$beta)),0)
CS_slack <- pmax(1-((cbind(train_data$X,1)%*%t(CS_fit$beta)*Y)%*%matrix(1, nrow = m, ncol = m) - cbind(train_data$X,1)%*%t(CS_fit$beta)),0)
Duchi_slack <- pmax(1-((cbind(train_data$X,1)%*%t(Duchi_fit$beta)*Y)%*%matrix(1, nrow = m, ncol = m) - cbind(train_data$X,1)%*%t(Duchi_fit$beta)),0)
MDuchi_slack <- pmax(1-((cbind(train_data$X,1)%*%t(MDuchi_fit$beta)*Y)%*%matrix(1, nrow = m, ncol = m) - cbind(train_data$X,1)%*%t(MDuchi_fit$beta)),0)


par(mfrow=c(2,2))
boxplot(WW_slack[y==1 & t==1, -1], main = paste("WW",round(mean(predict(WW_fit, oracle_data$X)[oracle_data$y==1 & t ==1]==1),3)), names=c("1vs2", "1vs3"))
boxplot(CS_slack[y==1& t==1,-1], main = paste("CS",round(mean(predict(CS_fit, oracle_data$X)[oracle_data$y==1& t ==1]==1),3)), names=c("1vs2", "1vs3"))
boxplot(Duchi_slack[y==1& t==1,-1], main = paste("Duchi", round(mean(predict(Duchi_fit, oracle_data$X)[oracle_data$y==1& t ==1]==1),3)), names=c("1vs2", "1vs3"))
boxplot(MDuchi_slack[y==1& t==1,-1], main = paste("Mduchi", round(mean(predict(MDuchi_fit, oracle_data$X)[oracle_data$y==1& t ==1]==1),3)), names=c("1vs2", "1vs3"))

boxplot(WW_slack[y==1 & t==2, -1], main = paste("WW",round(mean(predict(WW_fit, oracle_data$X)[oracle_data$y==1]==1),3)), names=c("1vs2", "1vs3"))
boxplot(CS_slack[y==1& t==2,-1], main = paste("CS",round(mean(predict(CS_fit, oracle_data$X)[oracle_data$y==1]==1),3)), names=c("1vs2", "1vs3"))
boxplot(Duchi_slack[y==1& t==2,-1], main = paste("Duchi", round(mean(predict(Duchi_fit, oracle_data$X)[oracle_data$y==1]==1),3)), names=c("1vs2", "1vs3"))
boxplot(MDuchi_slack[y==1& t==2,-1], main = paste("Mduchi", round(mean(predict(MDuchi_fit, oracle_data$X)[oracle_data$y==1]==1),3)), names=c("1vs2", "1vs3"))

boxplot(WW_slack[y==1 & t==3, -1], main = paste("WW",round(mean(predict(WW_fit, oracle_data$X)[oracle_data$y==1]==1),3)), names=c("1vs2", "1vs3"))
boxplot(CS_slack[y==1& t==3,-1], main = paste("CS",round(mean(predict(CS_fit, oracle_data$X)[oracle_data$y==1]==1),3)), names=c("1vs2", "1vs3"))
boxplot(Duchi_slack[y==1& t==3,-1], main = paste("Duchi", round(mean(predict(Duchi_fit, oracle_data$X)[oracle_data$y==1]==1),3)), names=c("1vs2", "1vs3"))
boxplot(MDuchi_slack[y==1& t==3,-1], main = paste("Mduchi", round(mean(predict(MDuchi_fit, oracle_data$X)[oracle_data$y==1]==1),3)), names=c("1vs2", "1vs3"))


par(mfrow=c(2,2))
boxplot(WW_slack[y==2,-2], main = paste("WW",round(mean(predict(WW_fit, oracle_data$X)[oracle_data$y==2]==2),3)), names=c("2vs1", "2vs3"))
boxplot(CS_slack[y==2,-2], main = paste("CS",round(mean(predict(CS_fit, oracle_data$X)[oracle_data$y==2]==2),3)), names=c("2vs1", "2vs3"))
boxplot(Duchi_slack[y==2,-2], main = paste("Duchi", round(mean(predict(Duchi_fit, oracle_data$X)[oracle_data$y==2]==2),3)), names=c("2vs1", "2vs1"))
boxplot(MDuchi_slack[y==2,-2], main = paste("Mduchi", round(mean(predict(MDuchi_fit, oracle_data$X)[oracle_data$y==2]==2),3)), names=c("2vs1", "2vs1"))

par(mfrow=c(2,2))
boxplot(WW_slack[y==3,-3], main = paste("WW",round(mean(predict(WW_fit, oracle_data$X)[oracle_data$y==3]==3),3)), names=c("3vs1", "3vs2"))
boxplot(CS_slack[y==3,-3], main = paste("CS",round(mean(predict(CS_fit, oracle_data$X)[oracle_data$y==3]==3),3)), names=c("3vs1", "3vs2"))
boxplot(Duchi_slack[y==3,-3], main = paste("Duchi", round(mean(predict(Duchi_fit, oracle_data$X)[oracle_data$y==3]==3),3)), names=c("3vs1", "3vs2"))
boxplot(MDuchi_slack[y==3,-3], main = paste("Mduchi", round(mean(predict(MDuchi_fit, oracle_data$X)[oracle_data$y==3]==3),3)), names=c("3vs1", "3vs2"))

