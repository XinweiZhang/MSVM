w2 <- c(sqrt(3)/2,-1,0)
w3 <- c(-sqrt(3)/2,-1,0)
w1 <- c(0,0,0) - w2 - w3
x1 <- seq(-3,3,length.out = 200)
x2 <- seq(-3,3,length.out = 200)
par(mfrow = c(2,2))
X <- lapply(x1, function(s){cbind(s,x2)})
X <- do.call(rbind,X)
X <- cbind(X,1)
head(X)
colnames(X) <- c("x1","x2","")
P1 <- X%*%cbind(w1,w2,0)

y <- apply(P1,1,which.max)

plot(X[,1:2], col = y)



P2 <- X%*%cbind(w1,w2,w3)
P2[,3] <- 1 - apply(as.matrix(cbind(P2[,1]+1,0)), MARGIN = 1, FUN = max) - apply(as.matrix(cbind(P2[,2]+1,0)), MARGIN = 1, FUN = max)


y2 <- apply(P2,1,which.max)

plot(X[,1:2], col = y2)

######################
w2 <- c(sqrt(3)/2,-1,0)
w1 <- c(-sqrt(3)/2,-1,0)
w3 <- c(0,0,0) - w1 - w2


P1 <- X%*%cbind(w1,w2,0)

y <- apply(P1,1,which.max)

plot(X[,1:2], col = y)



P2 <- X%*%cbind(w1,w2,w3)
P2[,3] <- 1 - apply(as.matrix(cbind(P2[,1]+1,0)), MARGIN = 1, FUN = max) - apply(as.matrix(cbind(P2[,2]+1,0)), MARGIN = 1, FUN = max)


y2 <- apply(P2,1,which.max)

plot(X[,1:2], col = y2)


###################

par(mfrow = c(1,1))
w2 <- c(sqrt(3)/2,-1,0)
w1 <- c(.2,.8,0)
w3 <- c(0,0,0) - w1 - w2

plot(1, type = "n", xlab="", ylab="", xlim = c(-1, 2), ylim = c(-1, 2), axes=T)
abline(a = -w2[3]/w2[2], b = - w2[1]/w2[2])
abline(a = -w1[3]/w1[2], b = - w1[1]/w1[2])
abline(a = -w3[3]/w3[2], b = - w3[1]/w3[2])
par(mfrow = c(2,2))

par(mfrow = c(1,2))

P1 <- X%*%cbind(w1,w2,w3)
y <- apply(P1,1,which.max)

plot(X[,1:2], col = y)


P2 <- X%*%cbind(w1,w2,w3)
P2[,3] <- 2 - apply(as.matrix(cbind(P2[,1]+1,0)), MARGIN = 1, FUN = max) - apply(as.matrix(cbind(P2[,2]+1,0)), MARGIN = 1, FUN = max)


y2 <- apply(P2,1,which.max)

plot(X[,1:2], col = y2)






################### 3 class
rm(list=ls())
x1 <- seq(-2,6,length.out = 200)
x2 <- seq(-2,6,length.out = 200)
X <- lapply(x1, function(s){cbind(s,x2)})
X <- do.call(rbind,X)
X <- cbind(X,1)
head(X)
colnames(X) <- c("x1","x2","")
par(mfrow = c(1,1))
w1 <- c(-1,1,-1)
w2 <- 1/2*c(1,-1,-1)
w3 <- c(0,0,0) - w2 - w1
par(mfrow = c(1,1))

plot(1, type = "n", xlab="", ylab="", xlim = c(-1, 2), ylim = c(-1, 2), axes=T)
abline(a = -w1[3]/w1[2], b = - w1[1]/w1[2])
abline(a = -w2[3]/w2[2], b = - w2[1]/w2[2])
abline(a = -w3[3]/w3[2], b = - w3[1]/w3[2])



par(mfrow = c(1,2))
P1 <- X%*%cbind(w1,w2,w3)
y <- apply(P1,1,which.max)
plot(X[,1:2], col = y)
# unique(y)

P2 <- X%*%cbind(w1,w2,w3) +1
P2[,3] <- 3 - apply(as.matrix(cbind(P2[,1],0)), MARGIN = 1, FUN = max) - apply(as.matrix(cbind(P2[,2],0)), MARGIN = 1, FUN = max) 

y2 <- apply(P2,1,which.max)

plot(X[,1:2], col = y2)




################### 4 class
rm(list=ls())
x1 <- seq(-2,4,length.out = 200)
x2 <- seq(-2,4,length.out = 200)
X <- lapply(x1, function(s){cbind(s,x2)})
X <- do.call(rbind,X)
X <- cbind(X,1)
head(X)
colnames(X) <- c("x1","x2","")
par(mfrow = c(1,1))
w1 <- c(-1,1,-1)
w2 <- c(1,-1,-1)
w3 <- c(1,1,-1)
w4 <- c(0,0,0) - w2 - w3 - w1
par(mfrow = c(1,1))

plot(1, type = "n", xlab="", ylab="", xlim = c(-1, 2), ylim = c(-1, 2), axes=T)
abline(a = -w1[3]/w1[2], b = - w1[1]/w1[2])
abline(a = -w2[3]/w2[2], b = - w2[1]/w2[2])
abline(a = -w3[3]/w3[2], b = - w3[1]/w3[2])
abline(a = -w4[3]/w4[2], b = - w4[1]/w4[2])
# 
# w12 <- w1-w2
# w23 <- w2-w3
# w34 <- w3-w4
# w14 <- w1-w4
# par(mfrow = c(1,1))
# 
# plot(1, type = "n", xlab="", ylab="", xlim = c(-1, 2), ylim = c(-1, 2), axes=T)
# abline(a = -w12[3]/w12[2], b = - w12[1]/w12[2])
# abline(a = -w23[3]/w23[2], b = - w23[1]/w2[2])
# abline(a = -w34[3]/w34[2], b = - w34[1]/w3[2])
# abline(a = -w14[3]/w14[2], b = - w14[1]/w4[2])



par(mfrow = c(1,2))
P1 <- X%*%cbind(w1,w2,w3,w4)
y <- apply(P1,1,which.max)
plot(X[,1:2], col = y)
unique(y)

P2 <- X%*%cbind(w1,w2,w3,w4) +1
P2[,4] <- 4 - apply(as.matrix(cbind(P2[,1],0)), MARGIN = 1, FUN = max) - apply(as.matrix(cbind(P2[,2],0)), MARGIN = 1, FUN = max) - apply(as.matrix(cbind(P2[,3],0)), MARGIN = 1, FUN = max)
  

y2 <- apply(P2,1,which.max)

plot(X[,1:2], col = y2)

par(mfrow = c(1,1))
plot(1, type = "n", xlab="", ylab="", xlim = c(-1, 2), ylim = c(-1, 2), axes=T)
