################ Functions#########
library(CVXR)
library(MASS)
library(R.matlab)
library(e1071)
library(kernlab)
library(ggplot2)
library(ggpubr)

rm(list=ls())
setwd("~/Desktop/Multiclass Classification/MSVM Code")
source("primary form functions.R")
source("kernel primary form functions.R")
set.seed(123)

plot_circle_boundary <- function(models, X,y, title = NULL, xlim =NULL, ylim = NULL, rule = "simple_max"){
  X_dat <- as.data.frame(X)
  colnames(X_dat) <- c("x","y")
  X_dat$class <- as.factor(y)
  if(is.null(xlim)){
    xlim = c( floor(min(X_dat$x,X_dat$y))-.1, ceiling(max(X_dat$x,X_dat$y))+.1)
  }
  if(is.null(ylim)){
    ylim = c( floor(min(X_dat$x,X_dat$y))-.1, ceiling(max(X_dat$x,X_dat$y))+.1)
  }
  
  nd.r1 = seq(from = 0.85, 0.9, length.out = 30)
  nd.r2 = seq(from = 0.7, 0.75, length.out = 30)
  nd.r3 = seq(from = 0.55, 0.6, length.out = 30)
  nd.r4 = seq(from = 0.4, 0.45, length.out = 30)
  nd.r5 = seq(from = 0.25, 0.3, length.out = 30)
  nd.theta = seq(from = 0, 2*pi, length.out = 500)
  nd.x = cbind(cos(nd.theta), sin(nd.theta))
  circ.theta = seq(from = 0, 2*pi, length.out = 2000)
  circ.x = cbind(cos(circ.theta), sin(circ.theta))
  circ.nd = as.data.frame(0.97 * circ.x)
  colnames(circ.nd) <- c("x","y")
  nd1 = as.data.frame(do.call(rbind,lapply(nd.r1, function(r){r*nd.x})))
  nd2 = as.data.frame(do.call(rbind,lapply(nd.r2, function(r){r*nd.x})))
  nd3 = as.data.frame(do.call(rbind,lapply(nd.r3, function(r){r*nd.x})))
  nd4 = as.data.frame(do.call(rbind,lapply(nd.r4, function(r){r*nd.x})))
  nd5 = as.data.frame(do.call(rbind,lapply(nd.r5, function(r){r*nd.x})))
  
  # nd1$class <- pred(as.matrix(nd1), beta[[1]], rule = rule)
  # nd2$class <- pred(as.matrix(nd2), beta[[2]], rule = rule)
  # nd3$class <- pred(as.matrix(nd3), beta[[3]], rule = rule)
  # nd4$class <- pred(as.matrix(nd4), beta[[4]], rule = rule)
  # nd5$class <- pred(as.matrix(nd5), beta[[5]], rule = rule)

  nd1$class <- rep(predict(models[[1]], nd.x, rule = rule),30)
  nd2$class <- rep(predict(models[[2]], nd.x, rule = rule),30)
  nd3$class <- rep(predict(models[[3]], nd.x, rule = rule),30)
  nd4$class <- rep(predict(models[[4]], nd.x, rule = rule),30)
  nd5$class <- rep(predict(models[[5]], nd.x, rule = rule),30)
  
  nd <-  rbind(nd1,nd2,nd3,nd4,nd5)
  colnames(nd) <- c("x","y","class")
  nd$class <- as.factor(nd$class)
  seg1 <- 0.97*data.frame(x1 = c(0,0,0), x2 = c(0,1,cos(11/10*pi)), y1 = c(0,0,0), y2 = c(1,0,sin(11/10*pi)))
  
  plt <- ggplot(nd, aes(x=x, y=y, color = class)) +
    geom_point( shape = 20) +
    geom_point(data = X_dat, size = 1) +
    geom_point(data = circ.nd, size = .5, colour = "black") +
    geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = seg1,  size = 1,linetype = 1, color="black") +
    scale_x_continuous(limits = xlim, expand=c(0,0)) +
    scale_y_continuous(limits = ylim, expand=c(0,0)) +
    ggtitle(title) +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                            panel.background = element_blank(), axis.line = element_blank(),
                            axis.ticks = element_blank(),
                            axis.text.x = element_blank(),
                            axis.text.y = element_blank(),
                            plot.title = element_text(hjust = 0.5)) + 
    scale_color_manual(values=c("1" ="blue","2"= "green", "red"))
  
  suppressWarnings(print(plt))
}

draw_X <- function(y){
  t <- sapply(y, function(z){
    if(z==1){
      runif(1, min = 0, max = 5)
    }else if(z==2){
      runif(1, min = 5, max = 11)
    }else if(z==3){
      runif(1, min = 11, max = 20)
    }
  })
  
  X = t(sapply(t, function(x){c(cos(x*pi/10), sin(x*pi/10))}))
  return(X)
}


plot_ridgeless_circle_boundary <- function(model, X,y, title = NULL, xlim =NULL, ylim = NULL, rule = "simple_max"){
  X_dat <- as.data.frame(X)
  colnames(X_dat) <- c("x","y")
  X_dat$class <- as.factor(y)
  if(is.null(xlim)){
    xlim = c( floor(min(X_dat$x,X_dat$y))-.1, ceiling(max(X_dat$x,X_dat$y))+.1)
  }
  if(is.null(ylim)){
    ylim = c( floor(min(X_dat$x,X_dat$y))-.1, ceiling(max(X_dat$x,X_dat$y))+.1)
  }
  
  nd.r1 = seq(from = 0.5, 0.9, length.out = 100)
  
  nd.theta = seq(from = 0, 2*pi, length.out = 500)
  nd.x = cbind(cos(nd.theta), sin(nd.theta))
  circ.theta = seq(from = 0, 2*pi, length.out = 2000)
  circ.x = cbind(cos(circ.theta), sin(circ.theta))
  circ.nd = as.data.frame(0.97 * circ.x)
  colnames(circ.nd) <- c("x","y")
  nd1 = as.data.frame(do.call(rbind,lapply(nd.r1, function(r){r*nd.x})))
  
  nd1$class <- rep(predict(model[[1]], nd.x, rule = rule),100)
  
  
  nd <-  nd1
  colnames(nd) <- c("x","y","class")
  nd$class <- as.factor(nd$class)
  seg1 <- 0.97*data.frame(x1 = c(0,0,0), x2 = c(0,1,cos(11/10*pi)), y1 = c(0,0,0), y2 = c(1,0,sin(11/10*pi)))
  
  plt <- ggplot(nd, aes(x=x, y=y, color = class)) +
    geom_point( shape = 20) +
    geom_point(data = X_dat, size = 1) +
    geom_point(data = circ.nd, size = .5, colour = "black") +
    geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = seg1,  size = 1,linetype = 1, color="black") +
    scale_x_continuous(limits = xlim, expand=c(0,0)) +
    scale_y_continuous(limits = ylim, expand=c(0,0)) +
    ggtitle(title) +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                            panel.background = element_blank(), axis.line = element_blank(),
                            axis.ticks = element_blank(),
                            axis.text.x = element_blank(),
                            axis.text.y = element_blank(),
                            plot.title = element_text(hjust = 0.5)) + 
    scale_color_manual(values=c("1" ="blue","2"= "green", "red"))
  
  suppressWarnings(print(plt))
}

generate_noise_y <- function(y, level){
  idx <- sample(length(y), size = length(y)*level, replace = F)
  y_fake <- sample(c(1,2,3), size = length(idx), replace = T)
  y[idx] <- y_fake
  return(y)
}
c(0.1 + 0.9*.4, 0.9*.35, .9*.25)
c(0.9*.35, 0.1 + 0.9*.35, .9*.3)

sum(c( 0.9*.23, .9*.37, 0.1 + 0.9*.4))

generate_asy_noise_y <- function(y, level){
  idx <- sample(length(y), size = length(y)*level, replace = F)
  y_fake <- sapply(idx, function(x){
    if(y[x]==1){
      sample(c(1,2,3), size =1, prob = c(0.4,.35,.25))
    }else if(y[x]==2){
      sample(c(1,2,3), size =1, prob = c(0.9*.35, 0.1 + 0.9*.35, .9*.3))
    }else if(y[x]==3){
      sample(c(1,2,3), size =1, prob = c( 0.9*.23, .9*.37, 0.1 + 0.9*.4))
    }
  })

  y[idx] <- y_fake
  return(y)
}


###########Nosieless case#########


y <- sample(c(1,2,3), size = 100, replace = T)
X <- draw_X(y)
C_list <- 10^seq(4,0)

models <- lapply(C_list, function(C){WW_pri_opt(X,y,C,intercept = F)})

plot_circle_boundary(lapply(C_list, function(C){WW_pri_opt(X,y,C,intercept = F)}), X, y, title = "WW")
plot_circle_boundary(lapply(C_list, function(C){CS_pri_opt(X,y,C,intercept = F)}), X, y, title = "CS")
plot_circle_boundary(lapply(C_list, function(C){Duchi_pri_opt(X,y,C,intercept = F)}), X, y, title = "Duchi")
plot_circle_boundary(lapply(C_list, function(C){MDuchi_pri_opt(X,y,C,intercept = F)}), X, y, title = "MDuchi")
plot_circle_boundary(lapply(C_list, function(C){New1_pri_opt(X,y,C,intercept = F)}), X, y, title = "New1", rule = "dagger")
plot_circle_boundary(lapply(C_list, function(C){New1_pri_opt(X,y,C,intercept = F)}), X, y, title = "New1", rule = "0_max")
plot_decision_boundary(X, y, MSVM_Weak_Hard_opt(X,y,"New1")[1:2,], rule = "0_max")
# New1_pri_opt(X,y,C = 10,intercept = T)
# MSVM_Weak_Hard_opt(X,y,"New1")[1:2,]

plot_decision_boundary(X, y, New3_pri_opt(X,y,C = 1000,intercept = F), rule = "0_max")
# New1_pri_opt(X,y,C = 1000,intercept = T)
# New3_pri_opt(X,y,C = 1000,intercept = T)

plot_decision_boundary(X, y, Duchi_pri_opt(X,y,C=1000,intercept = T))
plot_decision_boundary(X, y, MSVM_Weak_Hard_opt(X,y,"Duchi"))
plot_decision_boundary(X, y, WW_pri_opt(X,y,C=1000,intercept = T))

plot_decision_boundary(X, y, New1_pri_opt(X,y,C = 10,intercept = T), rule = "0_max")

plot_circle_boundary(lapply(C_list, function(C){New1_pri_opt(X,y,C,intercept = F)}), X, y, title = "New1", rule = "dagger")
plot_circle_boundary(lapply(C_list, function(C){New1_pri_opt(X,y,C,intercept = F)}), X, y, title = "New1", rule = "0_max")
# plot_circle_boundary(lapply(C_list, function(C){New3_pri_opt(X,y,C,intercept = F)}), X, y, title = "New3", rule = "dagger")
plot_circle_boundary(lapply(C_list, function(C){OVA_pri_opt(X,y,C,intercept = F)}), X, y, title = "OVA")
plot_circle_boundary(lapply(C_list, function(C){LLW_pri_opt(X,y,C,intercept = F)}), X, y, title = "LLW")
# plot_circle_boundary(lapply(C_list, function(C){MSVM7_pri_opt(X,y,C,intercept = F)}), X, y, title = "MSVM7", rule = "dagger")
plot_circle_boundary(lapply(C_list, function(C){MSVM8_pri_opt(X,y,C,intercept = F)}), X, y, title = "MSVM8")

##################Noise Circle Case##########################
set.seed(123)
n = 1000
y <- sample(c(1,2,3), size = n, replace = T)
X <- draw_X(y)


y_noise <- generate_noise_y(y, level=0.90)

C_list <- 10^seq(0,-4, length.out = 5)
lambda_list <- 0
plot_ridgeless_circle_boundary(lapply(lambda_list, function(lambda){New1_pri_opt(X,y_noise,C=1,lambda=lambda,intercept = T, base_class = 1)}), X, y_noise, title = "New1,base1,dagger", rule = "dagger")
New1_pri_opt(X,y_noise,C=1,lambda=0,intercept = F, base_class = 1)$beta
WW_pri_opt(X,y_noise,C=1,lambda=0,intercept = T)$beta
plot_ridgeless_circle_boundary(lapply(lambda_list, function(lambda){WW_pri_opt(X,y_noise,C=1,lambda=lambda,intercept = F)}), X, y_noise, title = "WW")
WW_pri_opt(X,y_noise,C=1,lambda=0,intercept = T)$beta
Duchi_pri_opt(X,y_noise,C=1,lambda=0,intercept = T)$beta

plot_ridgeless_circle_boundary(lapply(lambda_list, function(lambda){CS_pri_opt(X,y_noise,C=1,lambda=lambda,intercept = F)}), X, y_noise, title = "CS")

plot_circle_boundary(lapply(C_list, function(C){WW_pri_opt(X,y_noise,C,intercept = F)}), X, y_noise, title = "WW")
plot_circle_boundary(lapply(C_list, function(C){CS_pri_opt(X,y_noise,C,intercept = F)}), X, y_noise, title = "CS")
plot_circle_boundary(lapply(C_list, function(C){Duchi_pri_opt(X,y_noise,C,intercept = F)}), X, y_noise, title = "Duchi")
plot_circle_boundary(lapply(C_list, function(C){MDuchi_pri_opt(X,y_noise,C,intercept = F)}), X, y_noise, title = "MDuchi")
plot_circle_boundary(lapply(C_list, function(C){New1_pri_opt(X,y_noise,C,intercept = F)}), X, y_noise, title = "New1,0", rule = "0_max")


# g6 <- plot_circle_boundary(lapply(C_list, function(C){New1_pri_opt(X,y_noise,C,intercept = T)}), X, y_noise, title = "New1,dagger", rule = "dagger")

# Duchi_pri_opt(X,y_noise,1,intercept = F)
plot_decision_boundary(X, y, beta = New1_pri_opt(X, y_noise, 1, intercept = T))
New1_pri_opt(X, y_noise, 1, intercept = T)

plot_decision_boundary(X, y, beta = MSVM7_pri_opt(X, y_noise, 1, intercept = T), rule = "dagger_MSVM7", xlim = c(-100,100), ylim = c(-100,100))


debugonce(plot_decision_boundary)
plot_decision_boundary(X, y, beta = New1_pri_opt(X, y_noise, 1e-3, intercept = T), rule = "dagger", xlim = c(-100,100), ylim = c(-100,100))

New1_fit <- New1_pri_opt(X, y_noise, 1, intercept = F)

plot_decision_boundary(X, y, beta = New1_primary_beta , rule = "dagger", xlim = c(-100,100), ylim = c(-100,100))

g7 <- plot_circle_boundary(lapply(C_list, function(C){New3_pri_opt(X,y_noise,C,intercept = F)}), X, y_noise, title = "New3,0", rule = "0_max")
g8 <- plot_circle_boundary(lapply(C_list, function(C){New3_pri_opt(X,y_noise,C,intercept = F)}), X, y_noise, title = "New3,dagger", rule = "dagger")
g9 <- plot_circle_boundary(lapply(C_list, function(C){OVA_pri_opt(X,y_noise,C,intercept = F)}), X, y_noise, title = "OVA")
g10 <- plot_circle_boundary(lapply(C_list, function(C){LLW_pri_opt(X,y_noise,C,intercept = F)}), X, y_noise, title = "LLW")
g11 <- plot_circle_boundary(lapply(C_list, function(C){MSVM7_pri_opt(X,y_noise,C,intercept = F)}), X, y_noise, title = "MSVM7,dagger", rule = "dagger_MSVM7")
g12 <- plot_circle_boundary(lapply(C_list, function(C){MSVM8_pri_opt(X,y_noise,C,intercept = F)}), X, y_noise, title = "MSVM8")
g13 <- plot_circle_boundary(lapply(C_list, function(C){MSVM7_pri_opt(X,y_noise,C,intercept = F)}), X, y_noise, title = "MSVM7,0", rule = "0_max")

ggarrange(g1, g2, g3,g4, g5,g6, g7, g8, g9,g10,g11,g12, g13, ncol=4, nrow=4, common.legend = TRUE, legend="bottom")

# plot_decision_boundary(X,y,New1_pri_opt(X,y_noise,1e-1,intercept = T), rule = "0_max")




##################Noise Circle Case##########################
set.seed(123)
n = 1000
y <- sample(c(1,2,3), size = n, replace = T)
X <- draw_X(y)


y_noise <- generate_asy_noise_y(y, level=0.90)

C_list <- 10^seq(0,-4, length.out = 5)
lambda_list <- 0
table(y_noise)

C_list <- 10^seq(0,-4)


g1 <- plot_circle_boundary(lapply(C_list, function(C){WW_pri_opt(X,y_noise,C,intercept = F)}), X, y_noise, title = "WW")
g2 <- plot_circle_boundary(lapply(C_list, function(C){CS_pri_opt(X,y_noise,C,intercept = F)}), X, y_noise, title = "CS")
g3 <- plot_circle_boundary(lapply(C_list, function(C){Duchi_pri_opt(X,y_noise,C,intercept = F)}), X, y_noise, title = "Duchi")
g4 <- plot_circle_boundary(lapply(C_list, function(C){MDuchi_pri_opt(X,y_noise,C,intercept = F)}), X, y_noise, title = "MDuchi")

new1_models_b1 <- lapply(C_list, function(C){New1_pri_opt(X,y_noise,C,intercept = T, base_class = 1)})
new1_models_b2 <- lapply(C_list, function(C){New1_pri_opt(X,y_noise,C,intercept = T, base_class = 2)})
new1_models_b3 <- lapply(C_list, function(C){New1_pri_opt(X,y_noise,C,intercept = T, base_class = 3)})

g5.1 <- plot_circle_boundary(new1_models_b1, X, y_noise, title = "New1,base1,0", rule = "0_max")
g5.2 <- plot_circle_boundary(new1_models_b2, X, y_noise, title = "New1,base2,0", rule = "0_max")
g5.3 <- plot_circle_boundary(new1_models_b3, X, y_noise, title = "New1,base3,0", rule = "0_max")
g5.4 <- NULL
g6.1 <- plot_circle_boundary(new1_models_b1, X, y_noise, title = "New1,base1,dagger", rule = "dagger")
g6.2 <- plot_circle_boundary(new1_models_b2, X, y_noise, title = "New1,base2,dagger", rule = "dagger")
g6.3 <- plot_circle_boundary(new1_models_b3, X, y_noise, title = "New1,base3,dagger", rule = "dagger")
g6.4 <- NULL
new3_models_b1 <- lapply(C_list, function(C){New3_pri_opt(X,y_noise,C,intercept = T, base_class = 1)})
new3_models_b2 <- lapply(C_list, function(C){New3_pri_opt(X,y_noise,C,intercept = T, base_class = 2)})
new3_models_b3 <- lapply(C_list, function(C){New3_pri_opt(X,y_noise,C,intercept = T, base_class = 3)})

g7.1 <- plot_circle_boundary(new3_models_b1, X, y_noise, title = "New3,base1,0", rule = "0_max")
g7.2 <- plot_circle_boundary(new3_models_b2, X, y_noise, title = "New3,base2,0", rule = "0_max")
g7.3 <- plot_circle_boundary(new3_models_b3, X, y_noise, title = "New3,base3,0", rule = "0_max")
g7.4 <- NULL
g8.1 <- plot_circle_boundary(new3_models_b1, X, y_noise, title = "New3,base1,dagger", rule = "dagger")
g8.2 <- plot_circle_boundary(new3_models_b2, X, y_noise, title = "New3,base2,dagger", rule = "dagger")
g8.3 <- plot_circle_boundary(new3_models_b3, X, y_noise, title = "New3,base3,dagger", rule = "dagger")
g8.4 <- NULL

g9 <- plot_circle_boundary(lapply(C_list, function(C){OVA_pri_opt(X,y_noise,C,intercept = F)}), X, y_noise, title = "OVA")
g10 <- tryCatch({plot_circle_boundary(lapply(C_list, function(C){LLW_pri_opt(X,y_noise,C,intercept = F)}), X, y_noise, title = "LLW")}, error = function(e){NULL})
g11 <- plot_circle_boundary(lapply(C_list, function(C){MSVM8_pri_opt(X,y_noise,C,intercept = F)}), X, y_noise, title = "MSVM8")
g12 <- NULL
msvm_models_b1 <- lapply(C_list, function(C){MSVM7_pri_opt(X,y_noise,C,intercept = T, base_class = 1)})
msvm_models_b2 <- lapply(C_list, function(C){MSVM7_pri_opt(X,y_noise,C,intercept = T, base_class = 2)})
msvm_models_b3 <- lapply(C_list, function(C){MSVM7_pri_opt(X,y_noise,C,intercept = T, base_class = 3)})
g14.1 <- plot_circle_boundary(msvm_models_b1, X, y_noise, title = "MSVM7,base1,dagger", rule = "dagger")
g14.2 <- plot_circle_boundary(msvm_models_b2, X, y_noise, title = "MSVM7,base2,dagger", rule = "dagger")
g14.3 <- plot_circle_boundary(msvm_models_b3, X, y_noise, title = "MSVM7,base3,dagger", rule = "dagger")
g14.4 <- NULL

g13.1 <- plot_circle_boundary(msvm_models_b1, X, y_noise, title = "MSVM7,base1,0", rule = "0_max")
g13.2 <- plot_circle_boundary(msvm_models_b2, X, y_noise, title = "MSVM7,base2,0", rule = "0_max")
g13.3 <- plot_circle_boundary(msvm_models_b3, X, y_noise, title = "MSVM7,base3,0", rule = "0_max")
g13.4 <- NULL
ggarrange(g1, g2, g3,g4, g5.1,g5.2,g5.3,g5.4,g6.1,g6.2,g6.3,g6.4, g7.1, g7.2,g7.3,g7.4, g8.1,g8.2,g8.3,g8.4, g9,g10,g11,g12, g13.1,g13.2,g13.3,g13.4,g14.1,g14.2,g14.3,g14.4, ncol=4, nrow=8, common.legend = TRUE, legend="bottom")
