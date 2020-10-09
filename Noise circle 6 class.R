---
  title: "Four Class Circle Classification Series II"
author: "Xinwei"
date: "7/9/2020"
output: 
  html_document:
  toc: true
toc_depth: 3  
theme: united
number_sections: true
---
  
  ```{r setup, include=FALSE}
library(CVXR)
library(MASS)
library(R.matlab)
library(e1071)
library(kernlab)
library(ggpubr)
library(ggplot2)
rm(list=ls())
setwd("~/Desktop/Multiclass Classification/MSVM Code")
source("primary form functions.R")
source("kernel primary form functions.R")


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
  seg1 <- 0.97*data.frame(x1 = rep(0,6), x2 = cos(seq(0,15,by =3)/9*pi), y1 = rep(0,6), y2 = sin(seq(0,15,by =3)/9*pi))
  
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
                            plot.title = element_text(hjust = 0.5))   +
    scale_color_manual(values=c("1" ="blue","2"= "green", "red","yellow", "hotpink","orange"))
    # scale_color_manual(values=colors)
  
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
  seg1 <- 0.97*data.frame(x1 = rep(0,6), x2 = cos(seq(0,15,by =3)/9*pi), y1 = rep(0,6), y2 = sin(seq(0,15,by =3)/9*pi))
  
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
    scale_color_manual(values=c("1" ="blue","2"= "green", "red","yellow", "hotpink","orange"))
  # suppressWarnings(print(plt))
}



draw_X <- function(y){
  t <- sapply(y, function(z){
    if(z==1){
      runif(1, min = 0, max = 1)
    }else if(z==2){
      runif(1, min = 1, max = 3)
    }else if(z==3){
      runif(1, min = 3, max = 6)
    }else if(z==4){
      runif(1, min = 6, max = 10)
    }else if(z==5){
      runif(1, min = 10, max = 15)
    }else if(z==6){
      runif(1, min = 15, max = 20)
    }
  })
  
  X = t(sapply(t, function(x){c(cos(x*pi/10), sin(x*pi/10))}))
  return(X)
}

generate_noise_y <- function(y, level){
  idx <- sample(length(y), size = length(y)*level, replace = F)
  y_fake <- sample(seq(1,8), size = length(idx), replace = T)
  y[idx] <- y_fake
  return(y)
}
set.seed(1)
rprob <- round(rdirichlet(6, rep(1,6) ),2)
cprob <- t(sapply(seq(1,6), function(idx){
  max_idx <- which.max(rprob[idx,])
  rprob[idx,idx] <- rprob[idx,max_idx]
  rprob[idx,max_idx] <- 1 - sum(rprob[idx,-max_idx])
  return(rprob[idx,])}))

generate_asy_noise_y <- function(y, noise_level){
  
  
  idx <- sample(length(y), size = length(y)*noise_level, replace = F)
  y_fake <- sapply(idx, function(x){
    if(y[x]==1){
      sample(seq(1,6), size =1, prob = cprob[1,])
    }else if(y[x]==2){
      sample(seq(1,6), size =1, prob = cprob[2,])
    }else if(y[x]==3){
      sample(seq(1,6), size =1, prob = cprob[3,])
    }else if(y[x]==4){
      sample(seq(1,6), size =1, prob = cprob[4,])
    }else if(y[x]==5){
      sample(seq(1,6), size =1, prob = cprob[5,])
    }else if(y[x]==6){
      sample(seq(1,6), size =1, prob = cprob[6,])
    }
  })

  y[idx] <- y_fake
  return(y)
}


plot_ridge_batch <- function(X,y_noise,C_list){
  g1 <- plot_circle_boundary(lapply(C_list, function(C){WW_pri_opt(X,y_noise,C,intercept = F)}), X, y_noise, title = "WW")
  g2 <- plot_circle_boundary(lapply(C_list, function(C){CS_pri_opt(X,y_noise,C,intercept = F)}), X, y_noise, title = "CS")
  g3 <- plot_circle_boundary(lapply(C_list, function(C){Duchi_pri_opt(X,y_noise,C,intercept = F)}), X, y_noise, title = "Duchi")
  g4 <- plot_circle_boundary(lapply(C_list, function(C){MDuchi_pri_opt(X,y_noise,C,intercept = F)}), X, y_noise, title = "MDuchi")
  # new1_models_b1 <- lapply(C_list, function(C){New1_pri_opt(X,y_noise,C,intercept = F, base_class = 1)})
  # new1_models_b2 <- lapply(C_list, function(C){New1_pri_opt(X,y_noise,C,intercept = F, base_class = 2)})
  # new1_models_b3 <- lapply(C_list, function(C){New1_pri_opt(X,y_noise,C,intercept = F, base_class = 3)})
  # new1_models_b4 <- lapply(C_list, function(C){New1_pri_opt(X,y_noise,C,intercept = F, base_class = 4)})
  # g5.1 <- plot_circle_boundary(new1_models_b1, X, y_noise, title = "New1,base1,0", rule = "0_max")
  # g5.2 <- plot_circle_boundary(new1_models_b2, X, y_noise, title = "New1,base2,0", rule = "0_max")
  # g5.3 <- plot_circle_boundary(new1_models_b3, X, y_noise, title = "New1,base3,0", rule = "0_max")
  # g5.4 <- plot_circle_boundary(new1_models_b4, X, y_noise, title = "New1,base4,0", rule = "0_max")
  # g6.1 <- plot_circle_boundary(new1_models_b1, X, y_noise, title = "New1,base1,dagger", rule = "dagger")
  # g6.2 <- plot_circle_boundary(new1_models_b2, X, y_noise, title = "New1,base2,dagger", rule = "dagger")
  # g6.3 <- plot_circle_boundary(new1_models_b3, X, y_noise, title = "New1,base3,dagger", rule = "dagger")
  # g6.4 <- plot_circle_boundary(new1_models_b4, X, y_noise, title = "New1,base4,dagger", rule = "dagger")
  # new3_models_b1 <- lapply(C_list, function(C){New3_pri_opt(X,y_noise,C,intercept = F, base_class = 1)})
  # new3_models_b2 <- lapply(C_list, function(C){New3_pri_opt(X,y_noise,C,intercept = F, base_class = 2)})
  # new3_models_b3 <- lapply(C_list, function(C){New3_pri_opt(X,y_noise,C,intercept = F, base_class = 3)})
  # new3_models_b4 <- lapply(C_list, function(C){New3_pri_opt(X,y_noise,C,intercept = F, base_class = 4)})
  # g7.1 <- plot_circle_boundary(new3_models_b1, X, y_noise, title = "New3,base1,0", rule = "0_max")
  # g7.2 <- plot_circle_boundary(new3_models_b2, X, y_noise, title = "New3,base2,0", rule = "0_max")
  # g7.3 <- plot_circle_boundary(new3_models_b3, X, y_noise, title = "New3,base3,0", rule = "0_max")
  # g7.4 <- plot_circle_boundary(new3_models_b4, X, y_noise, title = "New3,base4,0", rule = "0_max")
  # g8.1 <- plot_circle_boundary(new3_models_b1, X, y_noise, title = "New3,base1,dagger", rule = "dagger")
  # g8.2 <- plot_circle_boundary(new3_models_b2, X, y_noise, title = "New3,base2,dagger", rule = "dagger")
  # g8.3 <- plot_circle_boundary(new3_models_b3, X, y_noise, title = "New3,base3,dagger", rule = "dagger")
  # g8.4 <- plot_circle_boundary(new3_models_b4, X, y_noise, title = "New3,base4,dagger", rule = "dagger")
  
  g9 <- plot_circle_boundary(lapply(C_list, function(C){OVA_pri_opt(X,y_noise,C,intercept = F)}), X, y_noise, title = "OVA")
  g10 <- tryCatch({plot_circle_boundary(lapply(C_list, function(C){LLW_pri_opt(X,y_noise,C,intercept = F)}), X, y_noise, title = "LLW")}, error = function(e){NULL})
  g11 <- plot_circle_boundary(lapply(C_list, function(C){MSVM8_pri_opt(X,y_noise,C,intercept = F)}), X, y_noise, title = "MSVM8")
  g12 <- NULL
  # msvm_models_b1 <- lapply(C_list, function(C){MSVM7_pri_opt(X,y_noise,C,intercept = F, base_class = 1)})
  # msvm_models_b2 <- lapply(C_list, function(C){MSVM7_pri_opt(X,y_noise,C,intercept = F, base_class = 2)})
  # msvm_models_b3 <- lapply(C_list, function(C){MSVM7_pri_opt(X,y_noise,C,intercept = F, base_class = 3)})
  # msvm_models_b4 <- lapply(C_list, function(C){MSVM7_pri_opt(X,y_noise,C,intercept = F, base_class = 4)})
  # g14.1 <- plot_circle_boundary(msvm_models_b1, X, y_noise, title = "MSVM7,base1,dagger", rule = "dagger")
  # g14.2 <- plot_circle_boundary(msvm_models_b2, X, y_noise, title = "MSVM7,base2,dagger", rule = "dagger")
  # g14.3 <- plot_circle_boundary(msvm_models_b3, X, y_noise, title = "MSVM7,base3,dagger", rule = "dagger")
  # g14.4 <- plot_circle_boundary(msvm_models_b4, X, y_noise, title = "MSVM7,base4,dagger", rule = "dagger")
  # 
  # g13.1 <- plot_circle_boundary(msvm_models_b1, X, y_noise, title = "MSVM7,base1,0", rule = "0_max")
  # g13.2 <- plot_circle_boundary(msvm_models_b2, X, y_noise, title = "MSVM7,base2,0", rule = "0_max")
  # g13.3 <- plot_circle_boundary(msvm_models_b3, X, y_noise, title = "MSVM7,base3,0", rule = "0_max")
  # g13.4 <- plot_circle_boundary(msvm_models_b4, X, y_noise, title = "MSVM7,base4,0", rule = "0_max")
  # ggarrange(g1, g2, g3,g4, g5.1,g5.2,g5.3,g5.4,g6.1,g6.2,g6.3,g6.4, g7.1, g7.2,g7.3,g7.4, g8.1,g8.2,g8.3,g8.4, g9,g10,g11,g12, g13.1,g13.2,g13.3,g13.4,g14.1,g14.2,g14.3,g14.4, ncol=4, nrow=8, common.legend = TRUE, legend="bottom")
  ggarrange(g1, g2, g3,g4, g9,g10,g11,g12,  ncol=4, nrow=2, common.legend = TRUE, legend="bottom")
  
}

plot_ridgeless_batch <- function(X,y_noise){
  lambda_list = 0
  
  g1 <- plot_ridgeless_circle_boundary(lapply(lambda_list, function(lambda){WW_pri_opt(X,y_noise,C=1,lambda=lambda,intercept = F)}), X, y_noise, title = "WW")
  g2 <- plot_ridgeless_circle_boundary(lapply(lambda_list, function(lambda){CS_pri_opt(X,y_noise,C=1,lambda=lambda,intercept = F)}), X, y_noise, title = "CS")
  g3 <- plot_ridgeless_circle_boundary(lapply(lambda_list, function(lambda){Duchi_pri_opt(X,y_noise,C=1,lambda=lambda,intercept = F)}), X, y_noise, title = "Duchi")
  g4 <- plot_ridgeless_circle_boundary(lapply(lambda_list, function(lambda){MDuchi_pri_opt(X,y_noise,C=1,lambda=lambda,intercept = F)}), X, y_noise, title = "MDuchi")
  # g5.1 <- plot_ridgeless_circle_boundary(lapply(lambda_list, function(lambda){New1_pri_opt(X,y_noise,C=1,lambda=lambda,intercept = F, base_class = 1)}), X, y_noise, title = "New1,base1,0", rule = "0_max")
  # g5.2 <- plot_ridgeless_circle_boundary(lapply(lambda_list, function(lambda){New1_pri_opt(X,y_noise,C=1,lambda=lambda,intercept = F, base_class = 2)}), X, y_noise, title = "New1,base2,0", rule = "0_max")
  # g5.3 <- plot_ridgeless_circle_boundary(lapply(lambda_list, function(lambda){New1_pri_opt(X,y_noise,C=1,lambda=lambda,intercept = F, base_class = 3)}), X, y_noise, title = "New1,base3,0", rule = "0_max")
  # g5.4 <- plot_ridgeless_circle_boundary(lapply(lambda_list, function(lambda){New1_pri_opt(X,y_noise,C=1,lambda=lambda,intercept = F, base_class = 4)}), X, y_noise, title = "New1,base4,0", rule = "0_max")
  # g6.1 <- plot_ridgeless_circle_boundary(lapply(lambda_list, function(lambda){New1_pri_opt(X,y_noise,C=1,lambda=lambda,intercept = F, base_class = 1)}), X, y_noise, title = "New1,base1,dagger", rule = "dagger")
  # g6.2 <- plot_ridgeless_circle_boundary(lapply(lambda_list, function(lambda){New1_pri_opt(X,y_noise,C=1,lambda=lambda,intercept = F, base_class = 2)}), X, y_noise, title = "New1,base2,dagger", rule = "dagger")
  # g6.3 <- plot_ridgeless_circle_boundary(lapply(lambda_list, function(lambda){New1_pri_opt(X,y_noise,C=1,lambda=lambda,intercept = F, base_class = 3)}), X, y_noise, title = "New1,base3,dagger", rule = "dagger")
  # g6.4<- plot_ridgeless_circle_boundary(lapply(lambda_list, function(lambda){New1_pri_opt(X,y_noise,C=1,lambda=lambda,intercept = F, base_class = 4)}), X, y_noise, title = "New1,base4,dagger", rule = "dagger")
  
  # g7.1 <- plot_ridgeless_circle_boundary(lapply(lambda_list, function(lambda){New3_pri_opt(X,y_noise,C=1,lambda=lambda,intercept = F, base_class = 1)}), X, y_noise, title = "New3,base1,0", rule = "0_max")
  # g7.2 <- plot_ridgeless_circle_boundary(lapply(lambda_list, function(lambda){New3_pri_opt(X,y_noise,C=1,lambda=lambda,intercept = F, base_class = 2)}), X, y_noise, title = "New3,base2,0", rule = "0_max")
  # g7.3 <- plot_ridgeless_circle_boundary(lapply(lambda_list, function(lambda){New3_pri_opt(X,y_noise,C=1,lambda=lambda,intercept = F, base_class = 3)}), X, y_noise, title = "New3,base3,0", rule = "0_max")
  # g7.4<- plot_ridgeless_circle_boundary(lapply(lambda_list, function(lambda){New3_pri_opt(X,y_noise,C=1,lambda=lambda,intercept = F, base_class = 4)}), X, y_noise, title = "New3,base4,0", rule = "0_max")
  # 
  # g8.1 <- plot_ridgeless_circle_boundary(lapply(lambda_list, function(lambda){New3_pri_opt(X,y_noise,C=1,lambda=lambda,intercept = F, base_class = 1)}), X, y_noise, title = "New3,base1,dagger", rule = "dagger")
  # g8.2 <- plot_ridgeless_circle_boundary(lapply(lambda_list, function(lambda){New3_pri_opt(X,y_noise,C=1,lambda=lambda,intercept = F, base_class = 2)}), X, y_noise, title = "New3,base2,dagger", rule = "dagger")
  # g8.3 <- plot_ridgeless_circle_boundary(lapply(lambda_list, function(lambda){New3_pri_opt(X,y_noise,C=1,lambda=lambda,intercept = F, base_class = 3)}), X, y_noise, title = "New3,base3,dagger", rule = "dagger")
  # g8.4 <- plot_ridgeless_circle_boundary(lapply(lambda_list, function(lambda){New3_pri_opt(X,y_noise,C=1,lambda=lambda,intercept = F, base_class = 4)}), X, y_noise, title = "New3,base4,dagger", rule = "dagger")
  
  g9 <- plot_ridgeless_circle_boundary(lapply(lambda_list, function(lambda){OVA_pri_opt(X,y_noise,C=1,lambda=lambda,intercept = F)}), X, y_noise, title = "OVA")
  g10 <-  tryCatch({plot_ridgeless_circle_boundary(lapply(lambda_list, function(lambda){LLW_pri_opt(X,y_noise,C=1,lambda=lambda,intercept = F)}), X, y_noise, title = "LLW")}, error = function(e){NULL})
  g11 <- plot_ridgeless_circle_boundary(lapply(lambda_list, function(lambda){MSVM8_pri_opt(X,y_noise,C=1,lambda=lambda,intercept = F)}), X, y_noise, title = "MSVM8")
  g12 <- NULL
  # g13.1 <- plot_ridgeless_circle_boundary(lapply(lambda_list, function(lambda){MSVM7_pri_opt(X,y_noise,C=1,lambda=lambda,intercept = F, base_class = 1)}), X, y_noise, title = "MSVM7,base1,dagger", rule = "dagger")
  # g13.2 <- plot_ridgeless_circle_boundary(lapply(lambda_list, function(lambda){MSVM7_pri_opt(X,y_noise,C=1,lambda=lambda,intercept = F, base_class = 2)}), X, y_noise, title = "MSVM7,base2,dagger", rule = "dagger")
  # g13.3 <- plot_ridgeless_circle_boundary(lapply(lambda_list, function(lambda){MSVM7_pri_opt(X,y_noise,C=1,lambda=lambda,intercept = F, base_class = 3)}), X, y_noise, title = "MSVM7,base3,dagger", rule = "dagger")
  # g13.4 <- plot_ridgeless_circle_boundary(lapply(lambda_list, function(lambda){MSVM7_pri_opt(X,y_noise,C=1,lambda=lambda,intercept = F, base_class = 4)}), X, y_noise, title = "MSVM7,base4,dagger", rule = "dagger")
  # g14.1 <- plot_ridgeless_circle_boundary(lapply(lambda_list, function(lambda){MSVM7_pri_opt(X,y_noise,C=1,lambda=lambda,intercept = F, base_class = 1)}), X, y_noise, title = "MSVM7,base1,0", rule = "0_max")
  # g14.2 <- plot_ridgeless_circle_boundary(lapply(lambda_list, function(lambda){MSVM7_pri_opt(X,y_noise,C=1,lambda=lambda,intercept = F, base_class = 2)}), X, y_noise, title = "MSVM7,base2,0", rule = "0_max")
  # g14.3 <- plot_ridgeless_circle_boundary(lapply(lambda_list, function(lambda){MSVM7_pri_opt(X,y_noise,C=1,lambda=lambda,intercept = F, base_class = 3)}), X, y_noise, title = "MSVM7,base3,0", rule = "0_max")
  # g14.4 <-  plot_ridgeless_circle_boundary(lapply(lambda_list, function(lambda){MSVM7_pri_opt(X,y_noise,C=1,lambda=lambda,intercept = F, base_class = 4)}), X, y_noise, title = "MSVM7,base4,0", rule = "0_max")
  # ggarrange(g1, g2, g3,g4, g5.1,g5.2,g5.3,g5.4,g6.1,g6.2,g6.3,g6.4, g7.1, g7.2,g7.3,g7.4, g8.1,g8.2,g8.3,g8.4, g9,g10,g11,g12, g13.1,g13.2,g13.3,g13.4,g14.1,g14.2,g14.3,g14.4, ncol=4, nrow=8, common.legend = TRUE, legend="bottom")
  # 
  ggarrange(g1, g2, g3,g4, g9,g10,g11,g12,  ncol=4, nrow=2, common.legend = TRUE, legend="bottom")
}


set.seed(123)
n <- 800
y <- sample(seq(1,6), size = n, replace = T)
X <- draw_X(y)

noise_level = 0.4

y_noise <- generate_asy_noise_y(y, noise_level)
table(y_noise)

# C_list <- 10^seq(4,0,length.out = 5)
# 
# plot_ridge_batch(X,y, C_list)

C_list <- 10^seq(0,-4,length.out = 5)
plot_ridge_batch(X,y_noise,C_list)

print(plot_circle_boundary(lapply(C_list, function(C){WW_pri_opt(X,y_noise,C,intercept = F)}), X, y_noise, title = "WW"))
