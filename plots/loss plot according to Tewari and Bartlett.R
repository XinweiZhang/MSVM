library(plotly)

rm(list=ls())

MSVM7 <- function(gamma){
  return(c(max(max(1-gamma[1] -  max(gamma[2],0),0)  +  max(gamma[2],0)),
           max(max(1-gamma[2] -  max(gamma[1],0),0)  +  max(gamma[1],0)),
           max(gamma[1],0) + max(gamma[2],0)))
}


MSVM7(c(0,0))
MSVM7(c(1,0))
MSVM7(c(0,1))
MSVM7(c(1,1))

gamma =matrix(runif(100000,-2,2),ncol=2)

L_MSVM7 <- apply(gamma,MARGIN=1,FUN=MSVM7)
L_MSVM7_data <- as.data.frame(t(L_MSVM7))

L_MSVM7_data$color <- rowSums(L_MSVM7_data)

fig <- plot_ly(L_MSVM7_data, x= ~V1, y= ~V2, z = ~V3, marker = list(color =~color, colorscale = c('#FFE1A1', '#683531'), showscale = T))
fig <- fig %>% layout(
  title="MSVM7", scene = list( xaxis = list(title = 'L1', range=c(0,2.5)),
                               yaxis = list(title = 'L2', range=c(0,2.5)),
                               zaxis = list(title = 'L3', range = c(0,2.5))))
fig


MSVM8 <- function(gamma){
  return(c(max(1-gamma[1],  max(gamma[2],0) + max(gamma[3],0)),
  max(1-gamma[2],  max(gamma[1],0) + max(gamma[3],0)),
  max(1-gamma[3],  max(gamma[1],0) + max(gamma[2],0))))
}


MSVM8(c(0,0,0))
MSVM8(c(1,0,0))
MSVM8(c(0,1,0))
MSVM8(c(0,0,1))
MSVM8(c(-1,0,1))
MSVM8(c(-2,0,1))

gamma =matrix(runif(90000,-2,2),ncol=3)

L_MSVM8 <- apply(gamma,MARGIN=1,FUN=MSVM8)
L_MSVM8_data <- as.data.frame(t(L_MSVM8))

L_MSVM8_data$color <- rowSums(L_MSVM8_data)

fig <- plot_ly(L_MSVM8_data, x= ~V1, y= ~V2, z = ~V3, marker = list(color =~color, colorscale = c('#FFE1A1', '#683531'), showscale = T))
fig <- fig %>% layout(
  title="MSVM8", scene = list( xaxis = list(title = 'L1', range=c(0,2.5)),
                                   yaxis = list(title = 'L2', range=c(0,2.5)),
                                   zaxis = list(title = 'L3', range = c(0,2.5))))
fig

RF <- function(gamma){
  a=3/4
  return(c(a*max(1-gamma[1],0)+ (1-a)*(max(gamma[2],0) + max(gamma[3],0)),
           a*max(1-gamma[2],0)+ (1-a)*(max(gamma[1],0) + max(gamma[3],0)),
           a*max(1-gamma[3],0)+ (1-a)*(max(gamma[1],0) + max(gamma[2],0))))
}


gamma =matrix(runif(120000,-3,3),ncol=3)
gamma <- t(apply(gamma,MARGIN=1, FUN = function(x){ y <- x + rep((1- sum(x))/3,3)}))
L_RF <- t(apply(gamma,MARGIN=1,FUN=RF))
L_RF_data <- as.data.frame(L_RF)
L_RF_data$color <- rowSums(L_RF_data)

fig <- plot_ly(L_RF_data, x= ~V1, y= ~V2, z = ~V3, marker = list(color =~color, colorscale = c('#FFE1A1', '#683531'), showscale = T))
fig <- fig %>% layout(
  title="RF(a=3/4)",scene = list(xaxis = list(title = 'L1', range=c(0,3)),
                           yaxis = list(title = 'L2', range=c(0,3)),
                           zaxis = list(title = 'L3', range = c(0,3))))
fig

RF <- function(gamma){
  a=1/2
  return(c(a*max(1-gamma[1],0)+ (1-a)*(max(gamma[2],0) + max(gamma[3],0)),
           a*max(1-gamma[2],0)+ (1-a)*(max(gamma[1],0) + max(gamma[3],0)),
           a*max(1-gamma[3],0)+ (1-a)*(max(gamma[1],0) + max(gamma[2],0))))
}


gamma =matrix(runif(120000,-3,3),ncol=3)
gamma <- t(apply(gamma,MARGIN=1, FUN = function(x){ y <- x + rep((1- sum(x))/3,3)}))
L_RF <- t(apply(gamma,MARGIN=1,FUN=RF))
L_RF_data <- as.data.frame(L_RF)
L_RF_data$color <- rowSums(L_RF_data)

fig <- plot_ly(L_RF_data, x= ~V1, y= ~V2, z = ~V3, marker = list(color =~color, colorscale = c('#FFE1A1', '#683531'), showscale = T))
fig <- fig %>% layout(
  title="RF(a=1/2)",scene = list(xaxis = list(title = 'L1', range=c(0,3)),
                                 yaxis = list(title = 'L2', range=c(0,3)),
                                 zaxis = list(title = 'L3', range = c(0,3))))
fig

RF <- function(gamma){
  a=1/4
  return(c(a*max(1-gamma[1],0)+ (1-a)*(max(gamma[2],0) + max(gamma[3],0)),
           a*max(1-gamma[2],0)+ (1-a)*(max(gamma[1],0) + max(gamma[3],0)),
           a*max(1-gamma[3],0)+ (1-a)*(max(gamma[1],0) + max(gamma[2],0))))
}


gamma =matrix(runif(120000,-3,3),ncol=3)
gamma <- t(apply(gamma,MARGIN=1, FUN = function(x){ y <- x + rep((1- sum(x))/3,3)}))
L_RF <- t(apply(gamma,MARGIN=1,FUN=RF))
L_RF_data <- as.data.frame(L_RF)
L_RF_data$color <- rowSums(L_RF_data)

fig <- plot_ly(L_RF_data, x= ~V1, y= ~V2, z = ~V3, marker = list(color =~color, colorscale = c('#FFE1A1', '#683531'), showscale = T))
fig <- fig %>% layout(
  title="RF(a=1/4)",scene = list(xaxis = list(title = 'L1', range=c(0,3)),
                                 yaxis = list(title = 'L2', range=c(0,3)),
                                 zaxis = list(title = 'L3', range = c(0,3))))
fig


LLW <- function(gamma){
  return(c(max(gamma[2],0) + max(gamma[3],0),
           max(gamma[1],0) + max(gamma[3],0),
           max(gamma[1],0) + max(gamma[2],0)))
}


gamma =matrix(runif(90000,-2,2),ncol=3)
gamma <- t(apply(gamma,MARGIN=1, FUN = function(x){ y <- x + rep((1- sum(x))/3,3)}))
L_LLW <- t(apply(gamma,MARGIN=1,FUN=LLW))
L_LLW_data <- as.data.frame(L_LLW)
L_LLW_data$color <- rowSums(L_LLW_data)

fig <- plot_ly(L_LLW_data, x= ~V1, y= ~V2, z = ~V3, marker = list(color =~color, colorscale = c('#FFE1A1', '#683531'), showscale = T))
fig <- fig %>% layout(
  title="LLW",scene = list(xaxis = list(title = 'L1', range=c(0,2.5)),
                                   yaxis = list(title = 'L2', range=c(0,2.5)),
                                   zaxis = list(title = 'L3', range = c(0,2.5))))
fig


Duchi <- function(gamma){
  s_gamma <- sort(gamma, decreasing = T)
  g <- max(s_gamma[1]-1, (s_gamma[1]+s_gamma[2]-1)/2,  (s_gamma[1]+s_gamma[2] + s_gamma[3]-1)/3)
  return(c(1-gamma[1]+g,
           1-gamma[2]+g,
           1-gamma[3]+g))
}

Duchi(c(0,0,0))
Duchi(c(1,0,0))
Duchi(c(0,1,0))
Duchi(c(0,0,1))

gamma =matrix(runif(60000,-2,2),ncol=3)
gamma <- t(apply(gamma,MARGIN=1, FUN = function(x){ y <- x + rep((1- sum(x))/3,3)}))
L_Duchi <- t(apply(gamma,MARGIN=1,FUN=Duchi))

L_Duchi_data <- as.data.frame(L_Duchi)
L_Duchi_data$color <- rowSums(L_Duchi_data)
# L_Duchi_data$color <- rowSums(gamma)

fig <- plot_ly(L_Duchi_data, x= ~V1, y= ~V2, z = ~V3, marker = list(color =~color, colorscale = c('#FFE1A1', '#683531'), showscale = T))

fig <- fig %>% layout(
  title="Duchi",
  scene = list(xaxis = list(title = 'L1', range=c(0,4)),
                                   yaxis = list(title = 'L2', range=c(0,4)),
                                   zaxis = list(title = 'L3', range = c(0,4))))
fig


MDuchi <- function(gamma){
  s_gamma1 <- sort(gamma[c(2,3)], decreasing = T)
  s_gamma2 <- sort(gamma[c(1,3)], decreasing = T)
  s_gamma3 <- sort(gamma[c(1,2)], decreasing = T)
   
  L1 <- max(0,(s_gamma1[1]-gamma[1]+1)/2, ((s_gamma1[1]-gamma[1]+1)/3+ (s_gamma1[2]-gamma[1]+1)/3))
  L2 <- max(0,(s_gamma2[1]-gamma[2]+1)/2, ((s_gamma2[1]-gamma[2]+1)/3+ (s_gamma2[2]-gamma[2]+1)/3))
  L3 <- max(0,(s_gamma3[1]-gamma[3]+1)/2, ((s_gamma3[1]-gamma[3]+1)/3+ (s_gamma3[2]-gamma[3]+1)/3))
  
  return(c(L1,L2,L3))
}

MDuchi(c(1,0,0))
MDuchi(c(0,1,0))
MDuchi(c(0,0,1))

gamma =matrix(runif(60000,-2,2),ncol=3)

gamma <- t(apply(gamma,MARGIN=1, FUN = function(x){ y <- x + rep((1- sum(x))/3,3)}))
L_MDuchi <- t(apply(gamma,MARGIN=1,FUN=MDuchi))
L_MDuchi_data <- as.data.frame(L_MDuchi)
L_MDuchi_data$color <- rowSums(L_MDuchi_data)
L_MDuchi_data$color <- rowSums(gamma)
fig <- plot_ly(L_MDuchi_data, x= ~V1, y= ~V2, z = ~V3, marker = list(color =~color, colorscale = c('#FFE1A1', '#683531'), showscale = T))

fig <- fig %>% layout(
  title="M_duchi",
  scene = list(xaxis = list(title = 'L1', range=c(0,2.5)),
                                   yaxis = list(title = 'L2', range=c(0,2.5)),
                                   zaxis = list(title = 'L3', range = c(0,2.5))))
fig


New1 <- function(gamma){
  s_gamma <- sort(gamma, decreasing = T)
  S <- max(0, s_gamma[1]-1, (s_gamma[1]+s_gamma[2]-1)/2)
  L1 <- 1-gamma[1] + S
  L2 <-  1-gamma[2] + S
  L3 <- max(gamma[1],0)  + max(gamma[2],0)
  return(c(L1,L2,L3))
}

New1(c(1,0))
New1(c(0,1))
New1(c(0,0))

gamma =matrix(runif(80000,-4,4),ncol=2)

L_New1 <- t(apply(gamma,MARGIN=1,FUN=New1))
L_New1_data <- as.data.frame(L_New1)
L_New1_data$color <- rowSums(L_New1_data)
L_New1_data$color <- rowSums(gamma)
fig <- plot_ly(L_New1_data, x= ~V1, y= ~V2, z = ~V3, marker = list(color =~color, colorscale = c('#FFE1A1', '#683531'), showscale = T))
fig <- fig %>% layout(
  title="New1",scene = list(xaxis = list(title = 'L1', range=c(0,4)),
                                   yaxis = list(title = 'L2', range=c(0,4)),
                                   zaxis = list(title = 'L3', range = c(0,4))))
fig


New3 <- function(gamma){
  L1 <- max(0, 1- gamma[1], (gamma[2]-gamma[1]+1)/2)
  L2 <- max(0, 1- gamma[2], (gamma[1]-gamma[2]+1)/2)
  L3 <- max(gamma[1],0)  + max(gamma[2],0)
  return(c(L1,L2,L3))
}

New3(c(1,0))
New3(c(0,1))
New3(c(0,0))
New3(c(1/2,1/2))

gamma =matrix(runif(80000,-4,4),ncol=2)

L_New3 <- t(apply(gamma,MARGIN=1,FUN=New3))
L_New3_data <- as.data.frame(L_New3)
L_New3_data$color <- rowSums(L_New3_data)
fig <- plot_ly(L_New3_data, x= ~V1, y= ~V2, z = ~V3, marker = list(color =~color, colorscale = c('#FFE1A1', '#683531'), showscale = T))
fig <- fig %>% layout(title="New3",scene = list(xaxis = list(title = 'L1', range=c(0,4)),
                                   yaxis = list(title = 'L2', range=c(0,4)),
                                   zaxis = list(title = 'L3', range = c(0,4))))
fig




OVA <- function(gamma){
  return(c(sum(1-gamma[1],  max(gamma[2],0) + max(gamma[3],0)),
           sum(1-gamma[2],  max(gamma[1],0) + max(gamma[3],0)),
           sum(1-gamma[3],  max(gamma[1],0) + max(gamma[2],0))))
}




gamma =matrix(runif(150000,-2,2),ncol=3)

L_OVA <- apply(gamma,MARGIN=1,FUN=OVA)
L_OVA_data <- as.data.frame(t(L_OVA))
fig <- plot_ly(L_OVA_data, x= ~V1, y= ~V2, z = ~V3, marker = list(color =~color, colorscale = c('#FFE1A1', '#683531'), showscale = T))
fig <- fig %>% layout(title="OVA", scene = list(xaxis = list(title = 'L1', range=c(0,2.5)),
                                   yaxis = list(title = 'L2', range=c(0,2.5)),
                                   zaxis = list(title = 'L3', range = c(0,2.5))))
fig

