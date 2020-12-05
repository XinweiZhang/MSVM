setwd("~/Desktop/Multiclass Classification/MSVM Code")
library(plotly)

rm(list=ls())

Top2 <- function(gamma){
  gamma_sort <- sort(gamma, decreasing = T)
  L1 <- ifelse(gamma_sort[2]> gamma[1], 1, 0)
  L2 <- ifelse(gamma_sort[2]> gamma[2], 1, 0)
  L3 <- ifelse(gamma_sort[2]> gamma[3], 1, 0)
  return(c(L1,L2,L3))
}



gamma =matrix(runif(120000,-2,2),ncol=3)
L_Top2 <- apply(gamma,MARGIN=1,FUN=Top2)
L_Top2_data <- as.data.frame(t(L_Top2))

L_Top2_data$color <- rowSums(L_Top2_data)

fig <- plot_ly(L_Top2_data, x= ~V1, y= ~V2, z = ~V3, marker = list(color =~color, colorscale = c('#FFE1A1', '#683531'), showscale = T))
fig <- fig %>% layout(
  title="Top2", scene = list( xaxis = list(title = 'L1', range=c(0,3)),
                                      yaxis = list(title = 'L2', range=c(0,3)),
                                      zaxis = list(title = 'L3', range = c(0,3))))
print(fig)





gamma =matrix(runif(120000,-2,2),ncol=3)
L_Top2 <- apply(gamma,MARGIN=1,FUN=Top2)
L_Top2_data <- as.data.frame(t(L_Top2))

L_Top2_data$color <- rowSums(L_Top2_data)

fig <- plot_ly(L_Top2_data, x= ~V1, y= ~V2, z = ~V3, marker = list(color =~color, colorscale = c('#FFE1A1', '#683531'), showscale = T))
fig <- fig %>% layout(
  title="Top2", scene = list( xaxis = list(title = 'L1', range=c(0,3)),
                              yaxis = list(title = 'L2', range=c(0,3)),
                              zaxis = list(title = 'L3', range = c(0,3))))
print(fig)


Top2_Duchi <- function(gamma){
  gamma_sort <- sort(gamma, decreasing = T)
  L1 <-  1-gamma[1] + max( gamma_sort[1]-1, (gamma_sort[1] + gamma_sort[2] +  gamma_sort[3] - 2)/3)
  L2 <-  1-gamma[2] + max( gamma_sort[1]-1, (gamma_sort[1] + gamma_sort[2] +  gamma_sort[3] - 2)/3)
  L3 <-  1-gamma[3] + max( gamma_sort[1]-1, (gamma_sort[1] + gamma_sort[2] +  gamma_sort[3] - 2)/3)
  return(c(L1,L2,L3))
}

# Top2_Duchi(c(1,1,0))

gamma =matrix(runif(60000,-1,1),ncol=3)
L_Top2_Duchi <- apply(gamma,MARGIN=1,FUN=Top2_Duchi)
L_Top2_Duchi_data <- as.data.frame(t(L_Top2_Duchi))

L_Top2_Duchi_data$color <- rowSums(L_Top2_Duchi_data)

fig <- plot_ly(L_Top2_Duchi_data, x= ~V1, y= ~V2, z = ~V3, marker = list(color =~color, colorscale = c('#FFE1A1', '#683531'), showscale = T))
fig <- fig %>% layout(
  title="Top2_Duchi", scene = list( xaxis = list(title = 'L1', range=c(0,2)),
                              yaxis = list(title = 'L2', range=c(0,2)),
                              zaxis = list(title = 'L3', range = c(0,2))))
print(fig)



Top2_MDuchi <- function(gamma){
  s_gamma1 <- sort(gamma[c(2,3)], decreasing = T)
  s_gamma2 <- sort(gamma[c(1,3)], decreasing = T)
  s_gamma3 <- sort(gamma[c(1,2)], decreasing = T)
  
  L1 <-  max(0, (s_gamma1[1] + s_gamma1[2] - 2*gamma[1] + 1)/3)
  L2 <-  max(0, (s_gamma2[1] + s_gamma2[2] - 2*gamma[2] + 1)/3)
  L3 <-  max(0, (s_gamma3[1] + s_gamma3[2] - 2*gamma[3] + 1)/3)
  return(c(L1,L2,L3))
}

# Top2_Duchi(c(1,1,0))

gamma =matrix(runif(60000,-2,2),ncol=3)
L_Top2_MDuchi <- apply(gamma,MARGIN=1,FUN=Top2_MDuchi)
L_Top2_MDuchi_data <- as.data.frame(t(L_Top2_MDuchi))

L_Top2_MDuchi_data$color <- rowSums(L_Top2_MDuchi_data)

fig <- plot_ly(L_Top2_MDuchi_data, x= ~V1, y= ~V2, z = ~V3, marker = list(color =~color, colorscale = c('#FFE1A1', '#683531'), showscale = T))
fig <- fig %>% layout(
  title="Top2_MDuchi", scene = list( xaxis = list(title = 'L1', range=c(0,2)),
                                    yaxis = list(title = 'L2', range=c(0,2)),
                                    zaxis = list(title = 'L3', range = c(0,2))))
print(fig)




#########Below the psi1,2,3,4,5 loss are adpoted from Yang and Koyejo 2020
Top2_psi1 <- function(gamma){
  s_gamma1 <- sort(gamma[c(2,3)], decreasing = T)
  s_gamma2 <- sort(gamma[c(1,3)], decreasing = T)
  s_gamma3 <- sort(gamma[c(1,2)], decreasing = T)
  
  L1 <-  max(0, 1 + s_gamma1[1] - gamma[1])
  L2 <-  max(0, 1 + s_gamma2[1] - gamma[2])
  L3 <-  max(0, 1 + s_gamma3[1] - gamma[3])
  return(c(L1,L2,L3))
}

# Top2_Duchi(c(1,1,0))

gamma =matrix(runif(60000,-2,2),ncol=3)
L_Top2_psi1 <- apply(gamma,MARGIN=1,FUN=Top2_psi1)
L_Top2_psi1_data <- as.data.frame(t(L_Top2_psi1))

L_Top2_psi1_data$color <- rowSums(L_Top2_psi1_data)

fig <- plot_ly(L_Top2_psi1_data, x= ~V1, y= ~V2, z = ~V3, marker = list(color =~color, colorscale = c('#FFE1A1', '#683531'), showscale = T))
fig <- fig %>% layout(
  title="Top2_psi1", scene = list( xaxis = list(title = 'L1', range=c(0,2)),
                                     yaxis = list(title = 'L2', range=c(0,2)),
                                     zaxis = list(title = 'L3', range = c(0,2))))
print(fig)


Top2_psi2 <- function(gamma){
  s_gamma1 <- sort(gamma + c(0,1,1), decreasing = T)
  s_gamma2 <- sort(gamma + c(1,0,1), decreasing = T)
  s_gamma3 <- sort(gamma + c(1,1,0), decreasing = T)
  
  L1 <-  max(0, (s_gamma1[1] + s_gamma1[2])/2 - gamma[1])
  L2 <-  max(0, (s_gamma2[1] + s_gamma2[2])/2 - gamma[2])
  L3 <-  max(0, (s_gamma3[1] + s_gamma3[2])/2 - gamma[3])
  return(c(L1,L2,L3))
}

# Top2_Duchi(c(1,1,0))

gamma =matrix(runif(90000,-3,3),ncol=3)
L_Top2_psi2 <- apply(gamma,MARGIN=1,FUN=Top2_psi2)
L_Top2_psi2_data <- as.data.frame(t(L_Top2_psi2))

L_Top2_psi2_data$color <- rowSums(L_Top2_psi2_data)

fig <- plot_ly(L_Top2_psi2_data, x= ~V1, y= ~V2, z = ~V3, marker = list(color =~color, colorscale = c('#FFE1A1', '#683531'), showscale = T))
fig <- fig %>% layout(
  title="Top2_psi2", scene = list( xaxis = list(title = 'L1', range=c(0,4)),
                                   yaxis = list(title = 'L2', range=c(0,4)),
                                   zaxis = list(title = 'L3', range = c(0,4))))
print(fig)


Top2_psi3 <- function(gamma){
  s_gamma1 <- sort(gamma + c(0,1,1), decreasing = T)
  s_gamma2 <- sort(gamma + c(1,0,1), decreasing = T)
  s_gamma3 <- sort(gamma + c(1,1,0), decreasing = T)
  
  L1 <-  1/2*(max(0, s_gamma1[1] - gamma[1]) +  max(0, s_gamma1[2] - gamma[1]))
  L2 <-  1/2*(max(0, s_gamma2[1] - gamma[2]) +  max(0, s_gamma2[2] - gamma[2]))
  L3 <-  1/2*(max(0, s_gamma3[1] - gamma[3]) +  max(0, s_gamma3[2] - gamma[3]))
  return(c(L1,L2,L3))
}


gamma =matrix(runif(90000,-3,3),ncol=3)
L_Top2_psi3 <- apply(gamma,MARGIN=1,FUN=Top2_psi3)
L_Top2_psi3_data <- as.data.frame(t(L_Top2_psi3))

L_Top2_psi3_data$color <- rowSums(L_Top2_psi3_data)

fig <- plot_ly(L_Top2_psi3_data, x= ~V1, y= ~V2, z = ~V3, marker = list(color =~color, colorscale = c('#FFE1A1', '#683531'), showscale = T))
fig <- fig %>% layout(
  title="Top2_psi3", scene = list( xaxis = list(title = 'L1', range=c(0,4)),
                                   yaxis = list(title = 'L2', range=c(0,4)),
                                   zaxis = list(title = 'L3', range = c(0,4))))
print(fig)

Top2_psi4 <- function(gamma){
  s_gamma1 <- sort(gamma[c(2,3)] + c(1,1), decreasing = T)
  s_gamma2 <- sort(gamma[c(1,3)] + c(1,1), decreasing = T)
  s_gamma3 <- sort(gamma[c(1,2)] + c(1,1), decreasing = T)
  
  L1 <-  max(0, 1/2*sum(s_gamma1) - gamma[1]) 
  L2 <-   max(0, 1/2*sum(s_gamma2) - gamma[2])
  L3 <-   max(0, 1/2*sum(s_gamma3) - gamma[3])
  return(c(L1,L2,L3))
}


gamma =matrix(runif(90000,-3,3),ncol=3)
L_Top2_psi4 <- apply(gamma,MARGIN=1,FUN=Top2_psi4)
L_Top2_psi4_data <- as.data.frame(t(L_Top2_psi4))

L_Top2_psi4_data$color <- rowSums(L_Top2_psi4_data)

fig <- plot_ly(L_Top2_psi4_data, x= ~V1, y= ~V2, z = ~V3, marker = list(color =~color, colorscale = c('#FFE1A1', '#683531'), showscale = T))
fig <- fig %>% layout(
  title="Top2_psi4", scene = list( xaxis = list(title = 'L1', range=c(0,4)),
                                   yaxis = list(title = 'L2', range=c(0,4)),
                                   zaxis = list(title = 'L3', range = c(0,4))))
print(fig)



Top2_psi5 <- function(gamma){
  gamma_sort <- sort(gamma, decreasing = T)
  
  L1 <-  max(0, (1 + gamma_sort[3] - gamma[1]))
  L2 <-  max(0, (1 + gamma_sort[3] - gamma[2]))
  L3 <-  max(0, (1 + gamma_sort[3] - gamma[3]))
  return(c(L1,L2,L3))
}

# Top2_Duchi(c(1,1,0))

gamma =matrix(runif(120000,-5,5),ncol=3)
L_Top2_psi5 <- apply(gamma,MARGIN=1,FUN=Top2_psi5)
L_Top2_psi5_data <- as.data.frame(t(L_Top2_psi5))

L_Top2_psi5_data$color <- rowSums(L_Top2_psi5_data)

fig <- plot_ly(L_Top2_psi5_data, x= ~V1, y= ~V2, z = ~V3, marker = list(color =~color, colorscale = c('#FFE1A1', '#683531'), showscale = T))
fig <- fig %>% layout(
  title="Top2_psi5", scene = list( xaxis = list(title = 'L1', range=c(0,4)),
                                   yaxis = list(title = 'L2', range=c(0,4)),
                                   zaxis = list(title = 'L3', range = c(0,4))))
print(fig)

Top2_psi5 <- function(gamma){
  gamma_sort <- sort(gamma, decreasing = T)
  
  L1 <-  max(0, (1 + gamma_sort[2] - gamma[1]))
  L2 <-  max(0, (1 + gamma_sort[2] - gamma[2]))
  L3 <-  max(0, (1 + gamma_sort[2] - gamma[3]))
  return(c(L1,L2,L3))
}

# Top2_Duchi(c(1,1,0))

gamma =matrix(runif(120000,-5,5),ncol=3)
L_Top2_psi5 <- apply(gamma,MARGIN=1,FUN=Top2_psi5)
L_Top2_psi5_data <- as.data.frame(t(L_Top2_psi5))

L_Top2_psi5_data$color <- rowSums(L_Top2_psi5_data)

fig <- plot_ly(L_Top2_psi5_data, x= ~V1, y= ~V2, z = ~V3, marker = list(color =~color, colorscale = c('#FFE1A1', '#683531'), showscale = T))
fig <- fig %>% layout(
  title="Top2_psi5", scene = list( xaxis = list(title = 'L1', range=c(0,5)),
                                   yaxis = list(title = 'L2', range=c(0,5)),
                                   zaxis = list(title = 'L4', range = c(0,5))))
print(fig)
