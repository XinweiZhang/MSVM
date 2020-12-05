
Duchi <- function(gamma){
  gamma <- c(gamma,1)
  s_gamma <- sort(gamma, decreasing = T)
  g <- max(s_gamma[1]-1, (s_gamma[1]+s_gamma[2]-1)/2, (sum(s_gamma)-1)/3)
  return(c(g))
}
# 
gamma =matrix(runif(60000,-5,5),ncol=2)

L_Duchi <- apply(gamma,MARGIN=1,FUN=Duchi)
L_Duchi[1:5]

L_Duchi_data <- as.data.frame(cbind(L_Duchi,gamma))

L_Duchi_data$color <- L_Duchi
colnames(L_Duchi_data)

fig <- plot_ly(L_Duchi_data, x= ~V2, y= ~V3, z = ~L_Duchi, marker = list(color =~color, colorscale = c('#FFE1A1', '#683531'), showscale = T))
fig <- fig %>% layout(
  title="Duchi",
  scene = list(xaxis = list(title = 'L1', range=c(-5,5)),
                                   yaxis = list(title = 'L2', range=c(-5,5)),
                                   zaxis = list(title = 'L3', range = c(-5,5))))
fig






MDuchi <- function(gamma){
  gamma <- c(gamma,0 - sum(gamma))
  s_gamma <- sort(gamma, decreasing = T)
  g <- max(0, (s_gamma+1)/2, (sum(s_gamma)+2)/3)
  return(g)
}

gamma =matrix(runif(60000,-5,5),ncol=2)


L_MDuchi <- apply(gamma,MARGIN=1,FUN=MDuchi)
L_MDuchi_data <- as.data.frame(cbind(L_MDuchi,gamma))

L_MDuchi_data$color <- L_MDuchi
colnames(L_MDuchi_data)

fig <- plot_ly(L_MDuchi_data, x= ~V2, y= ~V3, z = ~L_MDuchi, marker = list(color =~color, colorscale = c('#FFE1A1', '#683531'), showscale = T))
fig <- fig %>% layout(
  title="MDuchi",
  scene = list(xaxis = list(title = 'L1', range=c(-5,5)),
               yaxis = list(title = 'L2', range=c(-5,5)),
               zaxis = list(title = 'L3', range = c(-5,5))))
fig





MDuchi <- function(gamma){
  gamma <- gamma
  s_gamma <- sort(gamma, decreasing = T)
  g <- max(0, (s_gamma+1)/2)
  return(g)
}

gamma =matrix(runif(6000,-5,5),ncol=1)


L_MDuchi <- apply(gamma,MARGIN=1,FUN=MDuchi)
L_MDuchi_data <- as.data.frame(cbind(L_MDuchi,gamma))

plot(x=gamma, y=L_MDuchi, type = "p")

L_MDuchi_data$color <- L_MDuchi
colnames(L_MDuchi_data)

fig <- plot_ly(L_MDuchi_data, x= ~V2, y= ~V3, z = ~L_MDuchi, marker = list(color =~color, colorscale = c('#FFE1A1', '#683531'), showscale = T))
fig <- fig %>% layout(
  title="MDuchi",
  scene = list(xaxis = list(title = 'L1', range=c(-5,5)),
               yaxis = list(title = 'L2', range=c(-5,5)),
               zaxis = list(title = 'L3', range = c(-5,5))))
fig


