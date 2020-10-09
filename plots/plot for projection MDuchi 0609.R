library(plotly)
library(combinat)
rm(list=ls())

MDuchi <- function(gamma){
  s_gamma1 <- sort(gamma[c(2,3)], decreasing = T)
  s_gamma2 <- sort(gamma[c(1,3)], decreasing = T)
  s_gamma3 <- sort(gamma[c(1,2)], decreasing = T)

  L1 <- max(0,(s_gamma1[1]-gamma[1]+1)/2, ((s_gamma1[1]-gamma[1]+1)/3+ (s_gamma1[2]-gamma[1]+1)/3))
  L2 <- max(0,(s_gamma2[1]-gamma[2]+1)/2, ((s_gamma2[1]-gamma[2]+1)/3+ (s_gamma2[2]-gamma[2]+1)/3))
  L3 <- max(0,(s_gamma3[1]-gamma[3]+1)/2, ((s_gamma3[1]-gamma[3]+1)/3+ (s_gamma3[2]-gamma[3]+1)/3))

  return(c(L1,L2,L3))
}


# gamma =matrix(runif(6000,-3,3),ncol=3)
permn(c(1,2,3))

MDuchi(c(0,0,1))

gamma1 = rbind(c(0,0,1),c(-1, .5,1.5),c(-2, 1,2),c(-3, 1.5,2.5),c(-4, 2,3),c(-5, 2.5,3.5),c(-6, 3,4), 
               c(-.5,-.5,2),c(-1,-1,3),c(-1.5,-1.5,4),c(-2,-2,5),c(-2.5,-2.5,6),c(-3,-3,7))

gamma <- do.call(rbind,lapply(permn(c(1,2,3)), function(x){gamma1[,x]}))

gamma <- t(apply(gamma,MARGIN=1, FUN = function(x){ y <- x + rep((1- sum(x))/3,3)}))
xyz_gamma <- t(apply(gamma, MARGIN=1, FUN = function(x){ 1-x }))
L_MDuchi <- t(apply(gamma,MARGIN=1,FUN=MDuchi))
L_MDuchi
# L_MDuchi

L_MDuchi_data <- as.data.frame(rbind(L_MDuchi,xyz_gamma))

L_MDuchi_data$color <- rowSums(L_MDuchi_data)



fig <- plot_ly(L_MDuchi_data, x= ~V1, y= ~V2, z = ~V3, marker = list(color =~color, colorscale = c('#FFE1A1', '#683531'), showscale = T))

fig <- fig %>% 
  layout(
  title="Duchi",
  scene = list(xaxis = list(title = 'L1', range=c(0,3)),
               yaxis = list(title = 'L2', range=c(0,3)),
               zaxis = list(title = 'L3', range = c(0,3)))) 

fig

m <- list( symbol = 200, size = 5, line = list( color = toRGB("yellow"), width = 2 ) ) 
m1 <- list( symbol = "triangle", color="black", sizemod = "diameter", size = 5, line = list( color = "black", width = 2 ) )
xyz_gamma[10,]

fig <- plot_ly(x = L_MDuchi[,1], y = L_MDuchi[,2], z = L_MDuchi[,3]) %>%
  add_trace(marker = m1, type = "scatter3d", mode = "text+markers", 
            name ="original",   linetypes = NULL) %>%
  add_trace(x = xyz_gamma[,1], y = xyz_gamma[,2], z = xyz_gamma[,3], marker = m, type = "scatter3d", mode = "text+markers", 
            name = "projected",linetypes = NULL) %>% 
  add_trace(x = rbind(L_MDuchi[10,], xyz_gamma[10,])[,1], y = rbind(L_MDuchi[10,],xyz_gamma[10,])[,2], z = rbind(L_MDuchi[10,],xyz_gamma[10,])[,3],
                                                               type = "scatter3d", mode = "lines", name = "lines", showlegend = FALSE)


for(i in 1:nrow(xyz_gamma)){
  fig <- fig %>% add_trace(x = rbind(L_MDuchi[i,], xyz_gamma[i,])[,1], y = rbind(L_MDuchi[i,],xyz_gamma[i,])[,2], z = rbind(L_MDuchi[i,], xyz_gamma[i,])[,3],
                         type = "scatter3d", mode = "lines", name = "lines", showlegend = FALSE)
}

fig <- fig %>% layout(title="MDuchi", scene = list(xaxis = list(title = 'L1'),
                               yaxis = list(title = 'L2'),
                               zaxis = list(title = 'L3')))


fig
####################1
a <- c(7,0,1)-c(1,0,1)
b <- c(7,-3,-2)-c(1,0,1)

acos(sum(a*b)/(sqrt(sum(a^2))*sqrt(sum(b^2))))/pi*180

a <- c(1,0,1) - c(7,0,1)
b <- c(7,-3,-2)-c(7,0,1)

acos(sum(a*b)/(sqrt(sum(a^2))*sqrt(sum(b^2))))/pi*180


#################################2



a <- c(4.75,0,4.75)-c(1,0,1)
b <-  c(3.5,-5,3.5)-c(1,0,1)

c(4.75,0,4.75) - c(3.5,-5,3.5)
sum(a*b)/(sqrt(sum(a^2))*sqrt(sum(b^2)))
acos(sum(a*b)/(sqrt(sum(a^2))*sqrt(sum(b^2))))/pi*180

a <- c(1,0,1) - c(4.75,0,4.75)
b <-  c(3.5,-5,3.5)-c(4.75,0,4.75)
acos(sum(a*b)/(sqrt(sum(a^2))*sqrt(sum(b^2))))/pi*180



