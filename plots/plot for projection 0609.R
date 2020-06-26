library(plotly)
rm(list=ls())
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
# 
gamma =matrix(runif(6000,-2,2),ncol=3)

# gamma = rbind(c(1,2,3), c(1,2,4),
#               c(1,2,5), c(1,2,6),
#               c(1,2,7), c(1,2,8),
#               c(1,2,9), c(1,2,10),
#               c(1,2,11), c(1,2,12),
#               c(1,2,13), c(1,2,14),
#               c(1,2,3),c(1,3,3),
#               c(1,4,3),c(1,5,3),
#               c(1,6,3),c(1,7,3),
#               c(1,8,3),c(1,9,3),
#               c(1,10,3),c(1,11,3),
#               c(1,12,3),c(1,13,3),
#               c(1,2,3), c(2,2,3),
#               c(3,2,3), c(4,2,3),
#               c(5,2,3), c(6,2,3),
#               c(7,2,3), c(8,2,3),
#               c(9,2,3), c(10,2,3),
#               c(11,2,3), c(12,2,3),
#               c(13,2,3), c(14,2,3))-1

gamma <- t(apply(gamma,MARGIN=1, FUN = function(x){ y <- x + rep((1- sum(x))/3,3)}))
L_Duchi <- t(apply(gamma,MARGIN=1,FUN=Duchi))

proj <- function(L){
  x <- L[2] + L[3] - 2*L[1] - 2
  y <- L[1] + L[3] - 2*L[2] - 2
  z <- L[1] + L[2] - 2*L[3] - 2
  return(-c(x,y,z)/3)
}

xyz_gamma  <- t(apply(gamma, MARGIN=1, FUN = function(x){ 1-x }))

L_Duchi_proj <- t(apply(L_Duchi, MARGIN = 1, FUN = proj))

L_Duchi_data <- as.data.frame(rbind(L_Duchi,L_Duchi_proj))
  
L_Duchi_data$color <- rowSums(L_Duchi_data)



fig <- plot_ly(L_Duchi_data, x= ~V1, y= ~V2, z = ~V3, marker = list(color =~color, colorscale = c('#FFE1A1', '#683531'), showscale = T))

fig <- fig %>% layout(
  title="Duchi",
  scene = list(xaxis = list(title = 'L1', range=c(-10,20)),
               yaxis = list(title = 'L2', range=c(-10,20)),
               zaxis = list(title = 'L3', range = c(-10,20))))
fig

fig <- fig %>% layout(
  title="Duchi",
  scene = list(xaxis = list(title = 'L1', range=c(0,20)),
               yaxis = list(title = 'L2', range=c(0,20)),
               zaxis = list(title = 'L3', range = c(0,20))))
fig

c(12,0,10) - proj(c(12,0,10))

proj(c(2,0,3/4))
proj(c(2,1,0))


x<-  proj(c(2,1,0))[1]
y<-  proj(c(2,1,0))[2]
z<-  proj(c(2,1/4,3/4))[3]

(3*x-1)/2
(2*y + x -1)/2



m <- list( symbol = 200, size = 5, line = list( color = toRGB("yellow"), width = 2 ) ) 
m1 <- list( symbol = "triangle", color="black", sizemod = "diameter", size = 5, line = list( color = "black", width = 2 ) )

fig <- plot_ly(x = L_Duchi[,1], y = L_Duchi[,2], z = L_Duchi[,3]) %>%
  add_trace(marker = m1, type = "scatter3d", mode = "text+markers", 
            name ="original",   linetypes = NULL) %>%
  add_trace(x = xyz_gamma[,1], y = xyz_gamma[,2], z = xyz_gamma[,3], marker = m, type = "scatter3d", mode = "text+markers", 
            name = "projected",linetypes = NULL) %>% 
  add_trace(x = rbind(L_Duchi[10,], xyz_gamma[10,])[,1], y = rbind(L_Duchi[10,],xyz_gamma[10,])[,2], z = rbind(L_Duchi[10,],xyz_gamma[10,])[,3],
            type = "scatter3d", mode = "lines", name = "lines", showlegend = FALSE)


for(i in 1:nrow(xyz_gamma)){
  fig <- fig %>% add_trace(x = rbind(L_Duchi[i,], xyz_gamma[i,])[,1], y = rbind(L_Duchi[i,],xyz_gamma[i,])[,2], z = rbind(L_Duchi[i,], xyz_gamma[i,])[,3],
                           type = "scatter3d", mode = "lines", name = "lines", showlegend = FALSE)
}
fig <- fig %>% layout(
  title="Duchi")

fig
