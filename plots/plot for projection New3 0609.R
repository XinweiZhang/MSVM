library(plotly)
rm(list=ls())


New3 <- function(gamma){
  L1 <- max(0, 1- gamma[1], (gamma[2]-gamma[1]+1)/2)
  L2 <- max(0, 1- gamma[2], (gamma[1]-gamma[2]+1)/2)
  L3 <- max(gamma[1],0)  + max(gamma[2],0)
  return(c(L1,L2,L3))
}

# gamma =matrix(runif(4000,-3,3),ncol=2)
# 
gamma1 = rbind(c(0,0), c(-1, 0), c(-2, 0), c(-3, 0), c(-4, 0),  c(-5, 0),  c(-6, 0),
               c(0,-1), c(0,-2), c(0,-3), c(0,-4), c(0,-5),c(0,-6),
               c(1,0), c(1.5,.5),c(2,1),c(2.5,1.5),c(3,2),c(3.5,2.5),c(4,3),
               c(1,-1), c(1,-2),c(1,-3),c(1,-4),c(1,-5),c(1,-6),
               c(0,1), c(.5,1.5),c(1,2),c(1.5,2.5),c(2,3),c(2.5,3.5),c(3,4),
               c(-1,1), c(-2,1),c(-3,1),c(-4,1),c(-5,1),c(-6,1))

t(apply(gamma1[,1:2],MARGIN=1,FUN=New3))

gamma2 = rbind(c(0,0), c(-.5,-.5),c(-1,-1),c(-1.5,-1.5),c(-2,-2),c(-2.5,-2.5),c(-3,-3),c(-3.5,-3.5),
               c(1,0), c(2,-1), c(3,-2), c(4,-3), c(5,-4), c(6,-5),
               c(0,1), c(-1,2), c(-2,3), c(-3,4), c(-4,5), c(-5,6),
               c(0,2),c(0,3),c(0,4),c(0,5),c(0,6),c(0,7),
               c(2,0),c(3,0),c(4,0),c(5,0),c(6,0),c(7,0))
t(apply(gamma2[,1:2],MARGIN=1,FUN=New3))


gamma=rbind(gamma1,gamma2)


# gamma <- t(apply(gamma,MARGIN=1, FUN = function(x){ y <- x + rep((1- sum(x))/3,3)}))

xyz_gamma <- 1-cbind(gamma, 1-rowSums(gamma))

L_New3 <- t(apply(gamma,MARGIN=1,FUN=New3))

L_New3_data <- as.data.frame(rbind(L_New3,xyz_gamma))

L_New3_data$color <- rowSums(L_New3_data)

# 
# 
# fig <- plot_ly(L_New3_data, x= ~V1, y= ~V2, z = ~V3, marker = list(color =~color, colorscale = c('#FFE1A1', '#683531'), showscale = T))
# 
# fig <- fig %>% 
#   layout(
#     title="New3",
#     scene = list(xaxis = list(title = 'L1', range=c(-2,3)),
#                  yaxis = list(title = 'L2', range=c(-2,3)),
#                  zaxis = list(title = 'L3', range = c(-2,3)))) 
# 
# fig

m <- list( symbol = 200, size = 5, line = list( color = toRGB("yellow"), width = 2 ) ) 
m1 <- list( symbol = "triangle", color="black", sizemod = "diameter", size = 5, line = list( color = "black", width = 2 ) )

fig <- plot_ly(x = L_New3[,1], y = L_New3[,2], z = L_New3[,3]) %>%
  add_trace(marker = m1, type = "scatter3d", mode = "text+markers", 
            name ="original",   linetypes = NULL) %>%
  add_trace(x = xyz_gamma[,1], y = xyz_gamma[,2], z = xyz_gamma[,3], marker = m, type = "scatter3d", mode = "text+markers", 
            name = "projected",linetypes = NULL) %>% 
  add_trace(x = rbind(L_New3[10,], xyz_gamma[10,])[,1], y = rbind(L_New3[10,],xyz_gamma[10,])[,2], z = rbind(L_New3[10,],xyz_gamma[10,])[,3],
            type = "scatter3d", mode = "lines", name = "lines", showlegend = FALSE)


for(i in 1:nrow(xyz_gamma)){
  fig <- fig %>% add_trace(x = rbind(L_New3[i,], xyz_gamma[i,])[,1], y = rbind(L_New3[i,],xyz_gamma[i,])[,2], z = rbind(L_New3[i,], xyz_gamma[i,])[,3],
                           type = "scatter3d", mode = "lines", name = "lines", showlegend = FALSE)
}

fig <- fig %>% layout(title="New3", scene = list(xaxis = list(title = 'L1'),
                                                 yaxis = list(title = 'L2'),
                                                 zaxis = list(title = 'L3')))

fig


##################### 1 
a <- c(1,0,6)-c(1,0,1)
b <- c(-1.5,-2.5, 6)-c(1,0,1)
acos(sum(a*b)/(sqrt(sum(a^2))*sqrt(sum(b^2))))/pi*180

a <- c(1,0,6)-c(-1.5,-2.5, 6)
b <- c(1,0,1) - c(-1.5,-2.5, 6)
acos(sum(a*b)/(sqrt(sum(a^2))*sqrt(sum(b^2))))/pi*180


##########################2 

a <- c(6,0,6)-c(1,0,1)
b <-  c(6,-5,1)-c(1,0,1)

acos(sum(a*b)/(sqrt(sum(a^2))*sqrt(sum(b^2))))/pi*180

a <- c(6,0,6) -  c(6,-5,1)
b <- c(1,0,1) -  c(6,-5,1)

acos(sum(a*b)/(sqrt(sum(a^2))*sqrt(sum(b^2))))/pi*180


a <- c(6,-5,1)-c(1,0,1)
b <-  c(-1.5,-2.5,6)-c(1,0,1)

acos(sum(a*b)/(sqrt(sum(a^2))*sqrt(sum(b^2))))/pi*180

###################3


a <- c(6,1,-5)-c(1,1,0)
b <-  c(6,1,0)-c(1,1,0)

acos(sum(a*b)/(sqrt(sum(a^2))*sqrt(sum(b^2))))/pi*180

a <- c(6,1,0) - c(6,1,-5) 
b <-  c(1,1,0) - c(6,1,-5)

acos(sum(a*b)/(sqrt(sum(a^2))*sqrt(sum(b^2))))/pi*180


a <- c(6,1,-5)-c(1,1,0)
b <-  c(6,6,-10)-c(1,1,0)

acos(sum(a*b)/(sqrt(sum(a^2))*sqrt(sum(b^2))))/pi*180


###################4

a <- c(6,6,0)-c(1,1,0)
b <-  c(6,6,-10)-c(1,1,0)

acos(sum(a*b)/(sqrt(sum(a^2))*sqrt(sum(b^2))))/pi*180



a <- c(6,6,0)-c(1,1,0)
b <- c(6,6,0)- c(6,6,-10)
acos(sum(a*b)/(sqrt(sum(a^2))*sqrt(sum(b^2))))/pi*180


