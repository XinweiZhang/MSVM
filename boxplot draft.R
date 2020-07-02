
boxplot(t(res.array[1,1,,]), names = c("WW", "CS", "Duchi", "MDuchi", "New1,d","New1,0", "New3,d","New3,0", "OVA", "LLW", "MSVM7,d","MSVM7,0", "MSVM8"))

boxplot(t(res.array[1,2,,]), names = c("WW", "CS", "Duchi", "MDuchi", "New1,d","New1,0", "New3,d","New3,0", "OVA", "LLW", "MSVM7,d","MSVM7,0", "MSVM8"))
abline(h =  .6925, col= "red", lty = 2)



boxplot(t(res.array[1,1,-c(6,8,12),]), names = c("WW", "CS", "Duchi", "MDuchi", "New1", "New3", "OVA", "LLW", "MSVM7", "MSVM8"))
round(summary.res.arry[,1,]*100,2)
boxplot(t(res.array[1,2,-c(6,8,12),]), names = c("WW", "CS", "Duchi", "MDuchi", "New1", "New3", "OVA", "LLW", "MSVM7", "MSVM8"))
round(summary.res.arry[,2,]*100,2)
