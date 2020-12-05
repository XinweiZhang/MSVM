MSVMPack_interface <- function(type, X_train, y_train, X_test = NULL, y_test = NULL,  C, sigma, tmp.dir = ""){

  now <- Sys.time()
  train_file_name <- paste0("msvm",format(now, "%m%d%H%M%S"),".train")
  fileConn<-file(train_file_name)
  writeLines(paste(nrow(X_train),"\n",paste(ncol(X_train)), sep=""), fileConn)
  close(fileConn)
  write.table(cbind(X_train, y_train), file = train_file_name, append = T, quote = F, sep = " ", col.names = F, row.names = F)
  
  test_file_name <- paste0("msvm", format(now, "%m%d%H%M%S"),".test")
  fileConn<-file(test_file_name)
  writeLines(paste(nrow(X_test),"\n",paste(ncol(X_test)), sep=""), fileConn)
  close(fileConn)
  write.table(cbind(X_test,y_test), file = test_file_name, append = T, quote = F, sep = " ", col.names = F, row.names = F)
  
  pred_file_name <- paste0("msvm", format(now, "%m%d%H%M%S"),".outputs")
  
  train_file_name <- paste0(tmp.dir, train_file_name)
  test_file_name <- paste0(tmp.dir, test_file_name)
  pred_file_name <- paste0(tmp.dir, pred_file_name)
  

  if(type == WW){
    svm_options = paste("-m WW -t 1","-c", C, "-k", kernel, "-p", sigma, "-q")
  }else if(type == "CS"){
    svm_options = paste("-m CS -t 1","-c", C, "-k", kernel, "-p", sigma, "-q")
  }else if(type == "LLW"){
    svm_options = paste("-m LLW -t 1","-c", C, "-k", kernel, "-p", sigma, "-q")
  }
  
  system(paste("trainmsvm", train_file_name, test_file_name, pred_file_name, svm_options))
  pred <- read.table(pred_file_name, header = FALSE, sep = " ")
  pred_l <- pred[,ncol(pred)]
  
  file.remove(c(train_file_name, test_file_name, pred_file_name))
  
  return(pred_l)
}

