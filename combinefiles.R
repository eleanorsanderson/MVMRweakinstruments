


rm(list = ls(all=TRUE))

results2_all <- NULL
file.names2 <- dir(pattern ="output2exnc")


for (i in 1:length(file.names2)){
  
  load(file.names2[i])
  
  results2_all <- rbind(results2_all, output)
  
}

save(results2_all, file=("Results_2nc.Rda"))



rm(list = ls(all=TRUE))

results2_all <- NULL
file.names3 <- dir(pattern ="output3ex")


for (i in 1:length(file.names3)){
  
  load(file.names3[i])
  
  results3_all <- rbind(results3_all, output3)
  
}

save(results3_all, file=("Results_3.Rda"))



rm(list = ls(all=TRUE))

results2_all <- NULL
file.names3 <- dir(pattern ="outputjackknife")


for (i in 1:length(file.names3)){
  
  load(file.names3[i])
  
  resultsjk_all <- rbind(resultsjk_all, output)
  
}

save(resultsjk_all, file=("Results_jk.Rda"))



