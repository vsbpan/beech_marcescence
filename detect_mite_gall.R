library(tidyverse)
library(pliman)
file.dir <-"images/2022_May_isolate/"
file.names <-list.files(file.dir,pattern = ".png")


out.list <- list()
for (i in seq_along(file.names)){
  cur.file.path<-paste0(file.dir,file.names[i])
  px <- measure_mite_gall(img.path = cur.file.path,
                    qc.check.path = gsub(".png","_qc.png",paste0(file.dir,"qc/",file.names[i]))
  )
  out.list[[i]]<- data.frame("img_path" = cur.file.path,"px" = px)
  cat(i, "/", length(file.names),"               \r")
}


d<-do.call("rbind",out.list)


log(d$px+1) %>% hist()
