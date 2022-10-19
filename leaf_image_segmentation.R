library(tidyverse)
library(pliman)

cut_kmeans <- function (x, km = c(1L, 2L)) {
  if (km[1] >= km[2]) {
    stop("Chosen boundary must be less than the number of clusters")
  }
  if (km[2] < 2) {
    stop("Number of clusters must be greater than two")
  }
  x <- x[!is.na(x)]
  km_out <- kmeans(x, km[2])
  index <- which(km_out$centers == sort(km_out$centers, decreasing = TRUE)[km[2]])
  thr <- max(x[km_out$cluster == index])
  return(thr)
}



split_leaves <- function(dir.path, save.dir.path, index = NULL){
  file.names <- list.files(dir.path,pattern = ".jpeg")
  if(!is.null(index)){
    file.names <- file.names[index]
  }
  file.path <- paste0(dir.path, file.names)
  for (i in seq_along(file.path)){
    img <- pliman::image_import(file.path[i])
    img2.seg<-image_segment(img,index = "R",
                            invert = FALSE, 
                            show_image = FALSE, 
                            fill_hull = TRUE)
    cat(i,"/",length(file.path), "Finding optimal threshold...", "                        \r")
    km.thresh<- cut_kmeans(c(img2.seg$R$image),km = c(2,4))
    cat(i,"/",length(file.path), "Segmenting image...", "                                   \r")
    img.ids<- object_id(
      img2.seg$R$image,
      watershed = FALSE,
      fill_hull = TRUE, 
      threshold = km.thresh, 
      show_image = FALSE
    )
    ids <-sort(img.ids$results$id[order(img.ids$results$area,decreasing = TRUE)][1:3],
               decreasing = FALSE)
    if(nrow(img.ids$results) < 3){
      km.thresh <- "Otsu"
      
      img.ids<- object_id(
        img2.seg$R$image,
        watershed = FALSE,
        fill_hull = TRUE, 
        threshold = km.thresh, 
        show_image = FALSE
      )
      ids <-sort(img.ids$results$id[order(img.ids$results$area,decreasing = TRUE)][1:3],
                 decreasing = FALSE)
    }
    
    for (j in seq_len(3)){
      iso <- object_isolate(img2.seg$R$image,
                            watershed = FALSE,
                            fill_hull = TRUE, 
                            threshold = km.thresh,
                            id = ids[j]) 
      image_export(iso,name = paste0(save.dir.path,gsub(".jpeg","",file.names[i]),"_leaf",j,".png"))
      cat(i,"/",length(file.path), "Saving leaf image", j, "                                 \r")
    }
    cat(i,"/",length(file.path),"                                    \r") 
  }
}


split_leaves(dir.path = "images/marcescence_project_2022_May/",
             save.dir.path = "images/2022_May_isolate/",
             index = 14:64)

split_leaves(dir.path = "images/marcescence_project_2022_June/",
             save.dir.path = "images/2022_June_isolate/")

split_leaves(dir.path = "images/marcescence_project_2022_August/",
             save.dir.path = "images/2022_August_isolate/",
             index = 46:61)
