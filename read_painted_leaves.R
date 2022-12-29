library(tidyverse)
library(pliman)
library(herbivar)






read_painted_leaves <- function(file.dir){
  file.names <-list.files(file.dir,pattern = ".png")
  out.list <- list()
  
  for (i in seq_along(file.names)){
    cur.file.path<-paste0(file.dir,file.names[i])
    img <- load.image(cur.file.path) %>% 
      add_px_size(px.size = "dpi:600")
    
    crp_img <- crop_leaf(img)
    
    leaf_area_mm2 <- crp_img %>% 
      leaf_area() %>% 
      .["mm2"]
    
    prop_herb <- crp_img %>% 
      leaf_herb(type = "prop")
    
    gall_mm2 <- immask(img, !threshold2(color_index(img,index = "R",plot = FALSE)$R) | 
                         threshold2(color_index(img,index = "G",plot = FALSE)$G), 
                       background = "blue") %>% 
      crop_leaf() %>% 
      leaf_area() %>% 
      .["mm2"]
    
    out.list[[i]]<- data.frame("img_path" = cur.file.path,
                               "gall_mm2" = gall_mm2, 
                               "prop_herb" = prop_herb, 
                               "leaf_area_mm2" = leaf_area_mm2)
    
    cat(i, "/", length(file.names),"               \r")
  }
  
  out <- do.call("rbind", out.list)
  rownames(out) <- NULL
  return(out)
}



may_d <- read_painted_leaves(file.dir = "images/segmented_processed/2022_May_isolate/")
june_d <- read_painted_leaves(file.dir = "images/segmented_processed/2022_June_isolate/")
august_d <- read_painted_leaves(file.dir = "images/segmented_processed/2022_August_isolate/")


d <- rbind(may_d,june_d,august_d)
d$gall_mm2[d$gall_mm2/d$leaf_area_mm2 < 0.001] <- 0
d$month <- gsub(".*2022_|_isolate.*","",d$img_path)
d$leaf <- as.numeric(gsub(".*_leaf|.png","",d$img_path))
d$tree <- as.numeric(gsub(".*tree|_.*","",d$img_path))
d$exp.round <- with(d,case_when(
  month == "August" ~ 3, 
  month == "June" ~ 2, 
  month == "May" ~ 1
))
d$gall_mm2 <- ifelse(d$gall_mm2/d$leaf_area_mm2 < 0.001, 
       0, 
       d$gall_mm2)
d$prop_herb <- ifelse(d$prop_herb < 0.005, 
                     0, 
                     d$prop_herb)


d <- read_csv("cleaned_data/beech_leaf_scan_read.csv")
d2 <- read_csv("raw_data/data_2022.csv")
d3 <- read_csv("raw_data/tree coordinates.csv")

d.full <- d2 %>% 
  left_join(d, by = c("tree","leaf","exp.round")) %>% 
  left_join(d3, by = "tree")


d.full$B <- ifelse(is.na(d.full$B...18),0,
       d.full$B...18) + d.full$B...24


d.full <- d.full %>% 
  select(-c(m.leaves, `B...18`,`B...24`,img_path))


write_csv(d.full,"beech_marcescence_cleaned.csv")
write_csv(d,"beech_leaf_scan_read.csv")


# 
# d.disagree <- d.full %>% 
#   mutate(gall = ifelse(gall_mm2 > 0, 1, 0)) %>% 
#   filter(gall != eriophyid.pres) %>% 
#   select(field.date,tree,leaf,gall_mm2,leaf_area_mm2,eriophyid.pres)



