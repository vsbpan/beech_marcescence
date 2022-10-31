install.packages("YOUR PATH HERE/herbivar_0.1.0.tar.gz", repos = NULL, type = "source")

library(tidyverse)
library(herbivar)

img <- load.image("images/segmented_processed/2022_May_isolate/tree1_20220517_leaf3.png") 
plot(img)
img2 <- crop_leaf(img %>% # Crop the blue background out and set as NA
                    add_px_size("dpi:600") %>% # attach real units to pixel size
                    thin(10)) # Reduce image resolution to save computation time and RAM issue
leaf_length(img2) # The maximum distance between two leaf boarder pixels (leaf length)
leaf_area(img2) # Leaf area

# Load all processed images, crop the leaf, then measure length and area.
file_names <-list.files("images/segmented_processed/2022_May_isolate/")
file_paths <- paste0("images/segmented_processed/2022_May_isolate/",file_names)

out.list <- list()
for (i in 1:109){
  img <- crop_leaf(load.image(file_paths[i]) %>% 
                      add_px_size("dpi:600") %>% 
                      thin(10),
                   thr = "99")
  out.list[[i]] <- c("length" = leaf_length(img)[2], "area"=leaf_area(img)[2])
  cat(i,"    \r")
}
#Resulting data.frame
d <- do.call("rbind",out.list) %>% as.data.frame()


d %>% 
  ggplot(aes(x=log(length.mm2), y= log(area.mm2))) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  labs(subtitle = "Beech r2 = 0.84") + 
  theme_bw(base_size = 15)


m<-lm(log(area.mm2) ~ log(length.mm2), data = d) # simple linear model with log-log transformation
# m<-lm(area.mm2 ~ length.mm2, data = d) # simple linear model with pretty much the same performance
summary(m) # R2
mean((exp(predict(m)) - (d$area.mm2))^2)^0.5 #RMSE
mean((predict(m) - log(d$area.mm2))^2)^0.5 # RMSE
exp(0.1098775)


