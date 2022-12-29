library(tidyverse)
library(pliman)
library(herbivar)


img <- image_import("images/2022_June_isolate/tree58_20220701_leaf3.png")
img_disease <- image_import("images/2022_May_isolate/disease.png")
plot(img)

image_index(img, index="all")
image_index(img, index = c("S, BI, BIM, R, CI,NDRBI"))
imgf<-image_binary(img,index = "BI",threshold = 0.51,invert = TRUE) %>% 
  image_filter(size = 2)
image_binary(img,"all")
imgf$BI %>% 
  plot()


img.seg <- image_segment_iter(img, 
                   c("SHP","HI"),
                   nseg = 2,
                   invert = c(FALSE,TRUE),
                   threshold = c(0.1,
                                 "Otsu"))
plot

img.seg$images$seg3 %>% 
  image_binary(index = "all")

analyze_objects(
  img.seg$images$seg2,
  watershed = TRUE,
  object_size = "small",
  lower_size = 1000,
  upper_eccent = 1,
  lower_circ = 0.1,
  filter = 3
)

str(p)
plot(p[[1]])


analyze_objects(
  image_blur(img,4),
  index = "SCI",
  invert = F,
  filter = 3,
  upper_size = 100000000,
  upper_eccent = 0.99,
  watershed = FALSE)

analyze_objects(
  img,
  foreground = img_disease,
  index = c("GR,S"),
  invert = F,
  filter = 3,
  watershed = FALSE)

measure_disease(
  img,
  img_symptoms = p[[1]],
  img_background = p[[2]],
  img_healthy = p[[3]],
  watershed = TRUE,
  col_background = "blue",
  contour_col = "red",
)





debug(help_count)






list.files("images/marcescence_project_2022_August/") %>% 
          gsub("_.*","",.)





library(tidyverse)


paste0("tree",1:66) %>% 
  setdiff(list.files("images/marcescence_project_2022_August/") %>% 
            gsub("_.*","",.))















