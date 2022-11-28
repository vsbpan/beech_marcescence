library(tidyverse)
library(pliman)



img <- pliman::image_import("images/marcescence_project_2022_May/tree1_20220517.jpeg")



plot(img)
dim(img2)
img2 <- image_crop(img,width = 100:3150, height = 300:5100)

plot(img2)

image_index(img)

img2.seg<-image_segment(img2,index = "R",invert = F)

plot(img2.seg$R$image)

img2con<-object_contour(img2.seg$R$image,
               watershed = FALSE,fill_hull = TRUE, 
               threshold = cut_kmeans(c(img2.seg$R$image),km = c(1,3)))

img2an<-analyze_objects(img2.seg$R$image,lower_size = 10000,
                        watershed = FALSE,
                        threshold = cut_kmeans(c(img2.seg$R$image),
                                               km = c(1,3)))

img2an$results

dim(img2)
img3 <- img2
img3[cbind(as.matrix(img2con$`3`),1)] <- 0

plot(img3)


plot(img2con$`3`)

str(img2an)

img2con$`3`

image_seg

object_id(img2.seg$R$image,
          watershed = FALSE,
          fill_hull = TRUE, 
          threshold = cut_kmeans(c(img2.seg$R$image),km = c(1,3)))

img.iso1<-object_isolate(img2.seg$R$image,
               watershed = FALSE,
               fill_hull = TRUE, 
               threshold = cut_kmeans(c(img2.seg$R$image),km = c(1,3)),
               id = 3)
class(img.iso1)
image_export(img.iso1,name = "abc.png")


