library(tidyverse)
library(glmmTMB)
library(performance)
library(sjPlot)
library(vegan)



get_biplot <- function(x, choices = c(1,2), scaling = 2, 
                       display = c("sites", "species", "biplot", "centroids"),
                       group = NULL){
  display <- match.arg(display,several.ok = TRUE)
  s <- scores(x,choices = choices,scaling = scaling)
  name <- names(s)
  if(is.null(name)){
    name <- "sites"
    s <- list(s)
  }
  s <- lapply(seq_along(name), function(i,s,name){
    x <- as.data.frame(s[[i]])
    x <- cbind(x, "a" = rownames(x))
    names(x)[length(x)] <- name[i]
    x
  }, s = s, name = name)
  names(s) <- name
  dim1 <- names(s$sites)[1]
  dim2 <- names(s$sites)[2]
  s$dummy <- data.frame(1,2)
  names(s$dummy) <- c(dim1, dim2)
  
  g <- s$dummy %>% 
    ggplot(aes_string(paste0("x = ", dim1),paste0("y = ", dim2))) +
    geom_vline(aes(xintercept=0),linetype="dashed",color="grey",size=1) + 
    geom_hline(aes(yintercept=0),linetype="dashed",color="grey",size=1) + 
    theme_bw() + 
    labs(x = dim1, y= dim2)
  
  if("sites" %in% display){
    if(!is.null(group)){
      s$sites <- cbind(s$sites, "group" = group)
      g <- g + geom_point(data = s$sites, aes(color = group))
    } else {
      g <- g + geom_point(data = s$sites,color = "deepskyblue")
    }
  }
  
  if("species" %in% display && !is.null(s$species)){
    g <- g + geom_text(data = s$species, 
                       aes(label = species), color = "violetred")
  }
  
  if("biplot" %in% display && !is.null(s$biplot)){
    s$biplot$x <- 0
    s$biplot$y <- 0
    s$biplot$xend <- s$biplot[,dim1]
    s$biplot$yend <-s$biplot[,dim2]
    g <- g + geom_text(data = s$biplot, 
                       aes(label = biplot), color = "black") + 
      geom_segment(data = s$biplot, aes(x = x, 
                                        y = y, 
                                        xend = xend,
                                        yend = yend), 
                   arrow = arrow(length = unit(0.2, "cm")), 
                   color = "black", 
                   size = 1)
  }
  
  if("centroids" %in% display && !is.null(s$centroids)){
    g <- g + geom_text(data = s$centroids, 
                       aes(label = centroids), color = "darkolivegreen")
  }
  print(g)
  return(g)
}

inv.hsin <- function(x){
  log(sqrt(x^2+1) + x)
}



d <- read_csv("raw_data/data_2022.csv")
names(d)[grepl("B...",names(d))] <- c("B1", "B2")
d$B1[is.na(d$B1)] <- 0

all.equal(d$total.good.mites,
          d$RS + d$LB + d$OB + d$LL + d$LC + d$CLB + d$PR + d$BY + d$nymph+ d$B1 + d$B2 + d$SP + d$MISC)
d$exp.round <- factor(d$exp.round)
diff(unique(d$day.of.year))

d$B <- d$B1 + d$B2

d %>% 
  filter(tree == 1) %>% 
  View()






d %>% 
  ggplot(aes(x= treatment, y = dom, color = factor(day.of.year))) + 
  geom_boxplot() + 
  geom_point(position = position_jitterdodge())


d %>% 
  ggplot(aes(x= treatment, y = total.good.mites, color = factor(day.of.year))) + 
  geom_boxplot() + 
  geom_point(position = position_jitterdodge())

d %>% 
  ggplot(aes(x= treatment, y = eriophyid.pres, fill = factor(day.of.year))) + 
  geom_bar(stat = "summary", position = "dodge") +
  geom_errorbar(stat = "summary", position = "dodge") +
  geom_point(position = position_jitterdodge())

d %>% 
  ggplot(aes(x= treatment, y = fungal.count, color = factor(day.of.year))) + 
  geom_boxplot() + 
  geom_point(position = position_jitterdodge())


m1 <- glmmTMB(
  dom ~ treatment + factor(day.of.year) + (1|tree) + ar1(0 + exp.round | tree), 
  family = nbinom2(),
  data = d
);summary(m1)

m2 <- glmmTMB(
  eriophyid.pres ~ inv.hsin(total.good.mites) + treatment + factor(day.of.year) + 
    (1|tree)  + ar1(0 + exp.round | tree), 
  family = binomial(),
  data = d
); summary(m2)
check_overdispersion(m2)
check_model(m2)

m3 <- glmmTMB(
  fungal.count ~ treatment + factor(day.of.year) + (1|tree) + ar1(0 + exp.round | tree), 
  family = nbinom2(),
  data = d
);summary(m3)
check_model(m3)

m4 <- glmmTMB(
  total.good.mites ~ inv.hsin(dom) * treatment + factor(day.of.year) +  
    (1 |tree) + ar1(0 + exp.round | tree), 
  family = nbinom2(),
  data = d
);summary(m4)
check_model(m4)



m5 <- glmmTMB(
  fungal.count ~ inv.hsin(total.good.mites) + treatment + factor(day.of.year) + (1 |tree)  + 
    ar1(0 + exp.round | tree), 
  family = nbinom2(),
  data = d 
); summary(m5)


m6 <- glmmTMB(
  fungal.count ~ MDS1 + MDS2 + inv.hsin(total.good.mites) +  
    factor(day.of.year) + (1 |tree)  + 
    ar1(0 + exp.round | tree), 
  family = nbinom2(),
  data = d 
); summary(m6)


inv.hsin(d$total.good.mites)




names(d)

species_col <- c("RS", "LB", "nymph", "OB", "LL", "B", "LC", "CLB", "SP", "MISC", "BY", "PR")
d <- d %>% 
  mutate(MDS1 = scores(pca.out)[,1], 
         MDS2 = scores(pca.out)[,2])


pca.out <- dbrda((d[species_col]) ~ Condition(tree), data = d)

pca.out <- dbrda((d[species_col]) ~ inv.hsin(dom) + day.of.year * treatment + 
                 Condition(tree), data = d)
pca.out
anova(pca.out, by = "term",permutations = 1000)

undebug(get_biplot)

get_biplot(pca.out,choices = c(1,2), scaling = 2, 
           display =  c("sites", "species", "biplot"), 
           d$treatment)






