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
  return(g)
}

inv.hsin <- function(x){
  log(sqrt(x^2+1) + x)
}

hellinger_trans <- function(x){
  rs <- rowSums(x)
  rs[rs == 0] <- 1
  sqrt(x/rs)
}

d <- read_csv("raw_data/data_2022.csv") # Read data
names(d)[grepl("B...",names(d))] <- c("B1", "B2") # Rename problematic column
d$B1[is.na(d$B1)] <- 0 # Replace NA with 0

# Apparently these are the morpho species that were considered good mites
all.equal(d$total.good.mites,
          d$RS + d$LB + d$OB + d$LL + d$LC + d$CLB + d$PR + d$BY + d$nymph+ d$B1 + d$B2 + d$SP + d$MISC)

# Change to factor for AR1
d$exp.round <- factor(d$exp.round)
diff(unique(d$day.of.year)) # Time interval looks good
d$B <- d$B1 + d$B2 # Merge redundant column -- probably the same morphospecies

d$total.good.mites.adults <- d$total.good.mites - d$nymph

# Species considered in the composition analysis
species_col <- c("RS", "LB", "OB", "LL", "B", "LC", "CLB", "SP", "MISC", "BY", "PR")
community_mat <- d[species_col]
pca.out <- rda(hellinger_trans(community_mat) ~ 
                   Condition(tree), 
                 data = d) # Simple PCA, keeping track of the tree id 
screeplot(pca.out) # Looks like it is enough to grab the first two PCs
d <- d %>% 
  mutate(PC1 = scores(pca.out)$sites[,1], 
         PC2 = scores(pca.out)$sites[,2], 
         shannon.H = vegan::diversity(community_mat)) # Add them to the data as an index of community composition





# Some exploratory plots
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







# Domatia count is not affected by treatment 
m1 <- glmmTMB(
  dom ~ treatment + factor(day.of.year) + (1|tree) + ar1(0 + exp.round | tree), 
  family = nbinom2(),
  data = d
);summary(m1)


# eriophyid presence is affected total good mites and not by treatment
# inv.hsin is the inverse hyperbolic sine transformation. It is similar to the log transformation, but better than log because it can deal with 0 (Burbidge et al. 1988). Changing the transformation to log(x+1) does not change the result 
m2 <- glmmTMB(
  eriophyid.pres ~ treatment + factor(day.of.year) + 
    (1|tree)  + ar1(0 + exp.round | tree), 
  family = binomial(),
  data = d, 
  control=glmmTMBControl(optimizer=optim,
                         optArgs=list(method="BFGS"))
); summary(m2)

m2.1 <- glmmTMB(
  eriophyid.pres ~ inv.hsin(total.good.mites) + treatment + factor(day.of.year) + 
    (1|tree)  + ar1(0 + exp.round | tree), 
  family = binomial(),
  data = d
); summary(m2)
check_overdispersion(m2.1)
check_model(m2.1)

# Fungal count is not affected by treatment
# Fungal count is affected by good mite community composition, not so much mite abundance after you account for the community composition 
m3 <- glmmTMB(
  fungal.count ~ treatment + factor(day.of.year) + (1|tree) + ar1(0 + exp.round | tree), 
  family = nbinom2(),
  data = d
);summary(m3)
check_model(m3)

m3.1 <- glmmTMB(
  fungal.count ~ inv.hsin(total.good.mites) + factor(day.of.year) + (1|tree) + ar1(0 + exp.round | tree), 
  family = nbinom2(),
  data = d
);summary(m3.1)

m3.2 <- glmmTMB(
  fungal.count ~ PC1 + PC2 + inv.hsin(total.good.mites) + factor(day.of.year) + (1|tree) + ar1(0 + exp.round | tree), 
  family = nbinom2(),
  data = d
);summary(m3.2)

m3.3 <- glmmTMB(
  fungal.count ~ PC1 + PC2 + factor(day.of.year) + (1|tree) + ar1(0 + exp.round | tree), 
  family = nbinom2(),
  data = d
);summary(m3.3)

anova(m3.1,m3.2,m3.3)


# Total good mite is increased by domatia abundance and not treatment. There might be a domatia and treatment interaction, but it is weak. 
m4 <- glmmTMB(
  total.good.mites ~ inv.hsin(dom) * treatment + factor(day.of.year) +  
    (1|tree) + ar1(0 + exp.round | tree), 
  family = nbinom2(),
  data = d, 
  control=glmmTMBControl(optimizer=optim,
                         optArgs=list(method="BFGS"))
);summary(m4)

m4.1 <- glmmTMB(
  total.good.mites ~ inv.hsin(dom) + treatment + factor(day.of.year) +  
    (1|tree) + ar1(0 + exp.round | tree), 
  family = nbinom2(),
  data = d, 
  control=glmmTMBControl(optimizer=optim,
                         optArgs=list(method="BFGS"))
);summary(m4.1)

check_model(m4)
plot_model(m4, type = "pred", terms = c("dom", "treatment"))

# Mite diversity is increased by domatia count but not by treatment
m5 <- glmmTMB(
  shannon.H ~ inv.hsin(dom) + treatment + factor(day.of.year) +  
    (1|tree) + ar1(0 + exp.round | tree), 
  family = tweedie(),
  data = d,
  control=glmmTMBControl(optimizer=optim,
                         optArgs=list(method="BFGS"))
);summary(m5)



#Here, I fitted a tb-partial RDA to analyze the community position 
#the hellinger transformation is basically a square root transformation of the relative species abundance for each tree. 
rda.out <- rda(hellinger_trans(community_mat) ~ 
                   inv.hsin(dom) + day.of.year + treatment + 
                 Condition(tree), 
                 data = d)
rda.out
anova(rda.out, by = "term",permutations = 1000) # Permutation test shows that only domatia and day of year explained variance in the community composition. 
anova(rda.out, by = "axis",permutations = 1000) # Only the first constrain axis explained significant amount of variance. 

get_biplot(pca.out,choices = c(1,2), scaling = 2, 
           display =  c("sites", "species", "biplot"), 
           d$treatment) + labs(title = "PCA")

get_biplot(rda.out,choices = c(1,2), scaling = 2, 
           display =  c("sites", "species", "biplot"), 
           d$treatment) + labs(title = "tb-partial-RDA")

get_biplot(rda.out,choices = c(1,3), scaling = 2, 
           display =  c("sites", "species", "biplot"), 
           d$treatment) + labs(title = "tb-partial-RDA")














