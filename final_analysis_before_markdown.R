library(tidyverse)
library(glmmTMB)
library(performance)
library(sjPlot)
library(vegan)
library(herbivar)
library(ggpubr)





d <- read_csv("cleaned_data/beech_marcescence_cleaned.csv") # Read data

# Change to factor for AR1
d$exp.round <- factor(d$exp.round)

d$total.good.mites <- with(d, RS + LB + OB + LL + B + LC + CLB + SP + MISC + BY + PR + nymph)

d$total.good.mites.adults <- d$total.good.mites - d$nymph

# Species considered in the composition analysis
species_col <- c("RS", "LB", "OB", "LL", "B", "LC", "CLB", "SP", "MISC", "BY", "PR")
community_mat <- d[species_col]

# Simple PCA, keeping track of the tree id
# Hellinger transform the community matrix
pca.out <- rda(hellinger_trans(community_mat) ~ 
                 Condition(tree), 
               data = d)

screeplot(pca.out) # Looks like it is enough to grab the first two PCs
d <- d %>% 
  mutate(PC1 = scores(pca.out)$sites[,1], 
         PC2 = scores(pca.out)$sites[,2], 
         shannon.H = vegan::diversity(community_mat)) # Add them to the data as an index of community composition

d$dom_log <- log(d$dom)
d$leaf_area_mm2_log <- log(d$leaf_area_mm2)
d$herb_mm2 <- d$prop_herb * d$leaf_area_mm2






herb.null <- glmmTMB(
  herb_mm2 ~ 
    (1|tree)  + 
    ar1(0 + exp.round | tree), 
  family = tweedie(link = "log"),
  data = d
); summary(herb.null)

mcfadden_r2(herb.m1, herb.null)



mcfadden_r2 <- function(m1, m0, adjust = FALSE){
  m1ll <- logLik(m1)
  m0ll <- logLik(m0)
  k <- (attr(m1ll, "df") - attr(m0ll, "df"))
  stopifnot(k >= 0)
  if(adjust){
    as.numeric(1 - ((m1ll - k) / m0ll) )
  } else {
    as.numeric(1 - (m1ll / m0ll))
  }
}



m2 <- glmmTMB(
  log(leaf_area_mm2 )~ 
    treatment + 
    factor(day.of.year) +
    (1|tree)  + ar1(0 + exp.round | tree), 
  family = gaussian,
  data = d2, 
  control=glmmTMBControl(optimizer=optim,
                         optArgs=list(method="BFGS"))
); summary(m2)







library(brms)







psem.out

devtools::install_version("piecewiseSEM", version = "1.2.1", repos = "http://cran.us.r-project.org")




psem.out




psem.out2 <- piecewiseSEM::sem.fit(mod.list,
                                  data = d)
sem.coefs(modelList = mod.list, 
          data = d)



sem.plot(mod.list, data =d)






















