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




d$tot

ga.m3 <- glmmTMB(
  gall_mm2 ~ 
    treatment +
    scale(dom_log) +
    scale(leaf_area_mm2_log) + 
    pred_density_score + 
    PC1 +
    PC2 + 
    exp.round + 
    (1|tree)  + 
    ar1(0 + exp.round | tree), 
  family = ziGamma(link = "log"),
  ziformula = ~ .,
  data = d,
  control=glmmTMBControl(optimizer=optim,
                         optArgs=list(method="BFGS"))
); summary(ga.m3)




ga.m2 <- glmmTMB(
  gall_mm2 ~ 
    treatment +
    scale(dom_log) +
    scale(leaf_area_mm2_log) + 
    pred_density_score + 
    exp.round + 
    (1|tree)  + 
    ar1(0 + exp.round | tree), 
  family = ziGamma(link = "log"),
  ziformula = ~ .,
  data = d,
  control=glmmTMBControl(optimizer=optim,
                         optArgs=list(method="BFGS"))
); summary(ga.m2)



d$pred_density_score %>% hist()

herb.m2 <- glmmTMB(
  herb_mm2 ~ 
    treatment +
    dom_log.scale +
    leaf_area_mm2_log.scale + 
    exp.round + 
    (1|tree)  + 
    ar1(0 + exp.round | tree), 
  family = tweedie(link = "log"),
  data = d
); summary(herb.m2)



d$fungal.count %>% hist()

fg.m2 <- glmmTMB(
   fungal.count~ 
    treatment +
    dom_log.scale +
    leaf_area_mm2_log.scale + 
    exp.round + 
    (1|tree)  + 
    ar1(0 + exp.round | tree), 
  family = poisson(link = "log"),
  data = d
); summary(fg.m2)


names(d)










ga.m1 <- glmmTMB(
  gall_mm2 ~ 
    treatment +
    dom_log.scale +
    leaf_area_mm2_log + 
    exp.round + 
    (1|tree)  + 
    ar1(0 + exp.round | tree), 
  family = ziGamma(link = "log"),
  ziformula = ~ .,
  data = d %>% 
    mutate(dom_log.scale = scale(dom_log)),
  control=glmmTMBControl(optimizer=optim,
                         optArgs=list(method="BFGS"))
); summary(ga.m1)

ga.m2 <- glmmTMB(
  gall_mm2 ~ 
    treatment +
    scale(dom_log) +
    scale(leaf_area_mm2_log) + 
    exp.round + 
    (1|tree)  + 
    ar1(0 + exp.round | tree), 
  family = tweedie(link = "log"),
  data = d,
  control=glmmTMBControl(optimizer=optim,
                         optArgs=list(method="BFGS"))
); summary(ga.m2)

ga.m1 %>% summary()
AIC(ga.m1, ga.m2)

DHARMa::simulateResiduals(den.m1) %>% plot()





herb.m2 <- glmmTMB(
  herb_mm2 ~ 
    treatment +
    scale(dom_log) +
    scale(leaf_area_mm2_log) + 
    exp.round + 
    (1|tree)  + 
    ar1(0 + exp.round | tree), 
  family = tweedie(link = "log"),
  data = d,
  control=glmmTMBControl(optimizer=optim,
                         optArgs=list(method="BFGS"))
); summary(herb.m2)










library(brms)



check_model(m1)




d2 <- d %>% 
  mutate(gall = ifelse(gall_mm2 > 0, 1, 0), 
         prop_gall = gall_mm2 / leaf_area_mm2) %>% 
  mutate(prop_gall.adj = adjust_prop(prop_gall,
                                     nudge.method = "replace",
                                     nudge.size = "warton_min",
                                     trans = "identity"), 
         herb_mm2 = leaf_area_mm2 * prop_herb)

hist(log(d2$herb_mm2))

summary(d2$herb_mm2)

m2 <- glmmTMB(
  herb_mm2 ~ 
    treatment +
    factor(day.of.year) +
    log(leaf_area_mm2) + 
    (1|tree)  + ar1(0 + exp.round | tree), 
  family = ziGamma(link = "log"),
  ziformula = ~. ,
  data = d2, 
  control=glmmTMBControl(optimizer=optim,
                         optArgs=list(method="BFGS"))
); summary(m2)





d2$herb_mm2

performance(m2)
check_model(m2)

plot_model(m2,type = "pred", terms = c("treatment"))

DHARMa::simulateResiduals(m2) %>% plot()



?family_glmmTMB


m3 <- glmmTMB(
  gall_mm2 ~ 
    treatment +
    factor(day.of.year) +
    asinh(total.good.mites) + 
    PC1 + 
    PC2 +
    (1|tree)  + ar1(0 + exp.round | tree), 
  family = ziGamma(link = "log"),
  ziformula = ~. ,
  data = d2, 
  control=glmmTMBControl(optimizer=optim,
                         optArgs=list(method="BFGS"))
); summary(m3)

names(d2)


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

