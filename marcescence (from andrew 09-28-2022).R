
#going to have to make these changes permanent eventually -- see https://www.accelebrate.com/library/how-to-articles/r-rstudio-library
myPaths <- .libPaths()
myPaths <- c(myPaths, "C:/Users/Michaela TerAvest/OneDrive/MSU Weber Lab onedrive/Documents/R/win-library/4.1")
.libPaths(myPaths)
myPaths
myPaths <- .libPaths()   # get the paths
myPaths <- c(myPaths[3], myPaths[1], myPaths[2])  # switch them
.libPaths(myPaths)


setwd("C:/Users/Michaela TerAvest/OneDrive/MSU Weber Lab onedrive/Weber Lab/Marcescence paper/data")
#read in data
data <- read.csv(file="data_2022.csv", header = T, na.strings=c("","NA"))
colnames(data)[1] <- "tree" #for some reason this first column had a weird name, so changed it to tree
data$leaf<-as.factor(data$leaf) #change leaf to factor
data$tree<-as.factor(data$tree) #change leaf to factor
data$exp.round<-as.factor(data$exp.round) #change experimental round to factor

library(plyr)
library(ggthemes)
data.summary<-ddply(data, .(exp.round, treatment), summarize,
                      N=length(total.good.mites),
                      N.fungus=length(fungal.count[!is.na(fungal.count)]),
                      mean.good.mites=mean(total.good.mites),
                      mean.erineum=mean(eriophyid.pres),
                      mean.fungus=mean(fungal.count, na.rm=TRUE),
                      mean.eggs=mean(eggs),
                      sd.good.mites   = sd(total.good.mites),
                      sd.erineum   = sd(eriophyid.pres),
                      sd.fungus   = sd(fungal.count, na.rm=TRUE),
                      sd.eggs   = sd(eggs),
                      se.good.mites   = sd.good.mites / sqrt(N),
                      se.erineum   = sd.erineum / sqrt(N),
                      se.fungus   = sd.fungus / sqrt(N.fungus),
                      se.eggs   = sd.eggs / sqrt(N)
                    )
#make my colour palette
cols <- c("control" = "gold2", "treatment" = "gray")
library(ggplot2)
ggplot.good.mites<- ggplot(data.summary, 
                        aes(x=exp.round, y=mean.good.mites, shape=treatment, colour=treatment, fill=treatment))+
  geom_point()+
  geom_errorbar(aes(ymin=mean.good.mites-se.good.mites, ymax=mean.good.mites+se.good.mites), colour="black", width=.2, position="dodge")+
  scale_colour_manual(values=cols)+ 
  scale_fill_manual(values=cols)+
  geom_line(size=1)+
  theme_bw()+
  scale_x_continuous(expand = c(0, 0), limits = c(0.8, 3.1), breaks=c(1, 2, 3))+
  xlab("Sampling Round")+
  ylab("Mean Mites/Leaf")
ggsave('good.mites.png', width=4, height=2)
ggplot.good.mites

ggplot.eggs<- ggplot(data.summary, 
                           aes(x=exp.round, y=mean.eggs, shape=treatment, colour=treatment, fill=treatment))+
  geom_point()+
  geom_errorbar(aes(ymin=mean.eggs-se.eggs, ymax=mean.eggs+se.eggs), colour="black", width=.2, position="dodge")+
  scale_colour_manual(values=cols)+ 
  scale_fill_manual(values=cols)+
  geom_line(size=1)+
  theme_bw()+
  scale_x_continuous(expand = c(0, 0), limits = c(0.8, 3.1), breaks=c(1, 2, 3))+
  xlab("Sampling Round")+
  ylab("Mean Eggs/Leaf")
ggsave('eggs.png', width=4, height=2)
ggplot.eggs

ggplot.erineum<- ggplot(data.summary, 
                           aes(x=exp.round, y=mean.erineum, shape=treatment, colour=treatment, fill=treatment))+
  geom_point()+
  geom_errorbar(aes(ymin=mean.erineum-se.erineum, ymax=mean.erineum+se.erineum), colour="black", width=.2, position="dodge")+
  scale_colour_manual(values=cols)+ 
  scale_fill_manual(values=cols)+
  geom_line(size=1)+
  theme_bw()+
  scale_x_continuous(expand = c(0, 0), limits = c(0.8, 3.1), breaks=c(1, 2, 3))+
  xlab("Sampling Round")+
  ylab("Prop. w/ Erineum Galls")
ggsave('erineum.png', width=4, height=2)
ggplot.erineum
#models

ggplot.fungus<- ggplot(data.summary, 
                        aes(x=exp.round, y=mean.fungus, shape=treatment, colour=treatment, fill=treatment))+
  geom_point()+
  geom_errorbar(aes(ymin=mean.fungus-se.fungus, ymax=mean.fungus+se.fungus), colour="black", width=.2, position="dodge")+
  scale_colour_manual(values=cols)+ 
  scale_fill_manual(values=cols)+
  geom_line(size=1)+
  theme_bw()+
  scale_x_continuous(expand = c(0, 0), limits = c(0.8, 3.1), breaks=c(1, 2, 3))+
  xlab("Sampling Round")+
  ylab("Mean Fungal Ct/Leaf")
ggsave('fungus.png', width=4, height=2)
ggplot.fungus



#PSEM
#remove.packages("piecewiseSEM")
#devtools::install_version("piecewiseSEM", version = "1.2.1", repos = "http://cran.us.r-project.org")
library(piecewiseSEM)
#install.packages('TMB', type = 'source')
library(glmmTMB)

#make basis set
#this one includes all the data with an term for temporal autocorrelation
m1<-glmmTMB(total.good.mites ~ 
            treatment + dom + (1|tree) + ar1(0+exp.round|tree), 
            family=poisson, data=data)
summary(m1)

m2<-glmmTMB(eriophyid.pres ~ 
            treatment + total.good.mites + (1|tree) + ar1(0+exp.round|tree), 
            family=binomial, data=data)
summary(m2)

m3<-glmmTMB(fungal.count ~ 
            treatment + total.good.mites + eriophyid.pres + (1|tree) + ar1(0+exp.round|tree), 
            family=poisson, data=data)
summary(m3)

#one other thing I wanted to check. univariate model with good.mites ~ treatment
treatment.model<-glmmTMB(total.good.mites ~ 
                               treatment + (1|tree) + ar1(0+exp.round|tree), 
                             family=poisson, data=data)
summary(treatment.model)


#make model list
mod.list<-list(m1,m2,m3)
psem.out<-sem.fit(mod.list, data)
sem.coefs(mod.list,psem.out)
sem.plot(mod.list,psem.out)
dev.off()

#one other thing I wanted to check. univariate model with good.mites ~ treatment
treatment.model<-glmmTMB(total.good.mites ~ 
                           treatment + (1|tree) + ar1(0+exp.round|tree), 
                         family=poisson, data=data)
summary(treatment.model)

#ok now I want to do the pSEM again, but one for every individual exp round
m1.1<-glmmTMB(total.good.mites ~ 
              treatment + dom + (1|tree), 
            family=poisson, data=data[which(data$exp.round==1),])
summary(m1.1)

m2.1<-glmmTMB(eriophyid.pres ~ 
              treatment + total.good.mites + (1|tree), 
            family=binomial, data=data[which(data$exp.round==1),])
summary(m2.1)

#make model list
mod.list.1<-list(m1.1,m2.1)
psem.out.1<-sem.fit(mod.list.1, data[which(data$exp.round==1),])
sem.coefs(mod.list.1,psem.out.1)

#for the first round, no fungus model, because we didn't see any fungi on the leaves
#m3.1<-glmmTMB(fungal.count ~ 
#              treatment + total.good.mites + eriophyid.pres + (1|tree), 
#            family=poisson, data=data[which(data$exp.round==1),],)
#summary(m3.1)




m1.2<-glmmTMB(total.good.mites ~ 
                treatment + dom + (1|tree), 
              family=poisson, data=data[which(data$exp.round==2),])
summary(m1.2)

m2.2<-glmmTMB(eriophyid.pres ~ 
                treatment + total.good.mites + (1|tree), 
              family=binomial, data=data[which(data$exp.round==2),])
summary(m2.2)

m3.2<-glmmTMB(fungal.count ~ 
              treatment + total.good.mites + eriophyid.pres + (1|tree), 
            family=poisson, data=data[which(data$exp.round==2),],)
summary(m3.2)

#make model list
mod.list.2<-list(m1.2,m2.2, m2.3)
psem.out.2<-sem.fit(mod.list.2, data[which(data$exp.round==2),])
sem.coefs(mod.list.2,psem.out.2)


m1.3<-glmmTMB(total.good.mites ~ 
                treatment + dom + (1|tree), 
              family=poisson, data=data[which(data$exp.round==3),])
summary(m1.3)

m2.3<-glmmTMB(eriophyid.pres ~ 
                treatment + total.good.mites + (1|tree), 
              family=binomial, data=data[which(data$exp.round==3),])
summary(m2.3)

m3.3<-glmmTMB(fungal.count ~ 
                treatment + total.good.mites + eriophyid.pres + (1|tree), 
              family=poisson, data=data[which(data$exp.round==3),],)
summary(m3.3)

#make model list
mod.list.3<-list(m1.3,m2.3, m3.3)
psem.out.3<-sem.fit(mod.list.3, data[which(data$exp.round==3),])
sem.coefs(mod.list.3,psem.out.3)