##

#Analysis for Age Cohort manuscript

#Contact: cprouty@ufl.edu    https://scholar.google.com/citations?user=PpeDx78AAAAJ&hl=en

##

#load packages
library(lsmeans)
library(multcomp)
library(afex)
library(brglm2)
library(car)
##

#Import files
BeeCounts <- read.csv("BeeCounts.csv")
DyeVisuals <- read.csv("DyeVisuals.csv")
##

#Data organization
BeeCounts$Colony <- as.factor(BeeCounts$Colony)
BeeCounts$SampleDay <- as.factor(BeeCounts$SampleDay)
names(BeeCounts)[1] <- "Date"
names(DyeVisuals)[1] <- "Treatment"
DyeVisuals$Hive <- as.factor(DyeVisuals$Hive)
DyeVisuals$SampleDay <- as.factor(DyeVisuals$SampleDay)
###

#Bee Count Observations
BC <- subset(BeeCounts, Bees > -1)
BCs <- glm(Bees ~ Treatment + Date + Area + Cumulative, family=poisson(link="log"), data= BC, method="brglmFit")
Anova(BCs)
lsm<-lsmeans (BCs, list( ~ Treatment))
cld(lsm)

#AP23 treatment group analysis by date
BC1 <- subset(BC, Treatment == "treatment")
BCs1 <- mixed(Bees ~ SampleDay + (1|Date) + (1|Area) + (1|Colony) + (1|Cumulative), family=poisson(link="log"), data= BC1, method="LRT")

lsm<-lsmeans (BCs1, list( ~ SampleDay))
cld(lsm)
###

#Dye Visualized in Bee Guts
DV <- subset(DyeVisuals, Visual1>-1)
DVs <- glm(Visual1 ~ Treatment + Date, family=binomial(link="logit"), data= DV, method=brglmFit)
Anova(DVs)
lsm<-lsmeans (DVs, list( ~ Treatment))
cld(lsm)

#AP23 treatment group analysis by date
DV1 <- subset(DV, Treatment == "treatment")
DVs1 <- mixed(Visual1 ~ SampleDay + (1|Date) + (1|Hive), family=binomial(link="logit"), data= DV1, method="LRT")
DVs1

lsm<-lsmeans (DVs1, list( ~ SampleDay))
cld(lsm)
###
