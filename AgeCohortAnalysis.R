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
setwd("C:/Users/codyp/Desktop/WorkBackup/Research/HBREL/Emily/AgeCohort/")
#Import files
BeeCounts <- read.csv("BeeCounts.csv")
DyeVisuals <- read.csv("DyeVisuals.csv")
Weather1 <- read.csv("2020.csv")
##

#Data organization
BeeCounts$Colony <- as.factor(BeeCounts$Colony)
BeeCounts$SampleDay <- as.factor(BeeCounts$SampleDay)
names(BeeCounts)[1] <- "Date"
names(DyeVisuals)[1] <- "Treatment"
DyeVisuals$Hive <- as.factor(DyeVisuals$Hive)
DyeVisuals$SampleDay <- as.factor(DyeVisuals$SampleDay)
#Time coded this way makes it difficult to cut out unneeded times
Weather2 <- subset(Weather1, Time == c("9:00:00"))
Weather3 <- subset(Weather1, Time == c("10:00:00"))
Weather4 <- subset(Weather1, Time == c("11:00:00"))
Weather5 <- subset(Weather1, Time == c("12:00:00"))
Weather <- rbind(Weather2,Weather3, Weather4,Weather5)
#average weather variable by date over 4 hours bees were sampled.
SoilTemp <- aggregate(SoilTemp~Date, data=Weather, FUN=mean)
Temp60cm <- aggregate(Temp60cm~Date, data=Weather, FUN=mean)
Temp2m <- aggregate(Temp2m~Date, data=Weather, FUN=mean)
Temp10m <- aggregate(Tem10m~Date, data=Weather, FUN=mean)
Humidity <- aggregate(Humidity~Date, data=Weather, FUN=mean)
DPTemp <- aggregate(DPTemp~Date, data=Weather, FUN=mean)
Rainfall <- aggregate(Rainfall~Date, data=Weather, FUN=mean)
WindSp <- aggregate(WindSp~Date, data=Weather, FUN=mean)
WindDir <- aggregate(WindDir~Date, data=Weather, FUN=mean)
Radiation <- aggregate(Radiation~Date, data=Weather, FUN=mean)
BeeCounts <- merge(BeeCounts,SoilTemp, by="Date")
BeeCounts <- merge(BeeCounts,Temp60cm, by="Date")
BeeCounts <- merge(BeeCounts,Temp2m, by="Date")
BeeCounts <- merge(BeeCounts,Temp10m, by="Date")
BeeCounts <- merge(BeeCounts,Humidity, by="Date")
BeeCounts <- merge(BeeCounts,DPTemp, by="Date")
BeeCounts <- merge(BeeCounts,Rainfall, by="Date")
BeeCounts <- merge(BeeCounts,WindSp, by="Date")
BeeCounts <- merge(BeeCounts,WindDir, by="Date")
BeeCounts <- merge(BeeCounts,Radiation, by="Date")

DyeVisuals <- merge(DyeVisuals, Temp2m, by="Date")
DyeVisuals <- merge(DyeVisuals, Humidity, by="Date")
DyeVisuals <- merge(DyeVisuals, SoilTemp, by="Date")
DyeVisuals <- merge(DyeVisuals, Temp60cm, by="Date")
DyeVisuals <- merge(DyeVisuals, Temp10m, by="Date")
DyeVisuals <- merge(DyeVisuals, DPTemp, by="Date")
DyeVisuals <- merge(DyeVisuals, Rainfall, by="Date")
DyeVisuals <- merge(DyeVisuals, WindSp, by="Date")
DyeVisuals <- merge(DyeVisuals, WindDir, by="Date")
DyeVisuals <- merge(DyeVisuals, Radiation, by="Date")

###

#Bee Count Observations - PCA

#build PCA dataframe
PCABees <- data.frame(Date=BeeCounts$Date, SampleDay=BeeCounts$SampleDay, Colony=BeeCounts$Colony, Treatment=BeeCounts$Treatment, Area=BeeCounts$Area, Cumulative=BeeCounts$Cumulative, Bees=BeeCounts$Bees,
                      Temp2m=BeeCounts$Temp2m,Temp60cm=BeeCounts$Temp60cm, Radiation=BeeCounts$Radiation,Temp10m=BeeCounts$Tem10m, Humidity=BeeCounts$Humidity,DPTemp=BeeCounts$DPTemp, Rainfall=BeeCounts$Rainfall, WindDir=BeeCounts$WindDir)
#NA omit
PCABees <- na.omit(PCABees)

#Run PCA
DRBees <- prcomp(PCABees[8:15], center = TRUE, scale = TRUE)
#Those are my PCs


#Build new dataframe which includes PC 1 and 2
PCAscoresBees <- data.frame(Date = PCABees$Date, SampleDay=PCABees$SampleDay, Colony=PCABees$Colony, Treatment=PCABees$Treatment, Area=PCABees$Area, Cumulative=PCABees$Cumulative, Bees=PCABees$Bees, PC1 = DRBees$x[,1], PC2 = DRBees$x[,2])

#Calculate proportion of variation explained
PCA1POVBees <- DRBees$sdev^2/sum(DRBees$sdev^2)
###

BCs <- glm(Bees ~ Cumulative + Treatment + PC1 + PC2 + Area + Colony, family=poisson(link="log"), data= PCAscoresBees,method="brglmFit")
anova(BCs,test="Chisq")

BC1 <- subset(PCAscoresBees, Treatment != "CControl")

BCs <- mixed(Bees ~ Cumulative + Treatment*SampleDay + PC1 + PC2 + (1|Date) + (1|Area) + (1|Colony), family=poisson(link="log"), data= BC1,method="LRT")
BCs

BCs2 <- mixed(Bees ~ Cumulative + Treatment + PC1 + PC2 + (1|Date) + (1|Area) + (1|Colony), family=poisson(link="log"), data= BC1,method="LRT")
BCs2


#AP23 treatment group analysis by date
BC2 <- subset(BC1, Treatment == "treatment")
BCs1 <- mixed(Bees ~ Cumulative + SampleDay + PC1 + PC2+ (1|Date) + (1|Area) + (1|Colony), family=poisson(link="log"), data= BC2, method="LRT")
summary(BCs1)

lsm<-lsmeans (BCs1, list( ~ SampleDay))
cld(lsm)
###

#Dye Visualized in Bee Guts

#Same method as above
PCADye <- data.frame(SampleDay=DyeVisuals$SampleDay, Date=DyeVisuals$Date, Hive=DyeVisuals$Hive, Treatment=DyeVisuals$Treatment, Visual=DyeVisuals$Visual1,SoilTemp2=DyeVisuals$SoilTemp,
                     Temp2m2=DyeVisuals$Temp2m,Temp60cm2= DyeVisuals$Temp60cm,Temp10m2=DyeVisuals$Tem10m, Humidity2=DyeVisuals$Humidity,DPTemp2=DyeVisuals$DPTemp, Rainfall2= DyeVisuals$Rainfall, WindSp2=DyeVisuals$WindSp)

PCADye <- na.omit(PCADye)
DRDye <- prcomp(PCADye[6:13], center = TRUE, scale = TRUE)
PCAscoresDye <- data.frame(SampleDay = PCADye$SampleDay, Date=PCADye$Date, Hive=PCADye$Hive, Treatment=PCADye$Treatment, Visual=PCADye$Visual, PC1 = DRDye$x[,1], PC2 = DRDye$x[,2])
PCA1POVDye <- DRDye$sdev^2/sum(DRDye$sdev^2)

#Run overall model for significance
DVs <- glm(Visual ~ Treatment + PC1 + PC2 + Hive, family=binomial(link="logit"), data= PCAscoresDye, method=brglmFit)
anova(DVs, test="Chisq")
lsm<-lsmeans (DVs, list( ~ Treatment))
cld(lsm)

#AP23 treatment group analysis by date
DV <- subset(PCAscoresDye, Treatment == "treatment")
DVs1 <- mixed(Visual ~ SampleDay + PC1 + PC2 + (1|Date) + (1|Hive), family=binomial(link="logit"), data= DV, method="LRT")
DVs1

lsm<-lsmeans (DVs1, list( ~ SampleDay))
cld(lsm)
###
