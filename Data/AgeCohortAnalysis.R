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
library(ggplot2)
library(tidyverse)
library(glmmTMB)
library(MuMIn)
library(DHARMa)
library(emmeans)
library(scales) 
##getting manual color codes from viridis package               
##
setwd("C:/Users/codyp/Desktop/WorkBackup/Research/HBREL/AgeCohort/")

#Import files
BeeCounts <- read.csv("BeeCounts.csv")
Weather1 <- read.csv("2020.csv")
##

#Data organization
BeeCounts$Colony <- as.factor(BeeCounts$Colony)
names(BeeCounts)[1] <- "Date"
names(DyeVisuals)[1] <- "Treatment"
BeeCounts$Age <- BeeCounts$SampleDay +1
#Weather data will be summarized from sunrise to the average of our sampling time, 7:00 - 11:00
Weather2 <- subset(Weather1, Time == c("7:00:00"))
Weather3 <- subset(Weather1, Time == c("8:00:00"))
Weather4 <- subset(Weather1, Time == c("9:00:00"))
Weather5 <- subset(Weather1, Time == c("10:00:00"))
Weather6 <- subset(Weather1, Time == c("11:00:00"))
Weather <- rbind(Weather2,Weather3, Weather4,Weather5,Weather6)
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
#Scale all weather data
SoilTemp$SoilTemp <- scale(SoilTemp$SoilTemp)
Temp60cm["Temp60cm"] <- scale(Temp60cm["Temp60cm"])
Temp2m["Temp2m"] <- scale(Temp2m["Temp2m"])
Temp10m["Temp10m"] <- scale(Temp10m["Temp10m"])
Humidity["Humidity"] <- scale(Humidity["Humidity"])
DPTemp["DPTemp"] <- scale(DPTemp["DPTemp"])
Rainfall["Rainfall"] <- scale(Rainfall["Rainfall"])
WindSp["WindSp"] <- scale(WindSp["WindSp"])
WindDir["WindDir"] <- scale(WindDir["WindDir"])
Radiation["Radiation"] <- scale(Radiation["Radiation"])
#Merge weather data to bee counts
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

#Remove control treatments (we want to keep the main treatment and the undyed treatment)
BeeCounts <- subset(BeeCounts, Treatment != "CControl")
BeeCounts <- subset(BeeCounts, Treatment != "PControl")
BeeCounts <- na.omit(BeeCounts)

#plot correlations between variables
count_cor2 <- BeeCounts %>% dplyr::select(-Date, -Treatment, -Area, -Notes, -Colony, -Hive)
count_cor3 <- cor(as.matrix(count_cor2)) ##make correlation matrix
str(count_cor2)
corrplot(count_cor3) ## plots strength of correlation as color-coded circles

###

#####AICc######
##model building
##start with most complex full model

cp1 <- glmmTMB(Bees ~ Age + I(Age**2) + offset(log(BeesTotal)) +
                Temp60cm + Rainfall + Radiation + Humidity +
                ar1(0 + as.factor(Age)|Colony), data=BeeCounts, family="nbinom2") ##quadratic term and AR1 random effect

##AR1 correlation structure accounts for correlation of points in time closer together

cp2 <- glmmTMB(Bees ~ Age + I(Age**2)  + offset(log(BeesTotal)) + 
                Temp60cm + Rainfall + Radiation +
                ar1(0 + as.factor(Age)|Colony), data=BeeCounts, family="nbinom2") ##pairwise by weather variables

cp3 <- glmmTMB(Bees ~ Age + I(Age**2) + offset(log(BeesTotal)) + 
                Temp60cm + Rainfall + Humidity +
                ar1(0 + as.factor(Age)|Colony), data=BeeCounts, family="nbinom2")

cp4 <- glmmTMB(Bees ~ Age + I(Age**2) + offset(log(BeesTotal)) + 
                Rainfall + Radiation + Humidity +
                ar1(0 + as.factor(Age)|Colony), data=BeeCounts, family="poisson")


cp5 <- glmmTMB(Bees ~ Age + I(Age**2) + offset(log(BeesTotal)) + 
                Temp60cm + Radiation + Humidity +
                ar1(0 + as.factor(Age)|Colony), data=BeeCounts, family="poisson")

cp6 <- glmmTMB(Bees ~ Age + I(Age**2) + offset(log(BeesTotal)) + 
                Temp60cm + Rainfall +
                ar1(0 + as.factor(Age)|Colony), data=BeeCounts, family="nbinom2")

cp7 <- glmmTMB(Bees ~ Age + I(Age**2) + offset(log(BeesTotal)) + 
                Temp60cm + Humidity +
                ar1(0 + as.factor(Age)|Colony), data=BeeCounts, family="nbinom2")

cp8 <- glmmTMB(Bees ~ Age + I(Age**2) + offset(log(BeesTotal)) + 
                Temp60cm + Radiation +
                ar1(0 + as.factor(Age)|Colony), data=BeeCounts, family="nbinom2")

cp9 <- glmmTMB(Bees ~ Age + I(Age**2) + offset(log(BeesTotal)) + 
                Rainfall + Radiation +
                ar1(0 + as.factor(Age)|Colony), data=BeeCounts, family="nbinom2")

cp10 <- glmmTMB(Bees ~ Age + I(Age**2) + offset(log(BeesTotal)) + 
                 Rainfall + Humidity +
                 ar1(0 + as.factor(Age)|Colony), data=BeeCounts, family="nbinom2")

cp11 <- glmmTMB(Bees ~ Age + I(Age**2) + offset(log(BeesTotal)) + 
                 Radiation + Humidity +
                 ar1(0 + as.factor(Age)|Colony), data=BeeCounts, family="nbinom2")

cp12 <- glmmTMB(Bees ~ Age + I(Age**2) + offset(log(BeesTotal)) + 
                 Humidity +
                 ar1(0 + as.factor(Age)|Colony), data=BeeCounts, family="nbinom2")

cp13 <- glmmTMB(Bees ~ Age + I(Age**2) + offset(log(BeesTotal)) + 
                 Radiation +
                 ar1(0 + as.factor(Age)|Colony), data=BeeCounts, family="nbinom2")

cp14 <- glmmTMB(Bees ~ Age + I(Age**2) + offset(log(BeesTotal)) + 
                 Temp60cm +
                 ar1(0 + as.factor(Age)|Colony), data=BeeCounts, family="nbinom2")

cp15 <- glmmTMB(Bees ~ Age + I(Age**2) + offset(log(BeesTotal)) + 
                 Rainfall +
                 ar1(0 + as.factor(Age)|Colony), data=BeeCounts, family="nbinom2")

cp16 <- glmmTMB(Bees ~ Age + I(Age**2) + offset(log(BeesTotal)) + 
                 ar1(0 + as.factor(Age)|Colony), data=BeeCounts, family="nbinom2") ##no weather variables

cp17 <- glmmTMB(Bees ~ Age + I(Age**2) + 
                 ar1(0 + as.factor(Age)|Colony), data=BeeCounts, family="nbinom2") ##no mortality offset

cp18 <- glmmTMB(Bees ~ Age + 
                 ar1(0 + as.factor(Age)|Colony), data=BeeCounts, family="poisson")##no quadratic term

cp19 <- glmmTMB(Bees ~ Age 
                , data=BeeCounts, family="nbinom2") ##no AR1 correlation structure

cp20 <- glmmTMB(Bees ~ Temp60cm + Rainfall + Radiation + Humidity + offset(log(BeesTotal)) + 
                 ar1(0 + as.factor(Age)|Colony), data=BeeCounts) ##no age, only weather variables

cp21 <- glmmTMB(Bees ~ 1, data=BeeCounts, family="nbinom2") ##null


##generate AICc scores, weights, and delta
Model.sel <- model.sel(cp1, cp2, cp3, cp4, cp5, cp6, cp7, cp8, cp9, cp10, cp11, cp12, cp13, cp14, cp15, cp16, cp17, cp18, cp19, cp20, cp21)
write.csv(Model.sel, "Modelsel.csv")
##8 of 12 models within 6 AICc, will need to average
##all 8 have quadratic term and AR1 random effect

##average model
cpollen_models <- model.sel(cp1, cp2, cp3, cp4, cp5, cp6, cp7, cp8, cp9, cp10, cp11, cp12, cp13, cp14, cp15, cp16, cp17, cp18, cp19, cp20, cp21)

model.avg(cpollen_models, subset = delta < 6)

summary(model.avg(cpollen_models, subset = delta < 6)) ##conditional avg does not include zeros for absent terms
##in the averaged model, the effect sizes were small, so the environmental variables contributed very little
##justify use simplest (model 16) without env variables....most parsimonious
##report estimates and SE for completeness
##Could potentially use likelihood ratio test to provide p value working from most complex model to simplest. Not ideal, but if can't publish on AIC alone

ggplot(data=BeeCounts, aes(x=Humidity, y=Rainfall))+geom_point()+geom_smooth(method="glm", 
                                                                        formula = y~x+I(x^2))

##

#assess residuals of selected model

plot(simulateResiduals(cp1)) 
hist(simulateResiduals(cp1))

#random effect components

summary(cp1)

###significance of model terms

Anova(cp1)

#make a table of emmeans from our selected model for manuscript
A2 <- summary(emmeans(cp1, ~1, type='response', at=list(Age=2))) 
A4 <-summary(emmeans(cp1, ~1, type='response', at=list(Age=4))) 
A6 <-summary(emmeans(cp1, ~1, type='response', at=list(Age=6))) 
A8 <-summary(emmeans(cp1, ~1, type='response', at=list(Age=8))) 
A10 <-summary(emmeans(cp1, ~1, type='response', at=list(Age=10))) 
A12 <-summary(emmeans(cp1, ~1, type='response', at=list(Age=12))) 
A14 <-summary(emmeans(cp1, ~1, type='response', at=list(Age=14))) 
A16 <-summary(emmeans(cp1, ~1, type='response', at=list(Age=16))) 
A18 <-summary(emmeans(cp1, ~1, type='response', at=list(Age=18))) 
A20 <-summary(emmeans(cp1, ~1, type='response', at=list(Age=20))) 
A22 <-summary(emmeans(cp1, ~1, type='response', at=list(Age=22))) 
A24 <-summary(emmeans(cp1, ~1, type='response', at=list(Age=24))) 
A26 <-summary(emmeans(cp1, ~1, type='response', at=list(Age=26))) 
A28 <-summary(emmeans(cp1, ~1, type='response', at=list(Age=28))) 

means <- data.frame(Response = c(A2$response,A4$response,A6$response,A8$response,A10$response,A12$response,A14$response,
                        A16$response,A18$response,A20$response,A22$response,A24$response,A26$response,A28$response), 
           Lower = c(A2$lower.CL, A4$lower.CL, A6$lower.CL, A8$lower.CL,A10$lower.CL, A12$lower.CL, A14$lower.CL, A16$lower.CL,
                     A18$lower.CL, A20$lower.CL, A22$lower.CL, A24$lower.CL,A26$lower.CL, A28$lower.CL), 
           Upper = c(A2$upper.CL, A4$upper.CL, A6$upper.CL, A8$upper.CL,
             A10$upper.CL, A12$upper.CL, A14$upper.CL, A16$upper.CL, A18$upper.CL, A20$upper.CL, A22$upper.CL, A24$upper.CL,
             A26$upper.CL, A28$upper.CL))

write.csv(means, "means.csv")
 



####Plots#####
str(full_pollen2)
full_pollen2$Hive_number <- as.factor(full_pollen2$Hive_number)


show_col(viridis_pal()(9)) ##use from bottom right corner, right to left, bottom to top


##testing the proportion of bees/bees total. More accurate to model with offset
##this just for comparison with emmeans later on. Very similiar, meaning weather variables have no effect, just offset
BeeCounts$TotalRound <- round(BeeCounts$BeesTotal, digits=0)
BeeCounts$prop <- (BeeCounts$Bees/BeeCounts$TotalRound)
summary(BeeCounts$prop)

ggplot(BeeCounts, aes(x=Age, y=prop))+
  geom_smooth(BeeCounts, mapping=aes(x=Age, y=prop, linetype=Hive, color=Hive), method="glm", method.args=list(family="binomial"), 
              formula = y~x+I(x^2), se=FALSE, linetype="dashed", size=1) + theme_bw() + 
  geom_smooth(color="black", method="glm", method.args=list(family="binomial"), se=FALSE, formula = y ~ x+I(x^2), size=2)+ ##can change order of geom_smooths for whatver line(s) you want on top
  theme(strip.background = element_blank(), panel.grid.minor= element_blank(), panel.grid.major= element_blank(),
        text = element_text(size=16))+
  labs(x="Worker age (day)", y="Marked bees observed on patty")+
  scale_color_viridis(discrete = T, direction = -1, "Colony") +#, labels = c("A","B", "C", "D", "E", "F"))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 15)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))


####For averaged plot, panel A
##extract emmeans for the averaged model and each hive

call1 <- as.data.frame(emmeans(cp5, ~1, type='response', at=list(Age=1))) 
call2 <- as.data.frame(emmeans(cp5, ~1, type='response', at=list(Age=2))) 
call3 <- as.data.frame(emmeans(cp5, ~1, type='response', at=list(Age=3))) 
call4 <- as.data.frame(emmeans(cp5, ~1, type='response', at=list(Age=4))) 
call5 <- as.data.frame(emmeans(cp5, ~1, type='response', at=list(Age=5)))
call6 <- as.data.frame(emmeans(cp5, ~1, type='response', at=list(Age=6)))
call7 <- as.data.frame(emmeans(cp5, ~1, type='response', at=list(Age=7)))
call8 <- as.data.frame(emmeans(cp5, ~1, type='response', at=list(Age=8)))
call9 <- as.data.frame(emmeans(cp5, ~1, type='response', at=list(Age=9)))
call10 <-as.data.frame(emmeans(cp5, ~1, type='response', at=list(Age=10)))
call11 <- as.data.frame(emmeans(cp5, ~1, type='response', at=list(Age=11)))
call12 <-as.data.frame(emmeans(cp5, ~1, type='response', at=list(Age=12)))
call13 <- as.data.frame(emmeans(cp5, ~1, type='response', at=list(Age=13)))
call14 <-as.data.frame(emmeans(cp5, ~1, type='response', at=list(Age=14)))
call15 <- as.data.frame(emmeans(cp5, ~1, type='response', at=list(Age=15)))
call16 <-as.data.frame(emmeans(cp5, ~1, type='response', at=list(Age=16)))
call17 <- as.data.frame(emmeans(cp5, ~1, type='response', at=list(Age=17))) 
call18 <-as.data.frame(emmeans(cp5, ~1, type='response', at=list(Age=18))) 
call19 <-as.data.frame(emmeans(cp5, ~1, type='response', at=list(Age=19))) 
call20 <-as.data.frame(emmeans(cp5, ~1, type='response', at=list(Age=20))) 
call21 <-as.data.frame(emmeans(cp5, ~1, type='response', at=list(Age=21))) 
call22 <-as.data.frame(emmeans(cp5, ~1, type='response', at=list(Age=22)))
call23 <- as.data.frame(emmeans(cp5, ~1, type='response', at=list(Age=23))) 
call24 <-as.data.frame(emmeans(cp5, ~1, type='response', at=list(Age=24)))
call25 <- as.data.frame(emmeans(cp5, ~1, type='response', at=list(Age=25))) 
call26 <-as.data.frame(emmeans(cp5, ~1, type='response', at=list(Age=26)))
call27 <- as.data.frame(emmeans(cp5, ~1, type='response', at=list(Age=27))) 
call28 <-as.data.frame(emmeans(cp5, ~1, type='response', at=list(Age=28)))
call29 <- as.data.frame(emmeans(cp5, ~1, type='response', at=list(Age=29))) 

##binding all rows
call <- bind_rows(call1, call2, call3, call4, call5, call6, call7, call8, call9, call10, call11,
                 call12, call13, call14, call15, call16, call17, call18, call19, call20, call21,
                 call22, call23, call24, call25, call26, call27, call28, call29)

##add age and description for binding later
call$age <- c(1:29)
call$hive <- "All"


##two separate panels
##panel A, using emmeans above
ggplot()+
  geom_ribbon(data=call, aes(x= age, ymin=lower.CL, ymax=upper.CL), alpha=0.1)+
  geom_line(data=call, aes(x=age, y=rate), color="black", size=2.5)+
  theme_bw()+ theme(strip.background = element_blank(), panel.grid.minor= element_blank(), panel.grid.major= element_blank(),
                    text = element_text(size=16))+
  labs(x="Worker age (day)", y="# Marked workers on pollen patties")+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))

##pannel b
##using the predict function to calculate the model's equation for each point

BeeCounts$pred <- exp(predict(cp1))/BeeCounts$BeesTotal

ggplot()+
  geom_smooth(data=BeeCounts, aes(x=Age, y=pred, color=Colony, linetype= Colony), size=1.8, linetype="dashed", 
              method="glm", method.args=list(family="nbinom2"), se=FALSE, formula = y ~ x+I(x^2))+
  theme_bw()+ theme(strip.background = element_blank(), panel.grid.minor= element_blank(), panel.grid.major= element_blank(),
                    text = element_text(size=16))+
  labs(x="Worker age (day)", y="Proportion of marked workers on pollen patties")+
  scale_color_manual(values=c("#FDE725FF", "#AADC32FF","#5DC863FF","#27AD81FF","#21908CFF", "#2C728EFF", "#3B528BFF", '#472D7BFF', '#440154FF'), "Colony", labels = c("A","B", "C", "D", "E", "F", "G", "H", "I"))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))



####Weather plots####

##unscaled weather data for plots

BeeCounts_unscaled <- read.csv("BeeCounts.csv")
BeeCounts_unscaled$Colony <- as.factor(BeeCounts_unscaled$Colony)
#BeeCounts$SampleDay <- as.factor(BeeCounts$SampleDay)
names(BeeCounts_unscaled)[1] <- "Date"
BeeCounts_unscaled$Age <- BeeCounts_unscaled$SampleDay +1
BeeCounts_unscaled <- subset(BeeCounts_unscaled, Treatment != "CControl")
BeeCounts_unscaled <- subset(BeeCounts_unscaled, Treatment != "PControl")
BeeCounts_unscaled <- na.omit(BeeCounts_unscaled)

##unscaled weather data for plots
Weather1 <- read.csv("2020.csv")
##
Weather2 <- subset(Weather1, Time == c("7:00:00"))
Weather3 <- subset(Weather1, Time == c("8:00:00"))
Weather4 <- subset(Weather1, Time == c("9:00:00"))
Weather5 <- subset(Weather1, Time == c("10:00:00"))
Weather6 <- subset(Weather1, Time == c("11:00:00"))
Weather_unscaled <- rbind(Weather2,Weather3, Weather4,Weather5,Weather6)
#average weather variable by date over 4 hours bees were sampled.
SoilTemp <- aggregate(SoilTemp~Date, data=Weather_unscaled, FUN=mean)
Temp60cm <- aggregate(Temp60cm~Date, data=Weather_unscaled, FUN=mean)
Temp2m <- aggregate(Temp2m~Date, data=Weather_unscaled, FUN=mean)
Temp10m <- aggregate(Tem10m~Date, data=Weather_unscaled, FUN=mean)
Humidity <- aggregate(Humidity~Date, data=Weather_unscaled, FUN=mean)
DPTemp <- aggregate(DPTemp~Date, data=Weather_unscaled, FUN=mean)
Rainfall <- aggregate(Rainfall~Date, data=Weather_unscaled, FUN=mean)
WindSp <- aggregate(WindSp~Date, data=Weather_unscaled, FUN=mean)
WindDir <- aggregate(WindDir~Date, data=Weather_unscaled, FUN=mean)
Radiation <- aggregate(Radiation~Date, data=Weather_unscaled, FUN=mean)


BeeCounts_unscaled <- merge(BeeCounts_unscaled,SoilTemp, by="Date")
BeeCounts_unscaled <- merge(BeeCounts_unscaled,Temp60cm, by="Date")
BeeCounts_unscaled <- merge(BeeCounts_unscaled,Temp2m, by="Date")
BeeCounts_unscaled <- merge(BeeCounts_unscaled,Temp10m, by="Date")
BeeCounts_unscaled <- merge(BeeCounts_unscaled,Humidity, by="Date")
BeeCounts_unscaled <- merge(BeeCounts_unscaled,DPTemp, by="Date")
BeeCounts_unscaled <- merge(BeeCounts_unscaled,Rainfall, by="Date")
BeeCounts_unscaled <- merge(BeeCounts_unscaled,WindSp, by="Date")
BeeCounts_unscaled <- merge(BeeCounts_unscaled,WindDir, by="Date")
BeeCounts_unscaled <- merge(BeeCounts_unscaled,Radiation, by="Date")

BeeCounts_unscaled$Rainfall <- BeeCounts_unscaled$Rainfall*2.54 #convert inch to cm


##extract residuals from model without one weather variable and plot residuals against said weather variable

##temp
BeeCounts_unscaled$residT <-residuals(cp4)

png("Temperature.png", width=8, height=6, units="in", res=400)
ggplot(BeeCounts_unscaled, aes(x=Temp60cm, y=residT))+
  theme_bw() + geom_jitter(BeeCounts_unscaled, mapping=aes(x=Temp60cm, y=residT, color=Treatment), height=0.1, alpha=0.8, size=1.75)+
  geom_smooth(color="black", method="glm", size=2)+ 
  theme(strip.background = element_blank(), panel.grid.minor= element_blank(), panel.grid.major= element_blank(),
        text = element_text(size=16), legend.position="none")+ ylim(-1.1,1.1)+
  labs(x=expression(paste('Average Daily Temperature (',~degree,'C)',sep='')), y="Model Residuals")+
  scale_color_viridis(discrete = T, direction = -1, begin= .85, end=.85, alpha=0.08)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) 
dev.off()

##rain
BeeCounts_unscaled$residR <-residuals(cp5)

png("Rainfall.png", width=8, height=6, units="in", res=400)
ggplot(BeeCounts_unscaled, aes(x=Rainfall, y=residR))+
  theme_bw() + geom_jitter(BeeCounts_unscaled, mapping=aes(x=Rainfall, y=residR, color=Treatment), height=0.1, alpha=0.8, size=1.75)+
  geom_smooth(color="black", method="glm", size=2)+ 
  theme(strip.background = element_blank(), panel.grid.minor= element_blank(), panel.grid.major= element_blank(),
        text = element_text(size=16), legend.position="none")+ ylim(-1.1,1.1)+
  labs(x="Average Daily Rainfall (cm)", y="Model Residuals")+
  scale_color_viridis(discrete = T, direction = -1, begin= .55, end=.55, alpha=0.08)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) 
dev.off()

##radiation
BeeCounts_unscaled$residRd <-residuals(cp3)

png("Radiation.png", width=8, height=6, units="in", res=400)
ggplot(BeeCounts_unscaled, aes(x=Radiation, y=residRd))+
  theme_bw() + geom_jitter(BeeCounts_unscaled, mapping=aes(x=Radiation, y=residRd, color=Treatment), height=0.1, alpha=0.8, size=1.75)+
  geom_smooth(color="black", method="glm", size=2)+ 
  theme(strip.background = element_blank(), panel.grid.minor= element_blank(), panel.grid.major= element_blank(),
        text = element_text(size=16), legend.position="none")+ ylim(-1.1,1.1)+
  labs(x=expression(Average~Daily~Radiation~(W/m^{"2"})), y="Model Residuals")+ ##use expression for superscript
  scale_color_viridis(discrete = T, direction = -1, begin= .25, end=.25, alpha=0.08)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) 
dev.off()

#humidity
BeeCounts_unscaled$residH <-residuals(cp2)

png("Humidity.png", width=8, height=6, units="in", res=400)
ggplot(BeeCounts_unscaled, aes(x=Humidity, y=residH))+
  theme_bw() + geom_jitter(BeeCounts_unscaled, mapping=aes(x=Humidity, y=residH, color=Treatment), height=0.1, alpha=0.8, size=1.75)+
  geom_smooth(color="black", method="glm", size=2)+ 
  theme(strip.background = element_blank(), panel.grid.minor= element_blank(), panel.grid.major= element_blank(),
        text = element_text(size=16), legend.position="none")+ ylim(-1.1,1.1)+
  labs(x="Average Daily Humidity (%)", y="Model Residuals")+
  scale_color_viridis(discrete = T, direction = -1, begin= 0, end=0, alpha=0.08)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
dev.off()

