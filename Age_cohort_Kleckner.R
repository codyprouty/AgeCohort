#####Kleckner#####
##Experiments in ecology and agriculture - final project##
##Noordyke age cohort study##

##load in packages
library(tidyverse)
library(GGally)
library(glmmTMB)
library(DHARMa)
library(car)
library(corrplot)
library(MuMIn)
library(boot)

#####Bee-gut visuals####
pollen <- read.csv("Age_cohort_visuals.csv")
str(pollen)

##remove unnecessary columns and convert Y/N to 1/0
pollen1 <- pollen %>% dplyr::select(-X, -X.1, -X.2, -Date_Formatting) %>% ##used dplyr:: bc of problems with masking with the MASS package
  mutate(Visual1 = if_else(Visual1 == "Y", 1L, 0L),
         Visual2 = if_else(Visual2 == "Y", 1L , 0L))
str(pollen1)

##sampling day is one day less than worker age. 
##create column for worker age for easier interpretation of results
pollen1$age <- pollen1$Sample_day+1

##need to remove visual gradings that disagree
pollen1$visual_avg <- (pollen1$Visual1+pollen1$Visual2)/2

hist(pollen1$visual_avg) ##very few 0.5 values out of the whole data set, meaning very few instances the two visual graders disagreed

pollen2 <- filter(pollen1, visual_avg !="0.5") ##removed 0.5 values and NAs (marked bees no longer present at the end of the sampling period)

hist(pollen2$visual_avg) ##it worked!


##quick data exploration
ggplot(pollen2, aes(x=Sample_day, y=visual_avg, color=Treatment))+geom_jitter()
##Looks like all the controls work. No obvious pattern with the treatment

ggplot(pollen2, aes(x=Sample_day, y=visual_avg, color=Treatment))+geom_jitter(height=0.1, alpha=0.1)+
  facet_wrap(~Treatment)+geom_smooth(color="black", method="glm", method.args=list(family="binomial"), 
                                     formula = y ~ x+I(x^2)) ##quadratic term seems helpful, will include in models

##!##!##!##!##!##!##!##!##!##
#####environmental data#####

##True date is calendar date
##Date and period are for summarizing the previous 24 hours from sampling (12pm-11am)
pollen_env2 <- read.csv("Age_cohort_env2.csv") ##From UF IFAS FAWN Alachua Station 260, trimmed for sampling days in October 2020
str(pollen_env2)

#average the previous 24 hours from sampling using "Date"
pollen_env_avg <- pollen_env2 %>% dplyr::select(-Station.ID, -Station.Name, -Time, -True_Date) %>% 
  group_by(Period, Date) %>% summarise_each(funs(mean)) %>% rename_with( ~ paste("Avg", .x, sep = "_")) %>% 
  rename("Date"="Avg_Date") %>% filter(Avg_Period < 31)
head(pollen_env_avg)

##change rainfall from inches to cm
pollen_env_avg$Avg_Rainfall <- pollen_env_avg$Avg_Rainfall*2.54

##scale env variables BEFORE merging
pollen_env_avg[c("Avg_Temp60cm", "Avg_Rainfall", "Avg_Humidity", "Avg_Radiation")] <- scale(pollen_env_avg[c("Avg_Temp60cm", "Avg_Rainfall", "Avg_Humidity", "Avg_Radiation")])
##center is zero, a 1.00 unit increase represents one standard deviation from the original scale
##easier to interpret in summary later on

##cbind with visuals using dates
full_pollen <- left_join(pollen2, pollen_env_avg, by = "Date") 
head(full_pollen) ##wooho

##remove controls for modeling
full_pollen2 <- filter(full_pollen, Treatment =="treatment")

head(full_pollen2)
tail(full_pollen2)

##check correlations
pollen_cor2 <- full_pollen2 %>% dplyr::select(-Date, -Treatment)
pollen_cor3 <- cor(as.matrix(pollen_cor2)) ##make correlation matrix

corrplot(pollen_cor3, method = "circle", type = "upper") ## plots strength of correlation as color-coded circles
###all the different temperature variables are very correlated with each other (no surprise)
###soil temperature is negatively correlated with sample day, but not biologically relevant
###but no strong trends with the visual averages

##scale environmental variables of interest
#full_pollen2[c("Avg_Temp60cm", "Avg_Rainfall", "Avg_Humidity", "Avg_Radiation")] <- scale(full_pollen2[c("Avg_Temp60cm", "Avg_Rainfall", "Avg_Humidity", "Avg_Radiation")])
##center is zero, a 1.00 unit increase represents one standard deviation from the original scale
##easier to interpret in summary later on
head(full_pollen2)


#####AICc######
##model building
##start with most complex full model

p1 <- glmmTMB(visual_avg ~ age + I(age**2) + 
                Avg_Temp60cm + Avg_Rainfall + Avg_Radiation + Avg_Humidity +
                ar1(0 + as.factor(age)|Hive_number), data=full_pollen2, family="binomial") ##quadratic term and AR1 random effect
##AR1 correlation structure accounts for correlation of points in time closer togeher
##look for coffiecient later in summary

p2 <- glmmTMB(visual_avg ~ age + I(age**2) + 
                Avg_Temp60cm + Avg_Rainfall + Avg_Radiation + #no humidity
                ar1(0 + as.factor(age)|Hive_number), data=full_pollen2, family="binomial") ##pairwise by weather variables

p3 <- glmmTMB(visual_avg ~ age + I(age**2) + 
                Avg_Temp60cm + Avg_Rainfall +  Avg_Humidity + #no radiation
                ar1(0 + as.factor(age)|Hive_number), data=full_pollen2, family="binomial")

p4 <- glmmTMB(visual_avg ~ age + I(age**2) + 
                Avg_Rainfall  + Avg_Radiation + Avg_Humidity + #no temp
                ar1(0 + as.factor(age)|Hive_number), data=full_pollen2, family="binomial")

p5 <- glmmTMB(visual_avg ~ age + I(age**2) + 
                Avg_Temp60cm  + Avg_Radiation + Avg_Humidity + #no rainfall
                ar1(0 + as.factor(age)|Hive_number), data=full_pollen2, family="binomial")

p6 <- glmmTMB(visual_avg ~ age + I(age**2) + 
                Avg_Temp60cm + Avg_Rainfall + 
                ar1(0 + as.factor(age)|Hive_number), data=full_pollen2, family="binomial")

p7 <- glmmTMB(visual_avg ~ age + I(age**2) + 
                Avg_Temp60cm + Avg_Humidity +
                ar1(0 + as.factor(age)|Hive_number), data=full_pollen2, family="binomial")

p8 <- glmmTMB(visual_avg ~ age + I(age**2) + 
                Avg_Temp60cm + Avg_Radiation  +
                ar1(0 + as.factor(age)|Hive_number), data=full_pollen2, family="binomial")

p9 <- glmmTMB(visual_avg ~ age + I(age**2) + 
                Avg_Rainfall + Avg_Radiation + 
                ar1(0 + as.factor(age)|Hive_number), data=full_pollen2, family="binomial")

p10 <- glmmTMB(visual_avg ~ age + I(age**2) + 
                Avg_Rainfall + Avg_Humidity +
                ar1(0 + as.factor(age)|Hive_number), data=full_pollen2, family="binomial")

p11 <- glmmTMB(visual_avg ~ age + I(age**2) + 
                Avg_Radiation + Avg_Humidity +
                ar1(0 + as.factor(age)|Hive_number), data=full_pollen2, family="binomial")

p12 <- glmmTMB(visual_avg ~ age + I(age**2) + 
                Avg_Humidity + ar1(0 + as.factor(age)|Hive_number), data=full_pollen2, family="binomial")

p13 <- glmmTMB(visual_avg ~ age + I(age**2) + 
                Avg_Radiation + ar1(0 + as.factor(age)|Hive_number), data=full_pollen2, family="binomial")

p14 <- glmmTMB(visual_avg ~ age + I(age**2) + 
                Avg_Temp60cm + ar1(0 + as.factor(age)|Hive_number), data=full_pollen2, family="binomial")

p15 <- glmmTMB(visual_avg ~ age + I(age**2) + 
                Avg_Rainfall + ar1(0 + as.factor(age)|Hive_number), data=full_pollen2, family="binomial")

p16 <- glmmTMB(visual_avg ~ age + I(age**2) + ar1(0 + as.factor(age)|Hive_number), data=full_pollen2, family="binomial") ##no weather variables

p17 <- glmmTMB(visual_avg ~ age + ar1(0 + as.factor(age)|Hive_number), data=full_pollen2, family="binomial")##no quadratic term

p18 <- glmmTMB(visual_avg ~ age + (1|Hive_number), data=full_pollen2, family="binomial") ##no AR1 correlation structure

p19 <- glmmTMB(visual_avg ~ Avg_Temp60cm + Avg_Rainfall + Avg_Radiation + Avg_Humidity+
                ar1(0 + as.factor(age)|Hive_number), data=full_pollen2, family="binomial") ##no age, only weather variables

p20 <- glmmTMB(visual_avg ~ 1, data=full_pollen2, family="binomial") ##null

##generate AICc scores, weights, and delta
model.sel(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20)
##15 of 20 models within 6 AICc, will need to average
##all 15 have quadratic term and AR1 random effect

##average model
pollen_models <- model.sel(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20)

model.avg(pollen_models, subset = delta < 6)

summary(model.avg(pollen_models, subset = delta < 6)) ##conditional avg does not include zeros for absent terms
##in the averaged model, the effect sizes were small, so the environmental variables contributed very little
##justify use simplest (model 16) without env variables....most parsimonious
##report estimates and SE for completeness
##Could potentially use likelihood ratio test to provide p value working from most complex model to simplest. Not ideal, but if can't publish on AIC alone

##didn't achieve goal of model selection, going for delta 2
summary(model.avg(pollen_models, subset = delta < 2))


##going with top model
p7 <- glmmTMB(visual_avg ~ age + I(age**2) + 
                Avg_Temp60cm + Avg_Humidity +
                ar1(0 + as.factor(age)|Hive_number), data=full_pollen2, family="binomial")


#####assess residuals

plot(simulateResiduals(p7)) ##pretty good
hist(simulateResiduals(p7)) ##good

####random effect components

summary(p7) ##Variation= 0.3024 #Not too low, no singularity
##Correlation coefficient for temporal autocorrelation is 0.90 - glad we included that

###significance of model terms

Anova(p7) ##both sample day terms are significant, but bc we used AIC, don't mix p values

###estimates, effect sizes, compare means
##Throughout age for table

emmeans(p7, ~1, type='response', at=list(age=2)) #0.0904 prop w/ dye
emmeans(p7, ~1, type='response', at=list(age=4)) #0.156 prop w/ dye
emmeans(p7, ~1, type='response', at=list(age=6)) #0.242 prop w/ dye
emmeans(p7, ~1, type='response', at=list(age=8)) #0.339 prop w/ dye
emmeans(p7, ~1, type='response', at=list(age=10)) #0.433 prop w/ dye
emmeans(p7, ~1, type='response', at=list(age=12)) #0.512 prop w/ dye
emmeans(p7, ~1, type='response', at=list(age=14)) #0.573 prop w/ dye
emmeans(p7, ~1, type='response', at=list(age=16)) #0.614 prop w/ dye
emmeans(p7, ~1, type='response', at=list(age=18)) #0.636 prop w/ dye
emmeans(p7, ~1, type='response', at=list(age=19)) #0.640 prop w/ dye
emmeans(p7, ~1, type='response', at=list(age=20)) #0.640 prop w/ dye
emmeans(p7, ~1, type='response', at=list(age=21)) #0.635 prop w/ dye
emmeans(p7, ~1, type='response', at=list(age=22)) #0.626 prop w/ dye
emmeans(p7, ~1, type='response', at=list(age=24)) #0.595 prop w/ dye
emmeans(p7, ~1, type='response', at=list(age=26)) #0.543 prop w/ dye
emmeans(p7, ~1, type='response', at=list(age=28)) #0.472 prop w/ dye

####For plots####
##extract emmeans for the averaged model and each hive
##I know there has to be an easier way to do this, but I just couldn't figure it out
##could not get ggemeans or ggpredict to work with our model

all1 <- as.data.frame(emmeans(p7, ~1, type='response', at=list(age=1))) 
all2 <- as.data.frame(emmeans(p7, ~1, type='response', at=list(age=2))) 
all3 <- as.data.frame(emmeans(p7, ~1, type='response', at=list(age=3))) 
all4 <- as.data.frame(emmeans(p7, ~1, type='response', at=list(age=4))) 
all5 <- as.data.frame(emmeans(p7, ~1, type='response', at=list(age=5)))
all6 <- as.data.frame(emmeans(p7, ~1, type='response', at=list(age=6)))
all7 <- as.data.frame(emmeans(p7, ~1, type='response', at=list(age=7)))
all8 <- as.data.frame(emmeans(p7, ~1, type='response', at=list(age=8)))
all9 <- as.data.frame(emmeans(p7, ~1, type='response', at=list(age=9)))
all10 <-as.data.frame(emmeans(p7, ~1, type='response', at=list(age=10)))
all11 <- as.data.frame(emmeans(p7, ~1, type='response', at=list(age=11)))
all12 <-as.data.frame(emmeans(p7, ~1, type='response', at=list(age=12)))
all13 <- as.data.frame(emmeans(p7, ~1, type='response', at=list(age=13)))
all14 <-as.data.frame(emmeans(p7, ~1, type='response', at=list(age=14)))
all15 <- as.data.frame(emmeans(p7, ~1, type='response', at=list(age=15)))
all16 <-as.data.frame(emmeans(p7, ~1, type='response', at=list(age=16)))
all17 <- as.data.frame(emmeans(p7, ~1, type='response', at=list(age=17))) 
all18 <-as.data.frame(emmeans(p7, ~1, type='response', at=list(age=18))) 
all19 <-as.data.frame(emmeans(p7, ~1, type='response', at=list(age=19))) 
all20 <-as.data.frame(emmeans(p7, ~1, type='response', at=list(age=20))) 
all21 <-as.data.frame(emmeans(p7, ~1, type='response', at=list(age=21))) 
all22 <-as.data.frame(emmeans(p7, ~1, type='response', at=list(age=22)))
all23 <- as.data.frame(emmeans(p7, ~1, type='response', at=list(age=23))) 
all24 <-as.data.frame(emmeans(p7, ~1, type='response', at=list(age=24)))
all25 <- as.data.frame(emmeans(p7, ~1, type='response', at=list(age=25))) 
all26 <-as.data.frame(emmeans(p7, ~1, type='response', at=list(age=26)))
all27 <- as.data.frame(emmeans(p7, ~1, type='response', at=list(age=27))) 
all28 <-as.data.frame(emmeans(p7, ~1, type='response', at=list(age=28)))
all29 <- as.data.frame(emmeans(p7, ~1, type='response', at=list(age=29))) 

##binding all rows
all <- bind_rows(all1, all2, all3, all4, all5, all6, all7, all8, all9, all10, all11,
                 all12, all13, all14, all15, all16, all17, all18, all19, all20, all21,
                 all22, all23, all24, all25, all26, all27, all28, all29)

##add age and description for binding later
all$age <- c(1:29)
all$hive <- "All"

## ## ##
##hive 4
hive4 <- full_pollen2 %>% filter(Hive_number== 4)

h4 <- glmmTMB(visual_avg ~ age + I(age**2) + 
                Avg_Temp60cm + Avg_Humidity +
                ar1(0 + as.factor(age)|Hive_number), data=hive4, family="binomial")

plot(simulateResiduals(h4)) ##pretty good
hist(simulateResiduals(h4)) ##good

summary(h4)

h4_1 <- as.data.frame(emmeans(h4, ~1, type='response', at=list(age=1))) 
h4_2 <- as.data.frame(emmeans(h4, ~1, type='response', at=list(age=2))) 
h4_3 <- as.data.frame(emmeans(h4, ~1, type='response', at=list(age=3))) 
h4_4 <- as.data.frame(emmeans(h4, ~1, type='response', at=list(age=4))) 
h4_5 <- as.data.frame(emmeans(h4, ~1, type='response', at=list(age=5)))
h4_6 <- as.data.frame(emmeans(h4, ~1, type='response', at=list(age=6)))
h4_7 <- as.data.frame(emmeans(h4, ~1, type='response', at=list(age=7)))
h4_8 <- as.data.frame(emmeans(h4, ~1, type='response', at=list(age=8)))
h4_9 <- as.data.frame(emmeans(h4, ~1, type='response', at=list(age=9)))
h4_10 <-as.data.frame(emmeans(h4, ~1, type='response', at=list(age=10)))
h4_11 <- as.data.frame(emmeans(h4, ~1, type='response', at=list(age=11)))
h4_12 <-as.data.frame(emmeans(h4, ~1, type='response', at=list(age=12)))
h4_13 <- as.data.frame(emmeans(h4, ~1, type='response', at=list(age=13)))
h4_14 <-as.data.frame(emmeans(h4, ~1, type='response', at=list(age=14)))
h4_15 <- as.data.frame(emmeans(h4, ~1, type='response', at=list(age=15)))
h4_16 <-as.data.frame(emmeans(h4, ~1, type='response', at=list(age=16)))
h4_17 <- as.data.frame(emmeans(h4, ~1, type='response', at=list(age=17))) 
h4_18 <-as.data.frame(emmeans(h4, ~1, type='response', at=list(age=18))) 
h4_19 <-as.data.frame(emmeans(h4, ~1, type='response', at=list(age=19))) 
h4_20 <-as.data.frame(emmeans(h4, ~1, type='response', at=list(age=20))) 
h4_21 <-as.data.frame(emmeans(h4, ~1, type='response', at=list(age=21))) 
h4_22 <-as.data.frame(emmeans(h4, ~1, type='response', at=list(age=22)))
h4_23 <- as.data.frame(emmeans(h4, ~1, type='response', at=list(age=23))) 
h4_24 <-as.data.frame(emmeans(h4, ~1, type='response', at=list(age=24)))
h4_25 <- as.data.frame(emmeans(h4, ~1, type='response', at=list(age=25))) 
h4_26 <-as.data.frame(emmeans(h4, ~1, type='response', at=list(age=26)))
h4_27 <- as.data.frame(emmeans(h4, ~1, type='response', at=list(age=27))) 
h4_28 <-as.data.frame(emmeans(h4, ~1, type='response', at=list(age=28)))
h4_29 <- as.data.frame(emmeans(h4, ~1, type='response', at=list(age=29))) 

h4_emm <- bind_rows(h4_1, h4_2, h4_3, h4_4, h4_5, h4_6, h4_7, h4_8, h4_9, h4_10, h4_11,
                    h4_12, h4_13, h4_14, h4_15, h4_16, h4_17, h4_18, h4_19, h4_20, h4_21,
                    h4_22, h4_23, h4_24, h4_25, h4_26, h4_27, h4_28, h4_29)

h4_emm$age <- c(1:29)
h4_emm$hive <- "A"

## ## ##
##hive 5

hive5 <- full_pollen2 %>% filter(Hive_number== 5)

h5 <- glmmTMB(visual_avg ~ age + I(age**2) + 
                Avg_Temp60cm + Avg_Humidity +
                ar1(0 + as.factor(age)|Hive_number), data=hive5, family="binomial")

plot(simulateResiduals(h5)) ##pretty good
hist(simulateResiduals(h5)) ##good

summary(h5)

h5_1 <- as.data.frame(emmeans(h5, ~1, type='response', at=list(age=1))) 
h5_2 <- as.data.frame(emmeans(h5, ~1, type='response', at=list(age=2))) 
h5_3 <- as.data.frame(emmeans(h5, ~1, type='response', at=list(age=3))) 
h5_4 <- as.data.frame(emmeans(h5, ~1, type='response', at=list(age=4))) 
h5_5 <- as.data.frame(emmeans(h5, ~1, type='response', at=list(age=5)))
h5_6 <- as.data.frame(emmeans(h5, ~1, type='response', at=list(age=6)))
h5_7 <- as.data.frame(emmeans(h5, ~1, type='response', at=list(age=7)))
h5_8 <- as.data.frame(emmeans(h5, ~1, type='response', at=list(age=8)))
h5_9 <- as.data.frame(emmeans(h5, ~1, type='response', at=list(age=9)))
h5_10 <-as.data.frame(emmeans(h5, ~1, type='response', at=list(age=10)))
h5_11 <- as.data.frame(emmeans(h5, ~1, type='response', at=list(age=11)))
h5_12 <-as.data.frame(emmeans(h5, ~1, type='response', at=list(age=12)))
h5_13 <- as.data.frame(emmeans(h5, ~1, type='response', at=list(age=13)))
h5_14 <-as.data.frame(emmeans(h5, ~1, type='response', at=list(age=14)))
h5_15 <- as.data.frame(emmeans(h5, ~1, type='response', at=list(age=15)))
h5_16 <-as.data.frame(emmeans(h5, ~1, type='response', at=list(age=16)))
h5_17 <- as.data.frame(emmeans(h5, ~1, type='response', at=list(age=17))) 
h5_18 <-as.data.frame(emmeans(h5, ~1, type='response', at=list(age=18))) 
h5_19 <-as.data.frame(emmeans(h5, ~1, type='response', at=list(age=19))) 
h5_20 <-as.data.frame(emmeans(h5, ~1, type='response', at=list(age=20))) 
h5_21 <-as.data.frame(emmeans(h5, ~1, type='response', at=list(age=21))) 
h5_22 <-as.data.frame(emmeans(h5, ~1, type='response', at=list(age=22)))
h5_23 <- as.data.frame(emmeans(h5, ~1, type='response', at=list(age=23))) 
h5_24 <-as.data.frame(emmeans(h5, ~1, type='response', at=list(age=24)))
h5_25 <- as.data.frame(emmeans(h5, ~1, type='response', at=list(age=25))) 
h5_26 <-as.data.frame(emmeans(h5, ~1, type='response', at=list(age=26)))
h5_27 <- as.data.frame(emmeans(h5, ~1, type='response', at=list(age=27))) 
h5_28 <-as.data.frame(emmeans(h5, ~1, type='response', at=list(age=28)))
h5_29 <- as.data.frame(emmeans(h5, ~1, type='response', at=list(age=29))) 

h5_emm <- bind_rows(h5_1, h5_2, h5_3, h5_4, h5_5, h5_6, h5_7, h5_8, h5_9, h5_10, h5_11,
                    h5_12, h5_13, h5_14, h5_15, h5_16, h5_17, h5_18, h5_19, h5_20, h5_21,
                    h5_22, h5_23, h5_24, h5_25, h5_26, h5_27, h5_28, h5_29)

h5_emm$age <- c(1:29)
h5_emm$hive <- "B"

## ## ##
##hive 6

hive6 <- full_pollen2 %>% filter(Hive_number== 6)

h6 <- glmmTMB(visual_avg ~ age + I(age**2) + 
                Avg_Temp60cm + Avg_Humidity +
                ar1(0 + as.factor(age)|Hive_number), data=hive6, family="binomial")

plot(simulateResiduals(h6)) ##pretty good
hist(simulateResiduals(h6)) ##good

summary(h6)

h6_1 <- as.data.frame(emmeans(h6, ~1, type='response', at=list(age=1))) 
h6_2 <- as.data.frame(emmeans(h6, ~1, type='response', at=list(age=2))) 
h6_3 <- as.data.frame(emmeans(h6, ~1, type='response', at=list(age=3))) 
h6_4 <- as.data.frame(emmeans(h6, ~1, type='response', at=list(age=4))) 
h6_5 <- as.data.frame(emmeans(h6, ~1, type='response', at=list(age=5)))
h6_6 <- as.data.frame(emmeans(h6, ~1, type='response', at=list(age=6)))
h6_7 <- as.data.frame(emmeans(h6, ~1, type='response', at=list(age=7)))
h6_8 <- as.data.frame(emmeans(h6, ~1, type='response', at=list(age=8)))
h6_9 <- as.data.frame(emmeans(h6, ~1, type='response', at=list(age=9)))
h6_10 <-as.data.frame(emmeans(h6, ~1, type='response', at=list(age=10)))
h6_11 <- as.data.frame(emmeans(h6, ~1, type='response', at=list(age=11)))
h6_12 <-as.data.frame(emmeans(h6, ~1, type='response', at=list(age=12)))
h6_13 <- as.data.frame(emmeans(h6, ~1, type='response', at=list(age=13)))
h6_14 <-as.data.frame(emmeans(h6, ~1, type='response', at=list(age=14)))
h6_15 <- as.data.frame(emmeans(h6, ~1, type='response', at=list(age=15)))
h6_16 <-as.data.frame(emmeans(h6, ~1, type='response', at=list(age=16)))
h6_17 <- as.data.frame(emmeans(h6, ~1, type='response', at=list(age=17))) 
h6_18 <-as.data.frame(emmeans(h6, ~1, type='response', at=list(age=18))) 
h6_19 <-as.data.frame(emmeans(h6, ~1, type='response', at=list(age=19))) 
h6_20 <-as.data.frame(emmeans(h6, ~1, type='response', at=list(age=20))) 
h6_21 <-as.data.frame(emmeans(h6, ~1, type='response', at=list(age=21))) 
h6_22 <-as.data.frame(emmeans(h6, ~1, type='response', at=list(age=22)))
h6_23 <- as.data.frame(emmeans(h6, ~1, type='response', at=list(age=23))) 
h6_24 <-as.data.frame(emmeans(h6, ~1, type='response', at=list(age=24)))
h6_25 <- as.data.frame(emmeans(h6, ~1, type='response', at=list(age=25))) 
h6_26 <-as.data.frame(emmeans(h6, ~1, type='response', at=list(age=26)))
h6_27 <- as.data.frame(emmeans(h6, ~1, type='response', at=list(age=27))) 
h6_28 <-as.data.frame(emmeans(h6, ~1, type='response', at=list(age=28)))
h6_29 <- as.data.frame(emmeans(h6, ~1, type='response', at=list(age=29))) 

h6_emm <- bind_rows(h6_1, h6_2, h6_3, h6_4, h6_5, h6_6, h6_7, h6_8, h6_9, h6_10, h6_11,
                    h6_12, h6_13, h6_14, h6_15, h6_16, h6_17, h6_18, h6_19, h6_20, h6_21,
                    h6_22, h6_23, h6_24, h6_25, h6_26, h6_27, h6_28, h6_29)

h6_emm$age <- c(1:29)
h6_emm$hive <- "C"

## ## ##
##hive 7


hive7 <- full_pollen2 %>% filter(Hive_number== 7)

h7 <- glmmTMB(visual_avg ~ age + I(age**2) + 
                Avg_Temp60cm + Avg_Humidity +
                ar1(0 + as.factor(age)|Hive_number), data=hive7, family="binomial")

plot(simulateResiduals(h7)) ##pretty good
hist(simulateResiduals(h7)) ##good

summary(h7)

h7_1 <- as.data.frame(emmeans(h7, ~1, type='response', at=list(age=1))) 
h7_2 <- as.data.frame(emmeans(h7, ~1, type='response', at=list(age=2))) 
h7_3 <- as.data.frame(emmeans(h7, ~1, type='response', at=list(age=3))) 
h7_4 <- as.data.frame(emmeans(h7, ~1, type='response', at=list(age=4))) 
h7_5 <- as.data.frame(emmeans(h7, ~1, type='response', at=list(age=5)))
h7_6 <- as.data.frame(emmeans(h7, ~1, type='response', at=list(age=6)))
h7_7 <- as.data.frame(emmeans(h7, ~1, type='response', at=list(age=7)))
h7_8 <- as.data.frame(emmeans(h7, ~1, type='response', at=list(age=8)))
h7_9 <- as.data.frame(emmeans(h7, ~1, type='response', at=list(age=9)))
h7_10 <-as.data.frame(emmeans(h7, ~1, type='response', at=list(age=10)))
h7_11 <- as.data.frame(emmeans(h7, ~1, type='response', at=list(age=11)))
h7_12 <-as.data.frame(emmeans(h7, ~1, type='response', at=list(age=12)))
h7_13 <- as.data.frame(emmeans(h7, ~1, type='response', at=list(age=13)))
h7_14 <-as.data.frame(emmeans(h7, ~1, type='response', at=list(age=14)))
h7_15 <- as.data.frame(emmeans(h7, ~1, type='response', at=list(age=15)))
h7_16 <-as.data.frame(emmeans(h7, ~1, type='response', at=list(age=16)))
h7_17 <- as.data.frame(emmeans(h7, ~1, type='response', at=list(age=17))) 
h7_18 <-as.data.frame(emmeans(h7, ~1, type='response', at=list(age=18))) 
h7_19 <-as.data.frame(emmeans(h7, ~1, type='response', at=list(age=19))) 
h7_20 <-as.data.frame(emmeans(h7, ~1, type='response', at=list(age=20))) 
h7_21 <-as.data.frame(emmeans(h7, ~1, type='response', at=list(age=21))) 
h7_22 <-as.data.frame(emmeans(h7, ~1, type='response', at=list(age=22)))
h7_23 <- as.data.frame(emmeans(h7, ~1, type='response', at=list(age=23))) 
h7_24 <-as.data.frame(emmeans(h7, ~1, type='response', at=list(age=24)))
h7_25 <- as.data.frame(emmeans(h7, ~1, type='response', at=list(age=25))) 
h7_26 <-as.data.frame(emmeans(h7, ~1, type='response', at=list(age=26)))
h7_27 <- as.data.frame(emmeans(h7, ~1, type='response', at=list(age=27))) 
h7_28 <-as.data.frame(emmeans(h7, ~1, type='response', at=list(age=28)))
h7_29 <- as.data.frame(emmeans(h7, ~1, type='response', at=list(age=29))) 

h7_emm <- bind_rows(h7_1, h7_2, h7_3, h7_4, h7_5, h7_6, h7_7, h7_8, h7_9, h7_10, h7_11,
                    h7_12, h7_13, h7_14, h7_15, h7_16, h7_17, h7_18, h7_19, h7_20, h7_21,
                    h7_22, h7_23, h7_24, h7_25, h7_26, h7_27, h7_28, h7_29)

h7_emm$age <- c(1:29)
h7_emm$hive <- "D"

## ## ##
###hive 8

hive8 <- full_pollen2 %>% filter(Hive_number== 8)

h8 <- glmmTMB(visual_avg ~ age + I(age**2) + 
                Avg_Temp60cm + Avg_Humidity +
                ar1(0 + as.factor(age)|Hive_number), data=hive8, family="binomial")

plot(simulateResiduals(h8)) ##pretty good
hist(simulateResiduals(h8)) ##good

summary(h8)

h8_1 <- as.data.frame(emmeans(h8, ~1, type='response', at=list(age=1))) 
h8_2 <- as.data.frame(emmeans(h8, ~1, type='response', at=list(age=2))) 
h8_3 <- as.data.frame(emmeans(h8, ~1, type='response', at=list(age=3))) 
h8_4 <- as.data.frame(emmeans(h8, ~1, type='response', at=list(age=4))) 
h8_5 <- as.data.frame(emmeans(h8, ~1, type='response', at=list(age=5)))
h8_6 <- as.data.frame(emmeans(h8, ~1, type='response', at=list(age=6)))
h8_7 <- as.data.frame(emmeans(h8, ~1, type='response', at=list(age=7)))
h8_8 <- as.data.frame(emmeans(h8, ~1, type='response', at=list(age=8)))
h8_9 <- as.data.frame(emmeans(h8, ~1, type='response', at=list(age=9)))
h8_10 <-as.data.frame(emmeans(h8, ~1, type='response', at=list(age=10)))
h8_11 <- as.data.frame(emmeans(h8, ~1, type='response', at=list(age=11)))
h8_12 <-as.data.frame(emmeans(h8, ~1, type='response', at=list(age=12)))
h8_13 <- as.data.frame(emmeans(h8, ~1, type='response', at=list(age=13)))
h8_14 <-as.data.frame(emmeans(h8, ~1, type='response', at=list(age=14)))
h8_15 <- as.data.frame(emmeans(h8, ~1, type='response', at=list(age=15)))
h8_16 <-as.data.frame(emmeans(h8, ~1, type='response', at=list(age=16)))
h8_17 <- as.data.frame(emmeans(h8, ~1, type='response', at=list(age=17))) 
h8_18 <-as.data.frame(emmeans(h8, ~1, type='response', at=list(age=18))) 
h8_19 <-as.data.frame(emmeans(h8, ~1, type='response', at=list(age=19))) 
h8_20 <-as.data.frame(emmeans(h8, ~1, type='response', at=list(age=20))) 
h8_21 <-as.data.frame(emmeans(h8, ~1, type='response', at=list(age=21))) 
h8_22 <-as.data.frame(emmeans(h8, ~1, type='response', at=list(age=22)))
h8_23 <- as.data.frame(emmeans(h8, ~1, type='response', at=list(age=23))) 
h8_24 <-as.data.frame(emmeans(h8, ~1, type='response', at=list(age=24)))
h8_25 <- as.data.frame(emmeans(h8, ~1, type='response', at=list(age=25))) 
h8_26 <-as.data.frame(emmeans(h8, ~1, type='response', at=list(age=26)))
h8_27 <- as.data.frame(emmeans(h8, ~1, type='response', at=list(age=27))) 
h8_28 <-as.data.frame(emmeans(h8, ~1, type='response', at=list(age=28)))
h8_29 <- as.data.frame(emmeans(h8, ~1, type='response', at=list(age=29))) 

h8_emm <- bind_rows(h8_1, h8_2, h8_3, h8_4, h8_5, h8_6, h8_7, h8_8, h8_9, h8_10, h8_11,
                    h8_12, h8_13, h8_14, h8_15, h8_16, h8_17, h8_18, h8_19, h8_20, h8_21,
                    h8_22, h8_23, h8_24, h8_25, h8_26, h8_27, h8_28, h8_29)

h8_emm$age <- c(1:29)
h8_emm$hive <- "E"

## ## ##
##hive 9

hive9 <- full_pollen2 %>% filter(Hive_number== 9)

h9 <- glmmTMB(visual_avg ~ age + I(age**2) + 
                Avg_Temp60cm + Avg_Humidity +
                ar1(0 + as.factor(age)|Hive_number), data=hive9, family="binomial")

plot(simulateResiduals(h9)) ##pretty good
hist(simulateResiduals(h9)) ##good

summary(h9)

h9_1 <- as.data.frame(emmeans(h9, ~1, type='response', at=list(age=1))) 
h9_2 <- as.data.frame(emmeans(h9, ~1, type='response', at=list(age=2))) 
h9_3 <- as.data.frame(emmeans(h9, ~1, type='response', at=list(age=3))) 
h9_4 <- as.data.frame(emmeans(h9, ~1, type='response', at=list(age=4))) 
h9_5 <- as.data.frame(emmeans(h9, ~1, type='response', at=list(age=5)))
h9_6 <- as.data.frame(emmeans(h9, ~1, type='response', at=list(age=6)))
h9_7 <- as.data.frame(emmeans(h9, ~1, type='response', at=list(age=7)))
h9_8 <- as.data.frame(emmeans(h9, ~1, type='response', at=list(age=8)))
h9_9 <- as.data.frame(emmeans(h9, ~1, type='response', at=list(age=9)))
h9_10 <-as.data.frame(emmeans(h9, ~1, type='response', at=list(age=10)))
h9_11 <- as.data.frame(emmeans(h9, ~1, type='response', at=list(age=11)))
h9_12 <-as.data.frame(emmeans(h9, ~1, type='response', at=list(age=12)))
h9_13 <- as.data.frame(emmeans(h9, ~1, type='response', at=list(age=13)))
h9_14 <-as.data.frame(emmeans(h9, ~1, type='response', at=list(age=14)))
h9_15 <- as.data.frame(emmeans(h9, ~1, type='response', at=list(age=15)))
h9_16 <-as.data.frame(emmeans(h9, ~1, type='response', at=list(age=16)))
h9_17 <- as.data.frame(emmeans(h9, ~1, type='response', at=list(age=17))) 
h9_18 <-as.data.frame(emmeans(h9, ~1, type='response', at=list(age=18))) 
h9_19 <-as.data.frame(emmeans(h9, ~1, type='response', at=list(age=19))) 
h9_20 <-as.data.frame(emmeans(h9, ~1, type='response', at=list(age=20))) 
h9_21 <-as.data.frame(emmeans(h9, ~1, type='response', at=list(age=21))) 
h9_22 <-as.data.frame(emmeans(h9, ~1, type='response', at=list(age=22)))
h9_23 <- as.data.frame(emmeans(h9, ~1, type='response', at=list(age=23))) 
h9_24 <-as.data.frame(emmeans(h9, ~1, type='response', at=list(age=24)))
h9_25 <- as.data.frame(emmeans(h9, ~1, type='response', at=list(age=25))) 
h9_26 <-as.data.frame(emmeans(h9, ~1, type='response', at=list(age=26)))
h9_27 <- as.data.frame(emmeans(h9, ~1, type='response', at=list(age=27))) 
h9_28 <-as.data.frame(emmeans(h9, ~1, type='response', at=list(age=28)))
h9_29 <- as.data.frame(emmeans(h9, ~1, type='response', at=list(age=29))) 

h9_emm <- bind_rows(h9_1, h9_2, h9_3, h9_4, h9_5, h9_6, h9_7, h9_8, h9_9, h9_10, h9_11,
                    h9_12, h9_13, h9_14, h9_15, h9_16, h9_17, h9_18, h9_19, h9_20, h9_21,
                    h9_22, h9_23, h9_24, h9_25, h9_26, h9_27, h9_28, h9_29)

h9_emm$age <- c(1:29)
h9_emm$hive <- "F"

##bind individual hives together
manual_emm <- bind_rows(h4_emm, h5_emm, h6_emm, h7_emm, h8_emm, h9_emm)


####PLOTS#####

##plotting emmeans directly
ggplot()+
  geom_line(data=manual_emm, aes(x=age, y=prob, color=hive, linetype=hive), size=1.5, linetype="dashed")+
  geom_line(data=all, aes(x=age, y=prob), color="black", size=2)+ #geom_ribbon(data=all, aes(x= age, ymin=lower.CL, ymax=upper.CL), alpha=0.2)+
  theme_bw()+ theme(strip.background = element_blank(), panel.grid.minor= element_blank(), panel.grid.major= element_blank(),
        text = element_text(size=16))+
  labs(x="Worker age (day)", y="Dye presence")+
  scale_color_manual(values=c("#FDE725FF", "#AADC32FF","#5DC863FF","#27AD81FF","#21908CFF", "#2C728EFF"), "Colony", labels = c("A","B", "C", "D", "E", "F"))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 15)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))

library(scales) ##getting manual color codes from viridis package               
show_col(viridis_pal()(9)) ##use from bottom right corner, right to left, bottom to top

## ## ##
##plotting using ggplot's glm in geom_smooth
str(full_pollen2)
full_pollen2$Hive_number <- as.factor(full_pollen2$Hive_number)

##plain b&w with data points
##NOTE: sample day called, not age
ggplot(full_pollen2, aes(x=Sample_day, y=visual_avg))+geom_jitter(height=0.1, alpha=0.12, size=1.75)+
  geom_smooth(color="black", method="glm", method.args=list(family="binomial"), formula = y ~ x+I(x^2)) + theme_bw() + 
  theme(strip.background = element_blank(), panel.grid.minor= element_blank(), panel.grid.major= element_blank(),
        text = element_text(size=16))+
  labs(x="Sample Day", y="Dye Presence")


##Final plot
ggplot(full_pollen2, aes(x=age, y=visual_avg))+
  geom_smooth(full_pollen2, mapping=aes(x=age, y=visual_avg, linetype=Hive_number, color=Hive_number), method="glm", method.args=list(family="binomial"), 
              formula = y ~ x+I(x^2), se=FALSE, linetype="dashed", size=1) + theme_bw() + 
  geom_smooth(color="black", method="glm", method.args=list(family="binomial"), formula = y ~ x+I(x^2), size=2)+ ##can change order of geom_smooths for whatver line(s) you want on top
  theme(strip.background = element_blank(), panel.grid.minor= element_blank(), panel.grid.major= element_blank(),
        text = element_text(size=16))+
  labs(x="Worker age (day)", y="Dye presence")+
  scale_color_viridis(discrete = T, direction = -1, "Colony", labels = c("A","B", "C", "D", "E", "F"))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 15)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
  #changed hives 4-9 to A-F for lack of confusion

####ENV PLOTS####
##unscaled env data
pollen_env_avg_US <- pollen_env2 %>% dplyr::select(-Station.ID, -Station.Name, -Time, -True_Date) %>% 
  group_by(Period, Date) %>% summarise_each(funs(mean)) %>% rename_with( ~ paste("Avg", .x, sep = "_")) %>% 
  rename("Date"="Avg_Date") %>% filter(Avg_Period < 31)

unscaled_pollen <- left_join(pollen2, pollen_env_avg_US, by = "Date") 
head(unscaled_pollen) ##wooho

unscaled_pollen2 <- filter(unscaled_pollen, Treatment =="treatment")

unscaled_pollen2$Avg_Rainfall <- unscaled_pollen2$Avg_Rainfall*2.54



##extract residuals of model without weather variables and plot against each variable
p16 <- glmmTMB(visual_avg ~ age + I(age**2) + ar1(0 + as.factor(age)|Hive_number), data=full_pollen2, family="binomial") ##no weather variables

unscaled_pollen2$resid <-residuals(p16)

ggplot(unscaled_pollen2, aes(x=Avg_Temp60cm, y=resid))+
  theme_bw() + geom_jitter(unscaled_pollen2, mapping=aes(x=Avg_Temp60cm, y=resid, color=Treatment), height=0.1, size=1.75)+
  geom_smooth(color="black", method="glm", size=2)+ 
  theme(strip.background = element_blank(), panel.grid.minor= element_blank(), panel.grid.major= element_blank(),
        text = element_text(size=16), legend.position="none")+ ylim(-1.1,1.1)+
  labs(x=expression(paste('Average Daily Temperature (',~degree,'C)',sep='')), y="Model Residuals")+
  scale_color_viridis(discrete = T, direction = -1, begin= .85, end=.85, alpha=0.08)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) 

ggplot(unscaled_pollen2, aes(x=Avg_Rainfall, y=resid))+
  theme_bw() + geom_jitter(unscaled_pollen2, mapping=aes(x=Avg_Rainfall, y=resid, color=Treatment), height=0.1, alpha=0.1, size=1.75)+
  geom_smooth(color="black", method="glm", size=2)+ 
  theme(strip.background = element_blank(), panel.grid.minor= element_blank(), panel.grid.major= element_blank(),
        text = element_text(size=16), legend.position="none")+ ylim(-1.1,1.1)+
  labs(x="Average Daily Rainfall (cm)", y="Model Residuals")+
  scale_color_viridis(discrete = T, direction = -1, begin= .55, end=.55, alpha=0.08)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) 

ggplot(unscaled_pollen2, aes(x=Avg_Radiation, y=resid))+
  theme_bw() + geom_jitter(unscaled_pollen2, mapping=aes(x=Avg_Radiation, y=resid, color=Treatment), height=0.1, alpha=0.1, size=1.75)+
  geom_smooth(color="black", method="glm", size=2)+ 
  theme(strip.background = element_blank(), panel.grid.minor= element_blank(), panel.grid.major= element_blank(),
        text = element_text(size=16), legend.position="none")+ ylim(-1.1,1.1)+
  labs(x=expression(Average~Daily~Radiation~(W/m^{"2"})), y="Model Residuals")+ ##use expression for superscript
  scale_color_viridis(discrete = T, direction = -1, begin= .25, end=.25, alpha=0.08)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) 

ggplot(unscaled_pollen2, aes(x=Avg_Humidity, y=resid))+
  theme_bw() + geom_jitter(unscaled_pollen2, mapping=aes(x=Avg_Humidity, y=resid, color=Treatment), height=0.1, alpha=0.1, size=1.75)+
  geom_smooth(color="black", method="glm", size=2)+ 
  theme(strip.background = element_blank(), panel.grid.minor= element_blank(), panel.grid.major= element_blank(),
        text = element_text(size=16), legend.position="none")+ ylim(-1.1,1.1)+
  labs(x="Average Daily Humidity (%)", y="Model Residuals")+
  scale_color_viridis(discrete = T, direction = -1, begin= 0, end=0, alpha=0.08)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))

##effect of temp and humidity are so small, not biologically relevant
##can put these plots in supplemental?


##get citations for write up
citation(package="DHARMa") 

##for reproducing with same version of packages
sessionInfo()
