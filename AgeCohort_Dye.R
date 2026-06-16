#Kleckner et al. 2026 - DOI:10.1093/jee/toag155
#Script 1/2 - Dye in worker guts 

####1) Set up ####

##load in packages
library(tidyverse)
library(GGally)
library(glmmTMB)
library(DHARMa)
library(car)
library(corrplot)
library(MuMIn)
library(boot)
library(emmeans)
library(scales)
library(viridis)
library(corrplot)
library(performance)

#####1.1) Dye visuals ####
pollen <- read.csv("Age_cohort_visuals.csv")
str(pollen) #check structure

#metadata
  #Treatment - diet treatments: treatment (dyed PSP), negative (undyed PSP), consumption (modeling clay), and positive (dyed fondant)
  #Sample_day - count of sampling days, one less than worker age
  #Date - calendar date of sampling
  #Hive_number - hive ID number
  #Sample_number - count of unique bees sampled from one hive on one day
  #Visual1 and Visual 2 - Two blind graders of visual dye presence (Y) or absence (N) in sample

##convert Y/N to 1/0
pollen1 <- pollen %>%
  mutate(Visual1 = if_else(Visual1 == "Y", 1L, 0L),
         Visual2 = if_else(Visual2 == "Y", 1L , 0L))
str(pollen1)

##create column for worker age (one more than sampling day)
pollen1$age <- pollen1$Sample_day+1

##change hive to factor
pollen1$Hive_number <- as.factor(pollen1$Hive_number)

##remove visual gradings that disagree
pollen1$visual_avg <- (pollen1$Visual1+pollen1$Visual2)/2

hist(pollen1$visual_avg) ##very few 0.5 values out of the whole data set, meaning very few instances the two visual graders disagreed

pollen2 <- filter(pollen1, visual_avg !="0.5") ##removed 0.5 values and NAs (marked bees no longer present at the end of the sampling period)

hist(pollen2$visual_avg) ##verify removals


#####1.2) Environmental data ####

##From UF IFAS FAWN Alachua Station 260, trimmed for sampling days in October 2020
pollen_env <- read.csv("Age_cohort_enviro.csv") 
str(pollen_env) #check structure

#metadata
  #Station.ID - UF IFAS FAWN station ID
  #Station.Name - UF IFAS FAWN station name
  #True_Date - calendar date of weather event
  #Date - date corresponding with sampling dates (24 hours from previous sampling, 12pm-11am)
  #Time - hour of weather measurement
  #Period - 24 hour period corresponding to sampling date (12pm-11am)
  #Temp60cm - Temperature (C) from 60 cm above the ground
  #Humidity - % humidity
  #Rainfall - inches of rainfall
  #Radiation - W/m^2 of solar radiation


##remove unnecessary columns and average weather across sampling period (Date, not True_Date)
pollen_env_avg <- pollen_env %>%
  select(-Station.ID, -Station.Name, -Time, -True_Date) %>%
  group_by(Period, Date) %>%
  summarise(across(everything(), mean, .names = "Avg_{.col}"), .groups = "drop") %>%
  filter(Period < 31)

head(pollen_env_avg) #check new dataset

##change rainfall from inches to cm
pollen_env_avg$Avg_Rainfall <- pollen_env_avg$Avg_Rainfall*2.54

##scale env variables
pollen_env_avg[c("Avg_Temp60cm", "Avg_Rainfall", "Avg_Humidity", "Avg_Radiation")] <- scale(pollen_env_avg[c("Avg_Temp60cm", "Avg_Rainfall", "Avg_Humidity", "Avg_Radiation")])
##center is zero, a 1.00 unit increase represents one standard deviation from the original scale


##cbind with visuals using dates
full_pollen <- left_join(pollen2, pollen_env_avg, by = "Date") 
head(full_pollen)






####2) Model construction####

##remove controls for modeling
full_pollen2 <- filter(full_pollen, Treatment =="treatment")

head(full_pollen2)
tail(full_pollen2)

##check correlations
pollen_cor2 <- full_pollen2 %>% dplyr::select(-Date, -Treatment, -Hive_number) #remove non-numeric columns
pollen_cor3 <- cor(as.matrix(pollen_cor2)) ##make correlation matrix


corrplot(pollen_cor3, method = "circle", type = "upper") ## plots strength of correlation as color-coded circles
###some environmental variables (e.g. Humidity and Radiation) are correlated with each other, but none with the visual averages


##start with most complex full model

p1 <- glmmTMB(visual_avg ~ age + I(age**2) + 
                Avg_Temp60cm + Avg_Rainfall + Avg_Radiation + Avg_Humidity +
                ar1(0 + as.factor(age)|Hive_number), data=full_pollen2, family="binomial") ##quadratic term and AR1 random effect
                ##AR1 correlation structure accounts for correlation of points in time closer together

##pairwise by weather variables
p2 <- glmmTMB(visual_avg ~ age + I(age**2) + 
                Avg_Temp60cm + Avg_Rainfall + Avg_Radiation + #no humidity
                ar1(0 + as.factor(age)|Hive_number), data=full_pollen2, family="binomial") 

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


####3) Model selection####

##generate AICc scores, weights, and delta
pollen_models <- model.sel(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20)
pollen_models
  ##15 of 20 models within 6 AICc -> average with stricter cutoff
  ##all 15 have quadratic term and AR1 random effect

##average model
summary(model.avg(pollen_models, subset = delta < 2)) ##delta 6 cut off did not achieve goal of model selection, using 2 instead


##going with top model for visuals + rest of analysis
p7 <- glmmTMB(visual_avg ~ age + I(age**2) + 
                Avg_Temp60cm + Avg_Humidity +
                ar1(0 + as.factor(age)|Hive_number), data=full_pollen2, family="binomial")

##assess residuals

plot(simulateResiduals(p7)) ##pretty good
hist(simulateResiduals(p7)) ##good

####random effect components

summary(p7) ##Variation= 0.3024 #Not too low, no singularity
##Correlation coefficient for temporal autocorrelation is 0.90

###estimates, effect sizes, compare means across age
emmeans(p7, ~1, type='response', at=list(age=2)) #0.0904 prop w/ dye
emmeans(p7, ~1, type='response', at=list(age=4)) #0.156 prop w/ dye
emmeans(p7, ~1, type='response', at=list(age=6)) #0.242 prop w/ dye
emmeans(p7, ~1, type='response', at=list(age=8)) #0.339 prop w/ dye
emmeans(p7, ~1, type='response', at=list(age=10)) #0.433 prop w/ dye
emmeans(p7, ~1, type='response', at=list(age=12)) #0.512 prop w/ dye
emmeans(p7, ~1, type='response', at=list(age=14)) #0.573 prop w/ dye
emmeans(p7, ~1, type='response', at=list(age=16)) #0.614 prop w/ dye
emmeans(p7, ~1, type='response', at=list(age=18)) #0.636 prop w/ dye
emmeans(p7, ~1, type='response', at=list(age=19)) #0.640 prop w/ dye >>Peak
emmeans(p7, ~1, type='response', at=list(age=20)) #0.640 prop w/ dye >>Peak
emmeans(p7, ~1, type='response', at=list(age=21)) #0.635 prop w/ dye
emmeans(p7, ~1, type='response', at=list(age=22)) #0.626 prop w/ dye
emmeans(p7, ~1, type='response', at=list(age=24)) #0.595 prop w/ dye
emmeans(p7, ~1, type='response', at=list(age=26)) #0.543 prop w/ dye
emmeans(p7, ~1, type='response', at=list(age=28)) #0.472 prop w/ dye



####4) Treatment effects####

#visualize all treatments
ggplot(pollen2, aes(x=Sample_day, y=visual_avg, color=Treatment))+geom_jitter(height=0.1, alpha=0.1)+
  facet_wrap(~Treatment)+geom_smooth(color="black", method="glm", method.args=list(family="binomial"), 
                                     formula = y ~ x+I(x^2))

#no data for consumption or negative control, remove before stats
full_pollen3 <- full_pollen %>% 
  filter(Treatment== "treatment" | Treatment == "positive ")

#build glmm
vis_treatment1 <- glmmTMB(visual_avg ~ age * Treatment + I(age^2)*Treatment  #start with both interactions for linear and quadratic age
                          + ar1(0 + as.factor(age)|Hive_number), ##ar1 random effect for autocorrelation
                          data=full_pollen3, family = binomial())

check_singularity(vis_treatment1)#true, remove ar1 in random effect

#model 2
vis_treatment2 <- glmmTMB(visual_avg ~ age * Treatment + I(age^2)*Treatment + (1|Hive_number), data=full_pollen3, family = binomial())
check_singularity(vis_treatment2)#false, fixed

##assess residuals
plot(simulateResiduals(vis_treatment2))
hist(simulateResiduals(vis_treatment2)) 
hist(residuals(vis_treatment2)) 

Anova(vis_treatment2) #both interactions not significant, remove both

#model 3 - final
vis_treatment3 <- glmmTMB(visual_avg ~ age + Treatment + I(age^2) + (1|Hive_number), data=full_pollen3, family = binomial())

#assess residuals
plot(simulateResiduals(vis_treatment3))
hist(simulateResiduals(vis_treatment3)) 
hist(residuals(vis_treatment3)) 

Anova(vis_treatment3) #both treatment and age have effect

#compare emmeans of treatments
emmeans(vis_treatment3, pairwise~Treatment, type= "response")

#extract emmeans across age for figure
emm_vis_treatment <- emmip(vis_treatment3, Treatment~age, type='response', CIs = TRUE, 
                           at=list(age = c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29)), plotit = FALSE )


####5) Visualization####

#####5.1) Dye presence across worker age (Fig 3)####

##extract emmeans for plotting for the averaged model and each hive
  ##I know there has to be an easier way to do this, but I just couldn't figure it out
  ##could not get ggemeans, ggpredict, or emmip to work with our model

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


##getting manual color codes from viridis package               
show_col(viridis_pal()(9)) ##use from bottom right corner, right to left, bottom to top

##plot

##emmeans
ggplot()+
  geom_smooth(full_pollen2, mapping=aes(x=age, y=visual_avg, linetype=Hive_number, color=Hive_number), method="glm", method.args=list(family="binomial"), 
              formula = y ~ x+I(x^2), se=FALSE, linetype="dashed", size=1.5) + theme_bw() + 
  geom_line(data=all, aes(x=age, y=prob), color="black", size=2)+ geom_ribbon(data=all, aes(x= age, ymin=asymp.LCL, ymax=asymp.UCL), alpha=0.2)+
   # geom_smooth(color="black", method="glm", method.args=list(family="binomial"), formula = y ~ x+I(x^2), size=2)+ ##can change order of geom_smooths for whatver line(s) you want on top
  theme(strip.background = element_blank(), panel.grid.minor= element_blank(), panel.grid.major= element_blank(),
        text = element_text(size=14), legend.position = "none",
        axis.title.x = element_text(margin = margin(t = 8)), # Add space above x-axis title
        axis.title.y = element_text(margin = margin(r = 8)) )+
  labs(x="Adult worker age (day)", y="Proportion of workers with dye present")+
  scale_color_viridis(discrete = T, direction = -1)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))

##just geom_smooth
ggplot(full_pollen2, aes(x=age, y=visual_avg))+
  geom_smooth(full_pollen2, mapping=aes(x=age, y=visual_avg, linetype=Hive_number, color=Hive_number), method="glm", method.args=list(family="binomial"), 
              formula = y ~ x+I(x^2), se=FALSE, linetype="dashed", size=1.5) + theme_bw() + 
  geom_smooth(color="black", method="glm", method.args=list(family="binomial"), formula = y ~ x+I(x^2), size=2)+ ##can change order of geom_smooths for whatver line(s) you want on top
  theme(strip.background = element_blank(), panel.grid.minor= element_blank(), panel.grid.major= element_blank(),
        text = element_text(size=14), legend.position = "none",
        axis.title.x = element_text(margin = margin(t = 8)), # Add space above x-axis title
        axis.title.y = element_text(margin = margin(r = 8)) )+
  labs(x="Adult worker age (day)", y="Proportion of workers with dye present")+
  scale_color_viridis(discrete = T, direction = -1)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))


##NOTE: The second figure without emmeans accidentally ended up in the final manuscript, 
        #but does not match the caption (no emmean for the solid line).
        #The results are nearly identical, but the error ribbon is wider with emmeans

#####5.2) Environmental variables (supplemental)####

##using unscaled env data
pollen_env_avg_US <- pollen_env %>% select(-Station.ID, -Station.Name, -Time, -True_Date) %>%
  group_by(Period, Date) %>%
  summarise(across(everything(), mean, .names = "Avg_{.col}"), .groups = "drop") %>%
  filter(Period < 31)

unscaled_pollen <- left_join(pollen2, pollen_env_avg_US, by = "Date") 
head(unscaled_pollen)

unscaled_pollen2 <- filter(unscaled_pollen, Treatment =="treatment")

unscaled_pollen2$Avg_Rainfall <- unscaled_pollen2$Avg_Rainfall*2.54 #inches to cm



##extract residuals of model without weather variables and plot against each variable

##method one: residuals of model without any environmental variables
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

##method two: residuals of model without one environmental variables against that removed variable
#temp
no_temp <- glmmTMB(visual_avg ~ age + I(age**2) + 
                     Avg_Rainfall + Avg_Radiation + Avg_Humidity +
                     ar1(0 + as.factor(age)|Hive_number), data=full_pollen2, family="binomial")

unscaled_pollen2$no_temp <-residuals(no_temp)

ggplot(unscaled_pollen2, aes(x=Avg_Temp60cm, y=no_temp))+
  theme_bw() + geom_jitter(unscaled_pollen2, mapping=aes(x=Avg_Temp60cm, y=no_temp, color=Treatment), height=0.1, size=1.75)+
  geom_smooth(color="black", method="glm", size=2)+ 
  theme(strip.background = element_blank(), panel.grid.minor= element_blank(), panel.grid.major= element_blank(),
        text = element_text(size=16), legend.position="none")+ ylim(-1.1,1.1)+
  labs(x=expression(paste('Average Daily Temperature (',~degree,'C)',sep='')), y="Model Residuals")+
  scale_color_viridis(discrete = T, direction = -1, begin= .85, end=.85, alpha=0.08)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) 

#rainfall
no_rain <- glmmTMB(visual_avg ~ age + I(age**2) + 
                     Avg_Temp60cm + Avg_Radiation + Avg_Humidity +
                     ar1(0 + as.factor(age)|Hive_number), data=full_pollen2, family="binomial")

unscaled_pollen2$no_rain <-residuals(no_rain)

ggplot(unscaled_pollen2, aes(x=Avg_Rainfall, y=no_rain))+
  theme_bw() + geom_jitter(unscaled_pollen2, mapping=aes(x=Avg_Rainfall, y=no_rain, color=Treatment), height=0.1, alpha=0.1, size=1.75)+
  geom_smooth(color="black", method="glm", size=2)+ 
  theme(strip.background = element_blank(), panel.grid.minor= element_blank(), panel.grid.major= element_blank(),
        text = element_text(size=16), legend.position="none")+ ylim(-1.1,1.1)+
  labs(x="Average Daily Rainfall (cm)", y="Model Residuals")+
  scale_color_viridis(discrete = T, direction = -1, begin= .55, end=.55, alpha=0.08)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) 

#radiation
no_rad <- glmmTMB(visual_avg ~ age + I(age**2) + 
                    Avg_Temp60cm + Avg_Rainfall + Avg_Humidity +
                    ar1(0 + as.factor(age)|Hive_number), data=full_pollen2, family="binomial")

unscaled_pollen2$no_rad <-residuals(no_rad)

ggplot(unscaled_pollen2, aes(x=Avg_Radiation, y=no_rad))+
  theme_bw() + geom_jitter(unscaled_pollen2, mapping=aes(x=Avg_Radiation, y=no_rad, color=Treatment), height=0.1, alpha=0.1, size=1.75)+
  geom_smooth(color="black", method="glm", size=2)+ 
  theme(strip.background = element_blank(), panel.grid.minor= element_blank(), panel.grid.major= element_blank(),
        text = element_text(size=16), legend.position="none")+ ylim(-1.1,1.1)+
  labs(x=expression(Average~Daily~Radiation~(W/m^{"2"})), y="Model Residuals")+ ##use expression for superscript
  scale_color_viridis(discrete = T, direction = -1, begin= .25, end=.25, alpha=0.08)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) 

#humidity
no_hum <- glmmTMB(visual_avg ~ age + I(age**2) + 
                    Avg_Temp60cm + Avg_Rainfall + Avg_Radiation +
                    ar1(0 + as.factor(age)|Hive_number), data=full_pollen2, family="binomial")

unscaled_pollen2$no_hum <-residuals(no_hum)

ggplot(unscaled_pollen2, aes(x=Avg_Humidity, y=no_hum))+
  theme_bw() + geom_jitter(unscaled_pollen2, mapping=aes(x=Avg_Humidity, y=no_hum, color=Treatment), height=0.1, alpha=0.1, size=1.75)+
  geom_smooth(color="black", method="glm", size=2)+ 
  theme(strip.background = element_blank(), panel.grid.minor= element_blank(), panel.grid.major= element_blank(),
        text = element_text(size=16), legend.position="none")+ ylim(-1.1,1.1)+
  labs(x="Average Daily Humidity (%)", y="Model Residuals")+
  scale_color_viridis(discrete = T, direction = -1, begin= 0, end=0, alpha=0.08)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))



#####5.3) Treatment effects (Fig 4)####

treatment_counts <- full_pollen %>%
  group_by(Treatment) %>%
  summarise(count = n())

treatment_counts

new_labels <- c(
  "consumption " = "Dyed modeling clay",
  "negative " = "Undyed pollen substitute patty",
  "positive " = "Dyed fondant",
  "treatment" = "Dyed pollen substitute patty")


treatments_new <- ggplot(full_pollen, aes(x=age, y=visual_avg, color=Treatment))+geom_jitter(height=0.075, alpha=0.2)+
  facet_wrap(~Treatment, labeller = as_labeller(new_labels))+
  geom_ribbon(data=emm_vis_treatment, aes(x= age, ymin=LCL, ymax=UCL, group=Treatment), alpha=0.1, inherit.aes = FALSE)+
  geom_smooth(data=emm_vis_treatment, aes(x=age, y=yvar), se=FALSE, color = "black")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), panel.background = element_blank (), strip.background = element_blank(), panel.grid.minor= element_blank(), panel.grid.major= element_blank(),
        text = element_text(size=14), legend.position="none",
        axis.title.x = element_text(margin = margin(t = 8)), # Add space above x-axis title
        axis.title.y = element_text(margin = margin(r = 8)))+
  labs(x="Adult worker age (day)", y="Dye presence")+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + scale_y_continuous(breaks = scales::pretty_breaks(n = 1))+
  scale_color_manual(values=c("#AADC32FF","#27AD81FF", "#3B528BFF", "#440154FF"))

treatments_new
