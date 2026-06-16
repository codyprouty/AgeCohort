#Kleckner et al. 2026 - DOI:10.1093/jee/toag155
#Script 2/2 - Observations of marked workers

####1) Set up ####

##load in packAges
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

#####1.1) Bee counts ####
bee_counts <- read.csv("BeeCounts.csv")
str(bee_counts) #check structure

#metadata
  #Date - calendar date of observations
  #SampleDay - count of sampling day (one less than worker Age)
  #Treatment - diet treatments: treatment (dyed PSP), NControl (undyed PSP), CControl (modeling clay), and Pcontrol (dyed fondant)
  #Colony - unique colony number
  #Bees - count of marked worker bees interacting with diet treatment
  #Area - area left of diet treatment
  #Notes - comments on diet replacement, etc.
  #Removed - count of marked workers removed for dye sampling
  #Background - count of marked workers removed from background mortality
  #BeesTotal - estimate of total marked bee population in each colony

##create column for worker Age (one more than sampling day)
names(bee_counts)[1] <- "Date"
bee_counts$Age <- bee_counts$SampleDay +1

##make colony and treatment a factor
bee_counts$Colony <- as.factor(bee_counts$Colony)
bee_counts$Treatment <- as.factor(bee_counts$Treatment)

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


##summarize from sunrise to the averAge of our sampling time, 7:00 - 11:00
pollen_env1 <- rbind(subset(pollen_env, Time == c("7:00:00") | Time == c("8:00:00") | Time == c("9:00:00") | Time == c("10:00:00") | Time == c("11:00:00")))

##remove unnecessary columns and averAge weather across sampling window
pollen_env_avg2 <- pollen_env1 %>% dplyr::select(-Station.ID, -Station.Name, -Time, -True_Date) %>% 
  group_by(Date, Period) %>% 
  summarise(across(everything(), mean, .names = "Avg_{.col}"), .groups = "drop") %>%
  filter(Period < 31)

head(pollen_env_avg2)

##change rainfall from inches to cm
pollen_env_avg2$Avg_Rainfall <- pollen_env_avg2$Avg_Rainfall*2.54

##scale environmental variables
pollen_env_avg2[c("Avg_Temp60cm", "Avg_Rainfall", "Avg_Humidity", "Avg_Radiation")] <- scale(pollen_env_avg2[c("Avg_Temp60cm", "Avg_Rainfall", "Avg_Humidity", "Avg_Radiation")])

##cbind counts and environmental variables using dates
full_bee_counts <- left_join(bee_counts, pollen_env_avg2, by = "Date") 

head(full_bee_counts)



####2) Model construction####

##remove controls for modeling, keeping dyed and undyed PSPs
full_bee_counts2 <- subset(na.omit(full_bee_counts), Treatment != "CControl" & Treatment != "PControl")

head(full_bee_counts2)
tail(full_bee_counts2)

##check correlations
count_cor2 <- full_bee_counts2 %>% dplyr::select(-Date, -Treatment, -Area, -Notes, -Colony)  #remove non-numeric columns
count_cor3 <- cor(as.matrix(count_cor2)) ##make correlation matrix
corrplot(count_cor3, method = "circle", type = "upper") ## plots strength of correlation as color-coded circles)
###some environmental variables (e.g. Humidity and Radiation) are correlated with each other, but none with the visual averAges


##start with most complex full model
cp1 <- glmmTMB(Bees ~ Age + I(Age**2) + offset(log(BeesTotal)) +
                 Avg_Temp60cm + Avg_Rainfall + Avg_Radiation + Avg_Humidity +
                 ar1(0 + as.factor(Age)|Colony), data=full_bee_counts2, family="nbinom2") ##quadratic term and AR1 random effect

##AR1 correlation structure accounts for correlation of points in time closer together

cp2 <- glmmTMB(Bees ~ Age + I(Age**2)  + offset(log(BeesTotal)) + 
                 Avg_Temp60cm + Avg_Rainfall + Avg_Radiation +
                 ar1(0 + as.factor(Age)|Colony), data=full_bee_counts2, family="nbinom2") ##pairwise by weather variables

cp3 <- glmmTMB(Bees ~ Age + I(Age**2) + offset(log(BeesTotal)) + 
                 Avg_Temp60cm + Avg_Rainfall + Avg_Humidity +
                 ar1(0 + as.factor(Age)|Colony), data=full_bee_counts2, family="nbinom2")

cp4 <- glmmTMB(Bees ~ Age + I(Age**2) + offset(log(BeesTotal)) + 
                 Avg_Rainfall + Avg_Radiation + Avg_Humidity +
                 ar1(0 + as.factor(Age)|Colony), data=full_bee_counts2, family="poisson")


cp5 <- glmmTMB(Bees ~ Age + I(Age**2) + offset(log(BeesTotal)) + 
                 Avg_Temp60cm + Avg_Radiation + Avg_Humidity +
                 ar1(0 + as.factor(Age)|Colony), data=full_bee_counts2, family="poisson")

cp6 <- glmmTMB(Bees ~ Age + I(Age**2) + offset(log(BeesTotal)) + 
                 Avg_Temp60cm + Avg_Rainfall +
                 ar1(0 + as.factor(Age)|Colony), data=full_bee_counts2, family="nbinom2")

cp7 <- glmmTMB(Bees ~ Age + I(Age**2) + offset(log(BeesTotal)) + 
                 Avg_Temp60cm + Avg_Humidity +
                 ar1(0 + as.factor(Age)|Colony), data=full_bee_counts2, family="nbinom2")

cp8 <- glmmTMB(Bees ~ Age + I(Age**2) + offset(log(BeesTotal)) + 
                 Avg_Temp60cm + Avg_Radiation +
                 ar1(0 + as.factor(Age)|Colony), data=full_bee_counts2, family="nbinom2")

cp9 <- glmmTMB(Bees ~ Age + I(Age**2) + offset(log(BeesTotal)) + 
                 Avg_Rainfall + Avg_Radiation +
                 ar1(0 + as.factor(Age)|Colony), data=full_bee_counts2, family="nbinom2")

cp10 <- glmmTMB(Bees ~ Age + I(Age**2) + offset(log(BeesTotal)) + 
                  Avg_Rainfall + Avg_Humidity +
                  ar1(0 + as.factor(Age)|Colony), data=full_bee_counts2, family="nbinom2")

cp11 <- glmmTMB(Bees ~ Age + I(Age**2) + offset(log(BeesTotal)) + 
                  Avg_Radiation + Avg_Humidity +
                  ar1(0 + as.factor(Age)|Colony), data=full_bee_counts2, family="nbinom2")

cp12 <- glmmTMB(Bees ~ Age + I(Age**2) + offset(log(BeesTotal)) + 
                  Avg_Humidity +
                  ar1(0 + as.factor(Age)|Colony), data=full_bee_counts2, family="nbinom2")

cp13 <- glmmTMB(Bees ~ Age + I(Age**2) + offset(log(BeesTotal)) + 
                  Avg_Radiation +
                  ar1(0 + as.factor(Age)|Colony), data=full_bee_counts2, family="nbinom2")

cp14 <- glmmTMB(Bees ~ Age + I(Age**2) + offset(log(BeesTotal)) + 
                  Avg_Temp60cm +
                  ar1(0 + as.factor(Age)|Colony), data=full_bee_counts2, family="nbinom2")

cp15 <- glmmTMB(Bees ~ Age + I(Age**2) + offset(log(BeesTotal)) + 
                  Avg_Rainfall +
                  ar1(0 + as.factor(Age)|Colony), data=full_bee_counts2, family="nbinom2")

cp16 <- glmmTMB(Bees ~ Age + I(Age**2) + offset(log(BeesTotal)) + 
                  ar1(0 + as.factor(Age)|Colony), data=full_bee_counts2, family="nbinom2") ##no weather variables

cp17 <- glmmTMB(Bees ~ Age + I(Age**2) + 
                  ar1(0 + as.factor(Age)|Colony), data=full_bee_counts2, family="nbinom2") ##no mortality offset

cp18 <- glmmTMB(Bees ~ Age + 
                  ar1(0 + as.factor(Age)|Colony), data=full_bee_counts2, family="poisson")##no quadratic term

cp19 <- glmmTMB(Bees ~ Age + (1 | Colony)
                , data=full_bee_counts2, family="nbinom2") ##no AR1 correlation structure

cp20 <- glmmTMB(Bees ~ Avg_Temp60cm + Avg_Rainfall + Avg_Radiation + Avg_Humidity + offset(log(BeesTotal)) + 
                  ar1(0 + as.factor(Age)|Colony), data=full_bee_counts2) ##no Age, only weather variables

cp21 <- glmmTMB(Bees ~ 1, data=full_bee_counts2, family="nbinom2") ##null




####3) Model selection####

##generate AICc scores, weights, and delta
Model.sel <- model.sel(cp1, cp2, cp3, cp4, cp5, cp6, cp7, cp8, cp9, cp10, cp11, cp12, cp13, cp14, cp15, cp16, cp17, cp18, cp19, cp20, cp21)

Model.sel
##8 of 21 models within 6 AICc, will need to averAge
##all 8 have quadratic term and AR1 random effect

##averAge model
cpollen_models_subset <- model.sel(cp1, cp2, cp3, cp6) ##delaa threshold = 6

model.avg(cpollen_models_subset)

summary(model.avg(cpollen_models_subset)) ##conditional avg does not include zeros for absent terms
##in the averAged model, the effect sizes were small, so the environmental variables contributed very little


##use top model (no env variables) for rest of analyses/visualization > most parsimonious

#assess residuals of selected model
plot(simulateResiduals(cp1)) 
hist(simulateResiduals(cp1))

#random effect components
summary(cp1)

##emmeans
###estimates, effect sizes, compare means across Age
emmeans(cp1, ~1, type='response', at=list(Age=2)) #0.0823 marked bees
emmeans(cp1, ~1, type='response', at=list(Age=4)) #0.189 marked bees
emmeans(cp1, ~1, type='response', at=list(Age=6)) #0.391 marked bees
emmeans(cp1, ~1, type='response', at=list(Age=8)) #0.732 marked bees
emmeans(cp1, ~1, type='response', at=list(Age=10)) #1.24 marked bees
emmeans(cp1, ~1, type='response', at=list(Age=12)) #1.89 marked bees
emmeans(cp1, ~1, type='response', at=list(Age=14)) #2.61 marked bees
emmeans(cp1, ~1, type='response', at=list(Age=16)) #3.25 marked bees
emmeans(cp1, ~1, type='response', at=list(Age=18)) #3.65 marked bees
emmeans(cp1, ~1, type='response', at=list(Age=19)) #3.73 marked bees >>Peak
emmeans(cp1, ~1, type='response', at=list(Age=20)) #3.71 marked bees
emmeans(cp1, ~1, type='response', at=list(Age=21)) #3.6 marked bees
emmeans(cp1, ~1, type='response', at=list(Age=22)) #3.41 marked bees
emmeans(cp1, ~1, type='response', at=list(Age=24)) #2.83 marked bees
emmeans(cp1, ~1, type='response', at=list(Age=26)) #2.12 marked bees
emmeans(cp1, ~1, type='response', at=list(Age=28)) #1.43 marked bees


####4) Treatment effects####

#counts of observations per treatment
treatment_counts <- bee_counts %>%
  group_by(Treatment) %>%
  summarise(count = n())

treatment_counts


#visualize all treatments
ggplot(bee_counts, aes(x=SampleDay, y=Bees, color=Treatment))+geom_jitter(height=0.1, alpha=0.5)+
  facet_wrap(~Treatment)+geom_smooth(color="black", method="glm", formula = y ~ x+I(x^2))

##remove modeling clay, all zero data
bee_counts2 <- bee_counts %>% 
  filter(Treatment != "CControl")

##omit NAs
bee_counts2 <- na.omit(bee_counts2)

##round to whole bee
bee_counts2$TotalRound <- round(bee_counts2$BeesTotal, digits=0)


##build full model w/offset and AR1 random effect
obs_treatment <- glmmTMB(Bees ~ Age + I(Age**2) + Treatment + offset(log(BeesTotal)) + 
                            ar1(0 + as.factor(Age)|Colony), ##for autocorrelation
                          data=bee_counts2, family="poisson")

check_overdispersion(obs_treatment) #none
check_zeroinflation(obs_treatment) #ok
check_singularity(obs_treatment) #true, simplify random effect

#model 2 - no AR1
obs_treatment2 <- glmmTMB(Bees ~ Age + I(Age**2) + Treatment + offset(log(BeesTotal)) + 
                            (1|Colony), data=bee_counts2, family="poisson")

check_overdispersion(obs_treatment2) #now detected, switch to neg bin
check_singularity(obs_treatment2) #false, fixed

#model 3 - neg binomial
obs_treatment3 <- glmmTMB(Bees ~ Age*Treatment + I(Age**2) * Treatment + offset(log(BeesTotal)) + 
                            (1|Colony), data=bee_counts2, family="nbinom2")

check_overdispersion(obs_treatment3) #fixed
check_zeroinflation(obs_treatment3) #ok
check_singularity(obs_treatment3) #still false

##check residuals with DHarma
plot(simulateResiduals(obs_treatment3)) 
hist(simulateResiduals(obs_treatment3)) 
hist(residuals(obs_treatment3))

Anova(obs_treatment3)

#remove interaction btwn age and treatment
obs_treatment4 <- glmmTMB(Bees ~ Age + I(Age^2) * Treatment + offset(log(BeesTotal)) + 
                            (1|Colony), data=bee_counts2, family="nbinom2")

#residuals
plot(simulateResiduals(obs_treatment4))
hist(simulateResiduals(obs_treatment4)) 
hist(residuals(obs_treatment4)) 

#significance
Anova(obs_treatment4)

emmeans(obs_treatment4, pairwise~Treatment, type = "response", at=list(Age=22)) #p<0.05 at day 23 btwn dyed and undyed PSP
emmeans(obs_treatment4, pairwise~Treatment, type = "response", at=list(Age=23))
emmeans(obs_treatment4, pairwise~Treatment, type = "response", at=list(Age=24))
emmeans(obs_treatment4, pairwise~Treatment, type = "response", at=list(Age=25))
emmeans(obs_treatment4, pairwise~Treatment, type = "response", at=list(Age=26))
emmeans(obs_treatment4, pairwise~Treatment, type = "response", at=list(Age=27))
emmeans(obs_treatment4, pairwise~Treatment, type = "response", at=list(Age=28))
emmeans(obs_treatment4, pairwise~Treatment, type = "response", at=list(Age=29))
#p<0.05 at day 24-29 btwn undyed PSP and both dyed PSP and fondant


#for plotting
emm_obs_treatment <- emmip(obs_treatment4, Treatment~Age, type='response', CIs = TRUE, 
                           at=list(Age = c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29)), plotit = FALSE )


####5) Visualization####

#####5.1) Observations across worker age (Fig 1)####

##emmeans from 
call <- as.data.frame(emmeans(cp1, ~1, type='response', at=list(Age=1))) 
i<- 2
while(i<30){
  call[nrow(call) + 1,] = as.data.frame(emmeans(cp1, ~1, type='response', at=list(Age=i)))
  i = i+1
}

##add age and description for binding later
call$age <- c(1:29)
call$hive <- "All"


##two separate panels
##panel A, using emmeans above
fig2a <- ggplot()+
  geom_ribbon(data=call, aes(x= age, ymin=asymp.LCL, ymax=asymp.UCL), alpha=0.1)+
  geom_line(data=call, aes(x=age, y=response), color="black", size=2.5)+
  theme_bw()+ theme(strip.background = element_blank(), panel.grid.minor= element_blank(), panel.grid.major= element_blank(),
                    text = element_text(size=14),
                    axis.title.x = element_text(margin = margin(t = 8)), # Add space above x-axis title
                    axis.title.y = element_text(margin = margin(r = 8)))+
  labs(x="Adult worker age (day)", y="Marked workers observed on PSPs")+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  annotate("text", x = 1.5, y = 5, label = "A", color = "black", size = 7)

fig2a

ggsave("AgeCohort_Fig2a.png", plot = fig2a, width = 4, height = 4, units = "in", dpi = 600)


##panel b
##using the predict function to calculate the model's equation for each point

full_bee_counts2$pred <- exp(predict(cp1))/full_bee_counts2$BeesTotal

dyed <- full_bee_counts2 %>% filter(Treatment == "treatment")
undyed <- full_bee_counts2 %>% filter(Treatment == "NControl")

fig2b <- ggplot()+
  geom_smooth(data=dyed, aes(x=Age, y=pred, color=Colony, linetype= Colony), size=1, linetype="73", 
              method="glm", method.args=list(family="binomial"), se=FALSE, formula = y ~ x+I(x^2))+
  geom_smooth(data=undyed, aes(x=Age, y=pred, color=Colony, linetype= Colony), size=1, linetype="1343", 
              method="glm", method.args=list(family="binomial"), se=FALSE, formula = y ~ x+I(x^2))+
  theme_bw()+ theme(strip.background = element_blank(), panel.grid.minor= element_blank(), panel.grid.major= element_blank(),
                    text = element_text(size=14),
                    legend.position = "none",
                    axis.title.x = element_text(margin = margin(t = 8)), # Add space above x-axis title
                    axis.title.y = element_text(margin = margin(r = 8)) )+
  labs(x="Adult worker age (day)", y="Proportion of marked workers on PSPs")+
  scale_color_viridis(discrete = T, direction = -1)+
  #scale_color_manual(values=c( "#2C728EFF", "#5DC863FF","#3B528BFF", '#472D7BFF', "#27AD81FF","#FDE725FF", "#440154FF","#21908CFF", "#AADC32FF"))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  annotate("text", x = 2.35, y = 0.0295, label = "B", color = "black", size = 7)

fig2b

ggsave("AgeCohort_Fig2b.png", plot = fig2b, width = 4.3, height = 4, units = "in", dpi = 600)



#####5.2) Environmental variables (supplemental)####

##using unscaled env data
pollen_env_avg_US <- pollen_env %>% select(-Station.ID, -Station.Name, -Time, -True_Date) %>%
  group_by(Period, Date) %>%
  summarise(across(everything(), mean, .names = "Avg_{.col}"), .groups = "drop") %>%
  filter(Period < 31)

unscaled_counts <- left_join(bee_counts, pollen_env_avg_US, by = "Date") 
head(unscaled_counts)

#prep data all over, didn't save separate scaled version in initial set up
names(unscaled_counts)[1] <- "Date"
unscaled_counts$Age <- unscaled_counts$SampleDay +1
unscaled_counts <- subset(na.omit(unscaled_counts), Treatment != "CControl" & Treatment != "PControl")
unscaled_counts$Avg_Rainfall <- unscaled_counts$Avg_Rainfall*2.54 #inches to cm


##extract residuals from model without one weather variable and plot residuals against said weather variable

##temp
unscaled_counts$residT <-residuals(cp4)

png("Temperature.png", width=8, height=6, units="in", res=400)
ggplot(unscaled_counts, aes(x=Avg_Temp60cm, y=residT))+
  theme_bw() + geom_jitter(unscaled_counts, mapping=aes(x=Avg_Temp60cm, y=residT, color=Treatment), height=0.1, alpha=0.8, size=1.75)+
  geom_smooth(color="black", method="glm", size=2)+
  theme(strip.background = element_blank(), panel.grid.minor= element_blank(), panel.grid.major= element_blank(),
        text = element_text(size=16), legend.position="none")+ 
  labs(x=expression(paste('Average Daily Temperature (',~degree,'C)',sep='')), y="Model Residuals")+
  scale_color_viridis(discrete = T, direction = -1, begin= .85, end=.85, alpha=0.08)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) 
dev.off()

##rain
unscaled_counts$residR <-residuals(cp5)

png("Rainfall.png", width=8, height=6, units="in", res=400)
ggplot(unscaled_counts, aes(x=Avg_Rainfall, y=residR))+
  theme_bw() + geom_jitter(unscaled_counts, mapping=aes(x=Avg_Rainfall, y=residR, color=Treatment), height=0.1, alpha=0.8, size=1.75)+
  geom_smooth(color="black", method="glm", size=2)+ ylim(-2.5,5)+
  theme(strip.background = element_blank(), panel.grid.minor= element_blank(), panel.grid.major= element_blank(),
        text = element_text(size=16), legend.position="none")+
  labs(x="Average Daily Rainfall (cm)", y="Model Residuals")+
  scale_color_viridis(discrete = T, direction = -1, begin= .55, end=.55, alpha=0.08)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5))  #fix y scale in ppt 
dev.off()

##radiation
unscaled_counts$residRd <-residuals(cp3)

png("Radiation.png", width=8, height=6, units="in", res=400)
ggplot(unscaled_counts, aes(x=Avg_Radiation, y=residRd))+
  theme_bw() + geom_jitter(unscaled_counts, mapping=aes(x=Avg_Radiation, y=residRd, color=Treatment), height=0.1, alpha=0.8, size=1.75)+
  geom_smooth(color="black", method="glm", size=2)+ 
  theme(strip.background = element_blank(), panel.grid.minor= element_blank(), panel.grid.major= element_blank(),
        text = element_text(size=16), legend.position="none")+
  labs(x=expression(Average~Daily~Radiation~(W/m^{"2"})), y="Model Residuals")+ ##use expression for superscript
  scale_color_viridis(discrete = T, direction = -1, begin= .25, end=.25, alpha=0.08)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) 
dev.off()

#humidity
unscaled_counts$residH <-residuals(cp2)

png("Humidity.png", width=8, height=6, units="in", res=400)
ggplot(unscaled_counts, aes(x=Avg_Humidity, y=residH))+
  theme_bw() + geom_jitter(unscaled_counts, mapping=aes(x=Avg_Humidity, y=residH, color=Treatment), height=0.001, alpha=0.8, size=1.75)+
  geom_smooth(color="black", method="glm", size=2)+ 
  theme(strip.background = element_blank(), panel.grid.minor= element_blank(), panel.grid.major= element_blank(),
        text = element_text(size=16), legend.position="none")+ 
  labs(x="Average Daily Humidity (%)", y="Model Residuals")+
  scale_color_viridis(discrete = T, direction = -1, begin= 0, end=0, alpha=0.08)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5))
dev.off()



#####5.3) Treatment effects (Fig 2)####

##calculate proportion of marked bees interacting with diets/total marked bee population
  ##mimics model offset for figures
##round to whole bee
bee_counts$TotalRound <- round(bee_counts$BeesTotal, digits=0)
bee_counts$prop <- (bee_counts$Bees/bee_counts$TotalRound)

#change treatment names
bee_counts_fig <- bee_counts

unique(bee_counts_fig$Treatment)
bee_counts_fig$Treatment <- factor(bee_counts_fig$Treatment, levels=c("CControl", "NControl", "PControl", "treatment"),
                               labels= c("Dyed modeling clay", "Undyed pollen substitute patty", "Dyed fondant", "Dyed pollen substitute patty"))

unique(emm_obs_treatment$Treatment)
emm_obs_treatment$Treatment <- factor(emm_obs_treatment$Treatment, levels=c("CControl", "NControl", "PControl", "treatment"),
                               labels= c("Dyed modeling clay", "Undyed pollen substitute patty", "Dyed fondant", "Dyed pollen substitute patty"))


treatments_fig <- ggplot(bee_counts_fig, aes(x=Age, y=Bees, color=Treatment))+geom_jitter(height=0.005, alpha=0.8)+
  facet_wrap(~Treatment)+
  geom_ribbon(data=emm_obs_treatment, aes(x= Age, ymin=LCL, ymax=UCL, group=Treatment), alpha=0.1, inherit.aes = FALSE)+
  geom_smooth(data=emm_obs_treatment, aes(x=Age, y=yvar), se=FALSE, color = "black")+
  #geom_smooth(color="black", method="glm", method.args=list(family="binomial"), se=FALSE, formula = y ~ x+I(x^2))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), panel.background = element_blank (), strip.background = element_blank(), panel.grid.minor= element_blank(), panel.grid.major= element_blank(),
        text = element_text(size=14), legend.position="none",
        axis.title.x = element_text(margin = margin(t = 8)), # Add space above x-axis title
        axis.title.y = element_text(margin = margin(r = 8)))+
  labs(x="Adult worker age (day)", y="Marked workers observed\ninteracting with treatment")+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
  scale_color_manual(values=c("#AADC32FF","#27AD81FF", "#3B528BFF", "#440154FF"))

treatments_fig
