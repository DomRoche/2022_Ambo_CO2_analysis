

#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
install.packages("Rtools40")
# packages
library(lme4)
library(pander)
library(kableExtra)
library(ggplot2)
library(emmeans)
library(here)
library(lmerTest)
library(brms)
library(cmdstanr)
library(tidyverse)       
library(tidybayes)       
library(broom)          
library(extraDistr)       
library(ggdist)           
library(gghalves)        
library(patchwork) 
library(rptR)
library(Rtools40)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)

# loading functions
# a function to get prediction from the model
predict_dat <- function(model, dat){
  model %>% 
    emmeans(~ trial,
            by = "group",
            at = list(trial = seq(min(dat$trial),max(dat$trial),length.out = 100)),
            epred = TRUE,
            re_formula = NA) %>% as_tibble() -> data
  data
}

# a function to get posterior distributions

post_dat <- function(model, dat){
  new_mat <- expand_grid(
    group = c("control", "CO2"), 
    trial =  seq(min(dat$trial),max(dat$trial),length.out = 100), 
    SL = mean(dat$SL))
  
  epred_draws(model, newdata = new_mat,
              dpar = c("mu", "sigma"), 
              re_formula = NA
  ) %>% 
    group_by(trial, .draw) %>% summarise(trial = mean(trial),
                                         predicted = mean(.epred),
                                         mu = mean(mu),
                                         sigma = mean(sigma),
                                         group = group) 
}


# function to compare two rptR objects 

comp_rpt <- function(treat, cont){
  treat <- unlist(treat[["R_boot"]])
  cont <- unlist(cont[["R_boot"]])
  contrast <- treat - cont
  quant <- quantile(contrast, c(0.025, 0.5, 0.975))
  quant
}


# data
# data
full_dat <- read.csv(here("Dom_analysis", "Ambo_ind_CO2_data.csv"), header=T)
full_dat$fishID <- as.factor(full_dat$fishID)

# omitting missing data 
dat <- na.omit(full_dat)
dat$group <- as.factor(dat$group)
# looking at distributions

par(mfrow=c(3,2),mar=c(2.8,2,1,3),oma=c(2,2,3,0),las=1); hist(dat$activity); hist(dat$shelter); hist(dat$thigmotaxis); hist(dat$novel); hist(dat$emergT)

# looking at residuals after transformaiton

# activity
m_activity <- lmer(sqrt(activity) ~  1 + group*trial + SL +
                     (1 + trial|fishID),
                   data = dat)

#shelter
m_shelter <- lmer(shelter ~  1 + group*trial + SL +
                    (1 + trial|fishID),
                  data = dat)

#thigmotaxis
m_thigmotaxis <- lmer(sqrt(thigmotaxis) ~  1 + group*trial + SL +
                        (1 + trial|fishID),
                      data = dat)

#novel
m_novel <- lmer(log(novel+0.5) ~  1 + group*trial + SL +
                  (1 + trial|fishID),
                data = dat)

# emergT
m_emergT <- lmer(log(emergT) ~  1 + group*trial + SL +
                   (1 + trial|fishID),
                 data = dat)

par(mfrow=c(3,2),mar=c(2.8,2,1,3),oma=c(2,2,3,0),las=1); 
hist(residuals(m_activity)); hist(residuals(m_shelter)); hist(residuals(m_thigmotaxis));hist(residuals(m_novel)); hist(residuals(m_emergT))

#############
# modeling
##############

# activity
mod_activity <- bf(sqrt(activity) ~  1 + group*trial + SL +
                     (1 + I(trial-1)|p|fishID),
                   sigma ~ 1+ group*trial)

model_activity <- brm(mod_activity, data = dat,
                      chains = 2, cores = 2, iter = 10000, warmup = 7000,
                      backend = "cmdstanr")

summary(model_activity)

#saveRDS(model_activity, here("Dom_analysis","model_activity.Rds"))

model_activity <- readRDS(here("Dom_analysis", "model_activity.Rds"))

# shelter
mod_shelter <- bf(shelter ~  1 + group*trial + SL +
                    (1 + I(trial-1)|p|fishID),
                  sigma ~ group*trial)

model_shelter <- brm(mod_shelter, data = dat,
                     chains = 2, cores = 2, iter = 10000, warmup = 7000,
                     backend = "cmdstanr")

summary(model_shelter)

saveRDS(model_shelter, here("Dom_analysis","model_shelter.Rds"))

# thigmotaxis

mod_thigmotaxis <- bf(sqrt(thigmotaxis) ~  1 + group*trial + SL +
                        (1 + I(trial-1)|p|fishID),
                      sigma ~ group*trial)

model_thigmotaxis <- brm(mod_thigmotaxis, data = dat,
                         chains = 2, cores = 2, iter = 10000, warmup = 7000,
                         backend = "cmdstanr")

summary(model_thigmotaxis)

saveRDS(model_thigmotaxis, here("Dom_analysis","model_thigmotaxis.Rds"))

# novel
mod_novel <- bf(log(novel+0.5) ~  1 + group*trial + SL +
                  (1 + I(trial-1)|p|fishID),
                sigma ~ group*trial)

model_novel <- brm(mod_novel, data = dat,
                   chains = 2, cores = 2, iter = 10000, warmup = 7000,
                   backend = "cmdstanr")

summary(model_novel)

saveRDS(model_novel, here("Dom_analysis","model_novel.Rds"))

# emergT
mod_emergT <- bf(log(emergT) ~  1 + group*trial + SL +
                   (1 + I(trial-1)|p|fishID),
                 sigma ~ group*trial)


model_emergT <- brm(mod_emergT, data = dat,
                    chains = 2, cores = 2, iter = 10000, warmup = 7000,
                    backend = "cmdstanr", control = list(adapt_delta = 0.99))

summary(model_emergT)

saveRDS(model_emergT, here("Dom_analysis","model_emergT.Rds"))

#################
# plotting
#################
# Model data 
model_activity <- readRDS(here("Dom_analysis", "model_activity.Rds"))
model_shelter <- readRDS(here("Dom_analysis", "model_shelter.Rds"))
model_thigmotaxis <- readRDS(here("Dom_analysis", "model_thigmotaxis.Rds"))
model_novel <- readRDS(here("Dom_analysis", "model_novel.Rds"))
model_emergT <- readRDS(here("Dom_analysis", "model_emergT.Rds"))

# test
#model_shelter <- repair_stanfit_names(model_shelter)

# creating predictions
pred_activity <- predict_dat(model_activity, dat)
pred_shelter <- predict_dat(model_shelter, dat)
pred_thigmotaxis <- predict_dat(model_thigmotaxis, dat)
pred_novel <- predict_dat(model_novel, dat)
pred_emergT <- predict_dat(model_emergT, dat)

# creating posteriors
post_activity <- post_dat(model_activity, dat)
post_shelter <- post_dat(model_shelter, dat)
post_thigmotaxis <- post_dat(model_thigmotaxis, dat)
post_novel <- post_dat(model_novel, dat)
post_emergT <- post_dat(model_emergT, dat)

# 1: activity

p1 <- ggplot()+
  # confidence interval
  geom_smooth(data = pred_activity, aes(x = trial, y = lower.HPD), method =  "loess", formula = y~x, se = FALSE, lty =  "dotted", lwd = 1,  col = "black") +
  geom_smooth(data = pred_activity, aes(x = trial, y = upper.HPD), method =  "loess", formula = y~x, se = FALSE, lty = "dotted", lwd = 1,  col = "black") +
  # main line
  geom_smooth(data = pred_activity, aes(x = trial, y = emmean), method =  "loess", formula = y~x, se = FALSE, lwd = 1.5, col = "black") +
  facet_wrap(vars(group), ncol = 2) +
  # real data
  geom_point(data = dat, aes(x = trial, y = sqrt(activity), fill = fishID, color = fishID), shape = 21, alpha = 0.5) +
  geom_smooth(data = dat,  aes(x = trial, y = sqrt(activity), color = fishID), method = "lm", formula = y ~x, se = FALSE, lwd = 0.1) + 
  scale_fill_viridis_d(end = 0.8) +
  scale_colour_viridis_d(end = 0.8) +
  labs(y = "sqrt(activity)", x = "trial (time)") +
  # themes
  theme_bw(base_size = 12) +  theme(legend.position = "none")

s1 <- ggplot(post_activity, aes(x = trial, y = sigma, fill = group, colour = group)) +
  stat_lineribbon(aes(fill_ramp = stat(level))) +
  scale_fill_viridis_d(end = 0.8) +
  scale_color_viridis_d(end = 0.8) +
  scale_fill_ramp_discrete(range = c(0.2, 0.7)) +
  facet_wrap(vars(group), ncol = 2) +
  #labeller = labeller(quota = c(`TRUE` = "Quota",
  #                             `FALSE` = "No Quota"))) +
  labs(y = "log(sigma_activity)", x = "trial (time)",
       #fill = "group", color = "group",
       fill_ramp = "Credible interval") +
  #ylim(0, 0.30) + xlim(24.5, 32.5) +
  guides(fill = "none", color = "none") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

# 3: shelter

p2 <- ggplot()+
  # confidence interval
  geom_smooth(data = pred_shelter, aes(x = trial, y = lower.HPD), method =  "loess", formula = y~x, se = FALSE, lty =  "dotted", lwd = 1,  col = "black") +
  geom_smooth(data = pred_shelter, aes(x = trial, y = upper.HPD), method =  "loess", formula = y~x, se = FALSE, lty = "dotted", lwd = 1,  col = "black") +
  # main line
  geom_smooth(data = pred_shelter, aes(x = trial, y = emmean), method =  "loess", formula = y~x, se = FALSE, lwd = 1.5, col = "black") +
  facet_wrap(vars(group), ncol = 2) +
  # real data
  geom_point(data = dat, aes(x = trial, y = shelter, fill = fishID, color = fishID), shape = 21, alpha = 0.5) +
  geom_smooth(data = dat,  aes(x = trial, y = shelter, color = fishID), method = "lm", formula = y ~x, se = FALSE, lwd = 0.1) + 
  scale_fill_viridis_d(end = 0.8) +
  scale_colour_viridis_d(end = 0.8) +
  labs(y = "shelter", x = "trial (time)") +
  # themes
  theme_bw(base_size = 12) +  theme(legend.position = "none")

s2 <- ggplot(post_shelter, aes(x = trial, y = sigma, fill = group, colour = group)) +
  stat_lineribbon(aes(fill_ramp = stat(level))) +
  scale_fill_viridis_d(end = 0.8) +
  scale_color_viridis_d(end = 0.8) +
  scale_fill_ramp_discrete(range = c(0.2, 0.7)) +
  facet_wrap(vars(group), ncol = 2) +
  #labeller = labeller(quota = c(`TRUE` = "Quota",
  #                             `FALSE` = "No Quota"))) +
  labs(y = "log(sigma_shelter)", x = "trial (time)",
       #fill = "group", color = "group",
       fill_ramp = "Credible interval") +
  #ylim(0, 0.30) + xlim(24.5, 32.5) +
  guides(fill = "none", color = "none") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

# 3: thigmotaxis

p3 <- ggplot()+
  # confidence interval
  geom_smooth(data = pred_thigmotaxis, aes(x = trial, y = lower.HPD), method =  "loess", formula = y~x, se = FALSE, lty =  "dotted", lwd = 1,  col = "black") +
  geom_smooth(data = pred_thigmotaxis, aes(x = trial, y = upper.HPD), method =  "loess", formula = y~x, se = FALSE, lty = "dotted", lwd = 1,  col = "black") +
  # main line
  geom_smooth(data = pred_thigmotaxis, aes(x = trial, y = emmean), method =  "loess", formula = y~x, se = FALSE, lwd = 1.5, col = "black") +
  facet_wrap(vars(group), ncol = 2) +
  # real data
  geom_point(data = dat, aes(x = trial, y = sqrt(thigmotaxis), fill = fishID, color = fishID), shape = 21, alpha = 0.5) +
  geom_smooth(data = dat,  aes(x = trial, y = sqrt(thigmotaxis), color = fishID), method = "lm", formula = y ~x, se = FALSE, lwd = 0.1) + 
  scale_fill_viridis_d(end = 0.8) +
  scale_colour_viridis_d(end = 0.8) +
  labs(y = "sqrt(thigmotaxis)", x = "trial (time)") +
  # themes
  theme_bw(base_size = 12) +  theme(legend.position = "none")

s3 <- ggplot(post_thigmotaxis, aes(x = trial, y = sigma, fill = group, colour = group)) +
  stat_lineribbon(aes(fill_ramp = stat(level))) +
  scale_fill_viridis_d(end = 0.8) +
  scale_color_viridis_d(end = 0.8) +
  scale_fill_ramp_discrete(range = c(0.2, 0.7)) +
  facet_wrap(vars(group), ncol = 2) +
  #labeller = labeller(quota = c(`TRUE` = "Quota",
  #                             `FALSE` = "No Quota"))) +
  labs(y = "log(sigma_thigmotaxis)", x = "",
       #fill = "group", color = "group",
       fill_ramp = "Credible interval") +
  #ylim(0, 0.30) + xlim(24.5, 32.5) +
  guides(fill = "none", color = "none") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

# 4: novel

p4 <- ggplot()+
  # confidence interval
  geom_smooth(data = pred_novel, aes(x = trial, y = lower.HPD), method =  "loess", formula = y~x, se = FALSE, lty =  "dotted", lwd = 1,  col = "black") +
  geom_smooth(data = pred_novel, aes(x = trial, y = upper.HPD), method =  "loess", formula = y~x, se = FALSE, lty = "dotted", lwd = 1,  col = "black") +
  # main line
  geom_smooth(data = pred_novel, aes(x = trial, y = emmean), method =  "loess", formula = y~x, se = FALSE, lwd = 1.5, col = "black") +
  facet_wrap(vars(group), ncol = 2) +
  # real data
  geom_point(data = dat, aes(x = trial, y = log(novel+0.5), fill = fishID, color = fishID), shape = 21, alpha = 0.5) +
  geom_smooth(data = dat,  aes(x = trial, y = log(novel+0.5), color = fishID), method = "lm", formula = y ~x, se = FALSE, lwd = 0.1) + 
  scale_fill_viridis_d(end = 0.8) +
  scale_colour_viridis_d(end = 0.8) +
  labs(y = "log(novel+0.5)", x = "trial (time)") +
  # themes
  theme_bw(base_size = 12) +  theme(legend.position = "none")

s4 <- ggplot(post_novel, aes(x = trial, y = sigma, fill = group, colour = group)) +
  stat_lineribbon(aes(fill_ramp = stat(level))) +
  scale_fill_viridis_d(end = 0.8) +
  scale_color_viridis_d(end = 0.8) +
  scale_fill_ramp_discrete(range = c(0.2, 0.7)) +
  facet_wrap(vars(group), ncol = 2) +
  #labeller = labeller(quota = c(`TRUE` = "Quota",
  #                             `FALSE` = "No Quota"))) +
  labs(y = "log(sigma_novel)", x = "trial (time)",
       #fill = "group", color = "group",
       fill_ramp = "Credible interval") +
  #ylim(0, 0.30) + xlim(24.5, 32.5) +
  guides(fill = "none", color = "none") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

# 5: emergT

p5 <- ggplot()+
  # confidence interval
  geom_smooth(data = pred_emergT, aes(x = trial, y = lower.HPD), method =  "loess", formula = y~x, se = FALSE, lty =  "dotted", lwd = 1,  col = "black") +
  geom_smooth(data = pred_emergT, aes(x = trial, y = upper.HPD), method =  "loess", formula = y~x, se = FALSE, lty = "dotted", lwd = 1,  col = "black") +
  # main line
  geom_smooth(data = pred_emergT, aes(x = trial, y = emmean), method =  "loess", formula = y~x, se = FALSE, lwd = 1.5, col = "black") +
  facet_wrap(vars(group), ncol = 2) +
  # real data
  geom_point(data = dat, aes(x = trial, y = log(emergT), fill = fishID, color = fishID), shape = 21, alpha = 0.5) +
  geom_smooth(data = dat,  aes(x = trial, y = log(emergT), color = fishID), method = "lm", formula = y ~x, se = FALSE, lwd = 0.1) + 
  scale_fill_viridis_d(end = 0.8) +
  scale_colour_viridis_d(end = 0.8) +
  labs(y = "log(emergT)", x = "trial (time)") +
  # themes
  theme_bw(base_size = 12) +  theme(legend.position = "none")

s5 <- ggplot(post_emergT, aes(x = trial, y = sigma, fill = group, colour = group)) +
  stat_lineribbon(aes(fill_ramp = stat(level))) +
  scale_fill_viridis_d(end = 0.8) +
  scale_color_viridis_d(end = 0.8) +
  scale_fill_ramp_discrete(range = c(0.2, 0.7)) +
  facet_wrap(vars(group), ncol = 2) +
  #labeller = labeller(quota = c(`TRUE` = "Quota",
  #                             `FALSE` = "No Quota"))) +
  labs(y = "log(sigma_emergT)", x = "trial (time)",
       #fill = "group", color = "group",
       fill_ramp = "Credible interval") +
  #ylim(0, 0.30) + xlim(24.5, 32.5) +
  guides(fill = "none", color = "none") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# putting all together

lay_out <- "
AABB
CCDD
#EE#
"
p1 + p2 + p3 + p4 + p5 + plot_layout(design = lay_out) + plot_annotation(tag_levels = "A")
s1 + s2 + s3 + s4 + s5 + plot_layout(design = lay_out) + plot_annotation(tag_levels = "A")

##################
### repeatability
##################
# here we calculate repeatablity for both contol and CO2 and compare
# TODO - you may want to do nboot 10000
# activity
activity_con <- rpt(sqrt(activity) ~  1 + trial + SL + (1|fishID), grname = "fishID", 
                    data = subset(dat, group == "control"), datatype = "Gaussian", 
                    nboot = 1000, npermut = 1000)

activity_co2 <- rpt(sqrt(activity) ~  1 + trial + SL + (1|fishID), grname = "fishID", 
                    data = subset(dat, group == "CO2"), datatype = "Gaussian", 
                    nboot = 1000, npermut = 1000)

comp_rpt(activity_co2, activity_con)


#shelter
shelter_con <- rpt(shelter ~  1 + trial + SL + (1|fishID), grname = "fishID", 
                   data = subset(dat, group == "control"), datatype = "Gaussian", 
                   nboot = 1000, npermut = 1000)

shelter_co2 <- rpt(shelter ~  1 + trial + SL + (1|fishID), grname = "fishID", 
                   data = subset(dat, group == "CO2"), datatype = "Gaussian", 
                   nboot = 1000, npermut = 1000)

comp_rpt(shelter_co2, shelter_con)

#thigmotaxis
thigmotaxis_con <- rpt(sqrt(thigmotaxis) ~  1 + trial + SL + (1|fishID), grname = "fishID", 
                       data = subset(dat, group == "control"), datatype = "Gaussian", 
                       nboot = 1000, npermut = 1000)

thigmotaxis_co2 <- rpt(sqrt(thigmotaxis) ~  1 + trial + SL + (1|fishID), grname = "fishID", 
                       data = subset(dat, group == "CO2"), datatype = "Gaussian", 
                       nboot = 1000, npermut = 1000)

comp_rpt(thigmotaxis_co2, thigmotaxis_con)

#novel
novel_con <- rpt(log(novel+0.5) ~  1 + trial + SL + (1|fishID), grname = "fishID", 
                 data = subset(dat, group == "control"), datatype = "Gaussian", 
                 nboot = 1000, npermut = 1000)

novel_co2 <- rpt(log(novel+0.5) ~  1 + trial + SL + (1|fishID), grname = "fishID", 
                 data = subset(dat, group == "CO2"), datatype = "Gaussian", 
                 nboot = 1000, npermut = 1000)

comp_rpt(novel_co2, novel_con)


# emergT
emergT_con <- rpt(log(emergT) ~  1 + trial + SL + (1|fishID), grname = "fishID", 
                  data = subset(dat, group == "control"), datatype = "Gaussian", 
                  nboot = 1000, npermut = 1000)

emergT_co2 <- rpt(log(emergT) ~  1 + trial + SL + (1|fishID), grname = "fishID", 
                  data = subset(dat, group == "CO2"), datatype = "Gaussian", 
                  nboot = 1000, npermut = 1000)

comp_rpt(emergT_co2, emergT_con)