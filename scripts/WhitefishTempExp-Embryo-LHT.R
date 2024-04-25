#### CLEAR THE ENVIRONMENT FIRST ---------------------------------------------

rm(list = ls(all.names = TRUE))


#### LOAD PACKAGES -----------------------------------------------------------

library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)
library(afex)
library(buildmer)
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)
library(emmeans)


#### LOAD HATCHING DATA ------------------------------------------------------

hatch.FR <- read_excel("data/Coregonine-Temperature-Experiment-EuropeFrance-Hatch.xlsx", sheet = "hatching") %>% 
  select(population, species_form, family, male, female, female_tl_mm, female_fm_g, male_tl_mm, male_fm_g, block, no, temperature, plate, eye, hatch, dpf, ADD, include.incubation) %>% 
  mutate(female = factor(female),
         male = factor(male),
         family = factor(family),
         block = factor(block),
         # Create a variable with population and species combined
         group = factor(interaction(population, species_form)))

#hatch.FI <- read_excel("data/Coregonine-Temperature-Experiment-EuropeFinland-Hatch.xlsx", sheet = "hatching") %>% 
#  select(population, species_form = species, family, male, female, female_tl_mm, female_fm_g, male_tl_mm, male_fm_g, block, no, temperature, plate, eye, hatch, dpf, ADD, include.incubation) %>% 
#  filter(temperature %in% c(6.9, 8), species_form == "lavaretus") %>% 
#  mutate(female = factor(female),
#         male = factor(male),
#         family = factor(family),
#         block = factor(block),
#         # Create a variable with population and species combined
#         group = factor("konnevesi.littoral"))

hatch <- hatch.FR %>% 
  mutate(temperature = factor(temperature, ordered = TRUE, levels = c(6.9, 7, 9), labels = c(7, 7, 9)),
         group = factor(group, ordered = TRUE, levels = c("constance.littoral", "constance.pelagic", "leman.littoral", "bourget.littoral")),
         trans.dpf = dpf^(1/4),
         trans.ADD = ADD^(1/4)) %>% 
  filter(include.incubation == "y")


#### FILTER TO EACH TRAITS' DATASET --------------------------------------------------------------

## filter to only eyed embryos
hatch.survival <- hatch %>% filter(eye != 0)

## filter to only hatched embryos
hatch.survived <- hatch %>% filter(!is.na(dpf), hatch == 1)


#### STATISTICAL ANALYSIS - SURVIVAL -----------------------------------------------------

## backward elimination to select best model
hatch.survival.glm <- buildmer(hatch ~ temperature + group + temperature:group + 
                                 (1|family) + (1|male) + (1|female) + (1|block) + (1|plate),
                               data = hatch.survival, family = binomial, 
                               buildmerControl = buildmerControl(direction = 'backward', args = list(control = glmerControl(optimizer = 'bobyqa'))))
( hatch.survival.glm.formula <- formula(hatch.survival.glm@model))

## fit best model
hatch.survival.glm.final <- glmer(hatch.survival.glm.formula, data = hatch.survival, family = binomial)

## likelihood ratio test for fixed effects
mixed(hatch.survival.glm.formula, data = hatch.survival, method = "LRT")

## fit model without random effects for LRT
# family
hatch.survival.glm.family <- glmer(hatch ~ temperature + group + temperature:group + 
                                     (1|female) + (1|block) + (1|plate), data = hatch.survival, 
                                   family = binomial, control = glmerControl(optimizer = "bobyqa"))
# female
hatch.survival.glm.female <- glmer(hatch ~ temperature + group + temperature:group + 
                                     (1|family) + (1|block) + (1|plate), data = hatch.survival, 
                                   family = binomial, control = glmerControl(optimizer = "bobyqa"))
# block
hatch.survival.glm.block <- glmer(hatch ~ temperature + group + temperature:group + 
                                    (1|female) + (1|family) + (1|plate), data = hatch.survival, 
                                  family = binomial, control = glmerControl(optimizer = "bobyqa"))
# plate
hatch.survival.glm.plate <- glmer(hatch ~ temperature + group + temperature:group + 
                                    (1|female) + (1|family) + (1|block), data = hatch.survival, 
                                  family = binomial, control = glmerControl(optimizer = "bobyqa"))

## Compare full to reduced models (LRT)
# family
anova(hatch.survival.glm.family, hatch.survival.glm.final)
# female
anova(hatch.survival.glm.female, hatch.survival.glm.final)
# block
anova(hatch.survival.glm.block, hatch.survival.glm.final)
# plate
anova(hatch.survival.glm.plate, hatch.survival.glm.final)


#### STATISTICAL ANALYSIS - INCUBATION PERIOD (DPF) --------------------------------------

## fit full model
hatch.dpf.glm.full <- lmer(trans.dpf ~ temperature + group + temperature:group + 
                           (1|family) + (1|male) + (1|female) + (1|block) + (1|plate), 
                           REML = FALSE, data = hatch.survived)

## backward elimination to select best model
hatch.dpf.glm <- step(hatch.dpf.glm.full)
( hatch.dpf.glm.formula <- get_model(hatch.dpf.glm)@call[["formula"]])

## fit best model
hatch.dpf.glm.final <- lmer(hatch.dpf.glm.formula, data = hatch.survived, REML = FALSE)

## check residuals for normality
lattice::qqmath(hatch.dpf.glm.final, id = 0.1, idLabels = ~.obs)
hist(rstudent(hatch.dpf.glm.final))

## likelihood ratio test for fixed and random effects
mixed(hatch.dpf.glm.formula, data = hatch.survived, method = "LRT")
rand(hatch.dpf.glm.final)


#### STATISTICAL ANALYSIS - INCUBATION PERIOD (ADD) --------------------------------------

## fit full model
hatch.ADD.glm.full <- lmer(trans.ADD ~ 1 + temperature + group + temperature:group + 
                             (1|family) + (1|male) + (1|female) + (1|block) + (1|plate), 
                           data = hatch.survived)

## backward elimination to select best model
hatch.ADD.glm <- step(hatch.ADD.glm.full)
( hatch.ADD.glm.formula <- get_model(hatch.ADD.glm)@call[["formula"]])

## fit best model
hatch.ADD.glm.final <- lmer(hatch.ADD.glm.formula, data = hatch.survived)

## check residuals for normality
lattice::qqmath(hatch.ADD.glm.final, id = 0.1, idLabels = ~.obs)
hist(rstudent(hatch.ADD.glm.final))

## likelihood ratio test for fixed and random effects
mixed(hatch.ADD.glm.formula, data = hatch.survived, method = "LRT")
rand(hatch.ADD.glm.final)


#### CALCULATE MEAN AND SE FOR NA & FI POPULATIONS -----------------------------------------------

temp <- data.frame(group = c("constance.littoral", "constance.littoral",
                             "constance.pelagic", "constance.pelagic",
                             "leman.littoral", "leman.littoral",
                             "bourget.littoral", "bourget.littoral"),
                   temperature = factor(rep(c(7, 9), 4), ordered = TRUE, levels = c(7, 9)))

## Embryo Survival Overall
hatch.survival.summary <- hatch %>% filter(eye != 0) %>% 
  group_by(population, group, temperature) %>% 
  summarize(mean.trait = mean(hatch),
            se.trait = sd(hatch)/sqrt(n())) %>% 
  mutate(trait = "survival")

## Embryo Survival - Standardized Within Family
hatch.survival.summary.family <- hatch %>% filter(eye != 0) %>% 
  group_by(population, temperature, group, family) %>% 
  summarize(mean.hatch = mean(hatch)) %>% ungroup()

hatch.survival.stand <- hatch.survival.summary.family %>% filter(temperature == 7) %>% 
  select(group, family, local.survival = mean.hatch)

hatch.survival.summary.stand <- hatch.survival.summary.family %>% left_join(hatch.survival.stand) %>% 
  filter(group != "leman.littoral" | family != "F03M13") %>%  ## No data at 7C
  mutate(survival.diff = 100*((mean.hatch-local.survival)/local.survival)) %>% 
  filter(!is.nan(survival.diff), survival.diff != Inf) %>%
  group_by(population, temperature, group) %>% 
  summarize(mean.trait.stand = mean(survival.diff),
            se.trait.stand = sd(survival.diff)/sqrt(n())) %>% 
  left_join(temp) %>% 
  mutate(se.trait.stand = ifelse(se.trait.stand == 0, NA, se.trait.stand),
         group = factor(group, ordered = TRUE, levels = c("konnevesi.littoral", "constance.littoral", "constance.pelagic", "leman.littoral", "bourget.littoral")),
         trait = "survival")


## Days Post Fertilization
hatch.dpf.summary <- hatch %>% filter(!is.na(dpf), hatch == 1) %>% 
  group_by(population, group, temperature) %>% 
  summarize(mean.trait = mean(dpf),
            se.trait = sd(dpf)/sqrt(n())) %>% 
  mutate(trait = "dpf")

## Days Post Fertilization - Standardized Within Family
hatch.dpf.summary.family <- hatch %>% filter(!is.na(dpf), hatch == 1) %>% 
  group_by(population, temperature, group, family) %>% 
  summarize(mean.dpf = mean(dpf)) %>% ungroup()

hatch.dpf.stand <- hatch.dpf.summary.family %>% filter(temperature == 7) %>% 
  select(group, family, local.dpf = mean.dpf)

hatch.dpf.summary.stand <- hatch.dpf.summary.family %>% left_join(hatch.dpf.stand) %>% 
  filter(group != "leman.littoral" | family != "F03M13") %>%  ## No data at 7C
  mutate(dpf.diff = 100*((mean.dpf-local.dpf)/local.dpf)) %>%
  filter(!is.nan(dpf.diff), dpf.diff != Inf) %>%
  group_by(population, temperature, group) %>% 
  summarize(mean.trait.stand = mean(dpf.diff),
            se.trait.stand = sd(dpf.diff)/sqrt(n())) %>% 
  left_join(temp) %>% 
  mutate(se.trait.stand = ifelse(se.trait.stand == 0, NA, se.trait.stand),
         group = factor(group, ordered = TRUE, levels = c("konnevesi.littoral", "constance.littoral", "constance.pelagic", "leman.littoral", "bourget.littoral")),
         trait = "dpf")

## Accumulated Degree-Days
hatch.ADD.summary <- hatch %>% filter(!is.na(ADD), hatch == 1) %>% 
  group_by(population, group, temperature) %>% 
  summarize(mean.trait = mean(ADD),
            se.trait = sd(ADD)/sqrt(n())) %>% 
  mutate(trait = "ADD")

## Accumulated Degree-Days - Standardized Within Family
hatch.ADD.summary.family <- hatch %>% filter(!is.na(ADD), hatch == 1) %>% 
  group_by(population, temperature, group, family) %>% 
  summarize(mean.ADD = mean(ADD)) %>% ungroup()

hatch.ADD.stand <- hatch.ADD.summary.family %>% filter(temperature == 7) %>% 
  select(group, family, local.ADD = mean.ADD)

hatch.ADD.summary.stand <- hatch.ADD.summary.family %>% left_join(hatch.ADD.stand) %>% 
  filter(group != "leman.littoral" | family != "F03M13") %>%  ## No data at 7C
  mutate(ADD.diff = 100*((mean.ADD-local.ADD)/local.ADD)) %>%
  filter(!is.nan(ADD.diff), ADD.diff != Inf) %>%
  group_by(population, temperature, group) %>% 
  summarize(mean.trait.stand = mean(ADD.diff),
            se.trait.stand = sd(ADD.diff)/sqrt(n())) %>% 
  left_join(temp) %>% 
  mutate(se.trait.stand = ifelse(se.trait.stand == 0, NA, se.trait.stand),
         group = factor(group, ordered = TRUE, levels = c("konnevesi.littoral", "constance.littoral", "constance.pelagic", "leman.littoral", "bourget.littoral")),
         trait = "ADD")

## Combine traits
traitsOverall.all <- bind_rows(hatch.survival.summary, hatch.dpf.summary, hatch.ADD.summary) %>% 
  mutate(group = factor(group, ordered = TRUE, levels = c("konnevesi.littoral", "constance.littoral", "constance.pelagic", "leman.littoral", "bourget.littoral"),
                        labels = c("S. Konnevesi ", "Constance (Lit.) ", "Constance (Pel.) ", "Geneva ", "Bourget"))) %>% 
  group_by(temperature) %>% 
  mutate(errorbar_width = 0.04 * n())
traitsStand.all <- bind_rows(hatch.survival.summary.stand, hatch.dpf.summary.stand, hatch.ADD.summary.stand) %>% 
  mutate(group = factor(group, ordered = TRUE, levels = c("konnevesi.littoral", "constance.littoral", "constance.pelagic", "leman.littoral", "bourget.littoral"),
                        labels = c("S. Konnevesi", "Constance\n(Lit.)", "Constance\n(Pel.)", "Geneva", "Bourget")))


#### VISUALIZATIONS - MEANS ----------------------------------------------------------------------
## Embryo Survival
plot.survival.legend <- ggplot(filter(traitsOverall.all, trait == "survival"), 
                               aes(x = temperature, y = (mean.trait * 100), 
                                   group = group, color = group, shape = group, 
                                   linetype = group, width = errorbar_width)) + 
  geom_line(size = 0.6, position = position_dodge(0.17)) +
  geom_point(size = 1.9, position = position_dodge(0.17), stroke = 1) +
  geom_errorbar(aes(ymin = (mean.trait - se.trait) * 100, ymax = (mean.trait + se.trait) * 100), 
                position = position_dodge(0.17), size = 0.5, 
                linetype = "solid", show.legend = FALSE) +
  scale_x_discrete(expand = c(0, 0.2)) +
  scale_y_continuous(limits = c(20, 100), breaks = seq(20, 100, 20), expand = c(0, 0)) +
  scale_color_grey("combine", start = 0.0, end = 0.8) +
  scale_shape_manual("combine", values = c(2, 5, 1, 0, 6)) +
  scale_linetype_manual("combine", values = c("solid", "dashed", "dotted", "solid", "dotdash")) +
  labs(y = "Mean Embryo Survival (%)", x = "Mean Incubation Temperature (°C)") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color = "Black", size = 16, margin = margin(0, 8, 0, 0)),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.ticks.length = unit(2, 'mm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.key.width = unit(1.2, 'cm'),
        legend.position = "top")
plot.survival <- plot.survival.legend + theme(legend.position = "none")

## Days Post Fertilization
plot.dpf <- ggplot(filter(traitsOverall.all, trait == "dpf"), 
                   aes(x = temperature, y = mean.trait, 
                       group = group, color = group, shape = group, 
                       linetype = group, width = errorbar_width)) + 
  geom_line(size = 0.6, position = position_dodge(0.17)) +
  geom_point(size = 1.9, position = position_dodge(0.17), stroke = 1) +
  geom_errorbar(aes(ymin = mean.trait - se.trait, ymax = mean.trait + se.trait), 
                position = position_dodge(0.17), size = 0.5, 
                linetype = "solid", show.legend = FALSE) +
  scale_x_discrete(expand = c(0, 0.2)) +
  scale_y_continuous(limits = c(40, 62.5), breaks = seq(40, 60, 5), expand = c(0, 0)) +
  scale_color_grey("combine", start = 0.0, end = 0.8) +
  scale_shape_manual("combine", values = c(2, 5, 1, 0, 6)) +
  scale_linetype_manual("combine", values = c("solid", "dashed", "dotted", "solid", "dotdash")) +
  labs(y = "Mean DPF", x = "Mean Incuabtion Temperature (°C)") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color = "Black", size = 16, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.ticks.length = unit(2, 'mm'),
        legend.position = "none")
plot.dpf

## Accumulated Degree-Days
plot.add <- ggplot(filter(traitsOverall.all, trait == "ADD"), 
                   aes(x = temperature, y = mean.trait, 
                       group = group, color = group, shape = group, 
                       linetype = group, width = errorbar_width)) + 
  geom_line(size = 0.6, position = position_dodge(0.17)) +
  geom_point(size = 1.9, position = position_dodge(0.17), stroke = 1) +
  geom_errorbar(aes(ymin = mean.trait - se.trait, ymax = mean.trait + se.trait), 
                position = position_dodge(0.17), size = 0.5, 
                linetype = "solid", show.legend = FALSE) +
  scale_x_discrete(expand = c(0, 0.2)) +
  scale_y_continuous(limits = c(400, 450), breaks = seq(400, 450, 10), expand = c(0, 0)) +
  scale_color_grey("combine", start = 0.0, end = 0.8) +
  scale_shape_manual("combine", values = c(2, 5, 1, 0, 6)) +
  scale_linetype_manual("combine", values = c("solid", "dashed", "dotted", "solid", "dotdash")) +
  labs(y = "Mean ADD (°C)", x = "Mean Incuabtion Temperature (°C)") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color = "Black", size = 16, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.ticks.length = unit(2, 'mm'),
        legend.position = "none")
plot.add


#### VISUALIZATIONS - STANDARDIZED ---------------------------------------------------------------
plot.survival.stand <- ggplot(filter(traitsStand.all, trait == "survival",  temperature %in% c(8, 9)),
                              aes(x = group, y = mean.trait.stand, group = temperature)) + 
  geom_bar(stat = "identity", fill = "gray50", color = "gray20", width = 0.8) +
  geom_hline(yintercept = 0, color = "black", linetype = "solid", alpha = 0.5) +
  geom_errorbar(aes(ymin = (mean.trait.stand - se.trait.stand), ymax = (mean.trait.stand + se.trait.stand)), 
                position = position_dodge(0.6), size = 0.8, width = 0.4, show.legend = FALSE) +
  scale_y_continuous(limits = c(-75.0, 4), breaks = seq(-80.0, 4, 10), expand = c(0, 0)) +
  labs(y = "Standardized Embryo Survival (%)", x = "Population") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color = "Black", size = 16, margin = margin(0, 8, 0, 0)),
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.ticks.length = unit(2, 'mm'))
plot.survival.stand

plot.dpf.stand <- ggplot(filter(traitsStand.all, trait == "dpf", temperature %in% c(8, 9)),
                         aes(x = group, y = mean.trait.stand, group = temperature)) + 
  geom_bar(stat = "identity", fill = "gray50", color = "gray20", width = 0.8) +
  geom_hline(yintercept = 0, color = "black", linetype = "solid", alpha = 0.5) +
  geom_errorbar(aes(ymin = (mean.trait.stand - se.trait.stand), ymax = (mean.trait.stand + se.trait.stand)), 
                position = position_dodge(0.6), size = 0.8, width = 0.4, show.legend = FALSE) +
  scale_y_continuous(limits = c(-30.0, 1.75), breaks = seq(-30.0, 0, 5), expand = c(0, 0)) +
  labs(y = "Standardized DPF (%)", x = "Population") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color = "Black", size = 16, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.ticks.length = unit(2, 'mm'),
        legend.position = "none")
plot.dpf.stand

plot.add.stand <- ggplot(filter(traitsStand.all, trait == "ADD", temperature %in% c(8, 9)),
                         aes(x = group, y = mean.trait.stand, group = temperature)) + 
  geom_bar(stat = "identity", fill = "gray50", color = "gray20", width = 0.8) +
  geom_hline(yintercept = 0, color = "black", linetype = "solid", alpha = 0.5) +
  geom_errorbar(aes(ymin = (mean.trait.stand - se.trait.stand), ymax = (mean.trait.stand + se.trait.stand)), 
                position = position_dodge(0.6), size = 0.8, width = 0.4, show.legend = FALSE) +
  scale_y_continuous(limits = c(-7.5, 2.5), breaks = seq(-7.5, 2.5, 2.5), expand = c(0, 0)) +
  labs(y = "Standardized ADD (%)", x = "Population") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color = "Black", size = 16, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.ticks.length = unit(2, 'mm'),
        legend.position = "none")
plot.add.stand


## Combine all figures
plot.lht <- align_plots(plot.survival, plot.dpf, plot.add,
                        plot.survival.stand, plot.dpf.stand, plot.add.stand,
                        align = c("v"))

lht.title <- ggdraw() + draw_label("Mean Incubation Temperature (°C)", x = 0.525, y = 0.5, fontfamily = "Arial", size = 16)

plot.lht.legend <- get_legend(plot.survival.legend)
plot.lht.row <- plot_grid(ggdraw(plot.lht[[1]]),
                          ggdraw(plot.lht[[2]]),
                          ggdraw(plot.lht[[3]]),
                          ncol = 3)

top <- plot_grid(plot.lht.legend, 
                 plot.lht.row, 
                 lht.title, 
                 rel_heights = c(0.1, 1, 0.05),
                 ncol = 1)
top


lht.stand.title <- ggdraw() + draw_label("Populations", x = 0.525, y = 1, fontfamily = "Arial", size = 16)
plot.lht.stand.row <- plot_grid(ggdraw(plot.lht[[4]]),
                                ggdraw(plot.lht[[5]]),
                                ggdraw(plot.lht[[6]]),
                                ncol = 3)
bot <- plot_grid(NULL,
                 plot.lht.stand.row,
                 lht.stand.title,
                 ncol = 1,
                 rel_heights = c(0.05, 1, 0.05))
                 
plot.all <- plot_grid(top, bot, rel_heights = c(1, 1), ncol = 1)
plot.all

ggsave("figures/Fig1-LHT.png", plot = plot.all, width = 12, height = 10, dpi = 300)

