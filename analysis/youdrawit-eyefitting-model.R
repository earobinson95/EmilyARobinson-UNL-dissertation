# ------------------------------------------------------------------------------
# LOAD LIBRARIES ---------------------------------------------------------------
# ------------------------------------------------------------------------------
library(tidyverse)
library(readr)
library(openssl)
library(tictoc)
library(mgcv)
library(lme4)
library(lmerTest)
library(emmeans)
library(pls)
`%!in%` = function(x,y) !(x %in% y)
source("analysis/gamm-predict-function.R")
source("analysis/lmer-predict-function.R")

# ------------------------------------------------------------------------------
# IMPORT DATA ------------------------------------------------------------------
# ------------------------------------------------------------------------------
feedback_smooth  <- read_csv("data/youdrawit-feedback-smooth.csv")

eyefitting_data <- feedback_smooth %>%
  filter(parm_id %in% c("S", "F", "V", "N")) %>%
  select(participantID, nick_name, study_starttime, age, gender, academic_study, recruitment, plotID, start_time, end_time, parm_id, x, y, ydrawn, yloess, residualdrawn, residualloess) %>%
  arrange(participantID, x)

factorCols = c('participantID', 'nick_name', 'age', 'gender', 'academic_study', 'recruitment', 'plotID', 'parm_id')
eyefitting_data[,factorCols] <- lapply(eyefitting_data[,factorCols], factor)
summary(eyefitting_data)

# ------------------------------------------------------------------------------
# OBTAIN FIRST PRINCIPAL COMPONENT ---------------------------------------------
# ------------------------------------------------------------------------------
# trial <- eyefitting_data %>%
#   filter(participantID == "a47a513b05d2ae93a2778eefaae62831", parm_id == "F")
# trial.pcr <- pcr(ydrawn ~ x, data = trial)
# trial$ypcr <- predict(trial.pcr)
# trial %>%
#   ggplot(aes(x = x)) +
#   geom_line(aes(y = y)) +
#   geom_line(aes(y = ypcr), linetype = "dashed") +
#   geom_line(aes(y = yloess), color = "steelblue") +
#   facet_wrap(~parm_id) +
#   theme_bw() +
#   theme(aspect.ratio = 1)

# Fit PCR
pcrCalc <- function(data){
  pcr.mod <- pcr(ydrawn ~ x, data = data)
  predict(pcr.mod)
}

eyefitting_model_data <- eyefitting_data %>%
  tidyr::nest(-plotID) %>%
  dplyr::mutate(ypcr = purrr::map(data, pcrCalc)) %>%
  tidyr::unnest(cols = c(data, ypcr)) %>%
  mutate(residualpcr = ydrawn - ypcr,
         residualpcr.loess = yloess - ypcr) %>%
  rename(yols = y,
         residualols = residualdrawn,
         residualols.loess = residualloess) %>%
  dplyr::select(participantID, nick_name, study_starttime, age, gender, academic_study, recruitment, plotID, start_time, end_time, parm_id, x, yols, ypcr, ydrawn, yloess, residualols, residualols.loess, residualpcr, residualpcr.loess)

head(eyefitting_model_data)

# write.csv(feedback_smooth, file = "data/youdrawit-feedback-smooth.csv", row.names = F, na = "")

# ------------------------------------------------------------------------------
# PLOT RAW DATA ----------------------------------------------------------------
# ------------------------------------------------------------------------------

# yloess (OLS)
eyefitting_model_data %>%
  ggplot(aes(x = x)) +
  geom_line(aes(y = yloess, group = plotID), alpha = 0.5, color = "steelblue") +
  geom_line(alpha = 0.2, aes(y = yols, group = participantID)) +
  # geom_line(alpha = 0.2, aes(y = ypcr, group = participantID), linetype = "dashed") +
  facet_wrap(~ parm_id) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  scale_x_continuous(limits = c(0, 20)) +
  ggtitle("OLS")

# yloess (PCR)
eyefitting_model_data %>%
  ggplot(aes(x = x)) +
  geom_line(aes(y = yloess, group = plotID), alpha = 0.5, color = "steelblue") +
  # geom_line(alpha = 0.2, aes(y = yols, group = participantID)) +
  geom_line(alpha = 0.2, aes(y = ypcr, group = participantID)) +
  facet_wrap(~ parm_id) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  scale_x_continuous(limits = c(0, 20)) +
  ggtitle("PCR")

# Residual Loess
eyefitting_model_data %>%
  ggplot(aes(x = x)) +
  geom_line(aes(y = residualols.loess, group = plotID), color = "steelblue", alpha = 0.5) +
  geom_line(aes(y = residualpcr.loess, group = plotID), color = "orange", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~parm_id) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  scale_x_continuous(limits = c(0, 20))

# ------------------------------------------------------------------------------
# FIT GAMM MODEL  --------------------------------------------------------------
# ------------------------------------------------------------------------------

# Fit GAMM

# OLS
tic()
eyefitting.ols.gamm <- bam(residualols ~ -1 + parm_id + 
              s(x, by = parm_id) +
              s(participantID, bs = "re") +
              s(x,participantID, bs = "re"),
            method = "REML",
            data = eyefitting_model_data)
toc()
summary(eyefitting.ols.gamm)
anova(eyefitting.ols.gamm)

# PCR
tic()
eyefitting.pcr.gamm <- bam(residualpcr ~ -1 + parm_id + 
                             s(x, by = parm_id) +
                             s(participantID, bs = "re") +
                             s(x,participantID, bs = "re"),
                           method = "REML",
                           data = eyefitting_model_data)
toc()
summary(eyefitting.pcr.gamm)
anova(eyefitting.pcr.gamm)

# Obtain Predictions
eyefitting.grid.gamm <- expand_grid(parm_id = c("S", "V", "F", "N"),
                                    x = seq(0,20, 0.5),
                                    participantID = eyefitting_model_data$participantID[1])

# OLS
eyefitting.ols.preds <- predict_gamm(eyefitting.ols.gamm, newdata = eyefitting.grid.gamm, se = T, re_form = NA)
eyefitting.grid.gamm$ols.pred <- eyefitting.ols.preds$prediction
eyefitting.grid.gamm$ols.lower <- eyefitting.ols.preds$prediction - (1.96 * eyefitting.ols.preds$se)
eyefitting.grid.gamm$ols.upper <- eyefitting.ols.preds$prediction + (1.96 * eyefitting.ols.preds$se)

# PCR
eyefitting.pcr.preds <- predict_gamm(eyefitting.pcr.gamm, newdata = eyefitting.grid.gamm, se = T, re_form = NA)
eyefitting.grid.gamm$pcr.pred <- eyefitting.pcr.preds$prediction
eyefitting.grid.gamm$pcr.lower <- eyefitting.pcr.preds$prediction - (1.96 * eyefitting.pcr.preds$se)
eyefitting.grid.gamm$pcr.upper <- eyefitting.pcr.preds$prediction + (1.96 * eyefitting.pcr.preds$se)

# Plot Predictions
eyefitting.grid.gamm %>%
  filter((parm_id %in% c("F", "N", "S") | (x <= 16 & x >= 4))) %>%
  ggplot(aes(x = x)) +
  geom_line(data = eyefitting_model_data, aes(x = x, y = residualpcr, group = plotID, color = "PCR"), alpha = 0.1) +
  geom_line(data = eyefitting_model_data, aes(x = x, y = residualols, group = plotID, color = "OLS"), alpha = 0.1) +
  geom_ribbon(aes(ymin = ols.lower, ymax = ols.upper, fill = "OLS"), color = NA, alpha = 0.7) +
  geom_line(aes(y = ols.pred, color = "OLS")) +
  geom_ribbon(aes(ymin = pcr.lower, ymax = pcr.upper, fill = "PCR"), color = NA, alpha = 0.7) +
  geom_line(aes(y = pcr.pred, color = "PCR")) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  facet_wrap(~parm_id) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  scale_y_continuous("Residual") +
  scale_color_manual("Estimates", values = c("steelblue", "orange")) +
  scale_fill_manual("Estimates", values = c("steelblue", "orange")) +
  ggtitle("Eyefitting GAMM")

# ------------------------------------------------------------------------------
# FIT LMER MODEL (OLS) ---------------------------------------------------------
# ------------------------------------------------------------------------------

# Fit LMER

# OLS
tic()
eyefitting.ols.lmer <- lmer(residualols ~ -1 + parm_id + x:parm_id + (1|participantID),
                       data = eyefitting_model_data)
toc()
summary(eyefitting.ols.lmer)
anova(eyefitting.ols.lmer)

# PCR
tic()
eyefitting.pcr.lmer <- lmer(residualpcr ~ -1 + parm_id + x:parm_id + (1|participantID),
                            data = eyefitting_model_data)
toc()
summary(eyefitting.pcr.lmer)
anova(eyefitting.pcr.lmer)

# Obtain Predictions
eyefitting.grid.lmer <- expand_grid(parm_id = c("S", "V", "F", "N"),
                         x = seq(0,20, 0.5),
                         participantID = eyefitting_model_data$participantID[1])

eyefitting.ols.preds <- predict_lmer(eyefitting.ols.lmer, newdata=eyefitting.grid.lmer, re.form = NA, se.fit = TRUE, nsim = 50)
eyefitting.grid.lmer$ols.pred <- eyefitting.ols.preds$fit
eyefitting.grid.lmer$ols.lower <- eyefitting.ols.preds$ci.fit[1,]
eyefitting.grid.lmer$ols.upper <- eyefitting.ols.preds$ci.fit[2,]

eyefitting.pcr.preds <- predict_lmer(eyefitting.pcr.lmer, newdata=eyefitting.grid.lmer, re.form = NA, se.fit = TRUE, nsim = 50)
eyefitting.grid.lmer$pcr.pred <- eyefitting.pcr.preds$fit
eyefitting.grid.lmer$pcr.lower <- eyefitting.pcr.preds$ci.fit[1,]
eyefitting.grid.lmer$pcr.upper <- eyefitting.pcr.preds$ci.fit[2,]

# Plot Predictions
eyefitting.grid.lmer %>%
  filter((parm_id %in% c("F", "N", "S") | (x <= 16 & x >= 4))) %>%
  ggplot(aes(x = x)) +
  geom_line(data = eyefitting_model_data, aes(x = x, y = residualols, group = plotID, color = "OLS"), alpha = 0.1) +
  geom_line(data = eyefitting_model_data, aes(x = x, y = residualpcr, group = plotID, color = "PCR"), alpha = 0.1) +
  geom_ribbon(aes(ymin = ols.lower, ymax = ols.upper, fill = "OLS"), color = NA, alpha = 0.7) +
  geom_line(aes(y = ols.pred, color = "OLS")) +
  geom_ribbon(aes(ymin = pcr.lower, ymax = pcr.upper, fill = "PCR"), color = NA, alpha = 0.7) +
  geom_line(aes(y = pcr.pred, color = "PCR")) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  facet_wrap(~parm_id) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  scale_y_continuous("Residual") +
  scale_color_manual("Estimates", values = c("steelblue", "orange")) +
  scale_fill_manual("Estimates", values = c("steelblue", "orange")) +
  ggtitle("Eyefitting LMER")

# ------------------------------------------------------------------------------
# Compare Sum of Squares (OLS VS PCR) ------------------------------------------
# ------------------------------------------------------------------------------

ss_data <- eyefitting_model_data %>%
  group_by(participantID, nick_name, study_starttime, age, gender, academic_study, recruitment, plotID, parm_id, start_time, end_time) %>%
  summarise(olsSS = sum(residualols^2),
            pcrSS = sum(residualpcr^2),
            olsSS.loess = sum(residualols.loess^2),
            pcrSS.loess = sum(residualpcr.loess^2)) %>%
  ungroup() %>%
  pivot_longer(cols = c("olsSS", "pcrSS", "olsSS.loess", "pcrSS.loess"),
               names_to = "Fit",
               values_to = "SS")

ss.lmer <- lmer(log(SS) ~ parm_id*Fit + (1|participantID),
              data = ss_data %>% filter(Fit %in% c("olsSS", "pcrSS")))
anova(ss.lmer)
plot(ss.lmer)

ss.emmeans <- emmeans(ss.lmer, ~parm_id:Fit, type = "response") %>%
  as_tibble()
ss.slicediffs <- pairs(emmeans(ss.lmer, ~Fit | parm_id, type = "response"), infer = c(TRUE, TRUE)) %>% 
  as_tibble()

ss.emmeans %>%
  ggplot(aes(x = parm_id, y = response, group = Fit, fill = Fit)) +
  geom_bar(stat = "identity", position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL, width = 0.5), position = position_dodge(0.9)) +
  theme_bw() +
  theme(aspect.ratio = 0.75) +
  scale_fill_manual(values = c("steelblue", "orange"), labels = c("OLS", "PCR")) +
  scale_y_continuous("Sum of Squares") +
  scale_x_discrete("Data Set")

ss.slicediffs %>%
  ggplot(aes(x = ratio, y = parm_id)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), width = 0.5) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  theme_bw() +
  theme(aspect.ratio = 0.5) +
  scale_x_continuous("Sum of Squares Odds Ratio \n (OLS vs PCR)", limits = c(0,10), breaks = seq(0,10,2)) +
  scale_y_discrete("Data Set")
