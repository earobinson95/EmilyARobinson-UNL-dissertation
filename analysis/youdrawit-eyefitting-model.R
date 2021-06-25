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

# ------------------------------------------------------------------------------
# IMPORT DATA ------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FEEDBACK DATA
feedback_smooth  <- read_csv("data/youdrawit-feedback-smooth.csv")
eyefitting_data <- feedback_smooth %>%
  filter(parm_id %in% c("S", "F", "V", "N")) %>%
  select(participantID, nick_name, study_starttime, age, gender, academic_study, recruitment, plotID, start_time, end_time, parm_id, x, y, ydrawn, yloess, residualdrawn, residualloess) %>%
  arrange(participantID, x)
factorCols = c('participantID', 'nick_name', 'age', 'gender', 'academic_study', 'recruitment', 'plotID', 'parm_id')
eyefitting_data[,factorCols] <- lapply(eyefitting_data[,factorCols], factor)
# summary(eyefitting_data)

# SIMULATED DATA
simulated_data <- read_csv("data/youdrawit-simulated-data.csv") %>%
  filter(study_starttime > 1620152231)  %>%
  mutate(study_starttime = round(study_starttime)) %>%
  mutate(participantID = md5(paste(nick_name, study_starttime)),
         plotID = md5(paste(nick_name, study_starttime, parm_id)))

factorCols = c('participantID', 'nick_name', 'plotID', 'parm_id', 'dataset')
simulated_data[,factorCols] <- lapply(simulated_data[,factorCols], factor)
# summary(simulated_data)

# ------------------------------------------------------------------------------
# OBTAIN FIRST PRINCIPAL COMPONENT ---------------------------------------------
# ------------------------------------------------------------------------------

# Trial Data
participantIDs <- levels(eyefitting_data$participantID)
trial.sim <- simulated_data %>%
  filter(participantID == participantIDs[18], parm_id == "F", dataset == "point_data")
summary(trial.sim)
trial.feedback <- eyefitting_data %>%
  filter(participantID == participantIDs[18], parm_id == "F")
summary(trial.feedback)

# PCA
trial.pca <- prcomp(trial.sim[,c("x","y")])
trial.sim$PC1 <- predict(trial.pca)[,1]
trial.feedback$PC1 <- predict(trial.pca, newdata = trial.feedback)[,1]
pca.lm <- lm(y ~ PC1, data = trial.sim)
trial.feedback$ypca <- predict(pca.lm, newdata = trial.feedback)

trial.feedback %>%
  ggplot(aes(x = x)) +
  geom_line(aes(y = y)) +
  geom_line(aes(y = ypca), color = "black", linetype = "dashed") +
  geom_line(aes(y = yloess), color = "steelblue") +
  geom_point(data = trial.sim, aes(y = y)) +
  facet_wrap(~parm_id) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  scale_x_continuous(limits = c(0,20))

# Fit PCA
pcaCalc <- function(data){
  point.data <- simulated_data %>% 
    mutate(participantID = as.character(participantID),
           parm_id = as.character(parm_id)) %>%
    filter(dataset == "point_data", participantID == as.character(data[1,"participantID"]), parm_id == as.character(data[1, 'parm_id']))
  pca.mod <- prcomp(point.data[,c("x","y")])
  point.data$PC1 <- predict(pca.mod)[,1]
  data$PC1 <- predict(pca.mod, newdata = data)[,1]
  pca.lm.mod <- lm(y ~ PC1, data = point.data)
  predict(pca.lm.mod, newdata = data)
}

eyefitting_model_data <- eyefitting_data %>% 
  mutate(participantID = as.character(participantID),
         parm_id = as.character(parm_id)) %>%
  tidyr::nest(-plotID) %>%
  dplyr::mutate(ypca = purrr::map(data, pcaCalc)) %>%
  tidyr::unnest(cols = c(data, ypca)) %>%
  mutate(residualpca = ydrawn - ypca,
         residualpca.loess = yloess - ypca,
         participantID = factor(participantID),
         parm_id = factor(parm_id)) %>%
  rename(yols = y,
         residualols = residualdrawn,
         residualols.loess = residualloess) %>%
  dplyr::select(participantID, nick_name, study_starttime, age, gender, academic_study, recruitment, plotID, start_time, end_time, parm_id, x, yols, ypca, ydrawn, yloess, residualols, residualols.loess, residualpca, residualpca.loess)

head(eyefitting_model_data) 

# ------------------------------------------------------------------------------
# PLOT RAW DATA ----------------------------------------------------------------
# ------------------------------------------------------------------------------

# yloess (OLS)
eyefitting_model_data %>%
  ggplot(aes(x = x)) +
  geom_line(aes(y = yloess, group = plotID), alpha = 0.5, color = "steelblue") +
  geom_line(alpha = 0.2, aes(y = yols, group = participantID)) +
  # geom_line(alpha = 0.2, aes(y = ypca, group = participantID), linetype = "dashed") +
  facet_wrap(~ parm_id) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  scale_x_continuous(limits = c(0, 20)) +
  ggtitle("OLS")

# yloess (pca)
eyefitting_model_data %>%
  ggplot(aes(x = x)) +
  geom_line(aes(y = yloess, group = plotID), alpha = 0.5, color = "steelblue") +
  # geom_line(alpha = 0.2, aes(y = yols, group = participantID)) +
  geom_line(alpha = 0.2, aes(y = ypca, group = participantID)) +
  facet_wrap(~ parm_id) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  scale_x_continuous(limits = c(0, 20)) +
  ggtitle("PCA")

# Residual Loess
eyefitting_model_data %>%
  ggplot(aes(x = x)) +
  geom_line(aes(y = residualols.loess, group = plotID, color = "OLS"), alpha = 0.3) +
  geom_line(aes(y = residualpca.loess, group = plotID, color = "PCA"), alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~parm_id) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  scale_x_continuous(limits = c(0, 20)) +
  scale_color_manual("Estimate", values = c("steelblue", "orange"))

# ------------------------------------------------------------------------------
# FIT GAMM MODEL  --------------------------------------------------------------
# ------------------------------------------------------------------------------

# Fit GAMM

# OLS
# tic()
eyefitting.ols.gamm <- bam(residualols ~ -1 + parm_id + 
                             s(x, by = parm_id) +
                             s(participantID, bs = "re") +
                             s(x,participantID, bs = "re"),
                           method = "REML",
                           data = eyefitting_model_data)
# toc()
# summary(eyefitting.ols.gamm)
# anova(eyefitting.ols.gamm)

# pca
# tic()
eyefitting.pca.gamm <- bam(residualpca ~ -1 + parm_id + 
                             s(x, by = parm_id) +
                             s(participantID, bs = "re") +
                             s(x,participantID, bs = "re"),
                           method = "REML",
                           data = eyefitting_model_data)
# toc()
# summary(eyefitting.pca.gamm)
# anova(eyefitting.pca.gamm)

# Obtain Predictions
eyefitting.grid.gamm <- expand_grid(parm_id = c("S", "V", "F", "N"),
                                    x = seq(0,20, 0.5),
                                    participantID = eyefitting_model_data$participantID[1])

# OLS
eyefitting.ols.preds <- predict_gamm(eyefitting.ols.gamm, newdata = eyefitting.grid.gamm, se = T, re_form = NA)
eyefitting.grid.gamm$ols.pred <- eyefitting.ols.preds$prediction
eyefitting.grid.gamm$ols.lower <- eyefitting.ols.preds$prediction - (1.96 * eyefitting.ols.preds$se)
eyefitting.grid.gamm$ols.upper <- eyefitting.ols.preds$prediction + (1.96 * eyefitting.ols.preds$se)

# pca
eyefitting.pca.preds <- predict_gamm(eyefitting.pca.gamm, newdata = eyefitting.grid.gamm, se = T, re_form = NA)
eyefitting.grid.gamm$pca.pred <- eyefitting.pca.preds$prediction
eyefitting.grid.gamm$pca.lower <- eyefitting.pca.preds$prediction - (1.96 * eyefitting.pca.preds$se)
eyefitting.grid.gamm$pca.upper <- eyefitting.pca.preds$prediction + (1.96 * eyefitting.pca.preds$se)

# Plot Predictions
eyefitting.grid.gamm %>%
  filter((parm_id %in% c("F", "N", "S") | (x <= 16 & x >= 4))) %>%
  ggplot(aes(x = x)) +
  geom_line(data = eyefitting_model_data, aes(x = x, y = residualpca, group = plotID, color = "PCA"), alpha = 0.1) +
  geom_line(data = eyefitting_model_data, aes(x = x, y = residualols, group = plotID, color = "OLS"), alpha = 0.1) +
  geom_ribbon(aes(ymin = ols.lower, ymax = ols.upper, fill = "OLS"), color = NA, alpha = 0.7) +
  geom_line(aes(y = ols.pred, color = "OLS")) +
  geom_ribbon(aes(ymin = pca.lower, ymax = pca.upper, fill = "PCA"), color = NA, alpha = 0.7) +
  geom_line(aes(y = pca.pred, color = "PCA")) +
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
# tic()
eyefitting.ols.lmer <- lmer(residualols ~ -1 + parm_id + x:parm_id + (1|participantID),
                            data = eyefitting_model_data)
# toc()
# summary(eyefitting.ols.lmer)
# anova(eyefitting.ols.lmer)

# pca
# tic()
eyefitting.pca.lmer <- lmer(residualpca ~ -1 + parm_id + x:parm_id + (1|participantID),
                            data = eyefitting_model_data)
# toc()
# summary(eyefitting.pca.lmer)
# anova(eyefitting.pca.lmer)

# Obtain Predictions
eyefitting.ols.grid.lmer  <- ref_grid(eyefitting.ols.lmer, at = list(x = seq(1,20,0.5)))
eyefitting.ols.preds.lmer <- emmeans(eyefitting.ols.grid.lmer, ~ parm_id:x, cov.reduce = FALSE) %>% 
  as_tibble()

eyefitting.pca.grid.lmer  <- ref_grid(eyefitting.pca.lmer, at = list(x = seq(1,20,0.5)))
eyefitting.pca.preds.lmer <- emmeans(eyefitting.pca.grid.lmer, ~ parm_id:x, cov.reduce = FALSE) %>% 
  as_tibble()

eyefitting.preds.lmer <- eyefitting.ols.preds.lmer %>%
  full_join(eyefitting.pca.preds.lmer, by = c("x", "parm_id"), suffix = c(".ols", ".pca"))

# Plot Predictions
eyefitting.preds.lmer %>%
  filter((parm_id %in% c("F", "N", "S") | (x <= 16 & x >= 4))) %>%
  ggplot(aes(x = x)) +
  geom_line(data = eyefitting_model_data, aes(x = x, y = residualols, group = plotID, color = "OLS"), alpha = 0.1) +
  geom_line(data = eyefitting_model_data, aes(x = x, y = residualpca, group = plotID, color = "PCA"), alpha = 0.1) +
  geom_ribbon(aes(ymin = asymp.LCL.ols, ymax = asymp.UCL.ols, fill = "OLS"), color = NA, alpha = 0.7) +
  geom_line(aes(y = emmean.ols, color = "OLS")) +
  geom_ribbon(aes(ymin = asymp.LCL.pca, ymax = asymp.UCL.pca, fill = "PCA"), color = NA, alpha = 0.7) +
  geom_line(aes(y = emmean.pca, color = "PCA")) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  facet_wrap(~parm_id) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  scale_y_continuous("Residual") +
  scale_color_manual("Estimates", values = c("steelblue", "orange")) +
  scale_fill_manual("Estimates", values = c("steelblue", "orange")) +
  ggtitle("Eyefitting LMER")

# ------------------------------------------------------------------------------
# Compare Sum of Squares (OLS VS PCA) ------------------------------------------
# ------------------------------------------------------------------------------

ss_data <- eyefitting_model_data %>%
  group_by(participantID, nick_name, study_starttime, age, gender, academic_study, recruitment, plotID, parm_id, start_time, end_time) %>%
  summarise(olsSS = sum(residualols^2),
            pcaSS = sum(residualpca^2),
            olsSS.loess = sum(residualols.loess^2),
            pcaSS.loess = sum(residualpca.loess^2)) %>%
  ungroup() %>%
  pivot_longer(cols = c("olsSS", "pcaSS", "olsSS.loess", "pcaSS.loess"),
               names_to = "Fit",
               values_to = "SS")

ss.lmer <- lmer(log(SS) ~ parm_id*Fit + (1|participantID),
                data = ss_data %>% filter(Fit %in% c("olsSS", "pcaSS")))
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
  scale_fill_manual(values = c("steelblue", "orange"), labels = c("OLS", "PCA")) +
  scale_y_continuous("Sum of Squares") +
  scale_x_discrete("Data Set")

ss.slicediffs %>%
  ggplot(aes(x = ratio, y = parm_id)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), width = 0.5) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  theme_bw() +
  theme(aspect.ratio = 0.5) +
  scale_x_continuous("Sum of Squares Odds Ratio \n (OLS vs PCA)", limits = c(0,2), breaks = seq(0,2,0.5)) +
  scale_y_discrete("Data Set")
