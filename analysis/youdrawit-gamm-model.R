# ------------------------------------------------------------------------------
# LOAD LIBRARIES ---------------------------------------------------------------
# ------------------------------------------------------------------------------
library(tidyverse)
library(readr)
library(openssl)
library(mgcv)
library(lme4)
library(tictoc)
`%!in%` = function(x,y) !(x %in% y)
source("analysis/gamm-predict-function.R")

# ------------------------------------------------------------------------------
# IMPORT DATA ------------------------------------------------------------------
# ------------------------------------------------------------------------------
feedback_smooth  <- read_csv("data/youdrawit-feedback-smooth.csv")

exp_data <- feedback_smooth %>%
  filter(parm_id %!in% c("S", "F", "V", "N")) %>%
  separate(parm_id, into = c("beta", "points_end", "linear"), sep = "-") %>%
  mutate(beta = substr(beta, 5, 8),
         scale = as.factor(ifelse(linear == "true", "Linear", "Log"))) %>%
  select(participantID, nick_name, study_starttime, age, gender, academic_study, recruitment, plotID, start_time, end_time, beta, points_end, scale, x, y, ydrawn, yloess, residualdrawn, residualloess) %>%
  arrange(participantID, x)

factorCols = c('participantID', 'nick_name', 'age', 'gender', 'academic_study', 'recruitment', 'plotID', 'beta', 'points_end', 'scale')
exp_data[,factorCols] <- lapply(exp_data[,factorCols], factor)
summary(exp_data)

# ------------------------------------------------------------------------------
# PLOT RAW DATA ----------------------------------------------------------------
# ------------------------------------------------------------------------------

# yloess
exp_data %>%
  ggplot(aes(x = x)) +
  geom_line(aes(y = yloess, group = plotID, color = scale), alpha = 0.5) +
  geom_line(alpha = 0.1, aes(y = y, group = participantID)) +
  facet_grid(beta ~ points_end, scales = "free", labeller = labeller(beta = label_both, points_end = label_both)) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  scale_color_manual("Scale", values = c("steelblue", "orange")) +
  scale_x_continuous(limits = c(10, 20))

# residualloess
exp_data %>%
  ggplot(aes(x = x)) +
  geom_line(aes(y = residualloess, group = plotID, color = scale), alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(beta ~ points_end, scales = "free", labeller = labeller(beta = label_both, points_end = label_both)) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  scale_color_manual("Scale", values = c("steelblue", "orange")) +
  scale_x_continuous(limits = c(10, 20))

# ------------------------------------------------------------------------------
# FIT GAMM MODEL ---------------------------------------------------------------
# ------------------------------------------------------------------------------

# Fit GAMM
tic()
exp.gamm <- bam(residualdrawn ~ -1 + beta:points_end:scale + 
            s(x, by = beta:points_end:scale) +
            s(participantID, bs = "re") +
            s(x,participantID, bs = "re"),
            method = "REML",
            data = exp_data)
toc()
summary(exp.gamm)
anova(exp.gamm)

# Obtain Predictions
grid_data <- expand_grid(beta = c("0.1", "0.23"),
                         points_end = c("10", "15"),
                         scale = c("Linear", "Log"),
                         x = seq(10,20, 0.5),
                         participantID = exp_data$participantID[1])
preds <- predict_gamm(exp.gamm, newdata = grid_data, se = T, re_form = NA)
grid_data$estimate <- preds$prediction
grid_data$lower <- preds$prediction - (1.96 * preds$se)
grid_data$upper <- preds$prediction + (1.96 * preds$se)
head(grid_data)

# Plot Predictions
grid_data %>%
  ggplot(aes(x = x, y = estimate, group = scale, color = scale, fill = scale)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), color = NA, alpha = 0.4) +
  geom_line() +
  geom_line(data = exp_data, aes(x = x, y = residualdrawn, group = plotID), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(beta ~ points_end, scales = "free", labeller = labeller(beta = label_both, points_end = label_both)) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  scale_y_continuous("Residual \n (yloess - y)") +
  scale_color_manual("Scale", values = c("steelblue", "orange")) +
  scale_fill_manual("Scale", values = c("steelblue", "orange"))

