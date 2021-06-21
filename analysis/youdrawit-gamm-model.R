# ------------------------------------------------------------------------------
# LOAD LIBRARIES ---------------------------------------------------------------
# ------------------------------------------------------------------------------
library(tidyverse)
library(readr)
library(openssl)
library(mgcv)
library(lme4)
`%!in%` = function(x,y) !(x %in% y)
source("analysis/gamm-predict-function.R")

# ------------------------------------------------------------------------------
# IMPORT DATA ------------------------------------------------------------------
# ------------------------------------------------------------------------------
feedback_data  <- read_csv("data/youdrawit-feedback-data.csv") %>%
  filter(study_starttime > 1620152231) %>%
  mutate(study_starttime = round(study_starttime))
simulated_data <- read_csv("data/youdrawit-simulated-data.csv") %>%
  filter(study_starttime > 1620152231) %>%
  mutate(study_starttime = round(study_starttime))
users_data     <- read_csv("data/youdrawit-users-data.csv") %>%
  filter(study_starttime > 1620152231) %>%
  mutate(study_starttime = round(study_starttime))

researcher_nicknames <- unique(users_data$nick_name[users_data$recruitment == "I am the researcher"])

# ------------------------------------------------------------------------------
# CHECK ORDERS -----------------------------------------------------------------
# ------------------------------------------------------------------------------

feedback_check <- feedback_data %>%
  nest(feedback_vals = c("x", "y", "ydrawn")) %>%
  select(nick_name, study_starttime, parm_id, feedback_vals) %>%
  group_by(nick_name, study_starttime, parm_id) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(nick_name, study_starttime) %>%
  mutate(total = sum(count),
         maxcount = max(count)) %>%
  ungroup() %>%
  pivot_wider(id_cols = c("nick_name", "study_starttime", "total", "maxcount"),
              names_from = "parm_id",
              values_from = "count") %>%
  filter(nick_name %!in% researcher_nicknames) 
feedback_check

feedback_issues <- feedback_check %>%
  dplyr::select(nick_name, study_starttime, total, maxcount) %>%
  filter(maxcount > 1)

feedback_clean <- feedback_check %>%
  dplyr::select(nick_name, study_starttime, total, maxcount) %>%
  filter(maxcount == 1, total == 12)

simulated_check <- simulated_data %>%
  nest(simulated_vals = c("dataset", "x", "y")) %>%
  dplyr::select(nick_name, study_starttime, parm_id, simulated_vals) %>%
  group_by(nick_name, study_starttime, parm_id) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(nick_name, study_starttime) %>%
  mutate(total = sum(count),
         maxcount = max(count)) %>%
  ungroup() %>%
  pivot_wider(id_cols = c("nick_name", "study_starttime", "total", "maxcount"),
              names_from = "parm_id",
              values_from = "count") %>%
  filter(maxcount > 1 | total != 12)
simulated_check

# ------------------------------------------------------------------------------
# CHECK ORDERS -----------------------------------------------------------------
# ------------------------------------------------------------------------------

feedback_order <- feedback_data %>% 
  nest(feedbackVals = c("x", "y", "ydrawn")) %>%
  dplyr::select(nick_name, study_starttime, parm_id, feedbackVals) %>%
  rename(FeedbackParmID = parm_id) %>%
  group_by(nick_name, study_starttime) %>%
  mutate(order = 1:n()) %>%
  ungroup()

simulated_order <- simulated_data %>% 
  nest(simulatedVals = c("dataset", "x", "y")) %>%
  dplyr::select(nick_name, study_starttime, parm_id, simulatedVals) %>%
  rename(SimulatedParmID = parm_id) %>%
  group_by(nick_name, study_starttime) %>%
  mutate(order = 1:n()) %>%
  mutate(study_starttime = round(study_starttime)) %>%
  ungroup()

all_order <- feedback_order %>%
  full_join(simulated_order, by = c("nick_name", "study_starttime", "order")) %>%
  right_join(feedback_issues, by = c("nick_name", "study_starttime")) %>%
  mutate(match = ifelse(FeedbackParmID == SimulatedParmID, "1", "0")) %>%
  arrange(nick_name)

# ------------------------------------------------------------------------------
# COMBINE ALL DATA -------------------------------------------------------------
# ------------------------------------------------------------------------------
all_data <- simulated_data %>%
  filter(dataset == "line_data") %>%
  left_join(feedback_data, by = c("nick_name", "study_starttime", "ip_address", "parm_id", "x", "y")) %>%
  # right_join(feedback_clean, by = c("nick_name", "study_starttime")) %>%
  mutate(participantID = md5(paste(nick_name, study_starttime)),
         plotID = md5(paste(nick_name, study_starttime, parm_id))) %>%
  dplyr::select("participantID", "plotID", "nick_name", "study_starttime", "start_time", "end_time", "parm_id", "x", "y", "ydrawn") %>%
  arrange("nick_name", "study_starttime", "start_time", "end_time", "x")

exp_data <- all_data %>%
  filter(parm_id %!in% c("S", "F", "V", "N"))

# ------------------------------------------------------------------------------
# OBTAIN LOESS SMOOTHERS -------------------------------------------------------
# ------------------------------------------------------------------------------
# Fit Loess Smoother
loess.models <- exp_data %>%
  filter(!is.na(ydrawn)) %>%
  tidyr::nest(-plotID) %>%
  dplyr::mutate(
    # Perform loess calculation on each plotID
    loess.fit = purrr::map(data, loess,
                           formula = ydrawn ~ x),
    # Retrieve the fitted values from each model
    yloess = purrr::map(loess.fit, `[[`, "fitted")
  )

# Apply fitted y's as a new column
feedback_smooth <- loess.models %>%
  dplyr::select(-loess.fit) %>%
  tidyr::unnest(cols = c(data, yloess)) %>%
  mutate(residualdrawn = ydrawn - y,
         residualloess = yloess - y) %>%
  dplyr::select(participantID, plotID, parm_id, x, y, ydrawn, yloess, residualdrawn, residualloess)

# ------------------------------------------------------------------------------
# REMOVE OUTLIERS --------------------------------------------------------------
# ------------------------------------------------------------------------------

# simulated_pointdata <- simulated_data %>%
#   filter(dataset == "point_data") %>%
#   tidyr::nest(-plotID) %>%
#   dplyr::mutate(
#     # Perform exponential calculation on each plotID
#     nlme.fit = purrr::map(data, loess,
#                            formula = ydrawn ~ x),
#     # Retrieve the fitted values from each model
#     ymin = purrr::map(loess.fit, `[[`, "fitted"),
#     ymax = purrr::map(loess.fit, `[[`, "fitted")
#   )
# 
# 
# point_data <- simulated_data %>%
#   filter(dataset == "point_data") %>%
#   filter(nick_name == "a43e9f8c045e5762f6396aca7d3484f2", parm_id == "beta0.1-15-false")
# line_data <- simulated_data %>%
#   filter(dataset == "line_data") %>%
#   filter(nick_name == "a43e9f8c045e5762f6396aca7d3484f2", parm_id == "beta0.1-15-false")
# 
# predIntervalFunc <- function(point_data){
#   # Obtain starting value for beta
#   lm.fit <- lm(log(y) ~ x, data = point_data)
#   beta.0 <- coef(lm.fit)[1] %>% as.numeric()
#   # Use NLS to fit a better line to the data
#   start <- list(beta = beta.0)
#   nonlinear.fit <- nls(y ~ exp(x*beta),
#                        data = point_data,
#                        start = start)
#   yhat = predict(nonlinear.fit, newdata = line_data)
# }
# 
# library(investr)
# as_tibpredFit(nonlinear.fit, data = point_data, interval = "confidence", level = 0.95)
# predFit
# 
# library("propagate")
# predict(nonlinear.fit, newdata = line_data)
# predictNLS(nonlinear.fit, newdata = data.frame(x = line_data$x), interval="prediction", alpha=0.05, nsim=10000)$sim

# ------------------------------------------------------------------------------
# PLOT RAW DATA ----------------------------------------------------------------
# ------------------------------------------------------------------------------

model_data <- feedback_smooth %>%
  separate(parm_id, into = c("beta", "points_end", "linear"), sep = "-") %>%
  mutate(beta = substr(beta, 5, 8),
         scale = as.factor(ifelse(linear == "true", "Linear", "Log"))) %>%
  arrange(participantID, x)

factorCols = c('beta', 'points_end', 'scale', 'participantID', 'plotID')
model_data[,factorCols] <- lapply(model_data[,factorCols], factor)
summary(model_data)

# Spaghetti Plot
model_data %>%
  ggplot(aes(x = x)) +
  geom_line(aes(y = yloess, group = plotID, color = scale), alpha = 0.5) +
  geom_line(alpha = 0.1, aes(y = y, group = participantID)) +
  facet_grid(beta ~ points_end, scales = "free", labeller = labeller(beta = label_both, points_end = label_both)) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  scale_color_manual("Scale", values = c("steelblue", "orange")) +
  scale_x_continuous(limits = c(10, 20))

model_data %>%
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
gamm <- bam(residualloess ~ beta*points_end*scale + 
            s(x, by = beta:points_end:scale) +
            s(participantID, bs = "re") +
            s(x,participantID, bs = "re"),
            method = "REML",
            data = model_data)
summary(gamm)
anova(gamm)

# Obtain Predictions
grid_data <- expand_grid(beta = c("0.1", "0.23"),
                         points_end = c("10", "15"),
                         scale = c("Linear", "Log"),
                         x = seq(10,20, 0.5),
                         participantID = model_data$participantID[1])
preds <- predict_gamm(gamm, newdata = grid_data, se = T, re_form = NA)
grid_data$estimate <- preds$prediction
grid_data$lower <- preds$prediction - (2 * preds$se)
grid_data$upper <- preds$prediction + (2 * preds$se)
head(grid_data)

# Plot Predictions
grid_data %>%
  ggplot(aes(x = x, y = estimate, group = scale, color = scale, fill = scale)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), color = NA, alpha = 0.2) +
  geom_line() +
  geom_line(data = model_data, aes(x = x, y = residualloess, group = plotID), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(beta ~ points_end, scales = "free", labeller = labeller(beta = label_both, points_end = label_both)) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  scale_y_continuous("Residual \n (yloess - y)") +
  scale_color_manual("Scale", values = c("steelblue", "orange")) +
  scale_fill_manual("Scale", values = c("steelblue", "orange"))
