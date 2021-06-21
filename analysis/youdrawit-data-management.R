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
feedback_data  <- read_csv("data/youdrawit-feedback-data.csv") %>%
  filter(study_starttime > 1620152231) %>%
  mutate(study_starttime = round(study_starttime)) %>%
  mutate(participantID = md5(paste(nick_name, study_starttime)),
         plotID = md5(paste(nick_name, study_starttime, parm_id)))

simulated_data <- read_csv("data/youdrawit-simulated-data.csv") %>%
  filter(study_starttime > 1620152231)  %>%
  mutate(study_starttime = round(study_starttime)) %>%
  mutate(participantID = md5(paste(nick_name, study_starttime)),
         plotID = md5(paste(nick_name, study_starttime, parm_id)))

users_data     <- read_csv("data/youdrawit-users-data.csv") %>%
  filter(study_starttime > 1620152231) %>%
  mutate(study_starttime = round(study_starttime)) %>%
  mutate(participantID = md5(paste(nick_name, study_starttime)))

researcher_participantID <- unique(users_data$participantID[users_data$recruitment == "I am the researcher"])

# ------------------------------------------------------------------------------
# CHECK NUMBER OF COMBOS -------------------------------------------------------
# ------------------------------------------------------------------------------

feedback_check <- feedback_data %>%
  nest(feedback_vals = c("x", "y", "ydrawn")) %>%
  select(participantID, plotID, nick_name, study_starttime, parm_id, feedback_vals) %>%
  group_by(participantID, plotID, nick_name, study_starttime, parm_id) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(participantID, nick_name, study_starttime) %>%
  mutate(total = sum(count),
         maxcount = max(count)) %>%
  ungroup() %>%
  pivot_wider(id_cols = c("participantID", "total", "maxcount"),
              names_from = "parm_id",
              values_from = "count")
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
  left_join(feedback_data, by = c("participantID", "nick_name", "study_starttime", "ip_address", "parm_id", "x", "y")) %>%
  left_join(users_data, by = c("participantID", "nick_name", "study_starttime")) %>%
  #right_join(feedback_clean, by = c("nick_name", "study_starttime")) %>%
  mutate(participantID = md5(paste(nick_name, study_starttime)),
         plotID = md5(paste(nick_name, study_starttime, parm_id))) %>%
  dplyr::select("participantID", "nick_name", "study_starttime", "age", "gender", "academic_study", "recruitment", "plotID", "start_time", "end_time", "parm_id", "x", "y", "ydrawn") %>%
  arrange("participantID", "start_time", "end_time", "x") %>%
  filter(participantID %!in% researcher_participantID)

# ------------------------------------------------------------------------------
# OBTAIN LOESS SMOOTHERS -------------------------------------------------------
# ------------------------------------------------------------------------------
# Fit Loess Smoother
loess.models <- all_data %>%
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
  dplyr::select(participantID, nick_name, study_starttime, age, gender, academic_study, recruitment, plotID, start_time, end_time, parm_id, x, y, ydrawn, yloess, residualdrawn, residualloess)

head(feedback_smooth)

# write.csv(feedback_smooth, file = "data/youdrawit-feedback-smooth.csv", row.names = F, na = "")
