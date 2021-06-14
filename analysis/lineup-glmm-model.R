# --------------------------------------------------------------------------
# LOAD PACKAGES ------------------------------------------------------------
# --------------------------------------------------------------------------
library(readr)
library(tidyverse)
library(lme4)
library(emmeans)
library(cowplot)

# --------------------------------------------------------------------------
# IMPORT DATA --------------------------------------------------------------
# --------------------------------------------------------------------------

lineup_results_data <- read.csv("data/lineup-pilot-data-2020.11.30.csv")
summary(lineup_results_data)
model_data <- lineup_results_data %>%
  filter(participant_count > 5) %>%
  mutate(test_param = factor(test_param, levels = c("log", "linear")))

# --------------------------------------------------------------------------
# GENERATE IMAGES FOR ODDS RATIO PLOT --------------------------------------
# --------------------------------------------------------------------------

xMid_vals  <- c(14.5, 13, 11.5)
sigma_vals <- c(0.25, 0.12, 0.05)
yRange_vals = c(10,100)

# Obtain alphahat, betahat, and thetahat for different midpoints.
coefEst <- function(xMid, xRange = c(0,20), yRange = yRange_vals){
  
  # This creates the line y = -x (scaled to fit the x and y ranges)
  # |*            0
  # |  *
  # |    *
  # |      *
  # |        1
  # |          2
  # |0___________3
  #
  # where 1, 2, 3 represent different points used to determine the line curvature
  
  lineData   <- tibble(xLine = seq(xRange[1],xRange[2],0.1),
                       yLine = -(abs(diff(yRange))/abs(diff(xRange)))*(xLine-xRange[1])+yRange[2])
  pointsData <- tibble(xPoint = c(xRange[1], (xMid-0.1), (xMid+0.1), xRange[2]),
                       yPoint = c(yRange[1], lineData$yLine[lineData$xLine == xMid], lineData$yLine[lineData$xLine == xMid], yRange[2]))
  
  # Connecting the 0 points in the illustration above with the 3rd point that
  # determines curvature gives us a set of 3 points to use to fit an exponential
  # line to the data.
  
  # We fit a linear regression to the log-transformed data to get starting values
  lm.fit <- lm(log(yPoint) ~ xPoint, data = pointsData)
  
  alpha.0  <- exp(coef(lm.fit)[1]) %>% as.numeric()
  beta.0   <- coef(lm.fit)[2] %>% as.numeric()
  theta.0 <- min(pointsData$yPoint) * 0.5  # Why 0.5?
  
  # and then use NLS to fit a better line to the data
  start <- list(alpha = alpha.0, beta = beta.0, theta = theta.0)
  nonlinear.fit   <- nls(yPoint ~ alpha * exp(beta * xPoint) + theta ,
                         data = pointsData, start = start)
  
  coefficients <- tibble(alphahat = (coef(nonlinear.fit)[1] %>% as.numeric()),
                         betahat  = coef(nonlinear.fit)[2] %>% as.numeric(),
                         thetahat = coef(nonlinear.fit)[3] %>% as.numeric())
  
  return(coefficients)
}

expSim <- function(alphahat, betahat, thetahat, sigma, nReps = 1, N = 50, xRange = c(0,20), yRange = yRange_vals){
  
  alpha = alphahat/(exp(sigma^2/2))
  beta  = betahat
  theta = thetahat
  
  vals <- seq(xRange[1], xRange[2], length.out = N*3/4)
  xvals <- sample(vals, N, replace = T)
  xvals <- jitter(xvals)
  
  expData <- tibble(x = rep(xvals, nReps),
                    y = alpha*exp(beta*x + rnorm(N*nReps,0,sigma)) + theta)
  return(expData)
}

coefData <- tibble(xMid = xMid_vals) %>%
  mutate(coefficients = pmap(list(xMid),coefEst)) %>%
  unnest(coefficients)

set.seed(68505)
simData <- tibble(diff.num    = seq(1,3,1),
                   curvature   = c("E",  "M",  "H"),
                   variability = c("Lv", "Lv", "Lv"),
                   xMid        = c(rep(xMid_vals[1],1), rep(xMid_vals[2],1), rep(xMid_vals[3],1)),
                   sigma       = sigma_vals) %>%
  left_join(coefData, by = "xMid") %>%
  mutate(data = pmap(list(alphahat, betahat, thetahat,sigma),expSim)) %>%
  unnest(data)

simData  %>%
  mutate(curvature = factor(curvature, levels = c("E", "M", "H"))) %>%
  ggplot(aes(x = x, y = y, color = curvature)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio = 1) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none") +
  scale_color_manual(values = c("#004400", "#116611", "#55aa55"))

tE_nM <- simData  %>%
  filter(curvature %in% c("E", "M")) %>%
  mutate(curvature = factor(curvature, levels = c("E", "M", "H"))) %>%
  ggplot(aes(x = x, y = y, color = curvature)) +
  geom_point(size = 0.5) +
  theme_test() +
  theme(aspect.ratio = 1) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none") +
  # scale_color_manual(values = c("#004400", "#116611", "#55aa55")) +
  scale_color_manual(values = c("green3", "black"))
tE_nM
ggsave(filename = "images/tE_nM.png", plot = tE_nM, dpi = 50)

tE_nH <- simData  %>%
  filter(curvature %in% c("E", "H")) %>%
  mutate(curvature = factor(curvature, levels = c("E", "M", "H"))) %>%
  ggplot(aes(x = x, y = y, color = curvature)) +
  geom_point() +
  theme_test() +
  theme(aspect.ratio = 1) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none") +
  # scale_color_manual(values = c("#004400", "#116611", "#55aa55")) +
  scale_color_manual(values = c("green3", "black"))
tE_nH
ggsave(filename = "images/tE_nH.png", plot = tE_nH, dpi = 50)

tM_nE <- simData  %>%
  filter(curvature %in% c("E", "M")) %>%
  mutate(curvature = factor(curvature, levels = c("E", "M", "H"))) %>%
  ggplot(aes(x = x, y = y, color = curvature)) +
  geom_point() +
  theme_test() +
  theme(aspect.ratio = 1) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none") +
  # scale_color_manual(values = c("#004400", "#116611", "#55aa55")) +
  scale_color_manual(values = c("black", "green3"))
tM_nE
ggsave(filename = "images/tM_nE.png", plot = tM_nE, dpi = 50)

tM_nH <- simData  %>%
  filter(curvature %in% c("H", "M")) %>%
  mutate(curvature = factor(curvature, levels = c("E", "M", "H"))) %>%
  ggplot(aes(x = x, y = y, color = curvature)) +
  geom_point() +
  theme_test() +
  theme(aspect.ratio = 1) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none") +
  # scale_color_manual(values = c("#004400", "#116611", "#55aa55")) +
  scale_color_manual(values = c("green3", "black"))
tM_nH
ggsave(filename = "images/tM_nH.png", plot = tM_nH, dpi = 50)

tH_nE <- simData  %>%
  filter(curvature %in% c("E", "H")) %>%
  mutate(curvature = factor(curvature, levels = c("E", "M", "H"))) %>%
  ggplot(aes(x = x, y = y, color = curvature)) +
  geom_point() +
  theme_test() +
  theme(aspect.ratio = 1) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none") +
  # scale_color_manual(values = c("#004400", "#116611", "#55aa55")) +
  scale_color_manual(values = c("black", "green3"))
tH_nE
ggsave(filename = "images/tH_nE.png", plot = tH_nE, dpi = 50)

tH_nM <- simData  %>%
  filter(curvature %in% c("H", "M")) %>%
  mutate(curvature = factor(curvature, levels = c("E", "M", "H"))) %>%
  ggplot(aes(x = x, y = y, color = curvature)) +
  geom_point() +
  theme_test() +
  theme(aspect.ratio = 1) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none") +
  # scale_color_manual(values = c("#004400", "#116611", "#55aa55")) +
  scale_color_manual(values = c("black", "green3"))
tH_nM
ggsave(filename = "images/tH_nM.png", plot = tH_nM, dpi = 50)

# --------------------------------------------------------------------------
# Run GLMM Model -----------------------------------------------------------
# --------------------------------------------------------------------------

glmm_mod <- glmer(correct ~ curvature*test_param + 
                            (1 | run) +
                            (1 | data_name),
                  data = model_data,
                  family = binomial(link = "logit"))
summary(glmm_mod)
anova(glmm_mod)
lineup_anova_table <- car::Anova(glmm_mod) %>% 
                      as.data.frame() %>%
                      rownames_to_column("Effect") %>%
  mutate(Effect = c("Curvature", "Scale", "Curvature x Scale"),
         Chisq = round(Chisq, 2),
         `Pr(>Chisq)` = ifelse(`Pr(>Chisq)` < 0.0001, "<0.0001", round(`Pr(>Chisq)`, 5)))
write.csv(lineup_anova_table, file = "data/lineup-anova-table.csv", row.names = F, na = "")

# LSMEANS
lsmeans <- emmeans(glmm_mod, c("curvature", "test_param"), type = "response")
lsmeans

lsmeans %>%
  as_tibble() %>%
  mutate(curvature = factor(curvature, levels = c("t-E_n-H", "t-H_n-E", "t-E_n-M", "t-M_n-E", "t-M_n-H", "t-H_n-M"))) %>%
  ggplot(aes(x = curvature, y = prob, group = test_param)) +
  geom_bar(stat = "identity", fill = "lightgray") +
  # geom_text(aes(y = asymp.UCL + 0.05)) +
  facet_wrap(~test_param, ncol = 1) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.1) +
  theme_bw() +
  scale_y_continuous("Probability of target panel detected", limit = c(-0.05,1.1), breaks = seq(0,1,0.2)) +
  scale_x_discrete("Curvature")

# ODDS RATIOS
slices <- emmeans(glmm_mod, ~ test_param | curvature, type = "response")
slices
odds_ratios <- pairs(slices, infer = c(TRUE, TRUE)) %>% 
               as_tibble() %>%
               mutate(target = factor(substr(curvature,3,3), levels = c("E", "M", "H"), labels = c("Easy", "Medium", "Hard")),
                      null   = factor(substr(curvature,7,7), levels = c("E", "M", "H"), labels = c("Easy", "Medium", "Hard"))) %>%
  select(contrast, target, null, odds.ratio, SE, asymp.LCL, asymp.UCL, z.ratio, p.value)

write.csv(odds_ratios, file = "data/lineup-odds-ratios.csv", row.names = F, na = "")

# PLOT ODDS RATIOS

dodge <- position_dodge(width=0.9)
odds_ratio_plot <- odds_ratios %>%
  ggplot(aes(x = odds.ratio, y = null, color = target, shape = target)) + 
  geom_point(position = dodge, size = 3) + 
  geom_errorbar(aes(xmin = asymp.LCL, xmax = asymp.UCL), position = dodge, width = .1) +
  geom_vline(xintercept = 1) +
  theme_bw()  +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text  = element_text(size = 8),
        legend.key.size = unit(0.7, "line"),
        legend.position = "bottom"
  ) +
  scale_y_discrete("Null plot type", position = "left") +
  scale_x_continuous("Odds ratio (on log scale) \n (Log vs Linear)", trans = "log10") + 
  scale_color_manual("Target Plot Type", values = c("#004400", "#116611", "#55aa55")) + 
  scale_shape_discrete("Target Plot Type")
odds_ratio_plot

picsList <- c("images/tM_nH.png", "images/tE_nH.png", 
              "images/tH_nM.png", "images/tE_nM.png", 
              "images/tH_nE.png", "images/tM_nE.png"
)

pimage <- axis_canvas(odds_ratio_plot, axis = 'y') + 
  draw_image(picsList[1], y = 2.75, scale = 0.6) +
  draw_image(picsList[2], y = 2.25, scale = 0.6) +
  draw_image(picsList[3], y = 1.75, scale = 0.6) +
  draw_image(picsList[4], y = 1.25, scale = 0.6) +
  draw_image(picsList[5], y = 0.75, scale = 0.6) +
  draw_image(picsList[6], y = 0.25, scale = 0.6)

# insert the image strip into the plot
ggdraw(insert_yaxis_grob(odds_ratio_plot, pimage, position = "right", clip = "on"))


