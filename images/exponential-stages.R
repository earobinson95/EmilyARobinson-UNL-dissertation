library(tidyverse)
library(ggforce)

data <- tibble(x = seq(0,6,0.001),
               y = exp(x)) %>%
  mutate(group = ifelse(x <= 2, "G1", 
                        ifelse(x <= 4, "G2", "G3")))

data %>%
  ggplot(aes(x = x, y = y, color = group)) +
  geom_line(show.legend = T) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  # facet_zoom(xlim = c(0,2), shrink = F, zoom.size = 3) +
  # facet_zoom(xy = group == 'G1', zoom.size = 1/3) +
  geom_mark_rect(show.legend = T) +
  scale_color_brewer("Stage", palette = "Dark2", labels = c("Early", "Mid", "Late"))

data %>%
  filter(group == "G1") %>%
  ggplot(aes(x = x, y = y, color = group)) +
  geom_line(show.legend = F) +
  geom_line(show.legend = F) +
  theme_bw() +
  theme(aspect.ratio = 1)
  
data %>%
  filter(group == "G2") %>%
  ggplot(aes(x = x, y = y, color = group)) +
  geom_line(show.legend = F) +
  geom_line(show.legend = F) +
  theme_bw() +
  theme(aspect.ratio = 1) 


data %>%
  filter(group == "G3") %>%
  ggplot(aes(x = x, y = y, color = group)) +
  geom_line(show.legend = F) +
  geom_line(show.legend = F) +
  theme_bw() +
  theme(aspect.ratio = 1) 

