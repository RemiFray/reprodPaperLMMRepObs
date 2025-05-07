library(dplyr)
library(forcats)
library(ggplot2)

dta <- read.csv("results/otter_AllModels_Results.csv")

pointsSize <- 2
boxLinesWidth <- 0.6
plotLinesWidth <- 0.2
labelSize <- 11

myThemePdf <- 
  theme(panel.background = element_rect(fill = "white", colour = NA),  
        panel.grid = element_line(colour = "grey92"), 
        panel.grid.minor = element_line(linewidth = rel(0.5)),
        axis.line = element_line(colour = 'black', 
                                 linewidth=boxLinesWidth, linetype='solid'),
        strip.background = element_rect(fill = "white", colour = "black", 
                                        linewidth = rel(2)),
        axis.ticks = element_line(linewidth = boxLinesWidth),
        panel.grid.minor.x = element_line(linewidth = boxLinesWidth*0.8),
        panel.grid.minor.y = element_line(linewidth = boxLinesWidth*0.8)
  )

cols <- c("#999999", "#E69F00", "#009E73")

p <- dta %>% 
  mutate(N_max = N + N_sd,
         N_min = N - N_sd,
         year = case_when(Model == "L&B" ~ year - 0.2 + 2000,
                          Model == "Mt" ~ year - 0.1 + 2000,
                          Model == "Yoshizaki" ~ year + 0.1 + 2000,
                          TRUE ~ year + 2000),
         Model = fct_recode(Model, "M\u03bb\u03b1" = "LMM2")) %>% 
  mutate(Model = fct_relevel(Model, "Mt", "M\u03bb\u03b1", "Yoshizaki")) %>% 
  ggplot(aes(x = year, y = N, shape = Model, col = Model)) +
  geom_point() +
  geom_errorbar(aes(ymin=N_2.5, ymax=N_97.5,width = 0.08)) +
  ylab("Population size") +
  scale_x_continuous(minor_breaks = 2006:2012, 
                     breaks = 2006:2012) +
  scale_colour_manual(values=cols) +
  myThemePdf


p

cairo_pdf("./figures/otter_Nestim.pdf", family="DejaVu Sans")
plot(p)
dev.off()
