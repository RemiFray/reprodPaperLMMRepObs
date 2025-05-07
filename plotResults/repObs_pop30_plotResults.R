library(dplyr)
library(forcats)
library(ggplot2)

dta <- read.csv("./results/repObs_pop30.csv") %>% 
  mutate(in95interval = (N_97.5 > 30 & N_2.5 < 30) %>% as.factor)

pointsSize <- 2
boxLinesWidth <- 0.6
plotLinesWidth <- 0.2
labelSize <- 11

myThemePdf <- 
  theme(panel.background = element_rect(fill = "white", colour = NA),  
        panel.grid = element_line(colour = "grey92"), 
        panel.grid.minor.y = element_line(linewidth = boxLinesWidth*0.8),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line = element_line(colour = 'black', 
                                 linewidth=boxLinesWidth, linetype='solid'),
        strip.background = element_rect(fill = "white", colour = "black", 
                                        linewidth = rel(2)),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(linewidth = boxLinesWidth),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
  )

cols <- c("red", "black")

# ---- pop size ----
p1 <- dta %>%  
  ggplot(aes(x = iter, y = N_mean)) +
  geom_pointrange(aes(ymin = N_2.5, ymax = N_97.5, 
                      colour = in95interval)) +
  geom_hline(aes(yintercept = 30)) +
  ylab("Population size estimate") +
  scale_colour_manual(values=cols) +
  myThemePdf +
  theme(legend.position = "none")


p1

cairo_pdf("./figures/simulpop30_Nestim.pdf", family="DejaVu Sans")
plot(p)
dev.off()

# capture rates ----

cols <- c("purple", "red", "blue", "black")

p <- dta %>% 
  select(iter, matches("l[1-5]_[m29]"), in95interval) %>% 
  pivot_longer(cols = matches("l[1-5]"), 
               names_to = c("occas", ".value"),
               names_pattern = "l([1-5])_(.*)") %>% 
  rename(l_2.5 = "2.5", l_97.5 = "97.5") %>% 
  mutate(l95interval = (((l_97.5 > 0.5 & l_2.5 < 0.5)+2) * 
                          as.numeric(in95interval)) %>% as.factor) %>% 
  ggplot(aes(x = iter, y = mean)) +
  geom_pointrange(aes(ymin = l_2.5, ymax = l_97.5, 
                      colour = l95interval), 
                  position = "jitter") +
  geom_hline(aes(yintercept = 0.5)) +
  ylab("Mean number of capture (\u03bb) estimate") +
  scale_colour_manual(values=cols) +
  myThemePdf +
  theme(legend.position = "none")


p



cairo_pdf("./figures/simulpop30_Lestim.pdf",
          width = 12, heigh = 7, family="DejaVu Sans")
plot(p)
dev.off()

# ---- alpha ----
p2 <- dta %>%  
  ggplot(aes(x = iter, y = a_mean)) +
  geom_pointrange(aes(ymin = a_2.5, ymax = a_97.5, 
                      colour = in95interval)) +
  geom_hline(aes(yintercept = 0.95)) +
  ylab("Population size estimate") +
  scale_colour_manual(values=cols) +
  myThemePdf +
  theme(legend.position = "none")


p2

cairo_pdf("./figures/simulpop30_Aestim.pdf", family="DejaVu Sans")
plot(p2)
dev.off()
