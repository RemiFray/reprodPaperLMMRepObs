library(tidyverse)
library(scales)
library(gridExtra)
library(cowplot)

# ------ Data ------ ----

biasSummary <- read.csv("./results/repObs_AllModels_AllResults.csv") %>% 
  mutate(Convergence = as.factor((N_Rhat<= 1.1) * 1)) %>% 
  mutate(p = case_when(model == "Mt" ~ lambda-0.02,
                       model == "Yoshi" ~ lambda,
                       model == "LMM2" ~ lambda+0.02),
         model = as.factor(model),
         Convergence = fct_relevel(Convergence, "1", "0")) %>% 
  mutate(model = fct_recode(model, "Yoshizaki's" = "Yoshi",
                            "Mla" = "LMM2"))
biasSummaryMeans <- biasSummary %>% 
  group_by(model, p, alpha, S) %>% 
  summarise(N_estim = mean(N_mean),
            N_2.5 = mean(N_2.5),
            N_97.5 = mean(N_97.5)) %>% 
  ungroup()

minN <- min(floor(min(biasSummaryMeans$N_2.5)/10)*10,
            min(biasSummary$N_mean))
maxN <- max(500, 
            ceiling(max(biasSummaryMeans$N_97.5)/10)*10, 
            max(biasSummary$N_mean))
maxN <- ceiling(max(biasSummaryMeans$N_97.5)/10)*10

# ------ Theme ------ ----
S <- c(5, 7, 9)
models <- unique(biasSummary$model)
als <- unique(biasSummary$alpha)

# Themes and graphic parameters
pointsSize <- 2
boxLinesWidth <- 0.6
plotLinesWidth <- 0.2
labelSize <- 11

myThemePdf <- 
  # theme_classic() +
  theme(panel.background = element_rect(fill = "white", colour = NA), 
        # panel.border = element_rect(fill = NA, colour = "grey20"), 
        panel.grid = element_line(colour = "grey92"), 
        panel.grid.minor = element_line(size = rel(0.5)),
        axis.line = element_line(colour = 'black', 
                                 size=boxLinesWidth, linetype='solid'),
        strip.background = element_rect(fill = "white", colour = "black", 
                                        size = rel(2)),
        axis.ticks = element_line(size = boxLinesWidth),
        axis.title = element_blank(),
        panel.grid.minor.x = element_line(size = boxLinesWidth*0.8),
        panel.grid.minor.y = element_line(size = boxLinesWidth*0.8),
        legend.position = "none")


# ------ List of plots ------ ----
plts <- list("5" = list(), "7" = list(), "9" = list())
for(s in unique(biasSummary$S)){
  for(a in als){
    
    # subset of the data
    biasSummaryTmp <- biasSummary %>%
      filter(S == s, alpha == a)
    biasSummaryMeansTmp <- biasSummaryMeans %>% 
      filter(S == s, alpha == a)
    
    # boxplots of estimations of N 
    plts[[as.character(s)]][[as.character(a)]]  <- ggplot()+
      geom_hline(aes(yintercept = 500), lty = "dashed", col = "grey") +
      # all points
      geom_point(data = biasSummaryTmp,
                  aes(x=p, y=N_mean, shape = model, color = Convergence),
                  alpha = 0.4) +
      # means of 95% interval limits
      geom_errorbar(data = biasSummaryMeansTmp,
                    aes(x=p, ymin=N_2.5, ymax=N_97.5, group = model), 
                    col = "grey40", width = 0.008) +
      # means
      geom_line(data = biasSummaryMeansTmp,
                aes(x=p, y=N_estim, group = model), 
                size = plotLinesWidth) +
      geom_point(data = biasSummaryMeansTmp,
                 aes(x=p, y=N_estim, shape = model)) +
      # theme
      scale_color_manual(values = c("grey", "darkred")) +
      scale_shape_manual(values=c(15:18))+
      myThemePdf +
      scale_x_continuous(minor_breaks = c(0.11, 0.23, 0.36, 0.51), 
                         breaks = c(0.11, 0.23, 0.36, 0.51)) +
      scale_y_continuous(minor_breaks = seq(200, 1000, 150),
                         breaks = seq(200, 1000, 150),
                         limits = c(minN, maxN) ) 
    
    
  }
}


# ------ Arrange plots  ------ ----
# Legend
p <- biasSummaryMeansTmp %>% 
  ggplot(aes(x=p, y=N_estim))+
  geom_point(aes(shape = model)) +
  scale_shape_manual(values=c(15:18))+
  theme_bw()
legend <- plot_grid(get_legend(p))

# labels of the alpha
mod1 <- ggdraw() + draw_label("Alpha = 0.8", size = labelSize)
mod2 <- ggdraw() + draw_label("Alpha = 0.9", size = labelSize)
mod3 <- ggdraw() + draw_label("Alpha = 0.95", size = labelSize)
modGrid <- plot_grid(mod1, mod2, mod3, ncol = 3)

# label of the number of occasion
s1 <- ggdraw() + draw_label("5", size = labelSize)
s2 <- ggdraw() + draw_label("7", size = labelSize)
s3 <- ggdraw() + draw_label("9", size = labelSize)
sGrid <- plot_grid(s1, s2, s3, ncol = 1)

# label axes
labx <- ggdraw() + draw_label("Simulated lambda", size = 10)
laby <- ggdraw() + draw_label("Pop. size estimates", 
                              angle = 90, size = 10)

# grids of plots
plotGrid5 <- plot_grid(plotlist = plts[["5"]], ncol = 3)
plotGrid7 <- plot_grid(plotlist = plts[["7"]], ncol = 3)
plotGrid9 <- plot_grid(plotlist = plts[["9"]], ncol = 3)
plotGrid <- plot_grid(plotGrid5, plotGrid7, plotGrid9, ncol = 1)

# grid with labels and boxplots
pGrid <- plot_grid(NULL, NULL, modGrid, NULL,
                   sGrid, laby, plotGrid, legend,
                   NULL, NULL, labx, NULL,
                   rel_widths = c(0.05, 0.05, 1, 0.2), 
                   rel_heights = c(0.05, 1, 0.04))

# ------ Plot grid ------ ----

plot(pGrid)

pdf(file=paste("./figures/repObs_allModels_N500Estim.pdf", sep=""), 
    width = 7.8, height = 5.4,
    onefile = TRUE)
  plot(pGrid)
dev.off()

