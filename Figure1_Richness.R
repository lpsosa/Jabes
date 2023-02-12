rm(list = ls())
#   ____________________________________________________________________________
#  Plot richness per group                                                  ####

library(furrr)
library(ggplot2)
library(survival)
library(tidyverse)
library(tictoc)
library(reshape)
library(scales)

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Open multisession                                                       ####

future::plan(multisession)

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Read data                                                               ####

sepkoski00 = readRDS("sepkoski00.rds")

group_lab = c("Trilobita",
              "Brachiopoda",
              "Gastropoda",
              "Cephalopoda",
              "Anthozoa",
              "Bivalvia")

sepkoski00filtered = subset(sepkoski00, sepkoski00$grupoid < 7 &
                              sepkoski00$longevidad > 0)
sepkoski00filtered$FA = ceiling(sepkoski00filtered$FA)

##  ............................................................................
##  Mass extinctions, geological eras & periods                             ####

massextinct = read.csv("massextinct.csv")
eras = read.csv("eras_AVR.csv")
periods = read.csv("periods_AVR.csv")
periodsaux = data.frame(periods,
                        labelpos = periods$top + (periods$base - periods$top) /
                          2)
periodsaux$labelpos[12] = periodsaux$labelpos[12] - 5

nT = 520
minT = 1
maxT = 520
Tgrid = seq(from = minT,
            to = maxT,
            length.out = nT)

##  ............................................................................
##   Function to compute richness at geological time T                      ####

richness_T = function(T, groupID) {
  # groupID=0 means all data
  # groupID=1,2,3,4,5 means a specific group
  datafull = subset(sepkoski00filtered, FA > Tgrid[T] &
                      LA < Tgrid[T])
  if (groupID == 0) {
    data = datafull
  } else {
    data = subset(datafull, datafull$grupoid == groupID)
  }
  return(dim(data)[1])
}


##  ............................................................................
##  Function for richness plot all groups                                   ####

richness_all = function() {
  plot_title = ""
  legendtitle = "Group"
  
  plotdata = data.frame(matrix(NA, nrow = nT, ncol = 8))
  plotdata[, 1] = Tgrid
  names(plotdata) = c("T", group_lab, "all")
  
  for (groupID in 1:6) {
    richness = future_map_dbl(1:nT,
                              richness_T,
                              groupID =
                                groupID,
                              .progress = TRUE)
    plotdata[, groupID + 1] = richness
  }
  
  plotdata[, 8] = future_map_dbl(1:nT,
                                 richness_T,
                                 groupID =
                                   0,
                                 .progress = TRUE)
  
  cat("Maximum richness: ", max(plotdata[, 8]), "\n")
  
  plotdatalong = melt(plotdata,
                      id = "T",
                      measured = c(group_lab, "all"))
  names(plotdatalong) = c("T", "Group", "value")
  
  # ggplot production
  
  p = ggplot(data = plotdatalong,
             aes(
               x = T,
               y = value,
               col = Group,
               linetype = Group,
               size = Group
             )) +
    
    geom_rect(
      data = massextinct,
      mapping = aes(
        xmin = top,
        xmax = base,
        ymin = 0,
        ymax = Inf
      ),
      fill = "lightgray",
      alpha = .5,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    
    geom_line() +
    scale_x_reverse(breaks = seq(0, 500, 100), expand = c(0.01, 0)) +
    ylab("Richness (number of genera)") +
    xlab("Geologic Time (Ma)") +
    scale_y_continuous(limits = c(-200, 2566), expand = c(0, 0)) +
    scale_linetype_manual(values = rep("solid", 7)) +
    scale_size_manual(values = c(rep(.4, 6), 1)) +
    scale_color_manual(values = palette_groups) +
    guides(
      linetype = guide_legend(legendtitle),
      size = guide_legend(legendtitle),
      color = guide_legend(legendtitle)
    ) +
    ggtitle(plot_title) +
    theme_bw(base_size = 8) +
    theme(panel.grid.minor.x = element_blank(),
          panel.border = element_blank()) +
    
    # Geologic periods labels
  
  geom_rect(
    data = periods,
    mapping = aes(
      xmin = top,
      xmax = base,
      ymin = -200,
      ymax = -100
    ),
    fill = periods$color,
    color = "gray",
    alpha = 1,
    size = .2,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
    geom_text(
      data = periodsaux,
      aes(x = labelpos,
          y = -150,
          label = abbrev),
      inherit.aes = FALSE,
      show.legend = FALSE,
      size = 1.5
    ) +
    
    # Geologic eras labels
  
  geom_rect(
    data = eras,
    mapping = aes(
      xmin = top,
      xmax = base,
      ymin = -100,
      ymax = 0,
    ),
    fill = eras$color,
    color = "gray",
    alpha = 1,
    size = .2,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
    geom_text(
      data = eras,
      aes(
        x = top + (base - top) / 2,
        y = -50,
        label = name
      ),
      inherit.aes = FALSE,
      show.legend = FALSE,
      size = 1.5
    )
  
  print(p)
  ggsave(
    paste0("Figure1_Richness.png"),
    dpi = 600,
    width = 7,
    height = 4,
    scale = 1
  )
}

##  ............................................................................
##  Code for all groups plot                                                ####

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Color scales for groups                                                 ####

palette_groups = c(hue_pal()(6), "red")
names(palette_groups) = c(group_lab, "all")

richness_all()

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Close multisession                                                      ####

future::plan(sequential)
