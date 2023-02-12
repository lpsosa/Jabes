rm(list = ls())
#   ____________________________________________________________________________
#   Produce plots  of mean hazard ratios by group                           ####

library(reshape)
library(furrr)
library(ggplot2)
library(survival)
library(tidyverse)
library(tictoc)
library(scales)

##  ............................................................................
##  Open multisession                                                       ####

future::plan(multisession)

##  ............................................................................
##  Read data                                                               ####

source("Covariates_time.R")
source("cortes.R")
sepkoski00 = readRDS("sepkoski00.rds")

group_lab = c("Trilobita",
              "Brachiopoda",
              "Gastropoda",
              "Cephalopoda",
              "Anthozoa",
              "Bivalvia")

##  ............................................................................
##  Mass extinctions, geological eras & periods                             ####

massextinct = read.csv("massextinct.csv")
eras = read.csv("eras_AVR.csv")
periods = read.csv("periods_AVR.csv")
periodsaux = data.frame(periods,
                        labelpos = periods$top + (periods$base - periods$top) /
                          2)
periodsaux$labelpos[12] = periodsaux$labelpos[12] - 5

sepkoski00filtered = subset(sepkoski00, sepkoski00$grupoid < 7 &
                              sepkoski00$longevidad > 0)
sepkoski00filtered$FA = ceiling(sepkoski00filtered$FA)

##  ............................................................................
##  t grid for age                                                          ####

nt = 550
mint = 1
maxt = 550
tgrid = seq(from = mint,
            to = maxt,
            length.out = nt)

##  ............................................................................
##  T grid for geological time                                              ####

nT = 520
minT = 1
maxT = 520
Tgrid = seq(from = minT,
            to = maxT,
            length.out = nT)

##  ............................................................................
##  Read fitted model and coefficients                                      ####

# Read fitted Cox model

modelo = readRDS("modelo7.rds")
N_total_params = 6
N_non_group_params = 5
coeficientes = readRDS("coeficientes7.rds")

# Separate parameters by variables

coef_vars = coeficientes[1:N_non_group_params]
coef_grupos = matrix(coeficientes$grupos, nrow = 7, ncol = 6)
coef_grupos[7,] = 0

# Choose to print individual plots

indiv_flag = TRUE
joint_flag = TRUE

##  ............................................................................
##  Compute mean hazard at geological time T                                ####

hazard_T = function(T, modelo, coef_vars, coef_grupos, groupID) {
  # groupID=0 means all data
  # groupID=1,2,3,4,5,6 means specific group
  
  covarsT = covariables_temporales(T)
  datafull = subset(sepkoski00filtered,
                    sepkoski00filtered$FA > Tgrid[T] &
                      sepkoski00filtered$LA < Tgrid[T])
  if (groupID == 0) {
    data = datafull
  } else {
    data = subset(datafull, datafull$grupoid == groupID)
  }
  richness = dim(data)[1]
  
  cumul_hazard = 0
  if (richness > 0) {
    for (j in 1:richness) {
      age = data$FA[j] - T
      beta_vars_aux = vector("list", N_non_group_params)
      for (k in 1:N_non_group_params) {
        if (is.null(modelo[[k]]) == FALSE) {
          beta_vars_aux[[k]] = cortes(age, modelo[[k]], coef_vars[[k]])
        } else {
          beta_vars_aux[[k]] = coeficientes[[k]]
        }
      }
      beta_vars = do.call(cbind, beta_vars_aux)
      
      beta_grupos = cortes(age, modelo[[6]], coef_grupos[data$grupoid[j], ])
      linear_predictor =
        map(1:N_non_group_params, ~ beta_vars[1, .x] * covarsT[.x]) %>%
        cbind.data.frame %>%
        rowSums +
        beta_grupos
      cumul_hazard = cumul_hazard + exp(linear_predictor)
    }
  }
  
  ifelse(richness > 0, return(cumul_hazard / richness), return(NA))
  
}

##  ............................................................................
##  Compute global hazard at T                                              ####

global_hazard = future_map_dbl(
  1:nT,
  hazard_T,
  modelo,
  coef_vars,
  coef_grupos,
  groupID =
    0,
  .progress = TRUE
)


##  ............................................................................
##  Plot individual                                                         ####

grafica_tipo_2 = function(groupID, miny, maxy) {
  plot_title = paste0(group_lab[groupID])
  
  group_hazard =
    future_map_dbl(1:nT,
                   hazard_T,
                   modelo,
                   coef_vars,
                   coef_grupos,
                   groupID = groupID)
  
  ratio = group_hazard / global_hazard
  plotdata = data.frame(T = Tgrid, ratio = ratio)
  
  p = ggplot(data = plotdata, aes(x = T, y = ratio)) +
    
    geom_rect(
      data = massextinct,
      mapping = aes(
        xmin = top,
        xmax = base,
        ymin = -Inf,
        ymax = Inf
      ),
      fill = "gray",
      alpha = .5,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    
    geom_line() +
    
    scale_y_continuous(
      breaks = miny:maxy,
      labels = as.character(miny:maxy),
      limits = c(miny, maxy),
      minor_breaks = NULL
    ) +
    
    scale_x_reverse(breaks = seq(0, 500, 100)) +
    ylab("Ratio of average hazards") +
    xlab("My") +
    geom_hline(yintercept = 1, col = "red") +
    ggtitle(plot_title)
  
  print(p)
  ggsave(
    paste0("Figure6_AvgGrpHazard_", group_lab[groupID], ".png"),
    width = 7,
    height = 4,
    units = "in"
  )
}


##  ............................................................................
##  Plot all groups                                                         ####

grafica_tipo_2_integral = function(miny, maxy) {
  plot_title = ""
  legendtitle = "Group"
  plotdata = data.frame(matrix(NA, nrow = length(Tgrid), ncol = 7))
  names(plotdata) = c("T", group_lab)
  plotdata[, 1] = Tgrid
  for (groupID in 1:6) {
    group_hazard =
      future_map_dbl(1:nT,
                     hazard_T,
                     modelo,
                     coef_vars,
                     coef_grupos,
                     groupID = groupID)
    
    ratio = group_hazard / global_hazard
    plotdata[, groupID + 1] = ratio
  }
  
  plotdatalong = melt(plotdata,
                      id = "T",
                      measured = group_lab)
  names(plotdatalong) = c("T", "Group", "value")
  
  p = ggplot(data = plotdatalong, aes(x = T, y = value, col = Group)) +
    guides(linetype = guide_legend(legendtitle)) +
    geom_rect(
      data = massextinct,
      mapping = aes(
        xmin = top,
        xmax = base,
        ymin = -Inf,
        ymax = Inf
      ),
      fill = "gray",
      alpha = .5,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    
    geom_line() +
    
    scale_y_continuous(
      breaks = miny:maxy,
      labels = as.character(miny:maxy),
      limits = c(miny - 0.25, maxy),
      expand = c(0, 0),
      minor_breaks = NULL
    ) +
    scale_color_manual(values = palette_groups[1:6]) +
    scale_x_reverse(breaks = seq(0, 500, 100), expand = c(0.01, 0)) +
    ylab("Ratio of average hazards") +
    xlab("Geologic Time (Ma)") +
    geom_hline(yintercept = 1, col = palette_groups[7]) +
    ggtitle(plot_title) +
    
    # Geologic periods labels
    
    geom_rect(
      data = periods,
      mapping = aes(
        xmin = top,
        xmax = base,
        ymin = -.25,
        ymax = -.125
      ),
      fill = periods$color,
      color = "gray",
      alpha = 1,
      size = .25,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    geom_text(
      data = periodsaux,
      aes(
        x = labelpos,
        y = -(.125 + .25) / 2,
        label = abbrev
      ),
      inherit.aes = FALSE,
      show.legend = FALSE,
      size = 1
    ) +
    
    # Geologic eras labels
    
    geom_rect(
      data = eras,
      mapping = aes(
        xmin = top,
        xmax = base,
        ymin = -.125,
        ymax = -0
      ),
      fill = eras$color,
      color = "gray",
      alpha = 1,
      size = .25,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    geom_text(
      data = eras,
      aes(
        x = top + (base - top) / 2,
        y = -.125 / 2,
        label = name
      ),
      inherit.aes = FALSE,
      show.legend = FALSE,
      size = 1
    ) +
    
    theme_bw()
  
  print(p)
  ggsave(
    paste0("Figure6_AvgGrpHazard.png"),
    height = 4,
    width = 7,
    dpi = 600,
    scale = 1
  )
}

#   ____________________________________________________________________________
#   Code for producing individual plots                                    ####

if (indiv_flag) {
  miny = 0
  maxy = 4
  grafica_tipo_2(1, miny, maxy)
  grafica_tipo_2(2, miny, maxy)
  grafica_tipo_2(3, miny, maxy)
  grafica_tipo_2(4, miny, maxy)
  grafica_tipo_2(5, miny, maxy)
  grafica_tipo_2(6, miny, maxy)
}

#   ____________________________________________________________________________
#   Code for producing all groups plots                                     ####

# Color scales for groups
if (joint_flag) {
  palette_groups = c(hue_pal()(6), "red")
  names(palette_groups) = c(group_lab, "all")
  
  miny = 0
  maxy = 4
  grafica_tipo_2_integral(miny, maxy)
}

##  ............................................................................
##  Close multisession                                                      ####

future::plan(sequential)
