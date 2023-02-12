rm(list = ls())
#   ____________________________________________________________________________
#   Produce Type 1 plots: theoretical hazard per group                      ####

source("cortes.R")
source("Covariates_time.R")
require(tidyverse)
require(reshape)
require(RColorBrewer)
library(patchwork)

##  ............................................................................
##  Group labels                                                            ####

group_lab = c(
  "Trilobita",
  "Brachiopoda",
  "Gastropoda",
  "Cephalopoda",
  "Anthozoa",
  "Bivalvia",
  "otros"
)

##  ............................................................................
##  Mass extinctions                                                        ####

massextinct = read.csv("massextinct.csv")

##  ............................................................................
##  Age range                                                               ####

nt = 550
mint = 1
maxt = 80
tgrid = seq(from = mint,
            to = maxt,
            length.out = nt)

##  ............................................................................
##  Geologic time range                                                     ####

nT = 520
minT = 1
maxT = 520
Tgrid = seq(from = minT,
            to = maxT,
            length.out = nT)

##  ............................................................................
##  Choose model: 6=six groups, 7=seven groups                              ####

model_choice = 7
filename_suffix = ifelse(model_choice == 7, "_7grp", "_6grp")

##  ............................................................................
##   Read fitted Cox model                                                  ####

if (model_choice == 7) {
  modelo = readRDS("modelo7.rds")
  N_total_params = 6
  N_non_group_params = 5
  coeficientes = readRDS("coeficientes7.rds")
}

if (model_choice == 6) {
  modelo = readRDS("modelo6.rds")
  N_total_params = 6
  N_non_group_params = 5
  coeficientes = readRDS("coeficientes6.rds")
}

# Separate parameters by variables

coef_vars = coeficientes[1:N_non_group_params]
coef_grupos = matrix(coeficientes$grupos, nrow = model_choice, ncol = 6)[1:6, ]

# If model_choice==6 set bivalvia coeffs to zeros

if (model_choice == 6) {
  coef_grupos[6, ] = 0
}

# Evaluate regression parameters by intervals

beta_vars_aux = vector("list", N_non_group_params)
for (j in 1:N_non_group_params) {
  if (is.null(modelo[[j]]) == FALSE) {
    beta_vars_aux[[j]] = cortes(tgrid, modelo[[j]], coef_vars[[j]])
  } else {
    beta_vars_aux[[j]] = rep(coeficientes[[j]], nt)
  }
}
beta_vars = do.call(cbind, beta_vars_aux)

beta_grupos = matrix(NA, nrow = nt, ncol = 6)
for (j in 1:6) {
  beta_grupos[, j] = cortes(tgrid, modelo[[6]], coef_grupos[j, ])
}

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Choose to print individual plots                                        ####

indiv_flag = TRUE
joint_flag = TRUE

#   ____________________________________________________________________________
#   Main plotting function                                                  ####

minz = 9999999
maxz = -9999999

grafica_hazard = function(minz_share,
                          maxz_share,
                          n_zbreaks)
{
  plot_names = paste0("p", 1:6)
  
  for (groupID in 1:6) {
    # Breaks in contour plots
    
    mybreaks <-
      seq(minz_share, maxz_share, length.out = n_zbreaks + 1)
    
    # Function to return the desired number of colors
    
    mycolors <- function(x) {
      colors <- colorRampPalette(c("darkblue", "white"))(x)
      return(colors)
    }
    
    # Function to create labels for legend
    
    breaklabel <- function(x) {
      Labels <-
        paste0(
          formatC(mybreaks[1:n_zbreaks], digits = 1, format = "f"),
          ", ",
          formatC(mybreaks[2:(n_zbreaks + 1)], digits = 1, format = "f")
        )
      Labels[1:x]
    }
    
    plot_title = paste0(group_lab[groupID])
    
    zvals = matrix(NA, nrow = nT, ncol = nt)
    
    for (T in 1:nT) {
      coef_grupos = beta_grupos[, groupID]
      
      covarsT = covariables_temporales(T)
      
      linear_predictor =
        map(1:N_non_group_params, ~ beta_vars[, .x] * covarsT[.x]) %>%
        cbind.data.frame %>%
        rowSums
      linear_predictor = linear_predictor + coef_grupos
      
      zvals[T, ] = log(linear_predictor)
      
    }
    minz <<- min(c(minz, zvals), na.rm = TRUE)
    maxz <<- max(c(maxz, zvals), na.rm = TRUE)
    
    plotdata = data.frame(expand.grid(T = Tgrid, t = tgrid),
                          z = as.vector(zvals))
    
    # Theoretical contour plot
    
    p_tipo1 = ggplot(plotdata, aes(x = T, y = t, z = z)) +
      
      geom_contour_filled(breaks = mybreaks, show.legend = TRUE) +
      
      scale_fill_manual(
        palette = mycolors,
        labels = breaklabel(10),
        drop = FALSE,
        name = "log(predictor)"
      ) +
      
      scale_x_reverse(breaks = seq(0, 500, 50), expand = c(0, 0)) +
      scale_y_continuous(breaks = seq(0, 100, 10), expand = c(0, 0)) +
      
      geom_rect(
        data = massextinct,
        mapping = aes(
          xmin = top,
          xmax = base,
          ymin = 0,
          ymax = 80
        ),
        fill = "transparent",
        color = "red",
        linewidth = .25,
        inherit.aes = FALSE,
        show.legend = FALSE
      ) +
      ylab("Genera age (My)") +
      xlab("Geologic Time (Ma)") +
      ggtitle(plot_title) +
      theme_bw() +
      theme(text = element_text(size = 8))
    
    if (indiv_flag) {
      print(p_tipo1)
      ggsave(
        paste0("Figure4_Hazards_", group_lab[groupID],
               ".png"),
        width = 7,
        height = 4
      )
    }
    assign(plot_names[groupID], p_tipo1)
  }
  
  if (joint_flag) {
    png(
      filename = "Figure4_Hazards.png",
      width = 8,
      height = 10,
      units = "in",
      res = 600
    )
    print(p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 2))
    dev.off()
  }
}

#   ____________________________________________________________________________
#   Program calls for producing all plots                                   ####

n_breaks = 10
minz_manual = -8
maxz_manual = 2

grafica_hazard(minz_manual, maxz_manual, n_breaks)
