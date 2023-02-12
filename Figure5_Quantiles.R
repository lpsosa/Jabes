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
##   Read fitted Cox model                                                  ####

modelo = readRDS("modelo7.rds")
N_total_params = 6
N_non_group_params = 5
coeficientes = readRDS("coeficientes7.rds")

# Separate parameters by variables

coef_vars = coeficientes[1:N_non_group_params]
coef_grupos = matrix(coeficientes$grupos, nrow = 7, ncol = 6)[1:6, ]

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

# Prepare group-wise quantiles

quantiles_grp = c(.05, .25, .75, .95)
n_quantiles_grp = length(quantiles_grp)
source("Lifetime_quantiles.R")

# Prepare age slice plots per group

t_slices = readRDS("lifetime_quantiles.rds")
n_slices = dim(t_slices)[2]
ind_slices = matrix(findInterval(t_slices, tgrid),
                    nrow = 6,
                    ncol = n_slices)
slice_data = vector(mode = "list", length = 6)
age_colors = brewer.pal(n = n_slices, name = "RdBu")

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Choose to print individual plots                                        ####

indiv_flag = TRUE
joint_flag = TRUE

#   ____________________________________________________________________________
#   Main plotting function                                                  ####

minz = 9999999
maxz = -9999999

grafica_cuantiles = function(minz_share,
                             maxz_share,
                             n_zbreaks) {
  plot_names = paste0("p", 1:6)
  
  for (groupID in 1:6) {
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
    
    slice_data_aux = data.frame(zvals[, ind_slices[groupID, ]], T = Tgrid)
    names(slice_data_aux)[1:n_slices] = sprintf("%.1f", t_slices[groupID, ])
    slice_data[[groupID]] = melt(
      slice_data_aux,
      id = "T",
      measured = 1:6,
      variable_name = "age"
    )
    
    # Plot of log(predictor) by age by group
    
    p_tipo3 =
      ggplot(slice_data[[groupID]], aes(x = T, y = value, col = age)) +
      geom_line(size = 0.8, alpha = .5) +
      scale_x_reverse(breaks = seq(0, 500, 50)) +
      ylab("log(predictor)") +
      xlab("Geologic Time (Ma)") +
      scale_color_discrete(type = rev(age_colors),
                           name = "Genera \nage (My)") +
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
      ggtitle(plot_title) +
      theme_bw() +
      theme(text = element_text(size = 8))
    
    if (indiv_flag) {
      print(p_tipo3)
      ggsave(
        paste0("Figure5_Quartiles_", group_lab[groupID], ".png"),
        width = 7,
        height = 4
      )
    }
    assign(plot_names[groupID], p_tipo3)
  }
  
  if (joint_flag) {
    png(
      filename = paste0("Figure5_Quantiles.png"),
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

grafica_cuantiles(minz_manual, maxz_manual, n_breaks)
