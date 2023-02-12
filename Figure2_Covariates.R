rm(list = ls())
#   ____________________________________________________________________________
#   Plot of covariates used in model                                        ####


##  ............................................................................
##  Define and prepare covariates as functions of T                         ####

plotflag = FALSE
library(ggplot2)
library(reshape)

source("Covariates_time.R")

##  ............................................................................
##  Prepare matrix and function for covariates                              ####

COVARIABLES = cbind(
  nivelmar = nivelmar(1:536),
  sin1 = sinusoidal_1(1:536),
  sin2 = sinusoidal_2(1:536),
  evol = evolucion(1:536),
  co2 = CO2(1:536),
  sst = SST(1:536),
  sstdiff = abs(SSTdiff(1:536))
)

T = 1:536

COVARIABLESstandard = data.frame(T, scale(COVARIABLES))

names(COVARIABLESstandard) = c("T",
                               "Sea Level",
                               "sine 62",
                               "sine 140",
                               "fitness",
                               "CO2",
                               "SST",
                               "SSTdiff")


##  ............................................................................
##  Mass extinctions, geological eras & periods                             ####

massextinct = read.csv("massextinct.csv")
eras = read.csv("eras_AVR.csv")
periods = read.csv("periods_AVR.csv")

periodsaux = data.frame(periods,
                        labelpos = periods$top + (periods$base - periods$top) /
                          2)
periodsaux$labelpos[12] = periodsaux$labelpos[12] - 5


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Covariate names                                                         ####

sep = 8

var_labels = data.frame(
  x = rep(-40, 7),
  y = seq(sep, 7 * sep, sep),
  name = c(
    "Sea Level",
    "sine 62",
    "sine 140",
    "Fitness",
    "CO2",
    "SST",
    "SSTdiff"
  )
)

#   ____________________________________________________________________________
#   Main plot                                                               ####

covs = melt(
  COVARIABLESstandard,
  id = "T",
  measured = c("nivelmar", "sin1", "sin2", "evol", "co2", "sst", "sstdiff")
)
covs_trans = covs

levs = levels(covs$variable)

sep = 8

for (i in 1:7) {
  indics = which(covs[, 2] == levs[i])
  covs_trans[indics, 3] = covs[indics, 3] + sep * i
}

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### ggplot production                                                       ####

p2 = ggplot(data = covs_trans, aes(x = T, y = value, col = variable)) +
  theme(
    axis.text.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.background = element_blank(),
    axis.line.x.bottom = element_line(color = "blue"),
    panel.grid.major.x = element_line(color = "gray")
  ) +
  ylab("Standardized values") +
  xlab("Geologic Time (Ma)") +
  ylim(0, 8 * sep) +
  
  geom_rect(
    data = massextinct,
    mapping = aes(
      xmin = top,
      xmax = base,
      ymin = 0,
      ymax = Inf
    ),
    fill = "lightblue",
    alpha = 1,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  
  guides(color = guide_legend(reverse = TRUE)) +
  scale_x_reverse(breaks = seq(0, 500, 100), expand = c(0, 40)) +
  geom_hline(yintercept = seq(sep, 7 * sep, sep), color = "black") +
  geom_line(size = 1, show.legend = FALSE) +
  
  # Draw geological periods
  
  geom_rect(
    data = periods,
    mapping = aes(
      xmin = top,
      xmax = base,
      ymin = 0,
      ymax = 1
    ),
    fill = periods$color,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  geom_text(
    data = periodsaux,
    aes(x = labelpos,
        y = .5,
        label = abbrev),
    inherit.aes = FALSE,
    show.legend = FALSE,
    size = 2
  ) +
  
  # Draw geological eras
  
  geom_rect(
    data = eras,
    mapping = aes(
      xmin = top,
      xmax = base,
      ymin = 1,
      ymax = 2
    ),
    fill = eras$color,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  geom_text(
    data = eras,
    aes(
      x = top + (base - top) / 2,
      y = 1.5,
      label = name
    ),
    inherit.aes = FALSE,
    show.legend = FALSE,
    size = 2
  ) +
  geom_label(aes(x = x, y = y, label = name),
             data = var_labels,
             inherit.aes = FALSE)

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Print plot                                                              ####

print(p2)

ggsave(
  p2,
  scale = 1,
  width = 8,
  height = 10,
  dpi = 600,
  filename = "Figure2_Covariates.png"
)
