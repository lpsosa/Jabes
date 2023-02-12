rm(list = ls())
#   ____________________________________________________________________________
#   Kaplan-Meier contour plots                                              ####

require(ggplot2)
require(survival)


##  ............................................................................
##  Read data                                                               ####

sepkoski00 = readRDS("sepkoski00.RDS")

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Mass extinctions, periods & eras                                        ####

massextinct = read.csv("massextinct.csv")
eras = read.csv("eras_AVR.csv")
periods = read.csv("periods_AVR.csv")
periodsaux = data.frame(periods,
                        labelpos = periods$top + (periods$base - periods$top) /
                          2)
periodsaux$labelpos[12] = periodsaux$labelpos[12] - 5

# t grid for age

nt = 500
mint = 1
maxt = 550
tgrid = seq(from = mint,
            to = maxt,
            length.out = nt)

# T grid for geological time

nT = 536
minT = 1
maxT = 536
Tgrid = seq(from = minT,
            to = maxT,
            length.out = nT)

zvals = matrix(NA, nrow = nT, ncol = nt)
media = rep(NA, nT)
CI95L = rep(NA, nT)
CI95U = rep(NA, nT)

##  ............................................................................
##  Kaplan-Meier                                                            ####

for (i in 1:nT) {
  fit = survfit(
    Surv(longevidad, event = observed) ~ 1,
    data = sepkoski00,
    subset = (FA > Tgrid[i] & LA < Tgrid[i])
  )
  surv = summary(fit, times = tgrid, extend = T)
  zvals[i, ] = surv$surv
  mtemp = survival:::survmean(fit, rmean = 550)
  media[i] = mtemp[[1]]["*rmean"]
  CI95L[i] = media[i] - 1.96 * mtemp[[1]]["*se(rmean)"]
  CI95U[i] = media[i] + 1.96 * mtemp[[1]]["*se(rmean)"]
}

datos = data.frame(expand.grid(T = Tgrid, t = tgrid), z = as.vector(zvals))
datosmedia = data.frame(
  T = Tgrid,
  meansurv = media,
  L95 = CI95L,
  U95 = CI95U
)

#   ____________________________________________________________________________
#   ggplot definition                                                       ####

p = ggplot(datos, aes(x = T, y = t, z = z)) + geom_contour_filled() +
  scale_x_reverse(breaks = seq(0, 500, 50)) +
  scale_y_continuous(breaks = seq(0, 500, 50)) +
  geom_rect(
    data = massextinct,
    mapping = aes(
      xmin = top,
      xmax = base,
      ymin = 0,
      ymax = 550
    ),
    fill = "transparent",
    color = "white",
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  xlab("Geologic Time (Ma)") +
  ylab("Genera age (My)") +
  labs(fill = "Survival") +
  
  # Geologic periods labels
  
  geom_rect(
    data = periods,
    mapping = aes(
      xmin = top,
      xmax = base,
      ymin = -50,
      ymax = -25
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
    aes(x = labelpos,
        y = -37.5,
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
      ymin = -25,
      ymax = 0,
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
      y = -12.5,
      label = name
    ),
    inherit.aes = FALSE,
    show.legend = FALSE,
    size = 1.5
  ) +
  theme_bw(base_size = 8) +
  theme(panel.border = element_blank())

#   ____________________________________________________________________________
#   Print plot                                                              ####

print(p)
ggsave(
  "Figure3_KMcontours.png",
  scale = 1,
  width = 7,
  height = 4,
  units = "in",
  dpi = 600
)
