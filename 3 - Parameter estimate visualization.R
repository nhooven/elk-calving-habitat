# Title: Parturition site selection
# Subtitle: 3a - Parameter estimate plots (MCP)
# Author: Nathan D. Hooven
# Email: nathan.d.hooven@gmail.com
# Affiliation: Department of Forestry and Natural Resources, University of Kentucky
# Date began: 26 Apr 2021
# Date completed: 26 Apr 2021
# Date modified: 16 Aug 2022
# R version: 3.6.2

#_____________________________________________________________________________________________________________
# 1. Load in required packages ----
#_____________________________________________________________________________________________________________

library(tidyverse)
library(glmmTMB)     # modeling

#_____________________________________________________________________________________________________________
# 2. Read in .Rdata ----
#_____________________________________________________________________________________________________________

load("partsite_models.RData")

#_____________________________________________________________________________________________________________
# 3. Get data into correct format ----
#_____________________________________________________________________________________________________________

# top model 1
# extract estimates, SEs, and names
model.1.est <- fixef(final.models[[15]])$cond
model.1.se <- sqrt(diag(final.models[[15]]$s$cov.fixed))
model.1.names <- names(fixef(final.models[[15]])$cond)

# bind together
model.1.par <- data.frame(var = model.1.names,
                          est = model.1.est,
                          se = model.1.se,
                          hyp = c("INT",
                                  "FORCOV", "FORCOV", "FORCOV", 
                                  "LAND", "LAND", "LAND",
                                  "TOPO", "TOPO", "TOPO",
                                  "HUMAN"),
                          model = "1")

# subset out intercepts
model.1.par <- model.1.par %>% filter(var != "(Intercept)")

# reorder and relabel factor
model.1.par$var <- factor(model.1.par$var,
                          levels = rev(c("(Intercept)", 
                                         "EVI.2pix", "canopy", "canopy.Q",
                                         "IJI.250", "IJI.250.Q", "ED.500",
                                         "TRI", "TPI.250", "SRI.500.P",
                                         "dRoad.1000")))

# top model 2
# extract estimates, SEs, and names
model.2.est <- fixef(final.models[[13]])$cond
model.2.se <- sqrt(diag(final.models[[13]]$s$cov.fixed))
model.2.names <- names(fixef(final.models[[13]])$cond)

# bind together
model.2.par <- data.frame(var = model.2.names,
                          est = model.2.est,
                          se = model.2.se,
                          hyp = c("INT",
                                  "FORCOV", "FORCOV", "FORCOV", 
                                  "TOPO", "TOPO", "TOPO",
                                  "HUMAN"),
                          model = "2")

# subset out intercepts
model.2.par <- model.2.par %>% filter(var != "(Intercept)")

# reorder and relabel factor
model.2.par$var <- factor(model.2.par$var,
                          levels = rev(c("(Intercept)", 
                                         "EVI.2pix", "canopy", "canopy.Q",
                                         "TRI", "TPI.250", "SRI.500.P",
                                         "dRoad.1000")))

# bind together
both.models.par <- rbind(model.1.par, model.2.par)

# reorder "model" factor
both.models.par$model <- factor(both.models.par$model,
                          levels = rev(c("1", "2")))

#_____________________________________________________________________________________________________________
# 4. Horizontal plot ----
#_____________________________________________________________________________________________________________

ggplot(data = both.models.par, aes(x = est, 
                                   y = var, 
                                   group = model)) +
  geom_vline(xintercept = 0) +
  theme_bw() +
  geom_errorbar(aes(xmin = est - 1.645*se, 
                    xmax = est + 1.645*se,
                    color = hyp),
                width = 0,
                size = 1.25,
                position = position_dodge(width = 0.8)) +
  geom_point(aes(fill = hyp, 
                 shape = model),
             size = 2, 
             stroke = 1,
             position = position_dodge(width = 0.8)) +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        panel.grid = element_blank()) +
  xlab("Standardized selection coefficient") +
  scale_y_discrete(labels = rev(expression(EVI[1000], Canopy, Canopy^2,
                                           IJI[250], IJI[250]^2, ED[500],
                                           TRI, TPI[250], ln(SRI[500]),
                                           dRoad[1000]))) +
  scale_shape_manual(values = c(24, 21)) +
  scale_x_continuous(breaks = seq(-3, 3, 1)) +
  # hyp order: FORCOV, HUMAN, LAND, TOPO
  scale_color_manual(values = c("#66CC33", "#FF3300", "#3399FF", "#006633")) +
  scale_fill_manual(values = c("#66CC33", "#FF3300", "#3399FF", "#006633")) +
  geom_hline(yintercept = 7.5, linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = 4.5, linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = 1.5, linetype = "dashed", alpha = 0.5)

