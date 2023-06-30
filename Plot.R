library(ggplot2)
setwd("~/Desktop/Stages/Stage M2 CIRAD/Programme/community-generator/Results")

MEAN_PDI <- read.csv("mean_pdi.csv")
SD_PDI <- read.csv("sd_pdi.csv")
MEAN_SIGMA <- read.csv("mean_sigma.csv")
SD_SIGMA <- read.csv("sd_sigma.csv")

ggplot(MEAN_PDI, aes(Age, r, fill= Mean_PDI)) + 
  scale_x_continuous(breaks=seq(1, 15, 2)) +
  facet_grid(Guild ~ Model) +
  scale_fill_viridis_c(begin = 0, end = 1) +
  geom_tile()

ggplot(SD_PDI, aes(Age, r, fill= SD_PDI)) + 
  scale_x_continuous(breaks=seq(1, 15, 2)) +
  facet_grid(Guild ~ Model) +
  scale_fill_gradient(low="white", high="red", #colors in the scale
                      breaks=seq(0,max(SD_PDI$SD_PDI, na.rm = T),0.05), #breaks in the scale bar
                      limits=c(0, max(SD_PDI$SD_PDI, na.rm = T))) +
  geom_tile()

ggplot(MEAN_SIGMA, aes(Age, r, fill= Mean_Sigma)) + 
  scale_x_continuous(breaks=seq(1, 15, 2)) +
  facet_grid(Interaction ~ Model) +
  scale_fill_gradient2(low="red", mid="white", high="blue", #colors in the scale
                       midpoint=0,    #same midpoint for plots (mean of the range)
                       breaks=seq(-1,1,0.25), #breaks in the scale bar
                       limits=c(min(MEAN_SIGMA$Mean_Sigma), max(MEAN_SIGMA$Mean_Sigma))) +
  geom_tile()

ggplot(SD_SIGMA, aes(Age, r, fill= SD_Sigma)) + 
  scale_x_continuous(breaks=seq(1, 15, 2)) +
  facet_grid(Interaction ~ Model) +
  scale_fill_gradient(low="white", high="red", #colors in the scale
                      breaks=seq(0,max(SD_SIGMA$SD_Sigma),0.05), #breaks in the scale bar
                      limits=c(0, max(SD_SIGMA$SD_Sigma))) +
  geom_tile()
