### Fig 1a - PARAGON-I metaT manuscript ###
### map showing the hawaiin islands and the sampling location ###
### By: Isha Kalra ###
### Last updated: 07/21/2026 ###

library(tidyverse)
library(ggplot2)
library(maps)

## download Hawaii map data
hi_map <- map_data("world", region = "USA") %>%
  filter(subregion == "Hawaii")


# sampling sites
sampling_sites <- data.frame(
  lat = c(21.8319, 21.8006, 21.5790, 21.5054, 21.5488,
          21.7328, 21.9153, 22.0939, 22.2378, 22.3645,
          21.57293, 21.5620, 21.5604, 21.8276),
  lon = c(-155.3138, -155.2511, -155.4180, -155.6786, -155.9476,
          -156.2562, -156.3918, -156.4518, -156.4452, -156.3885,
          -155.3984, -155.3244, -155.91236, -156.2783),
  Sample = c(
    rep("Daily PIT", 10),
    rep("NetTrap", 4)
  ))

## plot the hawaii map with Paragon 1 coordinates and station ALOHA
fig_1a <- ggplot() +
  geom_polygon(data = hi_map, aes(x = long, y = lat, group = group), 
               fill = "grey80", color = "black") +
  geom_point(data = sampling_sites, aes(x = lon, y = lat, shape = Sample), #PARAGON-I sampling
             size = 3, color = "black", fill = "pink3") +
  scale_shape_manual(values = c("Daily PIT" = 21, "NetTrap" = 24))+
  annotate( "text", x = -155.50, y = 22.00, label = "PARAGON-I", 
            size = 3.5,fontface = "bold", color="pink3")+
  geom_point(aes(x = -158.00, y = 22.75), shape = 8, size = 3, color = "black" ) + #Station ALOHA
  annotate( "text", x = -158.00, y = 22.95, label = "Station ALOHA", 
            size = 3.5,fontface = "bold")+
  annotation_north_arrow(location = "tl", height = unit(1, "cm"), # direction arrow
                         width = unit(1, "cm")) +
  annotate("text", x = -155, y = 18.5, label = "Pacific Ocean", 
           fontface = "italic", color = "lightblue3", size = 5)+
  annotate("text", x = -158, y = 20, label = "HAWAII", 
           color = "darkgrey", size = 4)+
  geom_segment(aes(x = -160.5, xend = -159.54, y = 18.8, yend = 18.8), # km bar
               linewidth = 1.2) +
  geom_segment(aes(x = -160.5, xend = -160.5, y = 18.75, yend = 18.85),
               linewidth = 1.2) +
  geom_segment(aes(x = -159.54, xend = -159.54, y = 18.75, yend = 18.85),
               linewidth = 1.2) +
  annotate("text", x = -160.02, y = 18.95, label = "100 km", size = 3) + 
  coord_fixed(xlim = c(-161, -154), ylim = c(18.5, 23), ratio = 1.3) +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude")

#save map
ggsave("fig1a.pdf", plot=fig_1a, width = 6.5, height = 5.75, units = "in")
