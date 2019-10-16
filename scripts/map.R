

library(dplyr)
library(raster)
library(ggplot2)
library(ggthemes)
library(ggsn)
library(ggrepel) # to label points

mapdata <- getData("GADM", country = "usa", level = 1)
mymap <- fortify(mapdata)

lat <- c(37.662276, 37.970940, 38.096397 )
long <- c(-122.37413, -122.481423,-122.833672)
Name <- c("OP","LL", "TB")
mypoint <- data.frame(long = long, lat = lat, Name=Name)

g1 <- ggplot() +
      geom_blank(data = mymap, aes(x=long, y=lat)) +
      geom_map(data = mymap, map = mymap, 
               aes(group = group, map_id = id),
               fill = "#b2b2b2", color = "black", size = 0.3) +
      geom_point(data = mypoint, aes(x = long, y = lat, shape=Name),
                 color = "black", size = 7, bg="black") +
      scale_shape_manual(values=c( 21,22,24))+
      scale_x_continuous(limits = c(-123.2, -121.8), expand = c(0, 0)) +
      scale_y_continuous(limits = c(37.5, 38.3), expand = c(0, 0)) +
      theme_map()+
      ggsn::scalebar(location = "bottomleft", dist=10, dist_unit = "km",
               model = "WGS84", transform=TRUE,
               x.min = -123, x.max = -122,
               y.min = 37.55, y.max = 38.4, cex=0.2, st.size = rel(2.5)) +
      theme(legend.position = "none")

g1

states <- map_data("state")

ca_df <- states %>%
  filter(region == "california")

  g2 <- ggplotGrob(
          ggplot() +
          geom_polygon(data = ca_df,
                       aes(x = long, y = lat, group = group),
                       fill = "#b2b2b2", color = "black", size = 0.3) +
          geom_rect(aes(xmin=-123.1, xmax=-121.4, ymin=37.2, ymax=38.5), 
                color="black", fill=NA, lwd=0.7) +
          coord_map("polyconic") +
          theme_map() +
          theme(panel.background = element_rect(fill = NULL))
        )

  g3 <- g1 +
        annotation_custom(grob = g2, xmin = -122.1, xmax = -121.8,
                          ymin = 37.55, ymax = 37.9)
  g3


ggsave("~/map.pdf", g3, width=110, height=90, units="mm")


