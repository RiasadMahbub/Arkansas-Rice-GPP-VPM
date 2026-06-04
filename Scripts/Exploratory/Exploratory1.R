### VPM
###Load all the required libraries
library(raster)
library(tidyverse)
library(ggplot2)
library(tidyr)
library(dplyr)
library(readr)
library(tidyverse)
library(terra)
library(lubridate)
library(maptools)  ## For wrld_simpl
library(sp)
library(rgdal)
library(cowplot)
library(ggsn)
library(ggpubr)
library(sf)
library(gridExtra)
library(grid)
library(rasterVis)
library(scales)
library(viridis)  # better colors for everyone
library(ggthemes) # theme_map()
library(RColorBrewer)
library(magrittr)
library(dplyr)
library(conflicted)
library(exactextractr)
library(ggspatial)
library(osmdata)


library(conflicted)

# Check conflicts
conflict_scout()

#### all raster using rast function
#first import all files in a single folder as a list 
raster_dir <- "C:/Users/rbmahbub/Documents/Data/GeospatialData/CumulativeVPM/CumulativeVPMRasterPolygonCoverageFilter"
rastlist <- list.files(path = raster_dir, pattern='.tif$', 
                       all.files=TRUE, full.names=FALSE)
# Read the shapefile using st_read
ME2 <- st_read("C:/Users/rbmahbub/Documents/Data/GeospatialData/CumulativeVPM/ME_Shapefile/ME.shp")

setwd(raster_dir)
#import all raster files in folder using lapply
allrasters <- lapply(rastlist, raster)
stacked_raster <- stack(allrasters) ### all raster stack
stacked_raster_rast<-rast(stacked_raster)
#raster::calc
stacked_raster_rast_mean <- terra::app(stacked_raster_rast, "mean", na.rm=TRUE)
names(stacked_raster_rast_mean) <- "GPP"
### Create mean raster data
cumulativerasdf <- as.data.frame(stacked_raster_rast_mean,xy=TRUE)%>%drop_na()
## multiply by 8
cumulativerasdf$mean_8<-cumulativerasdf$GPP*8

# Output file to the exported files
output_file <- "C:/Users/rbmahbub/Documents/RProjects/VPMmodel/VPMmodel/Data/ExportedData/CSVCumulative2008_2020.csv"
write.csv(cumulativerasdf, file = output_file, row.names = FALSE)

# Calculate mean and standard deviation of mean_8
mean_mean_8 <- mean(cumulativerasdf$mean_8)
sd_mean_8 <- sd(cumulativerasdf$mean_8)

# Identify values that are 3 standard deviations away from the mean
outlier_threshold <- 3 * sd_mean_8
outliers <- which(abs(cumulativerasdf$mean_8 - mean_mean_8) > outlier_threshold)

# Remove outliers
cleaned_cumulativerasdf <- cumulativerasdf[-outliers, ]
nrow(cumulativerasdf)
nrow(cleaned_cumulativerasdf)

##########PLOTTING  CUMULATIVE ####################
## Plotting cumulative old technique
### cumulativerasdf normal GPP data
### cleaned_cumulativerasdf: 3 standard deviation filtered data
nameColor <- bquote(atop(Mean~Cumulative ~GPP~(2008-2020)~(g~C~m^-2~year^-1)~"  "))
my_breaks <- c(1100, 1500, 1800, 2100, 2500, 2800)
my_breaks <- seq(1000, 3000, by = 500)

##applydegree north with function with x axis label
nwbrks <- seq(31,36,1)
nwlbls <- unlist(lapply(nwbrks, function(x) paste(x, "°N")))

library(ggplot2)
library(sf)
library(ggspatial)  # For annotation_scale and annotation_north_arrow
library(cowplot)
library(viridis)

# Define my_breaks
my_breaks <- seq(1000, 3000, by = 500)

# Define nameColor
nameColor <- "Mean Cumulative GPP (g C m⁻² year⁻¹)"

# Plot
cumulativemap<-ggplot() +
  geom_sf(fill = 'transparent', data = ME2) +  # Add the shapefile
  geom_raster(aes(x = x, y = y, fill = mean_8), data = cleaned_cumulativerasdf) +  # Add the raster data
  scale_fill_viridis(option = "turbo", name = nameColor, direction = -1, 
                     breaks = my_breaks, limits = c(1000, 3000)) +  # Customize the fill scale
  labs(x = 'Longitude', y = 'Latitude') +  # Add axis labels
  scale_y_continuous(breaks = seq(34, 36, by = 1)) +  # Customize y-axis breaks
  cowplot::theme_cowplot(font_size = 24) +  # Apply a custom theme
  theme(legend.key.width = unit(4, "cm"), legend.spacing.x = unit(1, 'cm')) +  # Customize legend
  theme(axis.text = element_text(size = 25)) +  # Customize axis text size
  
  # Add scale bar
  annotation_scale(location = "bl", width_hint = 1, height = unit(0.5, "cm")) +  # Add a scale bar at the bottom-left
  
  # Add north arrow
  annotation_north_arrow(location = "tr", which_north = "true", 
                         style = north_arrow_fancy_orienteering,
                         height = unit(2, "cm")) +  # Add a north arrow at the top-right
  
  # Ensure the coordinate system is set to WGS84 (EPSG:4326)
  coord_sf(crs = 4326)

ggplot()+
  geom_sf(fill='transparent',data=ME2)+
  geom_raster(aes(x=x,y=y,fill=mean_8),data=cleaned_cumulativerasdf)+
  #scale_fill_viridis_c( limits = c(0, 300), option = "turbo", breaks = my_breaks, nameColor, direction = -1, oob = scales::squish)+
  scale_fill_viridis( option = "turbo", nameColor, direction = -1, 
                      breaks = my_breaks, 
                      limits = c(1000, 3000)) +
  labs(x='Longitude',y='Latitude', color = nameColor)+
  scale_y_continuous(breaks = seq(34, 36, by=1))+
  cowplot::theme_cowplot(font_size = 24)+
  theme(legend.key.width=unit(4,"cm"), legend.spacing.x = unit(1, 'cm'))+
  theme(axis.text = element_text(size = 25))  +
  
  ggsn::north(arkansasshp) +
  ggsn::scalebar(arkansasshp, dist = 50, dist_unit = "km",st.size=5, height=0.02,
                 transform = TRUE, model = "WGS84")

# Define the breaks for the y-axis to show 34.2, 34.4, 34.6, etc.
latitude_breaks <- seq(33, 36, by = 0.2)
latitude_breaksnwlbls <- unlist(lapply(latitude_breaks, function(x) paste(x, "°N")))
cumulativegam<-  ggplot(cleaned_cumulativerasdf, aes(x=y, y=mean_8)) +
  geom_smooth(level = 0.687)+labs(x='Latitude',y=expression(Mean~Cumulative ~GPP~(g~C~m^{-2}~year^{-1})))+
  scale_x_continuous(breaks = latitude_breaks, labels = latitude_breaksnwlbls, expand = c(0, 0)) +
  #scale_x_continuous(breaks = seq(33, 36, by=1))+
  
  coord_flip()+
  cowplot::theme_cowplot(font_size = 23)+
  theme(axis.text = element_text(size = 25))  

cumulativearranged<-ggarrange(cumulativemap, cumulativegam, labels = c("A", "B"),font.label = list(size = 35) , widths = c(1.6,1),
                              common.legend = TRUE, legend = "bottom")
cumulativearranged

cumulativemap
cumulativegam

ggsave("C:/Users/rbmahbub/Box/Research/ManuscriptFile/Optimum Air Temperature/Figure/VPMcumulativearrangedsameextent.png", plot=cumulativearranged, height=10, width=22, units="in", dpi=150)
ggsave("C:/Users/rbmahbub/Documents/RProjects/VPM_Spatial/Figure/VPMcumulativearranged.png", plot=cumulativearranged, height=10, width=22, units="in", dpi=150)
