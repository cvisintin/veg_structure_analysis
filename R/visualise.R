library(sf)
library(ggplot2)
library(ggimage)
library(grid)
library(gridExtra)
library(png)
library(dplyr)
library(tidyr)
library(viridis)
source('R/utils.R')

# Load background spatial data
plant_points <- read_sf("data/gis/Plant_Data.shp")
planting_polygons <- read_sf("data/gis/Planting_Data.shp")
building <- read_sf("data/gis/Building_Mass.shp")
site_boundary <- read_sf("data/gis/Site_Boundary.shp")

# Convert planting polygons to points and combine with existing plant points
plant_points <- convert_combine(point_locations = plant_points,
                                polygon_locations = planting_polygons)

# Load grids file
load(file = "output/years_1_30")

# Create spatial outputs
create_images(years_1_30, path = "figs/")

# Animate spatial outputs
system('convert -delay 50 -dispose previous -loop 0 figs/* figs/animations/veg_growth.gif')

### Density ###
density_plot <- plot_classes(plant_points,
                             variable_name = "density",
                             colour_palette = c("#397d53", "#b7e4c8", "#62a67c"),
                             image_path = "data/images/leaves_icon.png")

### Texture ###
texture_plot <- plot_classes(plant_points,
                             variable_name = "texture",
                             colour_palette = c("#4b6c90", "#afc6e0", "#6d8eb3"),
                             image_path = "data/images/texture_icon.png")

### Size ###
size_plot <- plot_classes(plant_points,
                          variable_name = "size",
                          colour_palette = c("#b8641d", "#ebc5a4", "#d9a171", "#c98449"),
                          image_path = "data/images/vegetation_icon.png")

### Endemism ###
endemism_plot <- plot_percent(plant_points,
                              variable_name = "endemism",
                              colour = "#a072a6",
                              label = "NATIVE")

### Type ###
type_plot <- plot_percent(plant_points,
                          variable_name = "type",
                          colour = "#a19a65",
                          label = "EVERGREEN")

### Species richness ###
richness_plot <- plot_circ_bar(spatial_points = plant_points,
                              variable_name = "richness",
                              colours = "#858383",
                              polar_rotation = 0.25)

### Phenology ###
phenology_plot <- plot_circ_bar(spatial_points = plant_points,
                               variable_name = "phenology",
                               colours = "#b0d4d6",
                               polar_rotation = 0.25)

### Coverage at 1m ###
coverage_plot <- plot_circ_bar(spatial_points = plant_points,
                                spatial_list = years_1_30,
                                variable_name = "coverage",
                                colours = "#c9837d")

### Connectivity ###
connectivity_plot <- plot_circ_bar(spatial_points = plant_points,
                               spatial_list = years_1_30,
                               variable_name = "connectivity",
                               stacked = TRUE)

png(paste0("figs/results.png"), height = 2400, width = 2400, pointsize = 4, res = 300)
grid.arrange(grobs = list(density_plot,
                          texture_plot,
                          size_plot,
                          endemism_plot,
                          richness_plot,
                          type_plot,
                          phenology_plot,
                          connectivity_plot,
                          coverage_plot
), ncol = 3)
dev.off()
