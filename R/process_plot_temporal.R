library(gridExtra)
library(png)
source("R/utils.R")

plant_points <- read_sf("data/gis/sample_temporal/plant_final_gis.shp")
boundary <- read_sf("data/gis/sample_temporal/boundary.shp")

years_1_30 <- process_structure(point_locations = plant_points,
                                boundary = boundary,
                                years = 1:30,
                                cellarea = 5)

save(years_1_30, file = "output/sample_temporal/years_1_30")

create_images(years_1_30, path = "figs/images/sample_temporal/", filename = "veg_structure")

m_plot <- readPNG('data/gis/sample_temporal/veg_plan_gis.png')
t1_plot <- readPNG('figs/images/sample_temporal/veg_structure_year_01.png')
t2_plot <- readPNG('figs/images/sample_temporal/veg_structure_year_03.png')
t3_plot <- readPNG('figs/images/sample_temporal/veg_structure_year_05.png')
t4_plot <- readPNG('figs/images/sample_temporal/veg_structure_year_10.png')

g_plot <- arrangeGrob(rasterGrob(m_plot),
                       rasterGrob(t1_plot),
                       rasterGrob(t2_plot),
                       rasterGrob(t3_plot),
                       rasterGrob(t4_plot),
                       ncol = 4,
                       nrow = 3,
                       layout_matrix = rbind(c(1,1,1,1), c(2,3,4,5), c(2,3,4,5)))

ggsave('figs/temporal_change.png', g_plot, height = 8, width = 10, pointsize = 8)