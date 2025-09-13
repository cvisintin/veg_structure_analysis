library(sf)
source('R/utils.R')

# Load background spatial data
home_ex_plant_points <- read_sf("data/gis/sample_landscape/Plant_Data_Home.shp")
home_ex_planting_polygons <- read_sf("data/gis/sample_landscape/Planting_Data_Home.shp")
#home_ex_buildings <- read_sf("data/gis/sample_landscape/Building_Mass_Home.shp")
home_ex_site_boundary <- read_sf("data/gis/sample_landscape/Site_Boundary_Home.shp")

# Convert planting polygons to points and combine with existing plant points
home_ex_plant_points <- convert_combine(point_locations = home_ex_plant_points,
                                        polygon_locations = home_ex_planting_polygons)

# Load spatial grids file
load(file = "output/sample_landscape/home_ex_years_1_30")

# Create export for visualisation
export_image_set(home_ex_years_1_30, path = "output/sample_landscape/", filename = "home_ex_veg_structure")

# Create spatial outputs
create_images(home_ex_years_1_30, path = "figs/images/sample_landscape/", filename = "home_ex_veg_structure")

# Animate spatial outputs
system('convert -delay 50 -dispose previous -loop 0 figs/images/sample_landscape/home_ex_veg* figs/animations/home_ex_veg_growth.gif')

# Load results file
load(file = "output/sample_landscape/home_ex_results")

# Create score sheet
create_score_sheet(analysis_results = home_ex_results,
                               path_directory = "figs/web/sample_landscape/",
                               title = "SAMPLE LANDSCAPE")

create_score_sheet(analysis_results = home_ex_results,
                   path_directory = "figs/pdf/sample_landscape/",
                   title = "SAMPLE LANDSCAPE",
                   web_based = FALSE)

# Plot changes in coverage and connectivity over time
png("figs/images/sample_landscape/cov_conn_change.png", height = 400, width = 600, pointsize = 4, res = 150)
plot_temporal_change(home_ex_results)
dev.off()

####################################################################

# Load background spatial data
park_tq_plant_points <- read_sf("data/gis/torquay/Plant_Data_TQ_Park.shp")
park_tq_site_boundary <- read_sf("data/gis/torquay/Site_Boundary_TQ_Park.shp")

# Load grids file
load(file = "output/torquay/park_tq_years_1_30")

# Create export for visualisation
export_image_set(park_tq_years_1_30, path = "output/torquay/", filename = "park_tq_veg_structure")

# Create spatial outputs
create_images(park_tq_years_1_30, path = "figs/images/torquay/", filename = "park_tq_veg_structure")

# Animate spatial outputs
system('convert -delay 50 -dispose previous -loop 0 figs/images/torquay/park_tq_veg* figs/animations/park_tq_veg_growth.gif')

# Load results file
load(file = "output/torquay/park_tq_results")

# Create score sheet
create_score_sheet(analysis_results = park_tq_results,
                               path_directory = "figs/web/torquay/",
                               title = "GREENSPACE")

# Plot changes in coverage and connectivity over time
png("figs/images/torquay/cov_conn_change.png", height = 400, width = 600, pointsize = 4, res = 150)
plot_temporal_change(park_tq_results, y_max = 30)
dev.off()

####################################################################

# Load background spatial data
park_av_plant_points <- read_sf("data/gis/averley/Plant_Data_AV_Park.shp")
park_av_planting_polygons <- read_sf("data/gis/averley/Planting_Data_AV_Park.shp")
park_av_site_boundary <- read_sf("data/gis/averley/Site_Boundary_AV_Park.shp")

# Convert planting polygons to points and combine with existing plant points
park_av_plant_points <- convert_combine(point_locations = park_av_plant_points,
                                        polygon_locations = park_av_planting_polygons)

# Load grids file
load(file = "output/averley/park_av_years_1_30")

# Create export for visualisation
export_image_set(park_av_years_1_30, path = "output/averley/", filename = "park_av_veg_structure")

# Create spatial outputs
create_images(park_av_years_1_30, path = "figs/images/averley/", filename = "park_av_veg_structure")

# Animate spatial outputs
system('convert -delay 50 -dispose previous -loop 0 figs/images/averley/park_av_veg* figs/animations/park_av_veg_growth.gif')

# Load results file
load(file = "output/averley/park_av_results")

# Create score sheet
create_score_sheet(analysis_results = park_av_results,
                               path_directory = "figs/web/averley/",
                               title = "NEIGHBOURHOOD PARK")

# Plot changes in coverage and connectivity over time
png("figs/images/averley/cov_conn_change.png", height = 400, width = 600, pointsize = 4, res = 150)
plot_temporal_change(park_av_results, y_max = 30)
dev.off()

####################################################################

# Load background spatial data
park_b_plant_points <- read_sf("data/gis/booyeembara/Plant_Data_B_Park.shp")
park_b_planting_polygons <- read_sf("data/gis/booyeembara/Planting_Data_B_Park.shp")
park_b_site_boundary <- read_sf("data/gis/booyeembara/Site_Boundary_B_Park.shp")

# Convert planting polygons to points and combine with existing plant points
park_b_plant_points <- convert_combine(point_locations = park_b_plant_points,
                                       polygon_locations = park_b_planting_polygons)

# Load grids file
load(file = "output/booyeembara/park_b_years_1_30")

# Create export for visualisation
export_image_set(park_b_years_1_30, path = "output/booyeembara/", filename = "park_b_veg_structure")

# Create spatial outputs
create_images(park_b_years_1_30, path = "figs/images/booyeembara/", filename = "park_b_veg_structure")

# Animate spatial outputs
system('convert -delay 50 -dispose previous -loop 0 figs/images/booyeembara/park_b_veg* figs/animations/park_b_veg_growth.gif')

# Load results file
load(file = "output/booyeembara/park_b_results")

# Create score sheet
create_score_sheet(analysis_results = park_b_results,
                               path_directory = "figs/web/booyeembara/",
                               title = "BUSHLAND")

# Plot changes in coverage and connectivity over time
png("figs/images/booyeembara/cov_conn_change.png", height = 400, width = 600, pointsize = 4, res = 150)
plot_temporal_change(park_b_results, y_max = 30)
dev.off()

#####################################################################

metrics <- c("connectivity", "coverage", "density", "endemism", "phenology", "richness", "size", "texture")
for (i in metrics) {
  combine_score_sheet_images(metric = i,
                             nrows = 1,
                             ncols = 3,
                             paths = c("figs/web/torquay/",
                                       "figs/web/averley/",
                                       "figs/web/booyeembara/"),
                             titles = c("Greenspace", "Neighborhood Park", "Bushland"))
}

# Plot changes in connectivity over time for all sites
png("figs/images/conn_change.png", height = 600, width = 800, pointsize = 4, res = 150)
plot_temporal_change(analysis_results_list = list(park_tq_results, park_av_results, park_b_results),
                     variable = "connectivity",
                     titles = c("Greenspace", "Neighborhood Park", "Bushland"))
dev.off()

