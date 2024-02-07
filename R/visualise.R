library(sf)
source('R/utils.R')

# Load background spatial data
home_ex_plant_points <- read_sf("data/gis/sample_landscape/Plant_Data_Home.shp")
home_ex_planting_polygons <- read_sf("data/gis/sample_landscape/Planting_Data_Home.shp")
#home_ex_buildings <- read_sf("data/gis/sample_landscape/Building_Mass_Home.shp")

# Convert planting polygons to points and combine with existing plant points
home_ex_plant_points <- convert_combine(point_locations = home_ex_plant_points,
                                        polygon_locations = home_ex_planting_polygons)

# Load grids file
load(file = "output/sample_landscape/home_ex_years_1_30")

# Create export for visualisation
export_image_set(home_ex_years_1_30, path = "output/sample_landscape/", filename = "home_ex_veg_structure")

# Create spatial outputs
create_images(home_ex_years_1_30, path = "figs/", filename = "home_ex_veg_structure")

# Animate spatial outputs
system('convert -delay 50 -dispose previous -loop 0 figs/home_ex_veg* figs/animations/home_ex_veg_growth.gif')

# Create score sheet
create_static_score_sheet(spatial_points = home_ex_plant_points,
                   spatial_list = home_ex_years_1_30,
                   path_filename = "figs/home_ex_results.png")

create_interactive_score_sheet(spatial_points = home_ex_plant_points,
                          spatial_list = home_ex_years_1_30,
                          path_directory = "figs/web/sample_landscape/")

####################################################################

# Load background spatial data
park_ex_plant_points <- read_sf("data/gis/Plant_Data_Park.shp")

# Load grids file
load(file = "output/park_ex_years_1_30")

# Create spatial outputs
create_images(park_ex_years_1_30, path = "figs/", filename = "park_ex_veg_structure")

# Animate spatial outputs
system('convert -delay 50 -dispose previous -loop 0 figs/park_ex_veg* figs/animations/park_ex_veg_growth.gif')

# Create score sheet
create_static_score_sheet(spatial_points = park_ex_plant_points,
                   spatial_list = park_ex_years_1_30,
                   path_filename = "figs/park_ex_results.png")

####################################################################

# Load background spatial data
park_av_plant_points <- read_sf("data/gis/Averley/Plant_Data_AV_Park.shp")
park_av_planting_polygons <- read_sf("data/gis/Averley/Planting_Data_AV_Park.shp")

# Convert planting polygons to points and combine with existing plant points
park_av_plant_points <- convert_combine(point_locations = park_av_plant_points,
                                        polygon_locations = park_av_planting_polygons)

# Load grids file
load(file = "output/park_av_years_1_30")

# Create spatial outputs
create_images(park_av_years_1_30, path = "figs/", filename = "park_av_veg_structure")

# Animate spatial outputs
system('convert -delay 50 -dispose previous -loop 0 figs/park_av_veg* figs/animations/park_av_veg_growth.gif')

# Create score sheet
create_static_score_sheet(spatial_points = park_av_plant_points,
                   spatial_list = park_av_years_1_30,
                   path_filename = "figs/park_av_results.png")

# Create score sheet
create_interactive_score_sheet(spatial_points = park_av_plant_points,
                               spatial_list = park_av_years_1_30,
                               path_directory = "figs/web/averley/")
