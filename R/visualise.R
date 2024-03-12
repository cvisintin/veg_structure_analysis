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

# Load grids file
load(file = "output/sample_landscape/home_ex_years_1_30")

# Create export for visualisation
export_image_set(home_ex_years_1_30, path = "output/sample_landscape/", filename = "home_ex_veg_structure")

# Create spatial outputs
create_images(home_ex_years_1_30, path = "figs/images/sample_landscape/", filename = "home_ex_veg_structure")

# Animate spatial outputs
system('convert -delay 50 -dispose previous -loop 0 figs/images/sample_landscape/home_ex_veg* figs/animations/home_ex_veg_growth.gif')

# Create score sheet
create_interactive_score_sheet(spatial_points = home_ex_plant_points,
                               boundary = home_ex_site_boundary,
                               spatial_list = home_ex_years_1_30,
                               path_directory = "figs/web/sample_landscape/")

####################################################################

# Load background spatial data
park_tq_plant_points <- read_sf("data/gis/torquay/Plant_Data_TQ_Park.shp")

# Load grids file
load(file = "output/torquay/park_tq_years_1_30")

# Create spatial outputs
create_images(park_tq_years_1_30, path = "figs/images/torquay/", filename = "park_tq_veg_structure")

# Animate spatial outputs
system('convert -delay 50 -dispose previous -loop 0 figs/images/torquay/park_tq_veg* figs/animations/park_tq_veg_growth.gif')

# Create score sheet
create_interactive_score_sheet(spatial_points = park_tq_plant_points,
                               spatial_list = park_tq_years_1_30,
                               path_directory = "figs/web/torquay/")

####################################################################

# Load background spatial data
park_av_plant_points <- read_sf("data/gis/averley/Plant_Data_AV_Park.shp")
park_av_planting_polygons <- read_sf("data/gis/averley/Planting_Data_AV_Park.shp")

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

####################################################################

# Load background spatial data
park_b_plant_points <- read_sf("data/gis/booyeembara/Plant_Data_B_Park.shp")
park_b_planting_polygons <- read_sf("data/gis/booyeembara/Planting_Data_B_Park.shp")

# Convert planting polygons to points and combine with existing plant points
park_b_plant_points <- convert_combine(point_locations = park_b_plant_points,
                                       polygon_locations = park_b_planting_polygons)

# Load grids file
load(file = "output/booyeembara/park_b_years_1_30")

# Create spatial outputs
create_images(park_b_years_1_30, path = "figs/", filename = "park_b_veg_structure")

# Animate spatial outputs
system('convert -delay 50 -dispose previous -loop 0 figs/park_b_veg* figs/animations/park_b_veg_growth.gif')

# Create score sheet
create_static_score_sheet(spatial_points = park_b_plant_points,
                          spatial_list = park_b_years_1_30,
                          path_filename = "figs/park_b_results.png")

# Create score sheet
create_interactive_score_sheet(spatial_points = park_b_plant_points,
                               spatial_list = park_b_years_1_30,
                               path_directory = "figs/web/booyeembara/")
