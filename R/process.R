library(sf)
source('R/utils.R')

##################### Example project ##############################

# Load background spatial data for residential project example
home_ex_plant_points <- read_sf("data/gis/sample_landscape/Plant_Data_Home.shp")
home_ex_planting_polygons <- read_sf("data/gis/sample_landscape/Planting_Data_Home.shp")
home_ex_building <- read_sf("data/gis/sample_landscape/Building_Mass_Home.shp")
home_ex_site_boundary <- read_sf("data/gis/sample_landscape/Site_Boundary_Home.shp")

# Convert planting polygons to points and combine with existing plant points
home_ex_plant_points <- convert_combine(point_locations = home_ex_plant_points,
                                polygon_locations = home_ex_planting_polygons)

# Process vegetation structure data for specified years 
home_ex_years_1_30 <- process_structure(point_locations = home_ex_plant_points,
                                obstructions = home_ex_building,
                                boundary = home_ex_site_boundary,
                                years = 1:30,
                                cellarea = 10)

# Save spatial grids file
save(home_ex_years_1_30, file = "output/sample_landscape/home_ex_years_1_30")

# Analyse spatial grids
home_ex_results <- analyse_spatial_data(spatial_points = home_ex_plant_points,
                               spatial_list = home_ex_years_1_30)

# Save analysis data
save(home_ex_results, file = "output/sample_landscape/home_ex_results")

####################################################################

####################### Torquay Greenspace #########################

# Load background spatial data for public greenspace example
park_tq_plant_points <- read_sf("data/gis/torquay/Plant_Data_TQ_Park.shp")
park_tq_site_boundary <- read_sf("data/gis/torquay/Site_Boundary_TQ_Park.shp")

# Process vegetation structure data for specified years 
park_tq_years_1_30 <- process_structure(point_locations = park_tq_plant_points,
                                        obstructions = NULL,
                                        boundary = park_tq_site_boundary,
                                        years = 1:30,
                                        cellarea = 10)

save(park_tq_years_1_30, file = "output/torquay/park_tq_years_1_30")

####################################################################

####################### Averley Park ###############################

# Load background spatial data for residential park example
park_av_plant_points <- read_sf("data/gis/averley/Plant_Data_AV_Park.shp")
park_av_planting_polygons <- read_sf("data/gis/averley/Planting_Data_AV_Park.shp")
park_av_site_boundary <- read_sf("data/gis/averley/Site_Boundary_AV_Park.shp")

# Convert planting polygons to points and combine with existing plant points
park_av_plant_points <- convert_combine(point_locations = park_av_plant_points,
                                        polygon_locations = park_av_planting_polygons)

# Process vegetation structure data for specified years
park_av_years_1_30 <- process_structure(point_locations = park_av_plant_points,
                                        boundary = park_av_site_boundary,
                                        years = 1:30,
                                        cellarea = 10)

save(park_av_years_1_30, file = "output/averley/park_av_years_1_30")

####################################################################

####################### Booyeembara Bushland #######################

# Load background spatial data for urban bushland example
park_b_plant_points <- read_sf("data/gis/booyeembara/Plant_Data_B_Park.shp")
park_b_planting_polygons <- read_sf("data/gis/booyeembara/Planting_Data_B_Park.shp")
park_b_site_boundary <- read_sf("data/gis/booyeembara/Site_Boundary_B_Park.shp")

# Convert planting polygons to points and combine with existing plant points
park_b_plant_points <- convert_combine(point_locations = park_b_plant_points,
                                        polygon_locations = park_b_planting_polygons)

# Process vegetation structure data for specified years
park_b_years_1_30 <- process_structure(point_locations = park_b_plant_points,
                                        boundary = park_b_site_boundary,
                                        years = 1:30,
                                        cellarea = 10)

save(park_b_years_1_30, file = "output/booyeembara/park_b_years_1_30")

####################################################################



####################### TESTING BEHAVIOUR #######################

# Load background spatial data for residential project example
conn_test_plant_points <- read_sf("data/gis/sample_landscape/Conn_Test_Data_Home.shp")
conn_test_site_boundary <- read_sf("data/gis/sample_landscape/Site_Boundary_Home.shp")

# Process vegetation structure data for specified years 
conn_test_years_1_30 <- process_structure(point_locations = conn_test_plant_points,
                                        boundary = conn_test_site_boundary,
                                        years = 1:30,
                                        cellarea = 10)

save(conn_test_years_1_30, file = "output/sample_landscape/conn_test_years_1_30")

# Load grids file
load(file = "output/sample_landscape/conn_test_years_1_30")

# Create score sheet
create_interactive_score_sheet(spatial_points = conn_test_plant_points,
                               boundary = conn_test_site_boundary,
                               spatial_list = conn_test_years_1_30,
                               path_directory = "figs/web/conn_test/",
                               title = "CONNECTIVITY TEST")
