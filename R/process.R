library(sf)
source('R/utils.R')

# Load background spatial data for residential project example
home_ex_plant_points <- read_sf("data/gis/Plant_Data_Home.shp")
home_ex_planting_polygons <- read_sf("data/gis/Planting_Data_Home.shp")
home_ex_building <- read_sf("data/gis/Building_Mass_Home.shp")
home_ex_site_boundary <- read_sf("data/gis/Site_Boundary_Home.shp")

# Convert planting polygons to points and combine with existing plant points
home_ex_plant_points <- convert_combine(point_locations = home_ex_plant_points,
                                polygon_locations = home_ex_planting_polygons)

# Process vegetation structure data for specified years 
home_ex_years_1_30 <- process_structure(point_locations = home_ex_plant_points,
                                obstructions = home_ex_building,
                                boundary = home_ex_site_boundary,
                                years = 1:30,
                                cellarea = 10)

save(home_ex_years_1_30, file = "output/home_ex_years_1_30")


# Load background spatial data for public park example
park_ex_plant_points <- read_sf("data/gis/Plant_Data_Park.shp")
park_ex_site_boundary <- read_sf("data/gis/Site_Boundary_Park.shp")

# Process vegetation structure data for specified years 
park_ex_years_1_30 <- process_structure(point_locations = park_ex_plant_points,
                                        obstructions = NULL,
                                        boundary = park_ex_site_boundary,
                                        years = 1:30,
                                        cellarea = 10)

save(park_ex_years_1_30, file = "output/park_ex_years_1_30")
