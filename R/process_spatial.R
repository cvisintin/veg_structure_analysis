library(sf)
source('R/utility_functions.R')

# Load background spatial data
plant_points <- read_sf("data/gis/Plant_Data.shp")
building <- read_sf("data/gis/Building_Mass.shp")
site_boundary <- read_sf("data/gis/Site_Boundary.shp")
# hex_grid <- read_sf("data/gis/Hex_Grid.shp")

# Process vegetation structure data for specified years 
all_data <- process_structure(locations = plant_points,
                              obstructions = building,
                              boundary = site_boundary,
                              years = c(1, 5, 10, 15, 20, 25, 30),
                              cellarea = 10)
