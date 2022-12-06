library(sf)
source('R/utility_functions.R')

# Load background spatial data
plant_points <- read_sf("data/gis/Plant_Data.shp")
building <- read_sf("data/gis/Building_Mass.shp")
hex_grid <- read_sf("data/gis/Hex_Grid.shp")


