library(sf)
source('R/utility_functions.R')

# Load background spatial data
plant_points <- read_sf("data/gis/Plant_Data.shp")
building <- read_sf("data/gis/Building_Mass.shp")
site_boundary <- read_sf("data/gis/Site_Boundary.shp")
hex_grid <- read_sf("data/gis/Hex_Grid.shp")

# Process vegetation structure data for specified years 
for(y in c(1, 5, 10, 15, 20)) {
  temp <- process_structure(cut_heights = c(0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0),
                            year = y,
                            locations = plant_points,
                            obstructions = building,
                            boundary = site_boundary,
                            hex_gridshape = TRUE,
                            cellarea = 10)
  assign(paste0("year_", sprintf("%02d", y)), temp)
}

