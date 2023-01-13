library(sf)
source('R/utils.R')

# Load background spatial data
plant_points <- read_sf("data/gis/Plant_Data.shp")
planting_polygons <- read_sf("data/gis/Planting_Data.shp")
building <- read_sf("data/gis/Building_Mass.shp")
site_boundary <- read_sf("data/gis/Site_Boundary.shp")

# Convert planting polygons to points and combine with existing plant points
plant_points <- convert_combine(point_locations = plant_points,
                                polygon_locations = planting_polygons)

# Process vegetation structure data for specified years 
years_1_30 <- process_structure(point_locations = plant_points,
                                obstructions = building,
                                boundary = site_boundary,
                                years = 1:30,
                                cellarea = 10)

save(years_1_30, file = "output/years_1_30")


