library(sf)
source('R/utility_functions.R')

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
                                polygon_locations = planting_polygons,
                                obstructions = building,
                                boundary = site_boundary,
                                years = 1:30,
                                cellarea = 10,
                                cut_heights = c(0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0))

save(years_1_30, file = "output/years_1_30")

