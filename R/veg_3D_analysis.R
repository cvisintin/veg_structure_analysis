library(sf)

plant_points <- read_sf("Plant_Data.shp")
building <- read_sf("Building_Mass.shp")
hex_grid <- read_sf("Hex_Grid.shp")

process_structure <- function(heights = c(0.1, 0.5, 1.0, 5.0, 10.0), # Input height datums in meters
                              year = 1, # Year of interest
                              locations = plant_points, # Spatial geometry and attributes of plants
                              obstructions = building, # Geometry that may obstruct plants/connectivity
                              grid = hex_grid # Grid of polygons to be used in the analysis
                              ) {

  # Calculate all plant heights and widths for the given year
  plant_heights <- locations$max_height / locations$year_max * pmin(year, locations$year_max)
  plant_widths <- locations$max_width / locations$year_max * pmin(year, locations$year_max)
  
  # Generate id lists for plants that intersect each height datum for the given year
  plant_idx <- lapply(heights, function(height) which(plant_heights >= height))
  
  map_data <- sapply(heights, function(height) {
    # Extract ids for plants that are within the height category
    ids <- unlist(plant_idx[which(heights == height)])
    
    if(length(ids) == 0) return(NULL)

    # Generate multipliers for width based on assumed geometry (note, this is moderated by height and type)
    mults <- sapply(ids, function(id) {
      ifelse(locations$type[id] == 'tree' & height < plant_heights[id] / 2, 0,
             ifelse(locations$type[id] == 'tree' & height >= plant_heights[id] / 2, sqrt((plant_heights[id] * 0.25)^2 - (height - (plant_heights[id] * 0.75))^2) / (plant_heights[id] * 0.25),
                    ifelse(locations$type[id] == 'shrub', sqrt((plant_heights[id] * 0.5)^2 - (height - (plant_heights[id] * 0.5))^2) / (plant_heights[id] * 0.5),
                           height / plant_heights[id])))
     })
    
    # Create a buffered region around all plants based on geometry type for the given year
    buffered_area <- st_union(st_buffer(locations[ids, ], (plant_widths[ids] / 2) * mults))
    
    st_write(buffered_area, paste0("year_", year, "_height_", height, ".shp"))
    
    # overlaps <- st_intersection(buffers, hex_grid)

  })
  
  
}

cbind(locations$type[ids], plant_heights[ids], mults)
