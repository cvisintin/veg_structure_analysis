library(sf)
library(dplyr)

process_structure <- function(heights = c(0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0), # Input height datums in meters
                              year = 1, # Year of interest
                              locations, # Spatial geometry and attributes of plants
                              obstructions, # Geometry that may obstruct plants/connectivity
                              grid # Grid of polygons to be used in the analysis
) {
  
  # Calculate all relative plant heights and widths for the given year (from reference height)
  plant_heights <- locations$max_height / locations$year_max * pmin(year, locations$year_max)
  plant_widths <- locations$max_width / locations$year_max * pmin(year, locations$year_max)
  
  # Generate id lists for plants that intersect each height datum for the given year
  plant_idx <- lapply(heights, function(height) which(plant_heights + locations$ref_height >= height | height >= locations$ref_height))
  
  proportion_grids <- lapply(heights, function(height) {
    
    # Create new empty grid
    out_grid <- grid[, -c(2:5)]
    out_grid$prop_ol <- NA
    
    # Extract ids for vegetation that are within the height category
    ids <- plant_idx[[which(heights == height)]]
    
    # Verify if any vegetation exists within height category 
    if(length(ids) > 0) {
      
      # Generate multipliers for width based on assumed geometry (note, this is moderated by height and type)
      mults <- sapply(ids, function(id) {
        ifelse(locations$type[id] == 'tree' & height < plant_heights[id]  / 2, 0,
               ifelse(locations$type[id] == 'tree' & height >= plant_heights[id] / 2, sqrt(pmax(0, (plant_heights[id] * 0.25)^2 - ((height - locations$ref_height[id]) - (plant_heights[id] * 0.75))^2)) / (plant_heights[id] * 0.25),
                      ifelse(locations$type[id] == 'shrub', sqrt(pmax(0, (plant_heights[id] * 0.5)^2 - ((height - locations$ref_height[id]) - (plant_heights[id] * 0.5))^2)) / (plant_heights[id] * 0.5),
                             (height - locations$ref_height[id]) / plant_heights[id])))
      })
      
      # Create a buffered region around all plants based on geometry type for the given year
      buffered_area <- st_union(st_buffer(locations[ids, ], (plant_widths[ids] / 2) * mults))
      
      # Remove areas with constructed buildings
      planted_area <- st_difference(buffered_area, obstructions)
      
      # Write out vegetation footprints for given year and selected heights
      # using centimeters to designate height (padded with zeros)
      st_write(planted_area, paste0("output/sp_vegetation/year_", sprintf("%02d", year), "_height_", sprintf("%04d", height * 100), "_veg.shp"), append = FALSE)
      
      # Identify the grid cells that have overlap with vegetation
      overlap_idx <- st_intersects(grid$geometry, st_union(planted_area), sparse = FALSE)
      
      # Determine the area of overlap in each qualified grid cell
      areas <- st_area(st_intersection(grid$geometry[overlap_idx], st_union(planted_area)))
      
      # Record proportions of overlap to eligible grids
      out_grid$prop_ol[overlap_idx] <- areas / st_area(grid$geometry[1])
      
    }
    
    # Write out grids for given year and selected heights using centimeters
    # to designate height (padded with zeros)
    st_write(out_grid, paste0("output/sp_grids/year_", sprintf("%02d", year), "_height_", sprintf("%04d", height * 100), "_grd.shp"), append = FALSE)
    
    out_grid
  })
  
  proportion_grids
}

#### Original function by Stefan Jünger
#### https://stefanjuenger.github.io/gesis-workshop-geospatial-techniques-R/slides/2_4_Advanced_Maps_II/2_4_Advanced_Maps_II.html#8
rotate_data <- function(data, x_add = 0, y_add = 0) {
  shear_matrix <- function () { 
    matrix(c(2, 1.2, 0, 1), 2, 2) 
  }
  rotate_matrix <- function(x) { 
    matrix(c(cos(x), sin(x), -sin(x), cos(x)), 2, 2) 
  }
  data %>% 
    dplyr::mutate(
      geometry = 
        .$geometry * shear_matrix() * rotate_matrix(pi/20) + c(x_add, y_add)
    )
}