library(sf)
library(dplyr)

process_structure <- function(cut_heights = c(0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0), # Input height datums in meters
                              year = 1, # Year of interest
                              locations, # Spatial geometry and attributes of plants
                              obstructions, # Geometry that may obstruct plant growth in horizontal plane
                              boundary, # Geometry that defines the boundary of the site
                              hex_gridshape = TRUE, # Use hexagonal-shaped grid cells in the analysis; default is yes
                              approx_cellarea # Target cell area in square meters; note, square-shaped cell areas will be slightly larger and hexagonal-shaped cell areas will be slightly smaller
) {
  
  # Create grid from boundary based on specified grid shape and approximate cell size
  if (hex_gridshape) {
    grid <- st_make_grid(boundary, n = as.numeric(floor(sqrt(ceiling(st_area(boundary)) / approx_cellarea))), square = FALSE)
    grid <- st_sf(data.frame("id" = 1:length(grid)), geometry = grid)
  } else {
    grid <- st_make_grid(boundary, n = as.numeric(floor(sqrt(ceiling(st_area(boundary)) / approx_cellarea))))
    grid <- st_sf(data.frame("id" = 1:length(grid)), geometry = grid)
  }
  
  # Calculate all relative plant heights and widths for the given year (from reference height)
  plant_heights <- locations$max_height * (pmin(year, locations$year_max) / locations$year_max)
  plant_widths <- locations$max_width * (pmin(year, locations$year_max) / locations$year_max)
  
  # Generate id lists for plants that intersect each height datum for the given year
  plant_idx <- lapply(cut_heights, function(cut_height) which(plant_heights + locations$ref_height >= cut_height & cut_height >= locations$ref_height))
  
  proportion_grids <- lapply(cut_heights, function(cut_height) {
    
    # Add new attribute to grid
    grid$prop_ol <- NA
    
    # Extract ids for vegetation that are within the height category
    ids <- plant_idx[[which(cut_heights == cut_height)]]
    
    # Verify if any vegetation exists within height category 
    if(length(ids) > 0) {
      
      # Generate multipliers for width based on assumed geometry (note, this is moderated by height and type)
      mults <- sapply(ids, function(id) {
        ifelse(locations$type[id] == 'tree' & cut_height < plant_heights[id]  / 2, 0,
               ifelse(locations$type[id] == 'tree' & cut_height >= plant_heights[id] / 2, sqrt(pmax(0, (plant_heights[id] * 0.25)^2 - ((cut_height - locations$ref_height[id]) - (plant_heights[id] * 0.75))^2)) / (plant_heights[id] * 0.25),
                      ifelse(locations$type[id] == 'shrub', sqrt(pmax(0, (plant_heights[id] * 0.5)^2 - ((cut_height - locations$ref_height[id]) - (plant_heights[id] * 0.5))^2)) / (plant_heights[id] * 0.5),
                             (cut_height - locations$ref_height[id]) / plant_heights[id])))
      })
      
      # Create a buffered region around all plants based on geometry type for the given year
      buffered_area <- st_union(st_buffer(locations[ids, ], (plant_widths[ids] / 2) * mults))
      
      # Determine which obstructions are above cut height
      obstruction_idx <- which(obstructions$height > cut_height)
      
      # Remove areas with obstructions (e.g. buildings)
      planted_area <- st_difference(buffered_area, obstructions$geometry[obstruction_idx])
      
      # Write out vegetation footprints for given year and selected heights
      # using centimeters to designate height (padded with zeros)
      st_write(planted_area, paste0("output/sp_vegetation/year_",
                                    sprintf("%02d", year),
                                    "_height_",
                                    sprintf("%04d", cut_height * 100),
                                    "_veg.shp"), append = FALSE)
      
      # Identify the grid cells that have overlap with vegetation
      overlap_idx <- st_intersects(grid$geometry, st_union(planted_area), sparse = FALSE)
      
      # Determine the area of overlap in each qualified grid cell
      areas <- st_area(st_intersection(grid$geometry[overlap_idx], st_union(planted_area)))
      
      # Record proportions of overlap to eligible grids
      grid$prop_ol[overlap_idx] <- areas / st_area(grid$geometry[1])
      
    }
    
    # Write out grids for given year and selected heights using centimeters
    # to designate height (padded with zeros)
    st_write(grid, paste0("output/sp_grids/year_",
                          sprintf("%02d", year),
                          "_height_",
                          sprintf("%04d", cut_height * 100),
                          "_grd_",
                          ifelse(hex_gridshape, "hex", "rect") ,
                          ".shp"), append = FALSE)
    
    grid
  })
  
  proportion_grids
}

prepare_conn_data <- function(grid, # Grid of polygons with proportions vegetation coverage to be used in the analysis
                              obstructions, # Geometry that may obstruct connectivity
                              cut_height, # Specified cut line height
                              grid_shape # Specify the shape of the grid cells
) {
  # Replace all missing values with zeros
  grid$prop_ol[is.na(grid$prop_ol)] <- 0
  
  # Determine which obstructions are at or above cut height
  obstruction_idx <- which(obstructions$height >= cut_height)
  
  # Verify if any obstructions exist 
  if(length(obstruction_idx) > 0) {
    
    # Identify the grid cells that have overlap with obstructions
    overlap_idx <- st_intersects(grid$geometry, st_union(obstructions$geometry[obstruction_idx]), sparse = FALSE)[ , 1]
    overlap_idx <- which(overlap_idx == TRUE)
    
    # Determine the proportion area of overlap in each qualified grid cell
    props <- as.numeric(st_area(st_intersection(grid$geometry[overlap_idx], st_union(obstructions)))) / as.numeric(st_area(grid$geometry[1]))
    
    # Identify the grid cells that have greater than 90% overlap with obstructions
    full_overlap_idx <- overlap_idx[which(props >= 0.90)]
    
    # Set all grid cells that have overlap with obstructions greater than 90% to NA
    grid$prop_ol[full_overlap_idx] <- NA
    
    # Identify grid cells that have some overlap with obstructions
    prop_overlap_idx <- setdiff(overlap_idx, full_overlap_idx)
    
    # Reduce proportion vegetation coverage by proportion overlap with obstructions
    grid$prop_ol[prop_overlap_idx] <- grid$prop_ol[prop_overlap_idx] * props[which(props < 0.90)]
  }
  
  # Invert all non-NA values
  grid$resist <- 1 - grid$prop_ol
  
  # Add attribute to identify grid cell shape type
  attr(grid, 'shape') <- grid_shape
    
  # Return grid with relative resistance values
  grid
}

# run_connectivity <- function(resistance_grid) {
# 
#   grid_shape <- attr(resistance_grid, 'shape')
#   
#   # Determine how many missing cells are in grid
#   n_missing <- sum(is.na(resistance_grid$resist))
#   
#   grid_noNA <- resistance_grid[!is.na(resistance_grid$resist), ]
#   
#   points <- st_centroid(grid_noNA)
#   
#   ext <- st_bbox(points)
#   n_points <- nrow(points)
#   idx <- points$id
#   
#   search_dist <- st_distance(points$geometry[1], points$geometry[2]) * 1.5
#   if (grid_shape == "rect") search_dist <- st_distance(points$geometry[1], points$geometry[2]) * 1.5
#   
#   nodes <- sapply(idx, function(p) {
#     ids <- st_intersection(points, st_buffer(points$geometry[points$id == p], search_dist))$id
#     ids <- setdiff(ids, p)
#     base_con <- points$resist[points$id == p]
#     means <- sapply(ids, function(i) mean(c(base_con, points$resist[which(points$id == i)])))
#     cbind(p, ids, means)
#   })
#   
#   node_list <- as.matrix(do.call(rbind, nodes))
#   
#   node_list <- node_list[!duplicated(apply(node_list[ , ], 1, function(row) paste(sort(row), collapse = ""))),]
#   
#   write.table(node_list, file = "data/circuitscape/network.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
#   
#   # write.table(points$id, file = "data/circuitscape/focal_nodes.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
#   
#   
#   ul_id <- points$id[st_nearest_feature(st_point(c(ext[1], ext[4])), points, check_crs = FALSE)]
#   ll_id <- points$id[st_nearest_feature(st_point(c(ext[1], ext[2])), points, check_crs = FALSE)]
#   ur_id <- points$id[st_nearest_feature(st_point(c(ext[3], ext[4])), points, check_crs = FALSE)]
#   lr_id <- points$id[st_nearest_feature(st_point(c(ext[3], ext[2])), points, check_crs = FALSE)]
#   write.table(c(ul_id, ll_id, ur_id, lr_id), file = "data/circuitscape/focal_nodes.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
#   
#   resistance_grid$current <- NA
#   sum_currents <- read.delim("data/circuitscape/connectivity_node_currents_cum.txt", header = FALSE)[idx , 2]
#   resistance_grid$current[idx] <- zero_to_one(sum_currents)
# }

zero_to_one <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

# create_grid <- function(boundary, # Boundary to create grid cells over
#                         gridshape = "hexagonal", # Shape of grid
#                         cellarea = 10 # Size of grid cells in square meters
# ) {
#   
#   # Extend the boundary based on 10% of the shortest dimension
#   b_box <- st_bbox(boundary)
#   x_length <- diff(b_box[c(1, 3)])
#   y_length <- diff(b_box[c(2, 4)])
#   offset <- min(x_length, y_length) * 0.1
#   new_extent <- c(b_box[1] - offset, b_box[2] - offset, b_box[3] + offset, b_box[4] + offset)
# 
#   if(gridshape == "hexagonal") {
#     
#     # Properties of hexagon based on specified area
#     side <- sqrt(cellarea / ((3 * sqrt(3)) / 2))
#     perimeter <- side * 6
#     apothem <- 2 * cellarea / perimeter
#     radius <- side
#     x_spacing <- 2 * radius
#     y_spacing <- 2 * apothem
#     
#     # Determine number of cells
#     n_x <- ceiling(x_length / x_diameter)
#     n_y <- ceiling(y_length / y_diameter)
#     
#     
#     
#   } else {
#     
#   }
#  
# }


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