library(sf)
library(dplyr)

convert_combine <- function(point_locations, # Spatial geometry and attributes of plants
                            polygon_locations # Spatial geometry and attributes of planted areas
) {
  n_polys <- nrow(polygon_locations)
  
  point_list <- lapply(seq_len(n_polys), function(i) {
    data <- st_drop_geometry(polygon_locations[i, ])
    cell_size <- 1 / sqrt(polygon_locations$spacing[i])
    grid <- st_make_grid(st_bbox(polygon_locations$geometry[i]), cellsize = cell_size, what = "centers")
    points <- st_intersection(grid, polygon_locations$geometry[i])
    n_points <- length(points)
    points <- points[sample(seq_len(n_points), floor(n_points * polygon_locations$coverage[i] / 100))]
    st_sf(data, geometry = points)
  })
  
  new_plant_points <- do.call(rbind, point_list)
  
  rbind(point_locations, new_plant_points)
}



process_structure <- function(point_locations, # Spatial geometry and attributes of plants
                              obstructions, # Geometry that may obstruct plant growth in horizontal plane
                              boundary, # Geometry that defines the boundary of the site
                              cut_heights = NULL, # Optional height datums in meters
                              years = 1, # Years of interest
                              hex_gridshape = TRUE, # Use hexagonal-shaped grid cells in the analysis; default is yes
                              cellarea, # Target cell area in square meters
                              veg_layers_path = NULL, # Path to folder where vegetation layers should be written
                              overlap_grids_path = NULL # Path to folder where grid layers should be written
) {
  
  # Create grid from boundary based on specified grid shape and approximate cell size
  master_grid <- create_grid(boundary = boundary, hex_gridshape = hex_gridshape, cellarea = cellarea)
  
  # Determine cell spacing distance
  cell_dist <- sqrt(2 * cellarea / (3 * sqrt(3))) * sqrt(3)
  if(hex_gridshape == FALSE) cell_dist <- sqrt(cellarea)
  
  # Determine maximum height of vegetation plus 5%
  max_veg_height <- max(point_locations$max_height) * 1.05
  
  # Number of vertical layers at cell spacing distance within maximum vegetation height
  n_vert_layers <- max_veg_height %/% cell_dist
  
  # If cut_heights are not specified, use regular spacing based on distance between cells
  if(is.null(cut_heights)) cut_heights <- (1:n_vert_layers) * cell_dist
  
  # If 1m is not specified in cut heights, add it
  if(!1 %in% cut_heights) cut_heights <- sort(c(1, cut_heights))
  
  out <- lapply(years, function(year) {  
    
    print(paste0("Processing year ", year))
    
    # Calculate all relative plant heights and widths for the given year (from reference height)
    plant_heights <- point_locations$max_height * (pmin(year, point_locations$year_max) / point_locations$year_max)
    plant_widths <- point_locations$max_width * (pmin(year, point_locations$year_max) / point_locations$year_max)
    
    # Generate id lists for plants that intersect each height datum for the given year
    plant_idx <- lapply(cut_heights, function(cut_height) which(plant_heights + point_locations$ref_height >= cut_height & cut_height >= point_locations$ref_height))
    
    proportion_grids <- lapply(cut_heights, function(cut_height) {
      
      # Copy master grid
      grid <- master_grid
      
      # Add new attribute to grid
      grid$prop_ol <- as.numeric(NA)
      
      # Extract ids for vegetation that are within the height category
      ids <- plant_idx[[which(cut_heights == cut_height)]]
      
      # Verify if any vegetation exists within height category 
      if(length(ids) > 0) {
        
        # Generate multipliers for width based on one of seven specified geometry forms (note, this is moderated by plant height)
        mults <- sapply(ids, function(id) {
          adj_cut_height <- cut_height - point_locations$ref_height[id]
          plant_height <- plant_heights[id]
          do.call(point_locations$form[id], list(adj_cut_height = adj_cut_height, plant_height = plant_height))
        })
        
        # Create a buffered region around all plants based on geometry forms for the given year
        buffered_area <- st_union(st_buffer(point_locations[ids, ], (plant_widths[ids] / 2) * mults))
        
        # Determine which obstructions are above cut height
        obstruction_idx <- which(obstructions$height > cut_height)
        
        # Remove areas with obstructions (e.g. buildings)
        planted_area <- st_difference(buffered_area, obstructions$geometry[obstruction_idx])
        
        # If path is provided, write out vegetation footprints for given year and selected heights
        # using centimeters to designate height (padded with zeros)
        if(!is.null(veg_layers_path)) {
          st_write(planted_area, paste0(veg_layers_path,
                                        "year_",
                                        sprintf("%02d", year),
                                        "_height_",
                                        sprintf("%04d", ceiling(cut_height * 100)),
                                        "_veg.shp"), append = FALSE)
        }
        
        # Identify the grid cells that have overlap with vegetation
        overlap_idx <- st_intersects(grid$geometry, st_union(planted_area), sparse = FALSE)
        
        # Determine the area of overlap in each qualified grid cell
        areas <- st_area(st_intersection(grid$geometry[overlap_idx], st_union(planted_area)))
        
        # Record proportions of overlap to eligible grids
        grid$prop_ol[overlap_idx] <- areas / st_area(grid$geometry[1])
        
      }
      
      # If path is provided, write out grids for given year and selected heights using centimeters
      # to designate height (padded with zeros)
      if(!is.null(overlap_grids_path)) {
        st_write(grid, paste0(overlap_grids_path,
                              "year_",
                              sprintf("%02d", year),
                              "_height_",
                              sprintf("%04d", ceiling(cut_height * 100)),
                              "_grd_",
                              ifelse(hex_gridshape, "hex", "rect") ,
                              ".shp"), append = FALSE)
      }
      
      grid
      
    })
    
    names(proportion_grids) <- paste0("cut_", sprintf("%04d", ceiling(cut_heights * 100)))
    
    proportion_grids
    
  })
  
  names(out) <- paste0("year_", sprintf("%02d", years))
  
  out
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

create_grid <- function(boundary, # Boundary to create grid cells over
                        hex_gridshape = TRUE, # Shape of grid; default is hexagonal
                        cellarea = 10 # Size of grid cells in square meters
) {
  
  # Determine horizontal and vertical spacing
  horizontal_spacing <- vertical_spacing <- sqrt(cellarea)
  
  # If hexagonal, calculate spatial properties and spacing for regular "flat-topped" cell
  if (hex_gridshape) {
    outer_radius <- sqrt(2 * cellarea / (3 * sqrt(3)))
    inner_radius <- (outer_radius * sqrt(3)) / 2
    horizontal_spacing <- 3 * outer_radius
    vertical_spacing <- inner_radius
  }
  
  # Increase boundary to ensure adequate coverage
  boundary_buff <- st_buffer(boundary, 2 * horizontal_spacing)
  
  # Determine extent and size of boundary
  extent <- st_bbox(boundary_buff)
  horizontal_dist <- diff(extent[c(1, 3)])
  vertical_dist <- diff(extent[c(2, 4)])
  
  # Determine total number of columns and rows
  n_columns <- horizontal_dist %/% horizontal_spacing
  n_rows <- vertical_dist %/% vertical_spacing
  
  # Create coordinates
  x_coords <- rep(c(extent[1], extent[1] + seq_len(n_columns - 1) * horizontal_spacing), n_rows)
  y_coords <- sort(rep(c(extent[4], extent[4] - seq_len(n_rows - 1) * vertical_spacing), n_columns), decreasing = TRUE)
  
  # Combine coordinates into a matrix of x and y values
  xy_coords <- cbind(x_coords, y_coords)
  
  # If hexagonal:
  if (hex_gridshape) {
    
    # identify indices of x coordinates to offset
    shift_idx <- rep(c(rep(FALSE, n_columns), rep(TRUE, n_columns)), sum(seq_len(n_rows) %% 2 == 0))
    if(n_rows %% 2 != 0) shift_idx <- c(shift_idx, c(rep(FALSE, n_columns)))
    
    # Offset coordinates for target indices
    xy_coords[shift_idx, 1] <- xy_coords[shift_idx, 1] + horizontal_spacing / 2
  }
  
  # Create lists of coordinates to bound hexagon cell for each centroid coordinate
  poly_coords <- lapply(seq_len(n_columns * n_rows), function(i) {
    
    x_coord <- xy_coords[i , 1]
    y_coord <- xy_coords[i , 2]
    
    ifelse(hex_gridshape,
           list(
             rbind(
               c(x_coord - (outer_radius / 2) , y_coord + inner_radius),
               c(x_coord + (outer_radius / 2) , y_coord + inner_radius),
               c(x_coord + outer_radius, y_coord),
               c(x_coord + (outer_radius / 2) , y_coord - inner_radius),
               c(x_coord - (outer_radius / 2) , y_coord - inner_radius),
               c(x_coord - outer_radius, y_coord),
               c(x_coord - (outer_radius / 2) , y_coord + inner_radius)
             )),
           list(
             rbind(
               c(x_coord - (horizontal_spacing / 2) , y_coord + (vertical_spacing / 2)),
               c(x_coord + (horizontal_spacing / 2) , y_coord + (vertical_spacing / 2)),
               c(x_coord + (horizontal_spacing / 2) , y_coord - (vertical_spacing / 2)),
               c(x_coord - (horizontal_spacing / 2) , y_coord - (vertical_spacing / 2)),
               c(x_coord - (horizontal_spacing / 2) , y_coord + (vertical_spacing / 2))
             ))
    )
  })
  
  # Create spatial polygons
  polys <- lapply(seq_len(n_columns * n_rows), function(i) st_sfc(st_polygon(poly_coords[[i]]), crs = st_crs(boundary)))
  polys <- do.call(c, polys)
  polys <- st_sf(data.frame("id" = 1:length(polys)), geometry = polys)
  
  # Crop to input boundary
  grid <- polys[st_intersection(polys, boundary)$id, ]
  
  # Reset numbering in ids
  grid$id <- seq_len(nrow(grid))
  
  grid
}

#### Geometric form functions

columnar <- function(adj_cut_height = 0, plant_height = 0) {
  ifelse(adj_cut_height < (plant_height * 0.33), 0, 1)
}

pyramidal <- function(adj_cut_height = 0, plant_height = 0) {
  ifelse(adj_cut_height < (plant_height * 0.33), 0, ((-1 * adj_cut_height) + plant_height) / (plant_height - (plant_height * 0.33)))
}

round <- function(adj_cut_height = 0, plant_height = 0) {
  sqrt(pmax(0, (plant_height * 0.33)^2 - (adj_cut_height - (plant_height * 0.66))^2)) / (plant_height * 0.33)
}

rounded <- function(adj_cut_height = 0, plant_height = 0) {
  sqrt(pmax(0, (plant_height * 0.5)^2 - (adj_cut_height - (plant_height * 0.5))^2)) / (plant_height * 0.5)
}

mounding <- function(adj_cut_height = 0, plant_height = 0) {
  1 - (adj_cut_height / plant_height)^2
}

vase <- function(adj_cut_height = 0, plant_height = 0) {
  (adj_cut_height / plant_height)
}

upright <- function(adj_cut_height = 0, plant_height = 0) {
  adj_cut_height / adj_cut_height
}

# function_names <- list("columnar", "pyramidal", "round", "rounded", "mounding", "vase", "upright")
# 
# par(mfrow = c(2, 4))
# 
# for (i in function_names) {
#   plot(seq(0, 2, 0.02), do.call(i, list(adj_cut_height = seq(0, 2, 0.02), plant_height = 2)),
#        type = 'l',
#        xlab = "",
#        ylab = "",
#        main = i)
# }

# Calculate the overall Shannon Diversity Evenness Index
shannon_evenness <- function(counts) {
  sp_props <- counts / sum(counts)
  sdi <- -sum(sp_props * log(sp_props))
  max_sdi <- log(length(counts)) 
  signif(sdi / max_sdi, 2)
}

extract_coverage <- function(spatial_list) {
  # Extract unique years from file list
  years <- names(spatial_list)
  n_year_groups <- length(years)
  
  # Extract unique heights from file list and convert to numeric and meter units
  heights <- names(spatial_list[[1]])
  heights <- as.numeric(substr(heights, (nchar(heights) + 1) - 4, nchar(heights))) / 100
  
  sapply(seq_len(n_year_groups), function(j) {
    
    spatial_year <- spatial_list[[j]]
    
    idx <- which(grepl("0100", names(spatial_year)))
    
    grid <- spatial_year[[idx]]
    
    grid$prop_ol[is.na(grid$prop_ol)] <- 0
    
    sum(grid$prop_ol) / length(grid$prop_ol)
    
  })
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