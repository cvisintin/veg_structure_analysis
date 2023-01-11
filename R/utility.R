library(sf)
library(dplyr)

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
  grid <- polys[which(sapply(st_intersects(polys, boundary), function(x) length(x)) == 1), ]
  
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

# Calculate the overall Shannon Diversity Evenness Index
shannon_evenness <- function(counts) {
  sp_props <- counts / sum(counts)
  sdi <- -sum(sp_props * log(sp_props))
  max_sdi <- log(length(counts)) 
  signif(sdi / max_sdi, 2)
}

# Extract the proportion coverage of vegetation per grid cell from spatial data at 1 meter height
extract_1m_coverage <- function(spatial_list) {
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