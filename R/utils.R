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
  
  # # Determine grid spacing
  # grid_spacing <- horizontal_spacing
  # if (hex_gridshape) {
  #   grid_spacing <- 2 * inner_radius
  # }
  
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
  
  # # Add attribute to identify grid spacing
  # attr(grid, 'gridspacing') <- grid_spacing
  
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
                              cut_heights = c(0.3, 1.0, 1.7, 2.4, 3.1), # Height datums in meters
                              years = 1, # Years of interest
                              hex_gridshape = TRUE, # Use hexagonal-shaped grid cells in the analysis; default is yes
                              cellarea, # Target cell area in square meters
                              veg_layers_path = NULL, # Path to folder where vegetation layers should be written
                              overlap_grids_path = NULL # Path to folder where grid layers should be written
) {
  
  # Create grid from boundary based on specified grid shape and approximate cell size
  master_grid <- create_grid(boundary = boundary, hex_gridshape = hex_gridshape, cellarea = cellarea)
  
  # # Get grid spacing
  # grid_spacing <- attr(master_grid, 'gridspacing')
  
  # Determine cell spacing distance
  cell_dist <- sqrt(2 * cellarea / (3 * sqrt(3))) * sqrt(3)
  if(hex_gridshape == FALSE) cell_dist <- sqrt(cellarea)
  
  # # Determine maximum height of vegetation plus 5%
  # max_veg_height <- max(point_locations$max_height) * 1.05
  # 
  # # Number of vertical layers at cell spacing distance within maximum vegetation height
  # n_vert_layers <- max_veg_height %/% cell_dist
  # 
  # # If cut_heights are not specified, use regular spacing based on distance between cells
  # if(is.null(cut_heights)) cut_heights <- (1:n_vert_layers) * cell_dist
  # 
  # # Set default switch for added 1m cut height
  # add_1m <- FALSE
  # 
  # # If 1m is not specified in cut heights, add it
  # if(!1 %in% cut_heights) {
  #   add_1m <- TRUE
  #   cut_heights <- sort(c(1, cut_heights))
  # }
  
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
  
  # # Add attribute to identify added cut
  # attr(out, 'cut_1m') <- FALSE
  # if(add_1m) attr(out, 'cut_1m') <- TRUE
  
  # Add attribute to identify gridshape
  attr(out, 'gridshape') <- "hex"
  if(hex_gridshape == FALSE) attr(out, 'gridshape') <- "rect"
  
  # Add attribute to identify grid spacing
  attr(out, 'gridspacing') <- cell_dist
  
  out
}

estimate_connectivity <- function(spatial_list,
                                  threshold = 0.3) {
  # Extract unique years from file list
  years <- names(spatial_list)
  n_year_groups <- length(years)
  
  # Extract unique heights from file list and convert to numeric and meter units
  heights <- names(spatial_list[[1]])
  heights <- as.numeric(substr(heights, (nchar(heights) + 1) - 4, nchar(heights))) / 100
  n_heights <- length(heights)
  
  master_grid <- spatial_list[[1]][[1]]
  
  grid_spacing <- attr(spatial_list, 'gridspacing')
  
  search_dist <- grid_spacing * 1.5
  
  buffers <- st_buffer(st_centroid(master_grid$geometry), search_dist)
  
  points <- st_centroid(master_grid$geometry)
  points <- st_sf(data.frame("id" = 1:length(points)), geometry = points)
  
  idx <- master_grid$id
  
  # Create master neighborhood list
  nb_list <- lapply(idx, function(p) {
    ids <- points$id[st_intersects(points, buffers[p], sparse = FALSE, prepared = FALSE)]
    ids <- setdiff(ids, p)
  })
  
  scores_df <- data.frame(matrix(NA, nrow = n_year_groups, ncol = n_heights))
  colnames(scores_df) <- names(spatial_list[[1]])
  
  for (i in seq_len(n_year_groups)) {
    
    spatial_year <- spatial_list[[i]]
    
    scores <- sapply(seq_len(n_heights), function(j) {
      
      grid <- spatial_year[[j]]
      
      grid$prop_ol[is.na(grid$prop_ol)] <- 0
      grid$prop_ol <- ifelse(grid$prop_ol > threshold, 1, 0)
      
      mean(sapply(idx, function(p) {
        sum(grid$prop_ol[nb_list[[p]]]) / length(sum(grid$prop_ol[nb_list[[p]]]))
      }))
      
    })
    
    scores_df[i, ] <- scores
    
  }
  
  scores_df
  
}

create_images <- function(spatial_list, path) {
  
  # Extract unique years from file list
  years <- names(spatial_list)
  n_year_groups <- length(years)
  
  # Extract unique heights from file list and convert to numeric and meter units
  cut_heights <- names(spatial_list[[1]])
  cut_heights <- as.numeric(substr(cut_heights, (nchar(cut_heights) + 1) - 4, nchar(cut_heights))) / 100
  
  sapply(seq_len(n_year_groups), function(j) {
    
    spatial_year <- spatial_list[[j]]
    
    plot_list <- list()
    
    for (i in 1:length(spatial_year)) {
      
      # Load processed spatial data
      spatial_data <- spatial_year[[i]]
      
      # Rotate spatial data
      spatial_data <- rotate_data(spatial_data)
      
      # Set annotation position
      x <- st_bbox(spatial_data)[3] * 0.999995
      y <- st_bbox(spatial_data)[2] * 1.000001 
      
      # Create plot of proportional overlap grid at currently specified height
      plot_list[[i]] <- ggplot() +
        geom_sf(data = spatial_data, aes(fill = prop_ol), color = 'gray40', size = 0.05, show.legend = FALSE) +
        annotate("text", label = paste0(cut_heights[i], " m"), x = x , y = y, hjust = 0, color = 'gray40', size = 2) +
        scale_fill_viridis_c(option = 'E', na.value = "transparent") +
        theme_void()
    } 
    
    # Write out PNG of overlap grids for all heights
    png(paste0(path, "veg_structure_", years[j], ".png"), height = 150 * length(spatial_year), width = 400, pointsize = 12, res = 300)
    grid.arrange(grobs = rev(plot_list), ncol = 1)
    dev.off()
    
    NULL # Because function is not meant to return anything
  })
}

plot_classes <- function(spatial_points, variable_name, colour_palette, image_path) {
  st_geometry(spatial_points) <- NULL
  spatial_points <- as.data.frame(spatial_points)
  
  height_ranges <- rbind(cbind(0, 0.75),
                         cbind(0.75, 1.5),
                         cbind(1.5, 3),
                         cbind(3, Inf))
  
  ifelse(variable_name == "size",
         data <- data.frame(
           category = c("0-0.75m", "0.75-1.5m", "1.5-3m", "+3m"),
           values = apply(height_ranges, 1, function(hts) length(which(spatial_points$max_height > hts[1] & spatial_points$max_height <= hts[2])))
         ),
         data <- data.frame(
           category = toupper(unique(spatial_points[ , variable_name])),
           values = sapply(unique(spatial_points[ , variable_name]), function(var) length(which(spatial_points[ , variable_name] == var)))
         ))
  
  # Compute percentages
  data$fraction = data$values / sum(data$values)
  
  # Compute the cumulative percentages (top of each rectangle)
  data$ymax = cumsum(data$fraction)
  
  # Compute the bottom of each rectangle
  data$ymin = c(0, head(data$ymax, n = -1))
  
  # Compute label position
  data$labelPosition <- (data$ymax + data$ymin) / 2
  
  # Make the plot
  ggplot(data, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = category)) +
    geom_rect() +
    geom_text(x = 3.5, aes(y = labelPosition, label = category), size = 4) +
    scale_fill_manual(values = colour_palette) +
    coord_polar(theta = "y") +
    xlim(c(0.8, 4)) +
    theme_void() +
    theme(legend.position = "none") +
    
    geom_image(data = data.frame(xx = 0.8, yy = 0, image = image_path), mapping = aes(xx, yy, image = image), size = .3, inherit.aes = FALSE)
}

plot_percent <- function(spatial_points, variable_name, colour, label) {
  st_geometry(spatial_points) <- NULL
  spatial_points <- as.data.frame(spatial_points)
  
  props <- table(spatial_points[ , variable_name])
  native_prop <- base::round((props[2] / sum(props)) * 100)
  exotic_prop <- base::round((props[1] / sum(props)) * 100)
  
  data <- data.frame(
    category = 1:100,
    values = c(rep(1, native_prop), rep(0, exotic_prop))
  )
  
  # Compute the cumulative percentages (top of each rectangle)
  data$ymax = 1:100
  
  # Compute the bottom of each rectangle
  data$ymin = c(0, head(data$ymax, n = -1))
  
  pr_palette <- c(rep(colour, native_prop),  rep("#e8e8e8", exotic_prop))
  
  # Make the plot
  ggplot(data, aes(ymax = ymax, ymin = ymin, xmax = 5, xmin = 1, fill = factor(category))) +
    geom_rect(color = "white", size = 1) +
    geom_point(aes(x = 0, y = -max(values) * 0.5), color = "white", size = 50) +
    #geom_text(x = 3.5, aes(y = labelPosition, label = toupper(category)), size = 4) +
    scale_fill_manual(values = pr_palette) +
    coord_polar(theta = "y") +
    xlim(c(0, 5)) +
    theme_void() +
    theme(legend.position = "none") +
    
    geom_text(data = data.frame(xx = 0, yy = -max(data$values) * 0.5, label = paste0(native_prop, "%\n", label)), mapping = aes(xx, yy, label = label), size = 4, inherit.aes = FALSE)
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