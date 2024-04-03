library(sf)
library(dplyr)
library(geos)
library(ggplot2)
library(ggimage)
library(grid)
library(gridExtra)
library(png)
library(tidyr)
library(viridis)
library(spatstat)

# Determine spatial patterning and assess whether there is uniform distribution or not
analyse_spatial_patterning <- function(spatial_points, boundary) {
  species <- unique(spatial_points$species)
  n_species <- length(species)
  scores <- rep(NA, n_species)
  for (i in seq_len(n_species)) {
    idx <- spatial_points$species == species[i]
    ppp <- ppp(x = st_coordinates(spatial_points)[idx, 1],
               y = st_coordinates(spatial_points)[idx, 2],
               window = as.owin(st_bbox(boundary)))
    scores[i] <- ifelse(ppp$n > 2,
                        {K <- envelope(ppp, silent = TRUE, verbose = FALSE)
                        auc_theo <- sum(diff(seq_len(length(K$theo))) * (head(K$theo, -1) + tail(K$theo, -1))) / 2
                        auc_obs <- sum(diff(seq_len(length(K$obs))) * (head(K$obs, -1) + tail(K$obs, -1))) / 2
                        auc_prop <- auc_theo / auc_obs
                        if(auc_prop <= 1) auc_score <- auc_prop
                        if(auc_prop > 1) auc_score <- 1 / auc_prop
                        auc_score
                        }, 0)
  }
  scores_idx <- which(scores > 0)
  final_score <- mean(scores[scores_idx])
  final_score
}

# Verify that the point and polygon spatial objects have the same properties and attributes
check_spatial_input <- function(point_locations, polygon_locations = NULL) {
  # Verify the geometries
  if(!all(st_geometry_type(point_locations) == "POINT")) stop("All specified geometry are not points, check input geometry type")
  
  if(!is.null(polygon_locations)) {
    if(!all(st_geometry_type(polygon_locations) == "POLYGON")) stop("All specified geometry are not polygons, check input geometry type")
  }
  
  # Verify the correct column names are used
  correct_names <- c("coverage",
                     "density",
                     "endemism",
                     "form",
                     "ini_height",
                     "max_height",
                     "max_width",
                     "phenology",
                     "ref_height",
                     "spacing",
                     "species",
                     "texture",
                     "year_max" )
  
  point_data <-as.data.frame(st_drop_geometry(point_locations))
  if(!identical(sort(names(point_data)), correct_names)) stop("Check column names in point data; they are not as required")
  
  if(!is.null(polygon_locations)) {
    polygon_data <-as.data.frame(st_drop_geometry(polygon_locations))
    if(!identical(sort(names(point_data)), correct_names)) stop("Check column names in polygon data; they are not as required")  }
  
  NULL 
}

# Create axonometric images of grids showing proportional coverage of vegetation
# at specified heights and for all specified timesteps 
create_images <- function(spatial_list, path, filename) {
  
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
      x <- st_bbox(spatial_data)[3] #* 0.999995
      y <- st_bbox(spatial_data)[2] * 1.000001
      x2 <- st_bbox(spatial_data)[1]
      y2 <- st_bbox(spatial_data)[4]
      
      # Create plot of proportional overlap grid at currently specified height
      plot_list[[i]] <- ggplot() +
        geom_sf(data = spatial_data, aes(fill = prop_ol), color = 'gray40', lwd = 0.03, show.legend = FALSE) +
        annotate("text",
                 label = paste0(cut_heights[i], " m"),
                 x = x ,
                 y = y,
                 hjust = "inward",
                 color = 'gray20',
                 size = 1,
                 fontface = "bold") +
        {if(i == length(spatial_year)) annotate("text",
                                                label = paste0("Year ", j),
                                                x = x2 ,
                                                y = y2,
                                                hjust = "inward",
                                                color = 'gray20',
                                                size = 1,
                                                fontface = "bold")} +
        scale_fill_viridis_c(option = 'E', na.value = "transparent") +
        theme_void()
      
    } 
    
    # Write out PNG of overlap grids for all heights
    png(paste0(path, filename, "_", years[j], ".png"), height = 150 * length(spatial_year), width = 400, pointsize = 12, res = 300)
    grid.arrange(grobs = rev(plot_list), ncol = 1)
    dev.off()
    
    NULL # Because function is not meant to return anything
  })
}

# Create a web-based, interactive score sheet to assess the outcome of a landscape design 
create_interactive_score_sheet <- function(spatial_points,
                                           boundary = NULL,
                                           spatial_list = NULL,
                                           path_directory) {
  
  ### Density ###
  print("Preparing density plot...")
  png(paste0(path_directory, "density.png"), height = 400, width = 400, pointsize = 4, res = 150)
  print(plot_classes(spatial_points = spatial_points,
                     variable_name = "density",
                     colour_palette = c("#397d53", "#b7e4c8", "#62a67c"),
                     image_path = "data/images/leaves_icon.png"))
  dev.off()
  
  ### Texture ###
  print("Preparing texture plot...")
  png(paste0(path_directory, "texture.png"), height = 400, width = 400, pointsize = 4, res = 150)
  print(plot_classes(spatial_points = spatial_points,
                     variable_name = "texture",
                     colour_palette = c("#4b6c90", "#afc6e0", "#6d8eb3"),
                     image_path = "data/images/texture_icon.png"))
  dev.off()
  
  ### Size ###
  print("Preparing sizes plot...")
  png(paste0(path_directory, "size.png"), height = 400, width = 400, pointsize = 4, res = 150)
  print(plot_classes(spatial_points = spatial_points,
                     variable_name = "size",
                     colour_palette = c("#b8641d", "#ebc5a4", "#d9a171", "#c98449"),
                     image_path = "data/images/vegetation_icon.png"))
  dev.off()
  
  ### Endemism ###
  print("Preparing endemism plot...")
  png(paste0(path_directory, "endemism.png"), height = 400, width = 400, pointsize = 4, res = 150)
  print(plot_percent(spatial_points = spatial_points,
                     variable_name = "endemism",
                     target_value = "native",
                     colour = "#a072a6",
                     label = "NATIVE"))
  dev.off()
  
  ### Patterning ###
  print("Preparing patterning plot...")
  png(paste0(path_directory, "patterning.png"), height = 400, width = 400, pointsize = 4, res = 150)
  print(plot_percent(spatial_points = spatial_points,
                     boundary = boundary,
                     variable_name = "patterning",
                     colour = "#a19a65",
                     label = "DISTRIBUTION\nSCORE"))
  dev.off()
  
  ### Species richness ###
  print("Preparing species plot...")
  png(paste0(path_directory, "richness.png"), height = 400, width = 400, pointsize = 4, res = 150)
  print(plot_circ_bar(spatial_points = spatial_points,
                      variable_name = "richness",
                      colours = "#858383",
                      polar_rotation = 0.25))
  dev.off()
  
  ### Phenology ###
  print("Preparing phenology plot...")
  png(paste0(path_directory, "phenology.png"), height = 400, width = 400, pointsize = 4, res = 150)
  print(plot_circ_bar(spatial_points = spatial_points,
                      variable_name = "phenology",
                      colours = "#b0d4d6",
                      polar_rotation = 0.25))
  dev.off()
  
  if(!is.null(spatial_list)) {
    ### Coverage at 1m ###
    print("Preparing coverage plot...")
    png(paste0(path_directory, "coverage.png"), height = 400, width = 400, pointsize = 4, res = 150)
    print(plot_circ_bar(spatial_points = spatial_points,
                        spatial_list = spatial_list,
                        variable_name = "coverage",
                        colours = "#c9837d"))
    dev.off()
    
    ### Connectivity ###
    print("Preparing connectivity plot...")
    png(paste0(path_directory, "connectivity.png"), height = 400, width = 400, pointsize = 4, res = 150)
    print(plot_circ_bar(spatial_points = spatial_points,
                        spatial_list = spatial_list,
                        #obstructions = obstructions,
                        variable_name = "connectivity",
                        stacked = TRUE))
    dev.off()
  }
  
  print("Creating score sheet.")
  file.copy("data/main.css", paste0(path_directory, "main.css"), overwrite = TRUE)
  file.copy("data/index.html", paste0(path_directory, "index.html"), overwrite = TRUE)
  
  if(!is.null(spatial_list)) {
    file.copy("data/index_full.html", paste0(path_directory, "index.html"), overwrite = TRUE)
  }
  
}

# Create a grid over site boundary to organise and analyse data
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
  
  # If hexagonal, provide offsets
  if (hex_gridshape) {
    
    # Identify indices of x coordinates to offset
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

# Convert polygons to points based on spacing and proportion coverage of plants
convert_combine <- function(point_locations, # Spatial geometry and attributes of plants
                            polygon_locations # Spatial geometry and attributes of planted areas
) {
  polygon_locations <- st_buffer(polygon_locations, 0)
  
  check_spatial_input(point_locations = point_locations, polygon_locations = polygon_locations)
  
  n_polys <- nrow(polygon_locations)
  
  point_list <- lapply(seq_len(n_polys), function(i) {
    data <- st_drop_geometry(polygon_locations[i, ])
    cell_size <- 1 / sqrt(polygon_locations$spacing[i])
    grid <- st_make_grid(st_bbox(polygon_locations$geometry[i]), cellsize = cell_size, what = "centers")
    points <- st_intersection(grid, polygon_locations$geometry[i])
    n_points <- length(points)
    points <- points[sample(seq_len(n_points), pmax(1, floor(n_points * polygon_locations$coverage[i] / 100)))]
    st_sf(data, geometry = points)
  })
  
  new_plant_points <- do.call(rbind, point_list)
  
  rbind(point_locations, new_plant_points)
}

# Estimate how well the site vegetation is connected using a threshold value
# for proportion of vegetation coverage in each cell
estimate_connectivity <- function(spatial_list, threshold = 0.3) {
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
      
      # Based on threshold value
      grid$prop_ol <- ifelse(grid$prop_ol > threshold, 1, 0)
      
      # # Based on Bernoulli probabilities
      # grid$prop_ol <- base::round(grid$prop_ol, 3)
      # grid$prop_ol <- rbinom(n = length(grid$prop_ol), size = 1, prob = grid$prop_ol)
      
      mean(sapply(idx, function(p) {
        sum(grid$prop_ol[nb_list[[p]]]) / length(nb_list[[p]])
      }))
      
    })
    
    scores_df[i, ] <- scores
    
  }
  
  scores_df
  
}

# Export a data file for use in an interactive Shiny application
export_image_set <- function(spatial_list, path, filename) {
  
  # Extract unique years from file list
  years <- names(spatial_list)
  n_year_groups <- length(years)
  
  # Extract unique heights from file list and convert to numeric and meter units
  cut_heights <- names(spatial_list[[1]])
  cut_heights <- as.numeric(substr(cut_heights, (nchar(cut_heights) + 1) - 4, nchar(cut_heights))) / 100
  
  plot_set <- sapply(seq_len(n_year_groups), function(j) {
    print(paste0("Processing year ", j))
    
    spatial_year <- spatial_list[[j]]
    
    plot_list <- list()
    
    for (i in 1:length(spatial_year)) {
      
      # Load processed spatial data
      spatial_data <- spatial_year[[i]]
      
      # Rotate spatial data
      spatial_data <- rotate_data(spatial_data)
      
      # Set annotation position
      x <- st_bbox(spatial_data)[3] #* 0.999995
      y <- st_bbox(spatial_data)[2] * 1.000001
      x2 <- st_bbox(spatial_data)[1]
      y2 <- st_bbox(spatial_data)[4]
      
      # Create plot of proportional overlap grid at currently specified height
      plot_list[[i]] <- ggplot() +
        geom_sf(data = spatial_data, aes(fill = prop_ol), color = 'gray20', lwd = 0.03, show.legend = FALSE) +
        annotate("text",
                 label = paste0(cut_heights[i], " m"),
                 x = x ,
                 y = y,
                 hjust = "inward",
                 color = 'gray20',
                 size = 5,
                 fontface = "bold") +
        {if(i == length(spatial_year)) annotate("text",
                                                label = paste0("Year ", j),
                                                x = x2 ,
                                                y = y2,
                                                hjust = "inward",
                                                color = 'gray20',
                                                size = 5,
                                                fontface = "bold")} +
        scale_fill_viridis_c(option = 'E', na.value = "transparent") +
        theme_void()
      
    } 
    
    arrangeGrob(grobs = rev(plot_list), ncol = 1)
  })
  
  save(plot_set, file = paste0(path, filename))
  
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

# Build the vegetation structure over time and horizontally cut at specified heights
process_structure <- function(point_locations, # Spatial geometry and attributes of plants
                              obstructions = NULL, # Geometry that may obstruct plant growth in horizontal plane
                              boundary, # Geometry that defines the boundary of the site
                              cut_heights = c(0.3, 1.0, 1.7, 2.4, 3.1), # Height datums in meters
                              years = 1, # Years of interest
                              hex_gridshape = TRUE, # Use hexagonal-shaped grid cells in the analysis; default is yes
                              cellarea, # Target cell area in square meters
                              veg_layers_path = NULL, # Path to folder where vegetation layers should be written
                              overlap_grids_path = NULL # Path to folder where grid layers should be written
) {
  # Verify input data
  check_spatial_input(point_locations = point_locations)
  
  # Create grid from boundary based on specified grid shape and approximate cell size
  master_grid <- create_grid(boundary = boundary, hex_gridshape = hex_gridshape, cellarea = cellarea)
  
  # Check for obstruction data
  obstruct <- !is.null(obstructions)
  
  # # Get grid spacing
  # grid_spacing <- attr(master_grid, 'gridspacing')
  
  # Determine cell spacing distance
  cell_dist <- sqrt(2 * cellarea / (3 * sqrt(3))) * sqrt(3)
  if(hex_gridshape == FALSE) cell_dist <- sqrt(cellarea)
  
  # Convert geometry to geos structure for faster processing
  point_locations_geom <- as_geos_geometry(point_locations)
  if (obstruct) obstructions_geom <- as_geos_geometry(obstructions)
  
  out <- lapply(years, function(year) {  
    
    print(paste0("Processing year ", year))
    
    # Calculate ages in years from initial heights
    plant_ages <- (point_locations$ini_height / point_locations$max_height) * point_locations$year_max
    
    # Calculate all relative plant heights and widths for the given year (from reference height)
    plant_heights <- point_locations$max_height * (pmin(year + plant_ages, point_locations$year_max) / point_locations$year_max)
    plant_widths <- point_locations$max_width * (pmin(year + plant_ages, point_locations$year_max) / point_locations$year_max)
    
    # Generate id lists for plants that intersect each height datum for the given year
    plant_idx <- lapply(cut_heights, function(cut_height) which(plant_heights + point_locations$ref_height >= cut_height & cut_height >= point_locations$ref_height))
    
    proportion_grids <- lapply(cut_heights, function(cut_height) {
      
      # Copy master grid
      grid <- master_grid
      
      # Convert geometry to geos structure for faster processing
      grid_geom <- as_geos_geometry(grid)
      
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
        planted_area <- geos_unary_union(geos_make_collection(geos_buffer(geom = point_locations_geom[ids], distance = (plant_widths[ids] / 2) * mults)))
        
        if (obstruct) {
          # Determine which obstructions are above cut height
          obstruction_idx <- obstructions$height > cut_height
          
          if (!is.null(obstruction_idx)) {
            # Remove areas with obstructions (e.g. buildings)
            planted_area <- geos_difference(planted_area, geos_unary_union(geos_make_collection(obstructions_geom[obstruction_idx])))
          }
        }
        
        # If path is provided, write out vegetation footprints for given year and selected heights
        # using centimeters to designate height (padded with zeros)
        if(!is.null(veg_layers_path)) {
          st_write(st_as_sfc(planted_area), paste0(veg_layers_path,
                                                   "year_",
                                                   sprintf("%02d", year),
                                                   "_height_",
                                                   sprintf("%04d", ceiling(cut_height * 100)),
                                                   "_veg.shp"), append = FALSE)
        }
        
        # Identify the grid cells that have overlap with vegetation
        overlap_idx <- geos_intersects_any(grid_geom, geos_strtree(planted_area))
        
        # Determine the area of overlap in each qualified grid cell
        areas <- geos_area(geos_intersection(grid_geom[overlap_idx], planted_area))
        
        # Record proportions of overlap to eligible grids
        grid$prop_ol[overlap_idx] <- areas / st_area(grid$geometry[1])
        
      }
      
      # If path is provided, write out grids for given year and selected heights using centimeters
      # to designate height (padded with zeros)
      if(!is.null(overlap_grids_path)) {
        st_write(st_as_sfc(grid), paste0(overlap_grids_path,
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

# Plot multidimensional data as stacked bars in a radial formation
plot_circ_bar <- function(spatial_points,
                          spatial_list = NULL,
                          #obstructions = NULL,
                          variable_name,
                          colours = NULL,
                          polar_rotation = NULL,
                          stacked = FALSE) {
  
  if(variable_name == "richness") {
    species <- sub("(\\w+\\s+\\w+).*", "\\1", spatial_points$species)
    
    data <- data.frame(
      id = factor(seq(1, length(unique(species)), 1)),
      value = sapply(unique(species), function(sp) length(which(species == sp)))
    )
  }
  
  if(variable_name == "phenology") {
    months <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")
    n_months <- length(months)
    
    data <- data.frame(
      id = factor(seq_len(n_months)),
      value = sapply(months, function(month) sum(grepl(month, spatial_points$phenology)))
    )
    
    data$value <- zero_to_one(data$value)
    
    img <- readPNG("data/images/flower_icon.png")
    g <- rasterGrob(img, interpolate = TRUE)
  }
  
  if(variable_name == "coverage" & !is.null(spatial_list)) {
    n_year_groups <- length(spatial_list)
    
    data <- data.frame(
      id = factor(seq_len(n_year_groups)),
      value = extract_1m_coverage(spatial_list),
      max = 1
    )
  }
  
  if(variable_name == "connectivity" ) {
    n_year_groups <- length(spatial_list)
    
    data_w <- estimate_connectivity(spatial_list)
    
    
    max_value <- ncol(data_w)
    
    score <- base::round(mean(as.matrix(data_w)), 2)
    
    data_w$remainder <- apply(data_w, 1, function(x) ncol(data_w) - sum(x))
    
    data <- gather(data_w, key = "group", value = "value")
    
    data$id <- factor(rep(1:n_year_groups, ncol(data_w)))
    data$group <- factor(data$group, levels = rev(unique(data$group)))
    
    colours <- c(viridis(ncol(data_w) - 1, end = 0.8), "#e1e1e1")
  }
  
  if(!stacked & variable_name != "coverage") p <- ggplot(data, aes(x = as.factor(id), y = value))
  if(!stacked & variable_name == "coverage") p <- ggplot(data, aes(x = as.factor(id), y = max))
  if(stacked) p <- ggplot(data, aes(x = id, y = value, fill = group))
  
  if(!stacked & variable_name != "coverage") p <- p + geom_bar(stat = "identity", fill = colours)
  if(!stacked & variable_name == "coverage") {
    p <- p + geom_bar(data, mapping = aes(x = as.factor(id), y = max), stat = "identity", fill = "#e1e1e1")
    p <- p + geom_bar(data, mapping = aes(x = as.factor(id), y = value), stat = "identity", fill = colours)
    }
  if(stacked) p <- p + geom_bar(position = "stack", stat = "identity")
  if(stacked) p <- p + scale_fill_manual(values = rev(colours))
  
  # Limits of the plot. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
  if(!stacked & variable_name != "coverage") p <- p + ylim(-max(data$value) * 1.1, max(data$value) * 1.5)
  if(!stacked & variable_name == "coverage") p <- p + ylim(-1.1, 1.15)
  if(stacked) p <- p + ylim(-max_value * 1.1, max_value * 1.3)
  
  # Custom the theme: no axis title and no cartesian grid
  p <- p + theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1, 4), "cm")
    )
  
  if(stacked) p <- p + theme(legend.position = "none")
  
  # This makes the coordinate polar instead of cartesian
  if(is.null(polar_rotation)) p <- p + coord_polar(start = -(2 / nrow(data)))
  if(!is.null(polar_rotation)) p <- p + coord_polar(start = polar_rotation)
  
  # Add additional annotation
  if(variable_name == "richness") {
    p <- p + geom_text(data = data.frame(xx = 0.5,
                                         yy = -max(data$value),
                                         label = paste0(shannon_evenness(data$value), "\nEVENNESS\nSCORE")),
                       mapping = aes(xx, yy, label = label),
                       size = 4,
                       inherit.aes = FALSE)
  }
  
  if(variable_name == "phenology") {
    p <- p + geom_text(data = data,
                       mapping = aes(x = id, y = value * 0.5, label = toupper(row.names(data))),
                       color = "black",
                       fontface = "bold",
                       alpha = 0.6,
                       size = 2.5) +
      
      geom_image(data = data.frame(xx = 0.5,
                                   yy = -max(data$value),
                                   image = "data/images/flower_icon.png"),
                 mapping = aes(xx, yy, image = image),
                 size = .2,
                 inherit.aes = FALSE)
  }
  
  if(variable_name == "coverage" & !is.null(spatial_list)) {
    p <- p + geom_text(data = data,
                       mapping = aes(x = id, y = 1, label = paste0("Y", toupper(row.names(data)))),
                       color = "black",
                       fontface = "bold",
                       alpha = 0.6,
                       size = 2.5) +
      
      geom_image(data = data.frame(xx = 1,
                                   yy = -1.1,
                                   image = "data/images/shrub_icon.png"),
                 mapping = aes(xx, yy, image = image),
                 size = .25,
                 inherit.aes = FALSE)
  }
  
  if(variable_name == "connectivity") {
    p <- p + geom_text(data = data[1:nrow(data_w), ],
                       mapping = aes(x = unique(id), y = max_value, label = paste0("Y", row.names(data_w))),
                       color = "black",
                       fontface = "bold",
                       alpha = 0.6,
                       size = 2.5) +
      
      geom_text(data = data.frame(xx = max_value,
                                  yy = -max_value * 1.1,
                                  label = paste0(score, "\nCONNECTIVITY\nSCORE")),
                mapping = aes(xx, yy, label = label),
                size = 4,
                inherit.aes = FALSE)
  }
  
  p
}

# Plot categories as bars on a ring
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
  
  # Check for any zero values
  idx <- which(data$values > 0)
  
  # Make the plot
  ggplot(data[idx, ], aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = category)) +
    geom_rect() +
    geom_text(x = 3.5, aes(y = labelPosition, label = category), size = 2.5) +
    scale_fill_manual(values = colour_palette) +
    coord_polar(theta = "y") +
    xlim(c(0.8, 4)) +
    theme_void() +
    theme(legend.position = "none") +
    
    geom_image(data = data.frame(xx = 0.8, yy = 0, image = image_path), mapping = aes(xx, yy, image = image), size = .3, inherit.aes = FALSE)
}

# Plot proportion of categories as ticks on a ring
plot_percent <- function(spatial_points,
                         boundary = NULL,
                         variable_name,
                         target_value,
                         colour,
                         label) {
  
  if(variable_name == "patterning") {
    score <- analyse_spatial_patterning(spatial_points, boundary)
    zero_prop <- base::round((1 - score) * 100)
    one_prop <- base::round(score * 100)
  }
  
  if(variable_name == "endemism") {
    st_geometry(spatial_points) <- NULL
    spatial_points <- as.data.frame(spatial_points)
    
    props <- table(spatial_points[ , variable_name])
    
    if(length(props) == 1 && names(props) %in% target_value) {
      props <- c(0, props)
    }
    
    zero_prop <- base::round((props[1] / sum(props)) * 100)
    one_prop <- base::round((props[2] / sum(props)) * 100)
    
  }
  
  data <- data.frame(
    category = 1:100,
    values = c(rep(0, zero_prop), rep(1, one_prop))
  )
  
  # Compute the cumulative percentages (top of each rectangle)
  data$ymax = 1:100
  
  # Compute the bottom of each rectangle
  data$ymin = c(0, head(data$ymax, n = -1))
  
  pr_palette <- c(rep(colour, one_prop),  rep("#e8e8e8", zero_prop))

  # Make the plot
  ggplot(data, aes(ymax = ymax, ymin = ymin, xmax = 5, xmin = 1, fill = factor(category))) +
    geom_rect(color = "white", linewidth = 1) +
    geom_point(aes(x = 0, y = -max(values) * 0.5), color = "white", size = 50) +
    #geom_text(x = 3.5, aes(y = labelPosition, label = toupper(category)), size = 4) +
    scale_fill_manual(values = pr_palette) +
    coord_polar(theta = "y") +
    xlim(c(0, 5)) +
    theme_void() +
    theme(legend.position = "none") +
    
    geom_text(data = data.frame(xx = 0,
                                yy = -max(data$values) * 0.5,
                                label = ifelse(variable_name == "patterning",
                                               paste0(base::round(one_prop / 100, digits = 2), "\n", label),
                                               paste0(one_prop, "%\n", label))),
              mapping = aes(xx, yy, label = label),
              size = 4,
              inherit.aes = FALSE)
}

# Original function by Stefan Jünger
# https://stefanjuenger.github.io/gesis-workshop-geospatial-techniques-R/slides/2_4_Advanced_Maps_II/2_4_Advanced_Maps_II.html#8
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

# Calculate the overall Shannon Diversity Evenness Index
shannon_evenness <- function(counts) {
  sp_props <- counts / sum(counts)
  sdi <- -sum(sp_props * log(sp_props))
  max_sdi <- log(length(counts)) 
  signif(sdi / max_sdi, 2)
}

# Rescale values to be between zero and one
zero_to_one <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

#### Geometric form functions ####
columnar <- function(adj_cut_height = 0, plant_height = 0) {
  ifelse(adj_cut_height < (plant_height * 0.33), 0, 1)
}

mounding <- function(adj_cut_height = 0, plant_height = 0) {
  sqrt(plant_height^2 - adj_cut_height^2)
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

upright <- function(adj_cut_height = 0, plant_height = 0) {
  adj_cut_height / adj_cut_height
}

vase <- function(adj_cut_height = 0, plant_height = 0) {
  (adj_cut_height / plant_height)
}
