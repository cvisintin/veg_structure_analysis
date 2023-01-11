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
                                obstructions = building,
                                boundary = site_boundary,
                                years = 1:30,
                                cellarea = 10)

save(years_1_30, file = "output/years_1_30")

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
  
  out
}
