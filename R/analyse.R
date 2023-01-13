library(sf)
source('R/utility.R')

# Load background spatial data
building <- read_sf("data/gis/Building_Mass.shp")

# Load grids file
load(file = "output/years_1_30")

# Extract unique years from file list
years <- names(years_1_30)
n_year_groups <- length(years)

# Extract unique heights from file list and convert to numeric and meter units
cut_heights <- names(years_1_30[[1]])
cut_heights <- as.numeric(substr(cut_heights, (nchar(cut_heights) + 1) - 4, nchar(cut_heights))) / 100

# Identify shape of grid cells
grid_shape <- attr(years_1_30, "gridshape")

# Convert data for connectivity analysis (note object names use centimeter units)
for(i in seq_len(n_year_groups)) {
  
  grid <- years_1_30[[i]]
  
  for (j in 1:length(cut_heights)) {
    assign(paste0(years[i], "_", cut_heights[j] * 100), prepare_conn_data(grid = read_sf(target_files[j]),
                                                                          obstructions = building,
                                                                          cut_height = cut_heights[j],
                                                                          grid_shape = grid_shape))
  }
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
  
  # Add attribute to identify grid cell shape
  attr(grid, 'gridshape') <- grid_shape
  
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
