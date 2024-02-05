library(sf)
source('R/utils.R')

# Load background spatial data for residential project example
home_ex_plant_points <- read_sf("data/gis/Plant_Data_Home.shp")
home_ex_planting_polygons <- read_sf("data/gis/Planting_Data_Home.shp")
home_ex_building <- read_sf("data/gis/Building_Mass_Home.shp")
home_ex_site_boundary <- read_sf("data/gis/Site_Boundary_Home.shp")

# Convert planting polygons to points and combine with existing plant points
home_ex_plant_points <- convert_combine(point_locations = home_ex_plant_points,
                                        polygon_locations = home_ex_planting_polygons)

# Process vegetation structure data for specified years 
spatial_list <- process_structure(point_locations = home_ex_plant_points,
                                  obstructions = home_ex_building,
                                  boundary = home_ex_site_boundary,
                                  years = 1:10,
                                  cellarea = 10)

n_year_groups <- length(spatial_list)

thresholds <- c(.001, seq(0.1, 0.9, 0.2), .999)

con_list <- lapply(thresholds, function(x) estimate_connectivity(spatial_list = spatial_list,
                                                                 threshold = x))

con_thresh_scores <- matrix(NA, nrow = length(con_list), ncol = 2)

con_thresh_scores[, 1] <- thresholds

for(i in 1:length(con_list)) {
  data_w <- con_list[[i]]
  
  score <- base::round(mean(as.matrix(data_w)), 2)
  
  con_thresh_scores[i, 2] <- score
  
  data_w$remainder <- apply(data_w, 1, function(x) ncol(data_w) - sum(x))
  
  data <- gather(data_w, key = "group", value = "value")
  
  data$id <- factor(rep(1:n_year_groups, ncol(data_w)))
  data$group <- factor(data$group, levels = rev(unique(data$group)))
  
  max_value <- max(data$value)
  colours <- c(viridis(ncol(data_w) - 1, end = 0.8), "#e1e1e1")
  
  p <- ggplot(data, aes(x = id, y = value, fill = group))
  p <- p + geom_bar(position = "stack", stat = "identity")
  p <- p + scale_fill_manual(values = rev(colours))
  p <- p + ylim(-max(data$value) * 1.1, max(data$value) * 1.3)
  p <- p + theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1, 4), "cm")
    )
  p <- p + theme(legend.position = "none")
  p <- p + coord_polar(start = -(2 / nrow(data)))
  p <- p + geom_text(data = data[1:nrow(data_w), ],
                     aes(x = unique(id), y = max_value, label = paste0("Y", row.names(data_w))),
                     color = "black",
                     fontface = "bold",
                     alpha = 0.6,
                     size = 2.5) +
    
    geom_text(data = data.frame(xx = max(data$value),
                                yy = -max(data$value) * 1.1,
                                label = paste0(score, "\nCONNECTIVITY\nSCORE")),
              mapping = aes(xx, yy, label = label),
              size = 4,
              inherit.aes = FALSE)
  print(p)
  
}

plot(con_thresh_scores[, 1], con_thresh_scores[, 2], type = 'l', xlab = "Threshold", ylab = "Connectivity Score")


estimate_connectivity <- function(spatial_list,
                                  proportion_target = NULL) {
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
  
  # Determine which points to exclude at each cut height
  #obstruction_idx <- which(obstructions$height >= cut_height)
  
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
      
      idx_0 <- which(grid$prop_ol == 0)
      
      n_grid <- length(grid$prop_ol)
      
      proportion_init <- (n_grid - length(idx_0)) / n_grid
      
      if(!is.null(proportion_target)) {
        proportion_change <- pmax(0, proportion_target - proportion_init)
        
        if(proportion_change > 0){
          sample_size <- pmax(1, base::round(proportion_change * n_grid, 0))
          grid$prop_ol[sample(idx_0, sample_size)] <- runif(sample_size, min = 0.1)
        }
      }
      
      # Based on Bernoulli probabilities
      grid$prop_ol <- base::round(grid$prop_ol, 3)
      grid$prop_ol <- rbinom(n = length(grid$prop_ol), size = 1, prob = grid$prop_ol)
      
      mean(sapply(idx, function(p) {
        sum(grid$prop_ol[nb_list[[p]]]) / length(nb_list[[p]])
      }))
      
    })
    
    scores_df[i, ] <- scores
    
  }
  
  scores_df
  
}

proportion_targets <- c(.001, seq(0.1, 0.9, 0.2), .999)

con_list <- lapply(proportion_targets, function(x) estimate_connectivity(spatial_list = spatial_list,
                                                                         proportion_target = x))

con_thresh_scores <- matrix(NA, nrow = length(con_list), ncol = 2)

con_thresh_scores[, 1] <- proportion_targets

for(i in 1:length(con_list)) {
  data_w <- con_list[[i]]
  
  score <- base::round(mean(as.matrix(data_w)), 2)
  
  con_thresh_scores[i, 2] <- score
  
  max_value <- ncol(data_w)
  
  data_w$remainder <- apply(data_w, 1, function(x) max_value - sum(x))
  
  data <- gather(data_w, key = "group", value = "value")
  
  data$id <- factor(rep(1:n_year_groups, ncol(data_w)))
  data$group <- factor(data$group, levels = rev(unique(data$group)))
  
  colours <- c(viridis(ncol(data_w) - 1, end = 0.8), "#e1e1e1")
  
  p <- ggplot(data, aes(x = id, y = value, fill = group))
  p <- p + geom_bar(position = "stack", stat = "identity")
  p <- p + scale_fill_manual(values = rev(colours))
  p <- p + ylim(-max_value * 1.1, max_value * 1.3)
  p <- p + theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1, 4), "cm")
    )
  p <- p + theme(legend.position = "none")
  p <- p + coord_polar(start = -(2 / nrow(data)))
  p <- p + geom_text(data = data[1:nrow(data_w), ],
                     aes(x = unique(id), y = max_value, label = paste0("Y", row.names(data_w))),
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
  print(p)
  
}

plot(con_thresh_scores[, 1], con_thresh_scores[, 2], type = 'l', xlab = "Threshold", ylab = "Connectivity Score")


############## Create example landscapes to test metrics ################

clumped_sparse_grid <- read_sf("data/gis/Sample_Grids/grid_sparse_clumped.shp")
dispersed_sparse_grid <- read_sf("data/gis/Sample_Grids/grid_sparse_dispersed.shp")
clumped_dense_grid <- read_sf("data/gis/Sample_Grids/grid_dense_clumped.shp")
dispersed_dense_grid <- read_sf("data/gis/Sample_Grids/grid_dense_dispersed.shp")

cell_idx <- sample(1:nrow(dispersed_dense_grid), 20)

near_full_grid <- dispersed_dense_grid
near_full_grid$prop_ol <- 1
near_full_grid$prop_ol[cell_idx] <- 0
near_zero_grid <- dispersed_dense_grid
near_zero_grid$prop_ol <- 0
near_zero_grid$prop_ol[cell_idx] <- 1

grid_rect <- create_grid(home_ex_site_boundary, cellarea = 10, hex_gridshape = FALSE)

clumped_sparse_grid_rect <- grid_rect
clumped_sparse_grid_rect$prop_ol <- 0
clumped_sparse_grid_rect$prop_ol <- clumped_sparse_grid$prop_ol[st_nearest_feature(clumped_sparse_grid_rect, clumped_sparse_grid)]

dispersed_sparse_grid_rect <- grid_rect
dispersed_sparse_grid_rect$prop_ol <- 0
dispersed_sparse_grid_rect$prop_ol <- dispersed_sparse_grid$prop_ol[st_nearest_feature(dispersed_sparse_grid_rect, dispersed_sparse_grid)]

clumped_dense_grid_rect <- grid_rect
clumped_dense_grid_rect$prop_ol <- 0
clumped_dense_grid_rect$prop_ol <- clumped_dense_grid$prop_ol[st_nearest_feature(clumped_dense_grid_rect, clumped_dense_grid)]

dispersed_dense_grid_rect <- grid_rect
dispersed_dense_grid_rect$prop_ol <- 0
dispersed_dense_grid_rect$prop_ol <- dispersed_dense_grid$prop_ol[st_nearest_feature(dispersed_dense_grid_rect, dispersed_dense_grid)]

near_full_grid_rect <- dispersed_dense_grid_rect
near_full_grid_rect$prop_ol <- 1
near_full_grid_rect$prop_ol[cell_idx] <- 0
near_zero_grid_rect <- dispersed_dense_grid_rect
near_zero_grid_rect$prop_ol <- 0
near_zero_grid_rect$prop_ol[cell_idx] <- 1


grid_list <- list(near_zero_grid, dispersed_sparse_grid, clumped_sparse_grid, dispersed_dense_grid, clumped_dense_grid, near_full_grid)
names(grid_list) <- c("near_empty", "sparse_disp", "sparse_clump", "dense_disp", "dense_clump", "near_full")

grid_list <- list(near_zero_grid_rect, dispersed_sparse_grid_rect, clumped_sparse_grid_rect, dispersed_dense_grid_rect, clumped_dense_grid_rect, near_full_grid_rect)
names(grid_list) <- c("near_empty", "sparse_disp", "sparse_clump", "dense_disp", "dense_clump", "near_full")


grid <- grid_list[[1]]

n_cell_sides <- nrow(st_coordinates(grid[1, ])) - 1

grid_area <- st_area(grid[1, ])
cell_dist <- ifelse(n_cell_sides == 6, sqrt(2 * grid_area / (3 * sqrt(3))) * sqrt(3), sqrt(grid_area))
search_dist <- cell_dist * 1.25
buffers <- st_buffer(st_centroid(grid$geometry), search_dist)

points <- st_centroid(grid$geometry)
points <- st_sf(data.frame("id" = 1:length(points)), geometry = points)

idx <- grid$id

# Create master neighborhood list
nb_list <- lapply(idx, function(p) {
  ids <- points$id[st_intersects(points, buffers[p], sparse = FALSE, prepared = FALSE)]
  ids <- setdiff(ids, p)
})

# total_perimeter_edges <- ifelse(n_cell_sides == 6,
#                                 sum(sapply(idx, function(p) 6 - length(nb_list[[p]]))),
#                                 sum(sapply(idx, function(p) 4 - length(nb_list[[p]]))))

scores <- sapply(grid_list, function(grid) {
  
  idx_1 <- which(grid$prop_ol == 1)
  idx_0 <- which(grid$prop_ol == 0)
  
  # cells_class <- tibble::tibble(class = c("ones", "zeros"),
  #                               value = c(length(idx_1), length(idx_0)))
  # cells_class$n <- trunc(sqrt(cells_class$value))
  # cells_class$m <- cells_class$value - cells_class$n ^ 2
  # cells_class$min_e <- ifelse(test = cells_class$m == 0,
  #                             yes = cells_class$n * n_cell_sides,
  #                             no = ifelse(test = cells_class$n ^ 2 < cells_class$value & cells_class$value <= cells_class$n * (1 + cells_class$n),
  #                                         yes = n_cell_sides * cells_class$n + 2,
  #                                         no = ifelse(test = cells_class$value > cells_class$n * (1 + cells_class$n),
  #                                                     yes = n_cell_sides * cells_class$n + 4,
  #                                                     no = NA)))
  
  min_e <- sqrt(48 * length(idx_1) - 12)
  
  proportion_area_occupied_1 <- sum(grid$prop_ol) / length(grid$prop_ol)
  proportion_area_occupied_0 <- 1 - proportion_area_occupied_1
  
  total_edges <- length(idx) * n_cell_sides
  
  n_neighbours <- sapply(idx, function(p) length(nb_list[[p]]))
  
  boundary_cells <- n_cell_sides - n_neighbours
  
  adjacency_values <- sapply(idx, function(p) grid$prop_ol[nb_list[[p]]])
  
  # adjacency_matrix <- matrix(NA, nrow = length(grid$prop_ol), ncol = length(grid$prop_ol))
  
  n_adjacencies <- sapply(idx, function(p) sum(grid$prop_ol[nb_list[[p]]]))
  n_like_adjacencies <- sapply(idx_1, function(p) sum(adjacency_values[[p]] == TRUE) + boundary_cells[p])
  n_nonlike_adjacencies <- c(sapply(idx_1, function(p) sum(adjacency_values[[p]] == FALSE)),
                             sapply(idx_0, function(p) sum(adjacency_values[[p]] == TRUE))) 
  
  adjacency_ratio <- n_like_adjacencies / (sum(n_nonlike_adjacencies) - min_e)
  
  adjacency_ratio <- pmin(sum(adjacency_ratio), 1)
  
  conn <- mean(n_adjacencies / n_neighbours)
  
  clump <- ifelse(adjacency_ratio >= proportion_area_occupied_1 | (adjacency_ratio < proportion_area_occupied_1 & proportion_area_occupied_1 >= 0.5),
                  (adjacency_ratio - proportion_area_occupied_1) / (1 - proportion_area_occupied_1),
                  (adjacency_ratio - proportion_area_occupied_1) / (proportion_area_occupied_1))
  
  # Pi <- proportion_area_occupied
  # Pk <- 1 - proportion_area_occupied
  # 
  # Gik <- sapply(idx_1, function(p) sum(adjacency_values[[p]] == FALSE))
  # Gki <- sapply(idx_0, function(p) sum(adjacency_values[[p]] == TRUE))
  # 
  # Gik_ratio <- Gik / (sum(Gik) + sum(Gki))
  # Gki_ratio <- Gki / (sum(Gki) + sum(Gki))
  # 
  # test <- sum((Pi * Gik_ratio) * log(Pi * Gik_ratio), na.rm = TRUE)
  # 
  # IK <- sum((Pi * Gik_ratio) * log(Pi * Gik_ratio), na.rm = TRUE)
  # KI <- sum((Pk * Gki_ratio) * log(Pk * Gki_ratio), na.rm = TRUE)
  # 
  # grand_sum <- IK + KI
  # 
  # cont <- (1 + (grand_sum / (2 * log(2)))) * 100
  
  c(conn, clump)
})

data.frame("connectivity" = scores[1 , ], "clump" = scores[2 , ])

png(paste0("working/grid_scores_", n_cell_sides, "sided.png"), width = 1200, height = 800, res = 150, pointsize = 10)
par(mfrow = c(2, 3))
for(i in 1:length(grid_list)) {
  plot(grid_list[[i]][ , "prop_ol"], main = paste0("clump: ", scores[2 , i]), key.pos = NULL, reset = FALSE, pal = c("ivory3", "tomato3"))
}
dev.off()


#######
hex_sides_table <- data.frame("x" = c(1:20), "y" = c(6, 10, 12, 14, 16, 18, 18, 20, 22, 24, 26, 24, 26, 26, 28, 28, 30, 30, 30, 32))
fit <- lm(y ~ sqrt(48 * x - 12), hex_sides_table)
summary(fit)



a <- seq(.01, .99, .1)
b <- seq(.01, .99, .1)
all <- expand.grid(a, b)
test <- ifelse(all[, 1] >= all[, 2] | (all[, 1] < all[, 2] & all[, 2] >= 0.5),
                (all[, 1] - all[, 2]) / (1 - all[, 2]),
                (all[, 1] - all[, 2]) / (all[, 2]))
