library(sf)
library(rgl)
options(rgl.printRglwidget = TRUE)
library(viridis)
library(JuliaCall)
julia_setup(installJulia = TRUE)
source('R/utility_functions.R')


grid <- create_grid(boundary = st_sfc(st_polygon(list(rbind(c(0,0),
                                               c(0, 10),
                                               c(10, 10),
                                               c(10, 0),
                                               c(0, 0)
                                               )))), hex_gridshape = FALSE, cellarea = 1)
grid$prop_ol <- 0

grid_dist <- st_distance(st_centroid(grid$geometry[1]), st_centroid(grid$geometry[2]))[1]

max_height <- 10

n_z_levels <- max_height %/% grid_dist

n_rows <- diff(st_bbox(grid)[c(2, 4)]) / grid_dist
n_columns <- diff(st_bbox(grid)[c(1, 3)]) / grid_dist

search_pattern <- expand.grid(seq(-1, 1, 1), seq(-1, 1, 1), seq(-1, 1, 1))
search_pattern <- search_pattern[which(apply(search_pattern, 1, function(x) all(x == 0)) == FALSE), ]

grid_array_prop_ol <- grid_array_ids <- array(1:prod(c(n_columns, n_rows, n_z_levels)), dim = c(n_columns, n_rows, n_z_levels))

for (i in seq_len(n_z_levels)) {
  grid_array_prop_ol[ , , i] <- runif(prod(c(n_columns, n_rows)))
}

set.seed(333)
na_coords <- cbind(sample(2:10, 8), sample(2:10, 8), sample(2:9, 8))
grid_array_prop_ol[na_coords] <- NA
na_idx <- which(is.na(grid_array_prop_ol))

coords <- do.call(rbind, lapply(1:prod(dim(grid_array_ids)), function(i) which(grid_array_ids == i, arr.ind = TRUE)))

nodes <- sapply(setdiff(1:nrow(coords), na_idx), function(id) {

  target_coords <- do.call(rbind, lapply(1:nrow(search_pattern), function(x) search_pattern[x, ] + c(coords[id, 1], coords[id, 2], coords[id, 3])))
  target_coords <- target_coords[target_coords[ , 1] %in% c(1:n_rows) & target_coords[ , 2] %in% c(1:n_columns) & target_coords[ , 3] %in% c(1:n_z_levels), ]

  target_coords_ids <- apply(target_coords, 1, function(t_id) {
    which(coords[ , 1] == t_id[1] & coords[ , 2] == t_id[2] & coords[ , 3] == t_id[3])
  })

  target_coords_ids <- setdiff(target_coords_ids, na_idx)
  
  means <- sapply(target_coords_ids, function(target_id) mean(c((1 - grid_array_prop_ol[id]), (1 - grid_array_prop_ol[target_id]))))
  
  cbind(id, target_coords_ids, means)
})

node_list <- as.matrix(do.call(rbind, nodes))

node_list <- node_list[!duplicated(apply(node_list[ , ], 1, function(row) paste(sort(row), collapse = ""))),]

# node_list[which(is.na(node_list))] <- "NODATA"

write.table(node_list, file = "data/circuitscape/network.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

# Identify grid array corners (8 in total)

array_extents <- expand.grid(c(1, n_rows),
                             c(1, n_columns),
                             c(1, n_z_levels))

focal_ids <- apply(array_extents, 1, function(t_id) {
  which(coords[ , 1] == t_id[1] & coords[ , 2] == t_id[2] & coords[ , 3] == t_id[3])
})

write.table(focal_ids, file = "data/circuitscape/focal_nodes.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)



sum_currents <- read.delim("data/circuitscape/connectivity_node_currents_cum.txt", header = FALSE)[-na_idx , 2]

grid_array_current <- grid_array_ids

grid_array_current[] <- NA
grid_array_current[-na_idx] <- zero_to_one(sum_currents)

color_idx <- sort(unique(round(zero_to_one(sum_currents) * 100)))


points3d(coords[-na_idx , 1], coords[-na_idx , 2], coords[-na_idx , 3])

# z_levels <- c(0, seq_len(n_z_levels)) * grid_dist
# 
# lapply(1:length(z_levels), function(level) {
# xyz <- t(sapply(seq_len(nrow(grid)), function(i) {
#   st_coordinates(st_zm(st_centroid(grid$geometry[i]), drop = FALSE, what = "Z"))
# }))
# 
# 
# })