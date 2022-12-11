library(sf)
library(dplyr)

test_boundary <- read_sf("data/gis/Test_Boundary.shp")

# Define cell area in square meters
cellarea <- 10

# Calculate spatial properties of hexagon cell and regular "flat-topped" spacing
outer_radius <- sqrt(2 * cellarea / (3 * sqrt(3)))
inner_radius <- (outer_radius * sqrt(3)) / 2
horizontal_spacing <- 3 * outer_radius
vertical_spacing <- inner_radius

# Increase boundary to ensure adequate coverage
boundary <- st_buffer(test_boundary, 2 * horizontal_spacing)

# Determine extent and size of boundary
extent <- st_bbox(boundary)
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

# Identify indices of x coordinates to offset
shift_idx <- rep(c(rep(FALSE, n_columns), rep(TRUE, n_columns)), sum(seq_len(n_rows) %% 2 == 0))
if(n_rows %% 2 != 0) shift_idx <- c(shift_idx, c(rep(FALSE, n_columns)))

# Offset coordinates for target indices
xy_coords[shift_idx, 1] <- xy_coords[shift_idx, 1] + horizontal_spacing / 2

# Create lists of coordinates to bound hexagon cell for each centroid coordinate
poly_coords <- lapply(seq_len(n_columns * n_rows), function(i) {
  x_coord <- xy_coords[i , 1]
  y_coord <- xy_coords[i , 2]
    list(
      rbind(
        c(x_coord - (outer_radius / 2) , y_coord + inner_radius),
        c(x_coord + (outer_radius / 2) , y_coord + inner_radius),
        c(x_coord + outer_radius, y_coord),
        c(x_coord + (outer_radius / 2) , y_coord - inner_radius),
        c(x_coord - (outer_radius / 2) , y_coord - inner_radius),
        c(x_coord - outer_radius, y_coord),
        c(x_coord - (outer_radius / 2) , y_coord + inner_radius)
      ))
})

# Create spatial polygons
polys <- lapply(seq_len(n_columns * n_rows), function(i) st_sfc(st_polygon(poly_coords[[i]]), crs = st_crs(test_boundary)))
polys <- do.call(c, polys)
polys <- st_sf(data.frame("id" = 1:length(polys)), geometry = polys)

# Crop to input boundary
grid <- polys[st_intersection(polys, test_boundary)$id, ]

# Reset numbering in ids
grid$id <- seq_len(nrow(grid))