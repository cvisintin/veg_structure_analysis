library(sf)
library(ggplot2)
library(gridExtra)
source('R/utility_functions.R')

# Extract unique years from file list
years <- unique(substr(list.files('output/sp_grids', pattern = '.shp'), 1, 7))
n_year_groups <- length(years)

# Extract unique heights from file list and convert to numeric and meter units
heights <- as.numeric(unique(gsub(".*height_(.+)_grd.*", "\\1", list.files('output/sp_grids', pattern = '.shp')))) / 100

# Generate lists of files grouped by years
sp_files <- lapply(years, function(y) grep(y, list.files('output/sp_grids', pattern = '.shp', full.names = TRUE), value = TRUE))

# Set parameters for spatial outputs
x <- 10636525
y <- 5583415
color <- 'gray40'

# Create spatial outputs
sapply(seq_len(n_year_groups), function(f) {

  target_files <- sp_files[[f]]
  
  plot_list <- list()
  
  for (i in 1:length(target_files)) {
    
    # Load processed spatial data
    spatial_data <- read_sf(target_files[i])
    
    # Rotate spatial data
    spatial_data <- rotate_data(spatial_data)
    
    # Create plot of proportional overlap grid at currently specified height
    plot_list[[i]] <- ggplot() +
      geom_sf(data = spatial_data, aes(fill = prop_ol), color = color, size = 0.05, show.legend = FALSE) +
      annotate("text", label = paste0(heights[i], " m"), x = x , y = y, hjust = 0, color = color, size = 2) +
      scale_fill_viridis_c(option = 'E', na.value = "transparent") +
      theme_void()
  } 
  
  # Write out PNG of overlap grids for all heights
  png(paste0("figs/veg_structure_", years[f], ".png"), height = 1200, width = 400, pointsize = 16, res = 300)
  grid.arrange(grobs = rev(plot_list), ncol = 1)
  dev.off()
  
  NULL # Because function is not meant to return anything
})
