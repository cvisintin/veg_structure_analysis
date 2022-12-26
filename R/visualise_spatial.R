library(sf)
library(ggplot2)
library(gridExtra)
source('R/utility_functions.R')

# Load grids file
load(file = "output/years_1_30")

# Extract unique years from file list
years <- names(years_1_30)
n_year_groups <- length(years)

# Extract unique heights from file list and convert to numeric and meter units
heights <- names(years_1_30[[1]])
heights <- as.numeric(substr(heights, (nchar(heights) + 1) - 4, nchar(heights))) / 100

# # Generate lists of files grouped by years
# sp_files <- lapply(years, function(y) grep(y, list.files('output/sp_grids', pattern = '.shp', full.names = TRUE), value = TRUE))

# Create spatial outputs
sapply(seq_len(n_year_groups), function(j) {

  spatial_year <- years_1_30[[j]]
  
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
      annotate("text", label = paste0(heights[i], " m"), x = x , y = y, hjust = 0, color = 'gray40', size = 2) +
      scale_fill_viridis_c(option = 'E', na.value = "transparent") +
      theme_void()
  } 
  
  # Write out PNG of overlap grids for all heights
  png(paste0("figs/veg_structure_", years[j], ".png"), height = 150 * length(spatial_year), width = 400, pointsize = 12, res = 300)
  grid.arrange(grobs = rev(plot_list), ncol = 1)
  dev.off()
  
  NULL # Because function is not meant to return anything
})

# Animate spatial outputs
system('convert -delay 50 -dispose previous -loop 0 figs/* figs/animations/veg_growth.gif')
