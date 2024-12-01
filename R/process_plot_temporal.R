library(gridExtra)
library(png)
source("R/utils.R")

plant_points <- read_sf("data/gis/sample_temporal/plant_final_gis.shp")
boundary <- read_sf("data/gis/sample_temporal/boundary.shp")

years_1_30 <- process_structure(point_locations = plant_points,
                                boundary = boundary,
                                years = 1:30,
                                cellarea = 5)

# save(years_1_30, file = "output/sample_temporal/years_1_30")
# create_images(years_1_30, path = "figs/images/sample_temporal/", filename = "veg_structure")

years <- names(years_1_30)[c(1, 3, 5, 10)]
n_year_groups <- length(years)

# Extract unique heights from file list and convert to numeric and meter units
cut_heights <- names(years_1_30[[1]])
cut_heights <- as.numeric(substr(cut_heights, (nchar(cut_heights) + 1) - 4, nchar(cut_heights))) / 100

y_plots <- lapply(seq_len(n_year_groups), function(j) {
  
  spatial_year <- years_1_30[[j]]
  
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
      geom_sf(data = spatial_data, aes(fill = prop_ol), color = 'gray60', lwd = 0.08, show.legend = FALSE) +
      # annotate("text",
      #          label = paste0(cut_heights[i], " m"),
      #          x = x ,
      #          y = y,
      #          hjust = "inward",
      #          color = 'gray20',
      #          size = 2,
      #          fontface = "bold") +
      # {if(i == length(spatial_year)) annotate("text",
      #                                         label = paste0("Year ", j),
      #                                         x = x2 ,
      #                                         y = y2,
      #                                         hjust = "inward",
      #                                         color = 'gray20',
      #                                         size = 2,
      #                                         fontface = "bold")} +
      scale_fill_viridis_c(option = 'E', na.value = "transparent") +
      theme_void()
    
  }
  grid.arrange(grobs = rev(plot_list), ncol = 1)
  
})


m_plot <- readPNG('data/gis/sample_temporal/veg_plan_gis.png')
# t1_plot <- readPNG('figs/images/sample_temporal/veg_structure_year_01.png')
# t2_plot <- readPNG('figs/images/sample_temporal/veg_structure_year_03.png')
# t3_plot <- readPNG('figs/images/sample_temporal/veg_structure_year_05.png')
# t4_plot <- readPNG('figs/images/sample_temporal/veg_structure_year_10.png')

g_plot <- arrangeGrob(rasterGrob(m_plot),
                      y_plots[[1]],
                      y_plots[[2]],
                      y_plots[[3]],
                      y_plots[[4]],
                      ncol = 4,
                      nrow = 3,
                      layout_matrix = rbind(c(1,1,1,1), c(2,3,4,5))
                      )

ggsave('figs/temporal_change.png',
       g_plot,
       height = 3000,
       width = 3000,
       device = 'png',
       dpi = 600,
       units = "px")
