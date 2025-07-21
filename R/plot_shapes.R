library(gridExtra)
source("R/utils.R")

shape_types <- c("rounded", "mounding", "vase", "upright", "columnar", "pyramidal", "round")

plant_height <- 1

adj_cut_height = seq(0, 1, 0.01)

plot_list <- vector("list", length = 7)

plot_color <- "palegreen4"

plot_layout <- rbind(c(1,1, 2,2, 3,3, 4,4),
                     c(NA, 5,5, 6,6, 7,7, NA))

for (i in shape_types) {
  plant_width <- do.call(i, list(adj_cut_height = adj_cut_height, plant_height = plant_height))
  df <- data.frame("x" = c(adj_cut_height, plant_height), "y" = c(plant_width, 0))
  # plot(df, type = 'l')
  
  plot_list[[which(i == shape_types)]] <- ggplot() + 
    geom_area(data = df, aes(x = x, y = y), linetype = 0, colour = plot_color, fill = plot_color) +
    geom_area(data = df, aes(x = x, y = -y), linetype = 0, colour = plot_color, fill = plot_color) +
    geom_hline(yintercept = 0, colour = plot_color) +
    geom_vline(xintercept = 0, colour = "black") +
    scale_x_continuous(name = "", 
                       expand = c(0, 0)) +
    scale_y_continuous(name = i, 
                       expand = c(0, 0)) +
    coord_flip() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
}

png("figs/geom_shapes.png", height = 1200, width = 1800, pointsize = 8, res = 300)
do.call("grid.arrange", c(plot_list,
                          list(ncol = 4,
                          layout_matrix = plot_layout)))
dev.off()
