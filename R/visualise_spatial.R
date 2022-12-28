library(sf)
library(ggplot2)
library(ggimage)
library(grid)
library(gridExtra)
library(png)
library(dplyr)
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

### Phenology ###
### From https://r-graph-gallery.com/295-basic-circular-barplot.html

months <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")

data <- data.frame(
  id = seq(1, 12, 1),
  value = sapply(months, function(month) sum(grepl(month, plant_points$phenology)))
)

img <- readPNG("data/images/flower.png")
g <- rasterGrob(img, interpolate = TRUE)

# label_data <- data
# 
# # calculate the ANGLE of the labels
# number_of_bar <- nrow(label_data)
# angle <-  75 - 360 * (label_data$id - 0.5) / number_of_bar     # I subtract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
# 
# # calculate the alignment of labels: right or left
# # If I am on the left part of the plot, my labels have currently an angle < -90
# label_data$hjust <- ifelse(angle < -90, 1, 0)
# 
# # flip angle BY to make them readable
# label_data$angle <- ifelse(angle < -90, angle + 180, angle)

# Make the plot
p <- ggplot(data, aes(x = as.factor(id), y = value)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  # This add the bars with a color
  geom_bar(stat = "identity", fill = alpha("red", 0.3)) +
  
  # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
  ylim(-max(data$value) * 0.5, max(data$value)) +
  
  # Custom the theme: no axis title and no cartesian grid
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1, 4), "cm")     # This remove unnecessary margin around plot
  ) +
  
  # This makes the coordinate polar instead of cartesian.
  coord_polar(start = 0.25) +
  
  # Add the labels, using the label_data dataframe that we have created before
  geom_text(data = data, aes(x = id, y = value * 0.5, label = toupper(row.names(data))),#, hjust = hjust),
            color = "black", fontface = "bold", alpha = 0.6, size = 2.5,
  ) +#angle = label_data$angle, inherit.aes = FALSE)
  
  geom_image(data = data.frame(xx = min(data$value), yy = -max(data$value) * 0.5, image = "data/images/flower_icon.png"), mapping = aes(xx, yy, image = image), size = .15, inherit.aes = FALSE)

p


### Species richness ###
data <- data.frame(
  id = seq(1, length(unique(plant_points$species)), 1),
  value = sapply(unique(plant_points$species), function(sp) sum(which(plant_points$species == sp)))
)

# label_data <- data
# 
# # calculate the ANGLE of the labels
# number_of_bar <- nrow(label_data)
# angle <-  90 - 360 * (label_data$id - 0.5) / number_of_bar     # I subtract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
# 
# # calculate the alignment of labels: right or left
# # If I am on the left part of the plot, my labels have currently an angle < -90
# label_data$hjust <- ifelse(angle < -90, 1, 0)
# 
# # flip angle BY to make them readable
# label_data$angle <- ifelse(angle < -90, angle + 180, angle)

# Make the plot
p <- ggplot(data, aes(x = as.factor(id), y = value)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  # This add the bars with a color
  geom_bar(stat = "identity", fill = alpha("green", 0.3)) +
  
  # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
  ylim(-max(data$value) * 0.5, max(data$value) * 1.05) +
  
  # Custom the theme: no axis title and no cartesian grid
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1, 4), "cm")     # This remove unnecessary margin around plot
  ) +
  
  # This makes the coordinate polar instead of cartesian.
  coord_polar(start = 0.25) +
  
  # Add the label
  geom_text(data = data.frame(xx = 0, yy = -max(data$value) * 0.5, label = shannon_evenness(data$value)), mapping = aes(xx, yy, label = label), size = 10, inherit.aes = FALSE)

# geom_text(data = label_data,
#           mapping = aes(x = id, y = value * 1.05, label = toupper(row.names(data)), hjust = hjust),
#           color = "black", fontface = "bold", alpha = 0.6, size = 2.5, angle = label_data$angle, inherit.aes = FALSE)

p


### Coverage at 1m ###
data <- data.frame(
  id = seq(1, n_year_groups, 1),
  value = extract_coverage(years_1_30)
)

# Make the plot
p <- ggplot(data, aes(x = as.factor(id), y = value)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  # This add the bars with a color
  geom_bar(stat = "identity", fill = alpha("green", 0.3)) +
  
  # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
  ylim(-max(data$value) * 0.5, max(data$value) * 1.05) +
  
  # Custom the theme: no axis title and no cartesian grid
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1, 4), "cm")     # This remove unnecessary margin around plot
  ) +
  
  # This makes the coordinate polar instead of cartesian.
  coord_polar(start = -(2 / nrow(data))) +
  
  # Add the labels, using the label_data dataframe that we have created before
  geom_text(data = data, aes(x = id, y = value * 0.5, label = paste0("Y", toupper(row.names(data)))),#, hjust = hjust),
            color = "black", fontface = "bold", alpha = 0.6, size = 2.5,
  ) +#angle = label_data$angle, inherit.aes = FALSE)
  
  geom_image(data = data.frame(xx = min(data$value), yy = -max(data$value) * 0.5, image = "data/images/shrub_icon.png"), mapping = aes(xx, yy, image = image), size = .17, inherit.aes = FALSE)

p


### Density ###


data <- data.frame(
  category = unique(plant_points$density),
  values = sapply(unique(plant_points$density), function(dens) sum(which(plant_points$density == dens))),
  palette = NA
)

# Compute percentages
data$fraction = data$values / sum(data$values)

# Compute the cumulative percentages (top of each rectangle)
data$ymax = cumsum(data$fraction)

# Compute the bottom of each rectangle
data$ymin = c(0, head(data$ymax, n = -1))

# Compute label position
data$labelPosition <- (data$ymax + data$ymin) / 2

# Compute a good label
data$label <- paste0(data$category, "\n value: ", data$values)


gr_palette <- c("#25603b", "#b7e4c8", "#62a67c")

# Make the plot
ggplot(data, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = category)) +
  geom_rect() +
  geom_text(x = 3.5, aes(y = labelPosition, label = toupper(category)), size = 4) +
  scale_fill_manual(values = gr_palette) +
  coord_polar(theta = "y") +
  xlim(c(0.8, 4)) +
  theme_void() +
  theme(legend.position = "none") +
  
  geom_image(data = data.frame(xx = 0.8, yy = 0, image = "data/images/leaves_icon.png"), mapping = aes(xx, yy, image = image), size = .3, inherit.aes = FALSE)



#### From https://pomvlad.blog/2018/05/03/gauges-ggplot2/
df <- data.frame(matrix(nrow = 5, ncol = 2))

names(df) <- c("variable", "percentage")
df$variable <- c("Carbohydrates", "Warming", "NGTnotPresent", "DrainNotPresent", "DrEaMing")
df$percentage <- c(0.67, 0.33, 0.86, 0.78, 0.58)

df <- mutate(df, group = ifelse(percentage < 0.6, "red",
                                ifelse(percentage >= 0.6 & percentage < 0.8, "orange", "green")),
             label = paste0(percentage * 100, "%"),
             title = dplyr::recode(variable, 'Carbohydrates' = "Preoperative\ncarbohydrate loading",
                                   'Warming' = "Intraoperative\nwarming",
                                   'NGTnotPresent' = "Patients without a\nnasogastric tube\non arrival in recovery",
                                   'DrainNotPresent' = "Patients without an\nabdominal drain\non arrival in recovery",
                                   'DrEaMing' = "Patients DrEaMing on\npostoperative day 1"))


ggplot(df, aes(fill = group, ymax = percentage, ymin = 0, xmax = 2, xmin = 1)) +
  geom_rect(aes(ymax = 1, ymin = 0, xmax = 2, xmin = 1), fill ="#ece8bd") +
  geom_rect() + 
  coord_polar(theta = "y", start = -pi / 2) + xlim(c(0, 2)) + ylim(c(0, 2)) +
  geom_text(aes(x = 0, y = 0, label = label, colour = group), size = 6.5, family = "Poppins SemiBold") +
  geom_text(aes(x = 1.5, y = 1.5, label = title), family = "Poppins Light", size = 4.2) + 
  facet_wrap(~title, ncol = 5) +
  theme_void() +
  scale_fill_manual(values = c("red" = "#C9146C", "orange" = "#DA9112", "green" = "#129188")) +
  scale_colour_manual(values = c("red" = "#C9146C", "orange" = "#DA9112", "green" = "#129188")) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) +
  guides(fill = "none") +
  guides(colour = "none")




# 