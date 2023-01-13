library(sf)
library(ggplot2)
library(ggimage)
library(grid)
library(gridExtra)
library(png)
library(dplyr)
library(tidyr)
library(viridis)
source('R/utils.R')

# Load background spatial data
plant_points <- read_sf("data/gis/Plant_Data.shp")
planting_polygons <- read_sf("data/gis/Planting_Data.shp")
building <- read_sf("data/gis/Building_Mass.shp")
site_boundary <- read_sf("data/gis/Site_Boundary.shp")

# Convert planting polygons to points and combine with existing plant points
plant_points <- convert_combine(point_locations = plant_points,
                                polygon_locations = planting_polygons)

# Load grids file
load(file = "output/years_1_30")

# Create spatial outputs
create_images(years_1_30, path = "figs/")

# Animate spatial outputs
system('convert -delay 50 -dispose previous -loop 0 figs/* figs/animations/veg_growth.gif')

### Density ###
density_plot <- plot_classes(plant_points,
                             variable_name = "density",
                             colour_palette = c("#397d53", "#b7e4c8", "#62a67c"),
                             image_path = "data/images/leaves_icon.png")

### Texture ###
texture_plot <- plot_classes(plant_points,
                             variable_name = "texture",
                             colour_palette = c("#4b6c90", "#afc6e0", "#6d8eb3"),
                             image_path = "data/images/texture_icon.png")

### Size ###
size_plot <- plot_classes(plant_points,
                          variable_name = "size",
                          colour_palette = c("#b8641d", "#ebc5a4", "#d9a171", "#c98449"),
                          image_path = "data/images/vegetation_icon.png")

### Endemism ###
endemism_plot <- plot_percent(plant_points,
                              variable_name = "endemism",
                              colour = "#a072a6",
                              label = "NATIVE")

### Type ###
type_plot <- plot_percent(plant_points,
                          variable_name = "type",
                          colour = "#a19a65",
                          label = "EVERGREEN")

### Species richness ###
richness_data <- data.frame(
  id = seq(1, length(unique(spatial_points$species)), 1),
  value = sapply(unique(spatial_points$species), function(sp) length(which(spatial_points$species == sp)))
)

# Make the plot
richness <- ggplot(richness_data, aes(x = as.factor(id), y = value)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  # This add the bars with a color
  geom_bar(stat = "identity", fill = "#858383") +
  
  # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
  ylim(-max(richness_data$value), max(richness_data$value) * 1.05) +
  
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
  geom_text(data = data.frame(xx = min(richness_data$value), yy = -max(richness_data$value), label = paste0(shannon_evenness(richness_data$value), "\nEVENNESS\nSCORE")), mapping = aes(xx, yy, label = label), size = 4, inherit.aes = FALSE)

# geom_text(data = label_data,
#           mapping = aes(x = id, y = value * 1.05, label = toupper(row.names(data)), hjust = hjust),
#           color = "black", fontface = "bold", alpha = 0.6, size = 2.5, angle = label_data$angle, inherit.aes = FALSE)

richness


### Phenology ###
### From https://r-graph-gallery.com/295-basic-circular-barplot.html

months <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")
n_months <- length(months)

phenology_data <- data.frame(
  id = seq_len(n_months),
  value = sapply(months, function(month) sum(grepl(month, spatial_points$phenology)))
)

img <- readPNG("data/images/flower_icon.png")
g <- rasterGrob(img, interpolate = TRUE)

# Make the plot
phenology <- ggplot(phenology_data, aes(x = as.factor(id), y = value)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  # This add the bars with a color
  geom_bar(stat = "identity", fill = "#b0d4d6") +
  
  # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
  ylim(-max(phenology_data$value), max(phenology_data$value)) +
  
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
  geom_text(data = phenology_data, aes(x = id, y = value * 0.5, label = toupper(row.names(phenology_data))),#, hjust = hjust),
            color = "black", fontface = "bold", alpha = 0.6, size = 2.5,
  ) +#angle = label_data$angle, inherit.aes = FALSE)
  
  geom_image(data = data.frame(xx = min(phenology_data$value), yy = -max(phenology_data$value), image = "data/images/flower_icon.png"), mapping = aes(xx, yy, image = image), size = .25, inherit.aes = FALSE)

phenology


### Coverage at 1m ###
coverage_data <- data.frame(
  id = factor(seq_len(n_year_groups)),
  value = extract_1m_coverage(spatial_list)
)

# Make the plot
coverage <- ggplot(coverage_data, aes(x = id, y = value)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  # This add the bars with a color
  geom_bar(stat = "identity", fill = "#c9837d") +
  
  # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
  ylim(-max(coverage_data$value) * 1.1, max(coverage_data$value) * 1.15) +
  
  # Custom the theme: no axis title and no cartesian grid
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1, 4), "cm")     # This remove unnecessary margin around plot
  ) +
  
  # This makes the coordinate polar instead of cartesian.
  coord_polar(start = -(2 / nrow(coverage_data))) +
  
  # Add the labels, using the label_data dataframe that we have created before
  geom_text(data = coverage_data, aes(x = id, y = value * 0.5, label = paste0("Y", toupper(row.names(coverage_data)))),#, hjust = hjust),
            color = "black", fontface = "bold", alpha = 0.6, size = 2.5,
  ) +#angle = label_data$angle, inherit.aes = FALSE)
  
  geom_image(data = data.frame(xx = 1, yy = -max(coverage_data$value) * 1.1, image = "data/images/shrub_icon.png"), mapping = aes(xx, yy, image = image), size = .25, inherit.aes = FALSE)

coverage


### Connectivity ###
connectivity_data_w <- estimate_connectivity(spatial_list)

connectivity_score <- base::round(mean(as.matrix(connectivity_data_w)), 2)

connectivity_data_w$remainder <- apply(connectivity_data_w, 1, function(x) ncol(connectivity_data_w) - sum(x))

connectivity_data <- gather(connectivity_data_w, key = "group", value = "value")

connectivity_data$id <- factor(rep(1:n_year_groups, ncol(connectivity_data_w)))
connectivity_data$group <- factor(connectivity_data$group, levels = rev(unique(connectivity_data$group)))

max_value <- max(connectivity_data$value)

conn_colours <- c(colorRampPalette(c("#9fbfa3", "#3e6343"))(ncol(connectivity_data_w) - 1), "#e1e1e1")
conn_colours <- c(viridis(ncol(connectivity_data_w) - 1, end = 0.8), "#e1e1e1")

# Make the plot
connectivity <- ggplot(connectivity_data, aes(x = id, y = value, fill = group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  # This add the bars with a color
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = rev(conn_colours)) +
  
  # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
  ylim(-max(connectivity_data$value) * 1.1, max(connectivity_data$value) * 1.3) +
  
  # Custom the theme: no axis title and no cartesian grid
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1, 4), "cm")     # This remove unnecessary margin around plot
  ) +
  
  # This makes the coordinate polar instead of cartesian.
  coord_polar(start = -(2 / nrow(connectivity_data_w))) +
  
  # Add the labels, using the label_data dataframe that we have created before
  geom_text(data = connectivity_data[1:nrow(connectivity_data_w), ], aes(x = unique(id), y = max_value, label = paste0("Y", row.names(connectivity_data_w))),#, hjust = hjust),
            color = "black", fontface = "bold", alpha = 0.6, size = 2.5,
  ) +#angle = label_data$angle, inherit.aes = FALSE)
  
  geom_text(data = data.frame(xx = max(connectivity_data$value), yy = -max(connectivity_data$value) * 1.1, label = paste0(connectivity_score, "\nCONNECTIVITY\nSCORE")), mapping = aes(xx, yy, label = label), size = 4, inherit.aes = FALSE)

connectivity


png(paste0("figs/results.png"), height = 2400, width = 2400, pointsize = 4, res = 300)
grid.arrange(grobs = list(density,
                          texture,
                          sizes,
                          endemism,
                          richness,
                          type,
                          phenology,
                          connectivity,
                          coverage
), ncol = 3)
dev.off()










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