library(ggplot2)

n <- 15
data <- data.frame(id = factor(seq(1, n)), value = round(runif(n, min = 500, max = 2000)))
colours = "#858383"
polar_rotation = 0.25

p <- ggplot(data, aes(x = as.factor(id), y = value))

p <- p + geom_bar(stat = "identity", fill = colours)

p <- p + ylim(-max(data$value) * 1.1, max(data$value) * 1.5)

p <- p + theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1, 4), "cm")
  )

p <- p + coord_polar(start = polar_rotation)

p <- p + geom_text(data = data.frame(xx = 0.5,
                                     yy = -max(data$value),
                                     label = paste0(0.98, "\nEVENNESS\nSCORE\n\n(", length(data$value), " SPECIES)")),
                   mapping = aes(xx, yy, label = label),
                   size = 3,
                   inherit.aes = FALSE)

p


#################################################################

data <- data.frame(id = factor(seq_len(30)),
                   value = sort(runif(30)),
                   max = 1)
colours = "#c9837d"

p <- ggplot(data, aes(x = as.factor(id), y = max))

p <- p + geom_bar(data, mapping = aes(x = as.factor(id), y = max), stat = "identity", fill = "#e1e1e1")

p <- p + geom_bar(data, mapping = aes(x = as.factor(id), y = value), stat = "identity", fill = colours)

p <- p + ylim(-1.1, 1.15)

p <- p + theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1, 4), "cm")
  )

p <- p + coord_polar(start = -(2 / nrow(data)))

p <- p + geom_text(data = data,
                   mapping = aes(x = id, y = 1, label = paste0("Y", toupper(row.names(data)))),
                   color = "black",
                   fontface = "bold",
                   alpha = 0.6,
                   size = 2.5) +
  
  geom_image(data = data.frame(xx = 1,
                               yy = -1.1,
                               image = "data/images/shrub_icon.png"),
             mapping = aes(xx, yy, image = image),
             size = .25,
             inherit.aes = FALSE)

p


#################################################################

months <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")
n_months <- length(months)

data <- data.frame(
  id = factor(seq_len(n_months)),
  value = runif(12)
)
row.names(data) <- months

img <- readPNG("data/images/flower_icon.png")
g <- rasterGrob(img, interpolate = TRUE)

colours = "#b0d4d6"
polar_rotation = 0.25

p <- ggplot(data, aes(x = as.factor(id), y = value))

p <- p + geom_bar(stat = "identity", fill = colours)

p <- p + ylim(-max(data$value) * 1.1, max(data$value) * 1.5)

p <- p + theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1, 4), "cm")
  )

p <- p + coord_polar(start = polar_rotation)

p <- p + geom_text(data = data,
                   mapping = aes(x = id, y = value * 0.5, label = toupper(row.names(data))),
                   color = "black",
                   fontface = "bold",
                   alpha = 0.6,
                   size = 2.5) +
  
  geom_image(data = data.frame(xx = 0.5,
                               yy = -max(data$value),
                               image = "data/images/flower_icon.png"),
             mapping = aes(xx, yy, image = image),
             size = .2,
             inherit.aes = FALSE)

p
