library(gdistance)
library(RColorBrewer)
library(JuliaCall)
julia_setup(installJulia = TRUE)

rex <- raster(matrix(1,4,4))
a <- rep(c(1.3333), times=5)
b <- c(-1.3333, -0.6666, 0, 0.6666, 1.3333)

x1 <- c(-a, b)
x2 <- c(a, b)
y1 <- c(b, -a)
y2 <- c(b, a)
x <- cbind(x1,x2)
y <- cbind(y1,y2)


par(mfrow=c(1,3), mar= c(2,2,2,2), oma = c(0,0,0,0) + 0.1, cex.main=1)

x4 <- transition(rex, mean, 4)
g4 <- graph.adjacency(transitionMatrix(x4), mode="undirected")
gridLayout <- xyFromCell(x4, 1:ncell(x4))
plot(g4,layout=gridLayout, edge.color="black", vertex.color="black", vertex.label=NA, main="4 neighbours")
for(i in 1:dim(x)[1]){lines(x[i,],y[i,], col="lightgray")}
plot(g4, layout=gridLayout, add=TRUE, edge.color="black", vertex.color="black", vertex.label=NA)

x8 <- transition(rex, mean, 8)
g8 <- graph.adjacency(transitionMatrix(x8), mode="undirected")
plot(g8,layout=gridLayout, edge.color="black", vertex.color="black", vertex.label=NA, main="8 neighbours")
for(i in 1:dim(x)[1]){lines(x[i,],y[i,], col="lightgray")}
plot(g8, layout=gridLayout, add=TRUE, edge.color="black", vertex.color="black", vertex.label=NA)

x16 <- transition(rex, mean, 16)
g16 <- graph.adjacency(transitionMatrix(x16), mode="undirected")
plot(g16, layout=gridLayout, edge.color="black", vertex.color="black", vertex.label=NA, main="16 neighbours")
for(i in 1:dim(x)[1]){lines(x[i,],y[i,], col="lightgray")}
plot(g16,layout=gridLayout, add=TRUE, edge.color="black", vertex.color="black", vertex.label=NA)


sample_grid <- year_10_100

n_missing <- sum(is.na(sample_grid$resist))

sample_grid_noNA <- sample_grid[!is.na(sample_grid$resist), ]

sample_points <- st_centroid(sample_grid_noNA)

ext <- st_bbox(sample_points)
n_points <- nrow(sample_points)
idx <- sample_points$id

search_dist <- st_distance(sample_points$geometry[1], sample_points$geometry[2]) * 1.5

sample_points$resist <- sample_points$resist

nodes <- sapply(idx, function(p) {
  ids <- st_intersection(sample_points, st_buffer(sample_points$geometry[sample_points$id == p], search_dist))$id
  ids <- setdiff(ids, p)
  base_con <- sample_points$resist[sample_points$id == p]
  means <- sapply(ids, function(i) mean(c(base_con, sample_points$resist[which(sample_points$id == i)])))
  cbind(p, ids, means)
})

node_list <- as.matrix(do.call(rbind, nodes))

node_list <- node_list[!duplicated(apply(node_list[ , ], 1, function(row) paste(sort(row), collapse = ""))),]

write.table(node_list, file = "data/circuitscape/network.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

# write.table(sample_points$id, file = "data/circuitscape/focal_nodes.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)


ul_id <- sample_points$id[st_nearest_feature(st_point(c(ext[1], ext[4])), sample_points, check_crs = FALSE)]
ll_id <- sample_points$id[st_nearest_feature(st_point(c(ext[1], ext[2])), sample_points, check_crs = FALSE)]
ur_id <- sample_points$id[st_nearest_feature(st_point(c(ext[3], ext[4])), sample_points, check_crs = FALSE)]
lr_id <- sample_points$id[st_nearest_feature(st_point(c(ext[3], ext[2])), sample_points, check_crs = FALSE)]
write.table(c(ul_id, ll_id, ur_id, lr_id), file = "data/circuitscape/focal_nodes.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

sample_grid$current <- NA
sum_currents <- read.delim("data/circuitscape/connectivity_node_currents_cum.txt", header = FALSE)[idx , 2]
sample_grid$current[idx] <- zero_to_one(sum_currents)

# focal_nodes <- expand.grid(c(ul_id, ll_id, ur_id, lr_id), c(ul_id, ll_id, ur_id, lr_id))
# stationary <- which(focal_nodes$Var1 == focal_nodes$Var2)
# focal_nodes <- focal_nodes[-stationary, ]

# write.table("mode include", file = "data/circuitscape/focal_nodes.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
# write.table(focal_nodes, file = "data/circuitscape/focal_nodes.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
# 
# filenames <- list.files("data/circuitscape/", pattern = "node_currents")

currents <- lapply(filenames, function(f) read.delim(paste0("data/circuitscape/", f), header = FALSE)[ , 2])
sum_currents <- apply(do.call(cbind, currents), 1, sum)
sample_grid$current <- NA
sample_grid$current[idx] <- sum_currents / max(sum_currents)
