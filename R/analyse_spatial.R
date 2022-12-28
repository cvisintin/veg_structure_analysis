library(sf)
source('R/utility_functions.R')

# Load background spatial data
building <- read_sf("data/gis/Building_Mass.shp")

# Extract unique years from file list
years <- unique(substr(list.files('output/sp_grids', pattern = '.shp'), 1, 7))
n_year_groups <- length(years)

# Extract unique cut heights from file list and convert to numeric and meter units
cut_heights <- as.numeric(unique(gsub(".*height_(.+)_grd.*", "\\1", list.files('output/sp_grids', pattern = '.shp')))) / 100

# Drop 1m cut if it is not part of the regular sequence


# Generate lists of files grouped by years
sp_files <- lapply(years, function(y) grep(y, list.files('output/sp_grids', pattern = '.shp', full.names = TRUE), value = TRUE))

# Read in files and convert data for connectivity analysis (note object names use centimeter units)
for(f in seq_len(n_year_groups)) {
  
  target_files <- sp_files[[f]]
  
  grid_shape <- gsub(".*grd_(.+).shp*", "\\1", target_files[1])
  
  for (i in 1:length(target_files)) {
    assign(paste0(years[f], "_", cut_heights[i] * 100), prepare_conn_data(grid = read_sf(target_files[i]),
                                                                obstructions = building,
                                                                cut_height = cut_heights[i],
                                                                grid_shape = grid_shape))
  }
}

# Run connectivity analysis on example grid (Year 10 @ 1m (100cm) height) - IN PROGRESS
# sample_connectivity <- run_connectivity(resistance_grid = year_10_100)
