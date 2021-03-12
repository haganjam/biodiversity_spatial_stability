
# Performance enhancing effect of biodiversity

# load relevant libraries
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(here)

# set script to call partition functions from
source(here("R_scripts/functions_BEF_partitions.R"))

# read in the data
ply_dat_raw <- read_delim( here("data/plymouth_rockpool_data.csv"), delim = ";" )

# duplicate this data
ply_dat <- ply_dat_raw

# check the variable structure
str(ply_dat)

# aggregate the data because some rows have multiple measurements

# set up a variable names vector
var_names <- names(ply_dat)

# set up a vector of id_vars
id_vars <- c("month", "shore", "pool", "composition", "removed")

# set up a vector of response aggregate variables
agg_vars <- var_names[!(names(ply_dat) %in% id_vars)]

# aggregate pools with multiple measurements by summing them
ply_dat <- 
  ply_dat %>%
  group_by_at(id_vars) %>%
  summarise_at(vars(agg_vars), ~sum(., na.rm = TRUE)) %>%
  ungroup()

# matches with Lars and Bagchi's script results

# examine which 'species' classes have notable cover

# analyse only the focal species
foc_spp <- 
  c("sargassum", "bifurcaria", "f_serratus", "l_digitata")

# create a few extra variables for this analysis

# create a species richness variable
# select only relevant columns
tot_spp <- 4

ply_dat <- 
  ply_dat %>%
  mutate(species_richness = (tot_spp - removed) ) %>%
  select( all_of(id_vars), grid_squares, species_richness, totcover, cover_nocrusts, all_of(foc_spp) )

# create a new mixture/monoculture column
# reorder the columns again
ply_dat <- 
  ply_dat %>%
  mutate(mono_mix = if_else(removed == "3", "monoculture",
                            if_else(removed == "4", "control", 
                                    if_else(removed == "0", "max_mixture", "mixture"))))


# remove the composition = "None" treatment as this is a general control
ply_dat <- 
  ply_dat %>% filter(composition != "None")

View(ply_dat)

# a column for total cover of the four focal species
ply_dat$total_focal_cover <- 
  ply_dat %>%
  select(all_of(foc_spp)) %>%
  rowSums(.)

# reorder columns
names(ply_dat)

ply_dat <- 
  ply_dat %>%
  select(shore, time = month, place = pool, grid_squares, 
         composition, mono_mix, total_focal_cover, cover_nocrusts)

# add a shore time variable
d <- mapply(function(x, y) {as.character(paste(x, y, sep = "_")) }, ply_dat$shore, ply_dat$time )
names(d) <- NULL

ply_dat$shore_time <- d


# get the mixtures only
ply_dat$composition

ply_mix <- 
  ply_dat %>%
  filter(composition == "FLB") %>%
  select(shore_time, place, grid_squares, Y = total_focal_cover) 

ply_mix <- split(ply_mix, ply_mix$shore_time)

# get the monocultures only
ply_mono <- 
  ply_dat %>%
  filter(composition %in% c("F", "L", "B")) %>%
  select(shore_time, place, grid_squares, composition, M = total_focal_cover) 

ply_mono <- split(ply_mono, ply_mono$shore_time)


adf_out <- 
  mapply(function(u, v){
  
  v_list <- split(v, v$composition)
  
  g <- 
    lapply(v_list, function(data) {
      
      v_size <- data$grid_squares
      
      i_out  <- 
        sapply(u$grid_squares, function(grid) {
          
          h <- abs(grid - v_size)
          
          j <- which( near(h, min(h) ) )
          
          if(length(j) > 1) {
            k <- sample(j, size = 1)
          } else { k <- j }
          
          return(k)
          
        })
      
      df <- data.frame(row_id = 1:nrow(u),
                       species = data$composition[1],
                       M = data[i_out, ]$M
      )
      
      return(df)
      
    })
  
  g <- bind_rows(g, .id = NULL)
  
  u$row_id <- 1:nrow(u)
  
  joined_dat <- full_join(u, g, by = "row_id")
  
  return(joined_dat)
  
}, ply_mix, ply_mono, SIMPLIFY = FALSE)


# bind this together
adf_out <- bind_rows(adf_out, .id = NULL)

adf_out <- 
  adf_out %>%
  select(shore_time, place, species, M, Y) %>%
  separate(shore_time, c("shore", "time"), "_")

# first try it with Challaborough
shore1 <- 
  adf_out %>%
  filter(shore == "Challaborough") %>%
  # filter(shore == "Kingsand", place != "k12") %>%
  select(-shore)


# add a sample column
n.sp <- length(unique(shore1$species))
contexts <- n_unique(shore1$time)*n_unique(shore1$place)

shore1$sample <- rep(1:contexts, each = n.sp)

shore1 <- 
  shore1 %>%
  select(sample, time, place, species, M, Y)

shore1$M <- sapply(shore1$M, function(x) ifelse(x < 0.1, 0.1, x))
shore1$Y <- sapply(shore1$Y, function(x) ifelse(x < 0.1, 0.1, x))

# why such high values? There is extremely high complementarity
neb <- isbell.2018.pt(adf = shore1, RY.exp =  c(0.33, 0.33, 0.33))
neb

# combined shores
comb_shores <- 
  adf_out %>%
  group_by(shore, time, species) %>%
  summarise(M = mean(M),
            Y = mean(Y), .groups = "drop")

comb_shores <- 
  comb_shores %>%
  rename(place = shore)

# add a sample column
n.sp <- length(unique(comb_shores$species))
contexts <- n_unique(comb_shores$time)*n_unique(comb_shores$place)

comb_shores$sample <- rep(1:contexts, each = n.sp)
comb_shores

comb_shores$M <- sapply(comb_shores$M, function(x) ifelse(x < 0.1, 0.1, x))
comb_shores$Y <- sapply(comb_shores$Y, function(x) ifelse(x < 0.1, 0.1, x))

comb_shores

neb <- isbell.2018.pt(adf = comb_shores, RY.exp =  c(0.33, 0.33, 0.33) )
neb





