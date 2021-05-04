
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
head(ply_dat)

# remove the composition = "None" treatment as this is a general control
ply_dat <- 
  ply_dat %>% 
  filter( composition %in% c("All", "F", "B", "S", "L") )

# add column for total cover of the four focal species
ply_dat$total_focal_cover <- 
  ply_dat %>%
  select(all_of(foc_spp)) %>%
  rowSums(.)

ply_dat

# reorder columns
names(ply_dat)

ply_sum <- 
  ply_dat %>%
  group_by(shore, month, mono_mix, composition) %>%
  summarise(sargassum = mean(sargassum, na.rm = TRUE),
            bifurcaria = mean(bifurcaria, na.rm = TRUE),
            f_serratus = mean(f_serratus, na.rm = TRUE),
            l_digitata = mean(l_digitata, na.rm = TRUE), .groups = "drop")

ply_mix <- 
  ply_sum %>%
  filter(composition == "All") %>%
  select(-composition, -mono_mix) %>%
  pivot_longer(cols = foc_spp,
               names_to = "species",
               values_to = "Y") %>%
  rename(place = shore, time = month) %>%
  arrange(place, time, species)


ply_mono <- 
  ply_sum %>% 
  filter(composition != "All") %>%
  pivot_longer(cols = foc_spp,
               names_to = "species",
               values_to = "M") %>%
  mutate(species_comp = toupper(substr(species, 1, 1)) ) %>%
  filter(composition == species_comp) %>%
  select(place = shore, time = month, species, M) %>%
  arrange(place, time, species)

# join the mixtures and the monocultures
ply_part <- full_join(ply_mono, ply_mix, by = c("time", "place", "species"))

# create a sample column
ply_part$sample <- paste(ply_part$place, ply_part$time, sep = "_")

# reorder the columns
ply_part <- 
  ply_part %>%
  select(sample, time, place, species, M, Y)

# make the monocultures non-zero
ply_part <- 
  ply_part %>%
  mutate(M = if_else(M < 0.1, 0.1, M))

# try the partition
isbell.2018.pt(adf = ply_part, RY.exp =  c(0.25, 0.25, 0.25, 0.25))

# write a function to generate random RY.exp
random_RY.exp <- function(n_sp = 4, by = 0.001) {
  x <- seq(0, 1, by)
  y <- sample(1:length(x), (n_sp-1) )
  
  z <- diff(c(0, x[sort(y)], 1))
  
  if( sum(z) != 1) {
    stop("error!")
  }
  
  return(z)
}

random_RY.exp(n_sp = 4, by = 0.001)

n <- 100
part_out <- vector("list", length = n)
for (i in 1:n) {
  
  z <- random_RY.exp(n_sp = 4, by = 0.001)
  part_out[[i]] <- isbell.2018.pt(adf = ply_part, RY.exp = z)
  
}

part_out <- 
  bind_rows(part_out, .id = "run") %>%
  group_by(biodiversity_effect) %>%
  summarise(effect_size_m = mean(effect_size),
            effect_size_sd = sd(effect_size), .groups = "drop")

part_out

df1 <- 
  part_out %>%
  filter(biodiversity_effect %in% c("local_complementarity_effect", "total_complementarity_effect",
                                    "local_selection_effect", "total_selection_effect"))

df2 <- 
  part_out %>%
  filter(biodiversity_effect %in% c("net_biodiversity_effect", "total_complementarity_effect",
                                    "non_random_overyielding", "total_insurance"))
  

df3 <- 
  part_out %>%
  filter(biodiversity_effect %in% c("spatial_insurance", "temporal_insurance", "spatio_temporal_insurance", "average_selection"))

ggplot(data = df1,
         mapping = aes(x = biodiversity_effect, y = effect_size_m)) +
  geom_point() +
  geom_errorbar(mapping = aes(ymin = effect_size_m - effect_size_sd,
                              ymax = effect_size_m + effect_size_sd),
                width = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

ggplot(data = df2,
       mapping = aes(x = biodiversity_effect, y = effect_size_m)) +
  geom_point() +
  geom_errorbar(mapping = aes(ymin = effect_size_m - effect_size_sd,
                              ymax = effect_size_m + effect_size_sd),
                width = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

ggplot(data = df3,
       mapping = aes(x = biodiversity_effect, y = effect_size_m)) +
  geom_point() +
  geom_errorbar(mapping = aes(ymin = effect_size_m - effect_size_sd,
                              ymax = effect_size_m + effect_size_sd),
                width = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))




# plot covariances between monoculture yields and mixture relative abundance
ply_part %>%
  group_by(sample) %>%
  mutate(Y_RA = Y/sum(Y) ) %>%
  ungroup() %>%
  group_by(time, species) %>%
  summarise(M = mean(M),
            Y_RA = mean(Y_RA)) %>%
  ggplot(data = .,
         mapping = aes(x = M, y = Y_RA, colour = species)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic()

ply_part %>%
  group_by(sample) %>%
  mutate(Y_RA = Y/sum(Y) ) %>%
  ungroup() %>%
  group_by(place, species) %>%
  summarise(M = mean(M),
            Y_RA = mean(Y_RA)) %>%
  ggplot(data = .,
         mapping = aes(x = M, y = Y_RA, colour = species)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic()


ply_part





### old code riffy raffy

ply_dat <- 
  ply_dat %>%
  select(shore, time = month, place = pool, grid_squares, 
         composition, mono_mix, total_focal_cover, cover_nocrusts)

# add a shore time variable
d <- mapply(function(x, y) {as.character(paste(x, y, sep = "_")) }, ply_dat$shore, ply_dat$time )
names(d) <- NULL

ply_dat$shore_time <- d

# summarise these data
ply_sum <- 
  ply_dat %>%
  group_by(shore, time, mono_mix, composition) %>%
  summarise(total_focal_cover = mean(total_focal_cover, na.rm = TRUE), .groups = "drop") 

ply_sum

# get the mixtures only
ply_dat$composition

ply_mix <- 
  ply_sum %>%
  filter(composition == "FLB") %>%
  select(shore, time, Y = total_focal_cover) 

ply_mix

# get the monocultures only
ply_mono <- 
  ply_dat %>%
  filter(composition %in% c("F", "L", "B")) %>%
  select(shore, time, M = total_focal_cover) 

ply_mono

# join the mixture and monoculture data


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



# function to make monocultures similar to mixtures
ply_mono <- split(ply_mono, ply_mono$shore_time)
ply_mix <- split(ply_mix, ply_mix$shore_time)
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

View(adf_out)

