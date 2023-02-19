###
# This function takes a focal species, finds its occurrence records (ORs) in GBIF
# then finds all other ORs (of any species eg. in a plant lineage) within a 
# specified distance of each of its ORs
###

# Author:  Mary Paz Mañé, Jason Pither
# Date: 2023-02-14

## Required Libraries
require(rgbif)
require(dplyr)
require(sf)
require(ggplot2)
require(bcdata)
require(bcmaps)
library(readr)

# set the focal species
focal.species <- "Bistorta vivipara"

# set proximity distance to look for co-occurring species (metres)
proximity.dist <- 100

# Getting the usageKey for your species of interest
backbone <- name_backbone(focal.species, verbose = TRUE) # Gives the 
# species with occurrences

# Checking the number of occurrences associated with that Key.
occ_count(taxonKey = backbone[1,1], georeferenced = TRUE)  

# NOTE: here you may wish to restrict to certain region

# Downloading the information of the records and reading the downloaded data as a
#       table in a new variable.
gbif_download <- occ_download(
  pred("taxonKey", backbone[1,1]), 
  pred("hasCoordinate", TRUE),
  pred("stateProvince", "British Columbia"),
  format = "SIMPLE_CSV")

gbif_focal <- occ_download_get(gbif_download, overwrite = TRUE)

gbif_results <- occ_download_import(gbif_focal)

# get rid of any NA lat/long
gbif_results <- gbif_results %>%
  filter(!is.na(decimalLongitude))%>%
  filter(!is.na(decimalLatitude))

# check for zero lat/longs... uncomment this if you wish to view

gbif_results %>% filter(decimalLatitude == 0 | decimalLongitude == 0)

# get rid of those with zeroes if necessary

gbif_results <- gbif_results %>% filter(decimalLatitude != 0)

# Changing data frame into sf object (coords become geometry = point)
# setting coordinate reference system to WGS84 (appropriate for GBIF data)
sf_focal <- st_as_sf(gbif_results, coords = c("decimalLongitude", "decimalLatitude"), 
                     crs = "EPSG:4326")

# get BC boundary
bc.boundary <- bc_bound()
# get BEC data 
bec.data <- bec()
bc.boundary <- st_transform(bc.boundary, crs = "EPSG:4326")

# the following takes a minute!
bec.data <- st_transform(bec.data, crs = "EPSG:4326")

# filter BEC to just ESSF zones
bec.essf <- bec.data %>% filter(ZONE == "ESSF")

# filter ORs to only those within ESSF

sf_focal.essf <- st_filter(sf_focal, bec.essf)
  
#
ggplot() + geom_sf(data = bc.boundary, fill = NA) +
   geom_sf(data = sf_focal.essf, colour = "red", shape = 1) +
   coord_sf(crs = st_crs(sf_focal.essf), xlim = c(-140, -110), ylim = c(48, 62)) +
   theme_bw()

# let's narrow in to below 52 degrees north
sf_focal.essf.south <- st_crop(sf_focal.essf, xmin = -140, ymin = 48, xmax = -110, ymax = 52)

# remap
ggplot() + geom_sf(data = bc.boundary, fill = NA) +
  geom_sf(data = sf_focal.essf.south, colour = "red", shape = 1) +
  coord_sf(crs = st_crs(sf_focal.essf.south), xlim = c(-130, -112), ylim = c(48, 53)) +
  theme_bw()

##########

# Adding a buffer to each point of our data frame (points become polygons)
sf_focal_buff <- st_buffer(sf_focal.essf.south, dist = proximity.dist)

  # Turning the data frame into a geometry object with a list for each record.
  focal_buff_geom <- st_geometry(sf_focal_buff)
  
  # Changing the geometry object into a wkt format for the occ_search function.
  focal_buff_wkt <- st_as_text(focal_buff_geom)

# Creating an empty list to store our desired co-occurring species
list_sps <- list()

# Getting the UsageKey for vascular plants 
backbone_vp <- name_backbone("Tracheophyta")

# For-loop checking each occurrence point for the focal species and giving a list 
#       of co-occurring vascular plants

## THIS WILL TAKE A WHILE

#for(i in 1:2){
  for(i in 1:length(focal_buff_wkt)){
  record <- gbif_wkt2bbox(focal_buff_wkt[i])    ## With this I stopped getting the
  # internal server error but nothing happens
  # now working, but **NOTE** the "occ_search" has default of 500 record limit
  # i'll change to 1000 here
  list_sps[[i]] <- occ_search(phylumKey = backbone_vp[1,1], geometry = record, 
                              hasCoordinate = TRUE, limit = 1000)$data
}

# So the result is a list, each element being a dataframe with rows as ORs and
# columns as data about the ORs.  for instance:

head(list_sps[[1]])

# We need to create a type of adjacency matrix where each row is an occurrence
# record for the focal species, and the columns are the various co-ocurring 
# species, and the cells are filled with numbers that represent the frequency
# of occurrence records of that co-occurring species at that focal location

# First get a full list of all species that occur at least once with focal

full.co.species.list <- sort(unique(unlist(lapply(list_sps, function(x){c(unique(x[,"acceptedScientificName"]))}))))

focal.co.long <- tibble(focal_gbifID = rep(st_drop_geometry(sf_focal_buff[1, "gbifID"]), 
                   length(rep(names(table(list_sps[[1]]$acceptedScientificName)), 
                   table(list_sps[[1]]$acceptedScientificName))),
                    species = rep(names(table(list_sps[[1]]$acceptedScientificName)), 
                        table(list_sps[[1]]$acceptedScientificName))))

# make dataframe from SF
sf_focal_buff.df <- st_drop_geometry(sf_focal_buff)

# initialize a tibble, to which we'll add from loop
# this is a long format, with the gbifID of the focal occurrence record,
# and the names of species for each co-occurring record, repeated as many
# times as it occurs

# don't mind the messy base code here. i'll "tidy" later. for now it works!

spec.co.long <- tibble(
  focal_gbifID = rep(as.character(sf_focal_buff.df$gbifID[1]),length(rep(names(table(list_sps[[1]]$acceptedScientificName)),
                                                                         table(list_sps[[1]]$acceptedScientificName)))),
  species = rep(names(table(list_sps[[1]]$acceptedScientificName)),table(list_sps[[1]]$acceptedScientificName)))

for (i in 2:length(list_sps)){
#for (i in 2:5){
  temp <- tibble(
    focal_gbifID = rep(as.character(sf_focal_buff.df$gbifID[i]),length(rep(names(table(list_sps[[i]]$acceptedScientificName)),
                                                                           table(list_sps[[i]]$acceptedScientificName)))),
    species = rep(names(table(list_sps[[i]]$acceptedScientificName)),table(list_sps[[i]]$acceptedScientificName)))
  
  # append
  spec.co.long <- dplyr::bind_rows(spec.co.long, temp)
  print(paste("done loop", i, sep = " "))
  
}

# now we can summarise a few things

# first a frequency table for co-occurring species, ordered most common to least
freq.species <- spec.co.long %>% 
  count(species) %>%
  arrange(desc(n))

# number of different species by focal record 
point.richness <- spec.co.long %>%
  group_by(focal_gbifID) %>%
  summarise(richness = n_distinct(species)) %>%
  arrange(desc(richness))

# filter by some other species like "Picea" and "Abies"

points.with.trees <- spec.co.long %>%
  filter(grepl('Abies|Picea', species)) %>%
  group_by(focal_gbifID) %>%
  count(species) %>%
  arrange(desc(n))

# map those points

points.of.interest <- sf_focal.essf.south %>%
  filter(as.character(gbifID) %in% points.with.trees$focal_gbifID)

ggplot() + geom_sf(data = bc.boundary, fill = NA) +
  geom_sf(data = points.of.interest, colour = "red", shape = 1) +
  coord_sf(crs = st_crs(points.of.interest), xlim = c(-130, -112), ylim = c(48, 53)) +
  theme_bw()

points.of.interest.csv <- st_drop_geometry(points.of.interest)

# write to csv

readr::write_csv(freq.species, "species_frequency.csv")
readr::write_csv(point.richness, "point_richness.csv")
readr::write_csv(points.of.interest.csv, "points_of_interest.csv")

# export KML for viewing in Google Earth

st_write(points.of.interest, "points_of_interest.kml", driver = "kml", delete_dsn = TRUE)

