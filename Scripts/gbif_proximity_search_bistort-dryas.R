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

## **NOTE** you need to wait for this to download completely
gbif_focal <- occ_download_get(gbif_download, overwrite = TRUE)

### ***CHECK** that the zip file has downloaded before proceeding

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

# get BC boundary **NOTE** it asks you for a response
bc.boundary <- bc_bound()

# get BEC data **NOTE** This takes a long time (10 min)
bec.data <- bec()

# **NOTE** you may wish to subset the BEC data at this stage, to 
# make things go faster...  e.g.:

bec.essf.alpine <- bec.data %>%
  filter(ZONE %in% c("BAFA", "ESSF", "IMA"))

# set CRS
bc.boundary <- st_transform(bc.boundary, crs = "EPSG:4326")

# **NOTE** the following takes a while
bec.essf.alpine <- st_transform(bec.essf.alpine, crs = "EPSG:4326")

# filter ORs to only those within originally selected BEC Zones

sf.focal.essf.alpine <- st_filter(sf_focal, bec.essf.alpine)

# **NOTE** This can take a while to render

#ggplot() + geom_sf(data = bc.boundary, fill = NA) +
#   geom_sf(data = sf_focal.essf, colour = "red", shape = 1) +
#   coord_sf(crs = st_crs(sf_focal.essf), xlim = c(-140, -110), ylim = c(48, 62)) +
#   theme_bw()

# let's narrow in to below 52 degrees north
sf.focal.essf.alpine.south <- st_crop(sf.focal.essf.alpine, xmin = -140, ymin = 48, xmax = -110, ymax = 52)

# remap
ggplot() + geom_sf(data = bc.boundary, fill = NA) +
  geom_sf(data = sf.focal.essf.alpine.south, colour = "red", shape = 1) +
  coord_sf(crs = st_crs(sf.focal.essf.alpine.south), xlim = c(-130, -112), ylim = c(48, 53)) +
  theme_bw()

########## BUFFER

# Adding a buffer to each point of our data frame (points become polygons)
sf_focal_buff <- st_buffer(sf.focal.essf.alpine.south, dist = proximity.dist)

# Turning the data frame into a geometry object with a list for each record.
focal_buff_geom <- st_geometry(sf_focal_buff)

# Changing the geometry object into a wkt format for the occ_search function.
focal_buff_wkt <- st_as_text(focal_buff_geom)

########

# Getting the UsageKey for vascular plants 
# change this as required
backbone_vp <- name_backbone("Tracheophyta")

# Creating an empty list to store our data about co-occurring species
list_sps <- list()

# For-loop that cycles through each focal species occurrence record (OR)
# and searches for other species' ORs whose lat/long coords fall within
# the buffer polygon

## THIS MAY TAKE A WHILE depending on the number of focal points

#for(i in 1:2){
for(i in 1:length(focal_buff_wkt)){
  record <- gbif_wkt2bbox(focal_buff_wkt[i])    
  # **NOTE** the "occ_search" has default of 500 record limit
  # we'll use 1000 here
  
  # first check that records were found, if not, return "NA", if yes
  # return the records
  
  if (is.null(occ_search(phylumKey = backbone_vp$usageKey, geometry = record, 
                         hasCoordinate = TRUE, limit = 1000)$data)
  ) { list_sps[[i]] <- NA } # value if no records found
  
  else(list_sps[[i]] <-occ_search(phylumKey = backbone_vp$usageKey, geometry = record, 
                                  hasCoordinate = TRUE, limit = 1000)$data
  )
  
  print(paste("done loop", i, sep = " ")) # comment out if you don't want
  # to keep track of progress
}

# So the result is a list, each element being a dataframe with rows as ORs and
# columns as data about the ORs.  for instance:

head(list_sps[[1]])

########### Create co-occurrence long format tibble

# First get a full list of all species that occur at least once with focal
# it first extracts the unique species from each focal buffer polygon
# then does this for all polygons, but only returns the final unique set
# sorted alphabetically

full.co.species.list <- sort( # sort alphabetically
  unique( # return only unique names
    # unlist simply coerces the output (here from lapply) to a vector rather
    # than a list
    unlist(
      # this "lapply" applies the function to each element in the list
      lapply(list_sps, function(x){
        # return a vector of unique species names from given focal polygon
        c(unique(x[,"acceptedScientificName"]))
      })
    )
  )
)

############# Initiate loop

# make dataframe from SF
sf_focal_buff.df <- st_drop_geometry(sf_focal_buff)

# create a long format tibble, with the following fields:
# "focal_gbifID", which is the gbif ID for the focal species occurrence record
# the values in this field will be repeated as many times as there are unique
# occurrence records of other species associated with that focal record 
# "species" is the list of unique species co-occurring at that focal record, 
# and each unique species name will be repeated the number of unique occurrence
# records that occur in that polygon; 
# we repeat the species name to keep track of the frequency of ORs

# Here we create the tibble for the first focal record, then we'll loop through
# subsequent records and append output to the tibble

focal.co.long <- tibble(
  # establish first field "focal_gbifID"
  # repeat the values as many times as there are unique species associated
  # that's what the combo of rep > length > names > table figures out
  focal_gbifID = rep(
    # the gbif ID should be kept character
    as.character(sf_focal_buff.df[1, "gbifID"]$gbifID), 
    length(
      rep(
        names(
          table(list_sps[[1]]$acceptedScientificName)
        ), 
        table(list_sps[[1]]$acceptedScientificName)
      )
    )
  ),
  # now get the list of unique species names, and repeat each name
  # the number of unique ORs for that species
  species = rep(
    names(
      table(list_sps[[1]]$acceptedScientificName)
    ), 
    table(list_sps[[1]]$acceptedScientificName)
  )
)

# now loop through remaining focal records, and append results each time

for (i in 2:length(list_sps)){
  #for (i in 2:5){
  temp <- tibble(
    # establish first field "focal_gbifID"
    # repeat the values as many times as there are unique species associated
    # that's what the combo of rep > length > names > table figures out
    focal_gbifID = rep(
      # the gbif ID should be kept character
      as.character(sf_focal_buff.df[i, "gbifID"]$gbifID), 
      length(
        rep(
          names(
            table(list_sps[[i]]$acceptedScientificName)
          ), 
          table(list_sps[[i]]$acceptedScientificName)
        )
      )
    ),
    # now get the list of unique species names, and repeat each name
    # the number of unique ORs for that species
    species = rep(
      names(
        table(list_sps[[i]]$acceptedScientificName)
      ), 
      table(list_sps[[i]]$acceptedScientificName)
    )
  )
  # append
  focal.co.long <- dplyr::bind_rows(focal.co.long, temp)
  print(paste("done loop", i, sep = " "))
  
}

# now we can summarise a few things

# first a frequency table for co-occurring species, ordered most common to least
freq.species <- focal.co.long %>% 
  count(species) %>%
  arrange(desc(n))

# number of different species by focal record 
point.richness <- focal.co.long %>%
  group_by(focal_gbifID) %>%
  summarise(richness = n_distinct(species)) %>%
  arrange(desc(richness))

# filter by some other species like "Picea" and "Abies"

points.with.dryas <- focal.co.long %>%
  filter(grepl('Dryas', species)) %>%
  group_by(focal_gbifID) %>%
  count(species) %>%
  arrange(desc(n))

# map those points

points.of.interest.dryas <- sf.focal.essf.alpine.south %>%
  filter(as.character(gbifID) %in% points.with.trees$focal_gbifID)

ggplot() + geom_sf(data = bc.boundary, fill = NA) +
  geom_sf(data = points.of.interest.dryas, colour = "red", shape = 1) +
  coord_sf(crs = st_crs(points.of.interest.dryas), xlim = c(-130, -112), ylim = c(48, 53)) +
  theme_bw()

points.of.interest.dryas.csv <- st_drop_geometry(points.of.interest.dryas)

# write to csv

readr::write_csv(freq.species, "species_frequency.csv")
readr::write_csv(point.richness, "point_richness.csv")
readr::write_csv(points.of.interest.dryas.csv, "points_of_interest.csv")

# export KML for viewing in Google Earth

st_write(points.of.interest.dryas, "points_of_interest_dryas.kml", driver = "kml", delete_dsn = TRUE)

