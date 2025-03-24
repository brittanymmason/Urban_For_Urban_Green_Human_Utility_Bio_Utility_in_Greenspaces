#############################################################################
# Summarize iNaturalist data for Broward County greenspaces
#############################################################################

library(sf)
library(dplyr)
library(readr)
library(spatialEco)

#read in shapefile of Broward county greenspaces
parks <- st_read("Data/shapefiles/urban_parks.shp") %>% st_set_crs(4326)

# read in iNaturalist data and set as a shapefile
inat <- readRDS("Data/iNaturalist_data.RDS") %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# calculate polygon area in hectares 
parks$poly_area <- as.numeric(st_area(parks)/10000)

# extract iNaturalist points that overlap the greenspace shapefile and then use the 
# st_join to combine attributes from parks to the iNaturalist data
bcp_inat <- inat[st_intersects(inat, parks) %>% lengths > 0,]
inat_in_parks <- st_join(bcp_inat, parks)

# read in broward_poly_metadata CSV file to provide meaningful names for each polygon id
poly_names <- read_csv("Data/broward_polygons_metadata.csv", locale=locale(encoding="latin1"))

# create poly_id, which is a unique name for each park
poly_names$poly_id <- paste(poly_names$municipality, poly_names$park_name,  sep="_")

# join the poly_names data frame with the iNaturalist data
inat_parks <- as.data.frame(left_join(inat_in_parks, poly_names, by="poly_id"))

# for the iNaturalist data, we are not interested in "casual" grade observations, so filter out this data
inat_obs <- inat_parks %>% dplyr::filter(quality_grade!="casual")

# count of observations by park for all iNaturalist data
unique_sp_park_inat <- inat_obs %>% 
  group_by(poly_id) %>% 
  summarise(park_name=first(park_name),
            park_size_ha=first(poly_area),
            total_observations=n(),
            total_observers=n_distinct(user_id),
            family_count_observations=n_distinct(taxon_family_name),
            order_count_observations=n_distinct(taxon_order_name),
            species_count_observations=n_distinct(taxon_species_name))

# remove parks without verified observations
gbif <- readRDS("Data/GBIF_data.RDS") %>%
  st_as_sf(., crs = 4326) %>%
  filter(st_intersects(., parks) %>% lengths > 0) %>% 
  st_join(parks)

unique_sp_park <- unique_sp_park_inat %>%
  filter(poly_id %in% gbif$poly_id)

# now add the all parks back into the data frame
alldata <- left_join(parks, unique_sp_park, by="poly_id")

# remove unnecessary columns 
alldata <- alldata %>% 
  as.data.frame() %>%
  dplyr::select(-park_name, -park_size_ha)

# join data frames
alldata <- as.data.frame(left_join(alldata, poly_names, by="poly_id"))

# clean data
alldata <- alldata %>% 
  select(-municipality, -park_name, -grns__2) %>%
  rename(park_size_ha=poly_area)

# there are two parks that have two polygons, we need to sum the area, but use the first of all other fields
alldata <- alldata %>% group_by(poly_id) %>%
  mutate(park_size_ha=sum(park_size_ha)) %>%
  distinct(poly_id, .keep_all=TRUE)

# save the final data frame
write_csv(alldata, "Data/broward_greenspace_summary/iNaturalist_park.csv")
