project_path <- "/Users/juancarlosvillasenorderbez/Library/CloudStorage/GoogleDrive-juancarlos@ucsb.edu/Shared drives/emlab/projects/current-projects/transferable-conservation"
data_path <- "/Users/juancarlosvillasenorderbez/Library/CloudStorage/GoogleDrive-juancarlos@ucsb.edu/Shared drives/emlab/data"

library(janitor)
library(data.table)
library(tidyverse)

area <-
  raster::raster(file.path(project_path,
                           "processed_data",
                           "suitability.tif")) %>% 
  raster::area() %>% 
  raster::as.data.frame(xy = T) %>% 
  rename(lon = x, lat = y, area = layer)

aquamaps_list <-
  read_csv(
    file.path(data_path, "aquamaps-v10-2019", "speciesoccursum.csv"),
    col_types = cols(
      SpeciesID = col_character(),
      SpecCode = col_double(),
      Genus = col_character(),
      Species = col_character(),
      FBname = col_character(),
      OccurCells = col_double(),
      Kingdom = col_character(),
      Phylum = col_character(),
      Class = col_character(),
      Order = col_character(),
      Family = col_character()
    )
  ) %>%
  clean_names() %>%
  filter(occur_cells >= 10) %>%
  mutate(sci_name = paste(genus, species)) %>%
  select(family, genus, sci_name, species_id)

## PROCESSING ##################################################################
# Create a vector of unique species
spec_codes <- aquamaps_list %>%
  pull(species_id) %>%
  unique()

# length(spec_codes) # Number of species contained in our analysis

# Load all aquamaps data
suitability_df <-
  fread(
    file = file.path(
      data_path, "aquamaps-v10-2019", "hcaf_species_native.csv"
    )
  ) %>%
  clean_names() %>%                                                             # Clean column names
  .[species_id %in% spec_codes] %>%                                             # Keep only species with 10 or more records
  .[probability >= 0.5] %>%                                                     # Keep only probability of existing
  as_tibble()

p_data <- 
  suitability_df %>%
  rename(lon = center_long, lat = center_lat) %>% 
  left_join(area, by = c("lon", "lat")) %>%
  mutate(spp_benefit = probability * area) %>%
  group_by(species_id) %>%
  mutate(spp_pct_area = area / sum(area),
         spp_pct_benefit = spp_benefit / sum(spp_benefit)) %>%
  ungroup() %>% 
  rename(spp_area = area)

saveRDS(object = p_data,
        file = here::here("data", "processed_data", "species_data.rds"))