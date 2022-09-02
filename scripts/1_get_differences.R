
library(data.table)
library(janitor)
library(furrr)
library(tidyverse)

project_path <- "/Users/juancarlosvillasenorderbez/Library/CloudStorage/GoogleDrive-juancarlos@ucsb.edu/Shared drives/emlab/projects/current-projects/transferable-conservation"
data_path <- "/Users/juancarlosvillasenorderbez/Library/CloudStorage/GoogleDrive-juancarlos@ucsb.edu/Shared drives/emlab/data"

spp_data <- readRDS(here::here("data", "processed_data", "species_data.rds"))


read <- function(bubble) {
  
  curves <-
    readRDS(
      file = file.path(
        project_path,
        "processed_data",
        "supply_curves",
        "no_mpas",
        paste0(bubble, "_eez_supply_curves_no_mpas.rds")
      ))
  
  return(curves)
}

protect <- function(r, data, spp_data) {
  # browser()
  protected <- data %>% 
    mutate(pixel_fraction = pmin(1 - ((ta - ((ta * r) / pct)) / area), 1)) %>% 
    filter(pixel_fraction >= 0) %>% 
    mutate(benefit = benefit * pixel_fraction,
           area = area * pixel_fraction,
           tb = cumsum(benefit),
           ta = cumsum(area))
  
  spp_benefits <- spp_data %>%
    filter(lat %in% protected$lat,
           lon %in% protected$lon) %>% 
    # inner_join(select(protected, lon, lat, pixel_fraction),
               # by = c("lon", "lat")) %>%
    # mutate(pixel_fraction = 1,
           # spp_pct_area = pixel_fraction * spp_pct_area,
           # spp_pct_benefit = pixel_fraction * spp_pct_benefit) %>%
    group_by(species_id) %>%
    summarize(pct_area_protected = sum(spp_pct_area, na.rm = T),
              pct_benefit_protected = sum(spp_pct_benefit, na.rm = T)) %>%
    ungroup()
  
  n_spp_area <- 1#sum(spp_benefits$pct_area_protected >= r)
  n_spp_benefit <- 2#sum(spp_benefits$pct_benefit_protected >= r)
    
    out <- protected %>% 
      tail(1) %>% 
      mutate(n_spp_area = n_spp_area,
             n_spp_benefit = n_spp_benefit)
    
    return(out)
}

process <- function(data, bubble, R = seq(0.01, 1, by = 0.01), spp_data){
  # browser()
  
  a <- data %>% 
    select(iso3, {{bubble}}, lat, lon, area, benefit, cost) %>% 
    mutate(mc = cost / area) %>% 
    group_by_at(c("iso3", bubble)) %>% 
    arrange(mc) %>% 
    mutate(ta = cumsum(area),
           tb = cumsum(benefit),
           pct = ta / max(ta)) %>% 
    ungroup()
  
  b <- data %>% 
    select(iso3, {{bubble}}, lat, lon, area, benefit, cost) %>% 
    mutate(mc = cost / benefit) %>% 
    group_by_at(c("iso3", bubble)) %>%  
    arrange(mc) %>% 
    mutate(ta = cumsum(area),
           tb = cumsum(benefit),
           pct = ta / max(ta)) %>% 
    ungroup()
  
  g <- data %>% 
    select(iso3, {{bubble}}, lat, lon, area, benefit, cost) %>% 
    group_by_at(c("iso3", bubble)) %>%  
    arrange(desc(benefit)) %>% 
    mutate(ta = cumsum(area),
           tb = cumsum(benefit),
           pct = ta / max(ta)) %>% 
    ungroup()
  
  A <- a %>% 
    group_by_at(c("iso3", bubble)) %>%
    nest() %>% 
    expand_grid(r = R) %>% 
    mutate(res = future_map2(r, data, protect, spp_data = spp_data)) %>% 
    unnest(res) %>% 
    group_by(r) %>% 
    summarize(tb = sum(tb),
              ta = sum(ta),
              n_spp_area = unique(n_spp_area),
              n_spp_benefit = unique(n_spp_benefit)) %>% 
    ungroup() %>% 
    mutate(approach = "area_target")
  
  B <- b %>% 
    group_by_at(c("iso3", bubble)) %>% 
    nest() %>% 
    expand_grid(r = R) %>% 
    mutate(res = future_map2(r, data, protect, spp_data = spp_data)) %>% 
    unnest(res) %>% 
    group_by(r) %>% 
    summarize(tb = sum(tb),
              ta = sum(ta),
              n_spp_area = unique(n_spp_area),
              n_spp_benefit = unique(n_spp_benefit)) %>% 
    ungroup() %>% 
    mutate(approach = "benefit_target")
  
  G <- g %>%
    group_by_at(c("iso3", bubble)) %>% 
    nest() %>% 
    expand_grid(r = R) %>% 
    mutate(res = future_map2(r, data, protect, spp_data = spp_data)) %>% 
    unnest(res) %>% 
    group_by(r) %>% 
    summarize(tb = sum(tb),
              ta = sum(ta),
              n_spp_area = unique(n_spp_area),
              n_spp_benefit = unique(n_spp_benefit)) %>% 
    ungroup() %>% 
    mutate(approach = "griddy")
  
  data <- rbind(A, B, G)
  
  return(data)
}



# PROCESSS

plan(multisession)

data <- tibble(bubble = c("global", "realm", "province", "ecoregion")) %>% 
  filter(bubble == "global") %>% 
  mutate(data = map(bubble, read)) %>% 
  mutate(data = map2(data, bubble, process, R = 0.1, spp_data = spp_data))

# saveRDS(object = data, 
#         file = "data.rds")