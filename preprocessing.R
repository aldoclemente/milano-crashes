## PREPROCESSING 

if(!require(pacman)) install.packages("pacman")
pacman::p_load("rstudioapi")

# setting working directory
setwd(dirname(getActiveDocumentContext()$path))

source("packages.R")
source("load_data.R")
source("LN-simulations/utils/utils.R")

if(!dir.exists("data/")) {
  dir.create("data/")
}

# 1. building road network -----------------------------------------------------

district <- "Lambrate"
foldername <- paste0("data/", district, "/")

if(!dir.exists(foldername)){
  dir.create(foldername)
}

# mapview()

its_bbox = st_bbox(c(xmin = 9.217 , ymin = 45.472 , 
                     xmax = 9.2421 , ymax = 45.486), crs = 4326) %>% st_as_sfc()

my_vectortranslate = c(
  # SQL-like query where we select only the following fields
  "-select", "highway", 
  
  "-where", "highway IN ('primary', 'secondary', 'tertiary','unclassified',
                        'primary_link', 'tertiary_link','residential')"
)

#'residential', 'living_street', 'path',
# 'service', 'tertiary_link', 'steps', 'trunk', 'trunk_link', 'path',
# 'unclassified','pedestrian', 'road', 'motorway'

osm_data = oe_get("Milan", vectortranslate_options = my_vectortranslate,
                  boundary= its_bbox)

mapview(osm_data) 
save(osm_data, file=paste0(foldername, district,"_osm.RData"))

# 2. building mesh relying on sfneworks package --------------------------------
sf_network <- as_sfnetwork(osm_data, directed = FALSE)

sf_network2 <- sf_network %>% 
  convert(to_spatial_subdivision, .clean = TRUE)

sf_network3 <- sf_network2 %>% 
  convert(to_components, .clean = TRUE, .select = 1L)

sf_network2
sf_network3

p1 = st_point(c(9.21906, 45.47212))
p2 = st_point(c(9.21856, 45.48433))
p3 = st_point(c(9.23521, 45.48569))
p4 = st_point(c(9.23707, 45.47238))

poly = st_multipoint(c(p1, p2, p3, p4)) %>%
  st_cast("POLYGON") %>%
  st_sfc(crs = 4326)

filtered = sf_network3 %>%
  activate("edges") %>%
  filter(edge_intersects(poly)) %>%
  activate("nodes") %>%
  filter(!node_is_isolated())

filtered_sf <- filtered %>%
  activate("edges") %>%
  st_as_sf()

mesh <- as.mesh.1.5D(filtered) # LN-simulations/utils/utlis.R

save(mesh, file=paste0(foldername, paste0("mesh.RData")))

# 3. building dataset ----------------------------------------------------------
folder.img <- paste0(foldername,"imgs/")
if(!dir.exists(folder.img))
  dir.create(folder.img)

incidenti_bbox <- data.frame(x.min = (min(mesh$nodes[,1]) + 0.1 * sd(mesh$nodes[,1])), 
                             x.max = (max(mesh$nodes[,1]) - 0.1 * sd(mesh$nodes[,1])),
                             y.min = (min(mesh$nodes[,2]) + 0.1 * sd(mesh$nodes[,2])), 
                             y.max = (max(mesh$nodes[,2]) - 0.1 * sd(mesh$nodes[,2])))

crashes_MI <- load_data()

mesi <- c("Jen", "Feb", "Mar", "Apr", "May", "June", "July",
          "Aug", "Sep", "Oct", "Nov", "Dec")

num_data <- matrix(0, nrow=12,ncol=1)
LOCATIONS <- list()
IDXS <- list()
for(i in 1:12){
  idxs =  intersect(which(crashes_MI$MESE_INCIDENTE == i),
                    intersect(
                      intersect(
                        which(crashes_MI$longitudine >= incidenti_bbox$x.min),
                        which(crashes_MI$longitudine <= incidenti_bbox$x.max)), 
                      intersect(
                        which(crashes_MI$latitudine >= incidenti_bbox$y.min), 
                        which(crashes_MI$latitudine <= incidenti_bbox$y.max))
                    ))
  IDXS[[i]] <- idxs
  
  locations <-  data.frame(lon = crashes_MI$longitudine[idxs], 
                           lat = crashes_MI$latitudine[idxs])
  
  LOCATIONS[[i]] <- locations
  
  locations <- sf::st_as_sf(locations, 
                            coords = c("lon", "lat"),
                            crs = 4269) 
  num_data[i] <- nrow(locations)
  {  
    milano_incidenti_plot <- ggplot() +
      geom_sf(data = st_as_sf(filtered,"edges"), 
              inherit.aes = FALSE, 
              color = "black", 
              size = .0005,
              alpha = .6) +
      geom_sf(data = locations, 
              inherit.aes = FALSE, 
              color = "red3", 
              size = 0.75,
              alpha = 1.) +
      geom_rect( aes(xmin=incidenti_bbox$x.min, xmax=incidenti_bbox$x.max, 
                     ymin=incidenti_bbox$y.min, ymax=incidenti_bbox$y.max), 
                 color="red", fill="transparent")+
      theme_bw() +
      labs(title = paste("car crashes (", mesi[i], " 2018)",sep=""),
           x = "Longitudine",
           y = "Latitudine")+
      theme(plot.title = element_text(hjust = 0.5))
  }
  
  ggsave(filename = paste0(folder.img, district,"_crashes_", i,"_.pdf"), 
         plot=milano_incidenti_plot,
         width=6, height=6, units = "in")
  
}

locs <- matrix(nrow=0,ncol=2)
idxs <- rep(0,length=0)
crashes <- matrix(nrow=0, ncol=ncol(crashes_MI))
for(i in 1:12){
  idxs <- c(idxs, IDXS[[i]])
  crashes <- rbind(crashes, crashes_MI[IDXS[[i]],])
  locs <- rbind(locs, LOCATIONS[[i]])
}

crashes$ID <- idxs
save(LOCATIONS, IDXS, crashes, file=paste0(foldername,"data_raw.RData"))

# ------------------------------------------------------------------------------
# projecting 
locs <- fdaPDE::projection.points.1.5D(mesh, locs)
crashes$longitudine <- locs[,1]
crashes$latitudine <- locs[,2]

LOCS <- data.frame(lon = locs[,1], lat = locs[,2])
LOCS <- sf::st_as_sf(LOCS, coords = c("lon", "lat"), crs = 4269) 
{  
  incidenti_plot <- ggplot() +
    geom_sf(data = st_as_sf(filtered,"edges"), 
            inherit.aes = FALSE, 
            color = "black", 
            size = .0005,
            alpha = .6) +
    geom_sf(data = LOCS, 
            inherit.aes = FALSE, 
            color = "red3", 
            size = 0.75,
            alpha = 1.) +
    geom_rect( aes(xmin=incidenti_bbox$x.min, xmax=incidenti_bbox$x.max, 
                   ymin=incidenti_bbox$y.min, ymax=incidenti_bbox$y.max), 
               color="red", fill="transparent")+
    theme_bw() +
    labs(title = paste0(district," crashes (2018)",sep=""),
         x = "Longitudine",
         y = "Latitudine")+
    theme(plot.title = element_text(hjust = 0.5))
}

ggsave(filename = paste0(folder.img, district, "_crashes_2018.pdf"), 
       plot=incidenti_plot,
       width=6, height=6, units = "in")