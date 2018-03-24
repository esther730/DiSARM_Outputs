# Overarching function that receives polygon and 
# gives back buildings with predictions

library(osmdata)
library(leaflet)
require(tcltk)
require(sp)
require(rworldmap)

# input longtitudes and latitudes
long1 <- readline("What is the first number of longitude?") 
long2 <- readline("What is the second number of longitude?")
lat1 <- readline("What is the first number of latitude?") 
lat2 <- readline("What is the second number of latitude?")
long1 = min(as.numeric(long1),as.numeric(long2))
long2 = max(as.numeric(long1),as.numeric(long2))
lat1 = min(as.numeric(lat1),as.numeric(lat2))
lat2 = max(as.numeric(lat1),as.numeric(lat2))
# test 
long1 =  34.5
long2 = 34.0
lat1 = -2.35
lat2 = -2.30

# checking latitude and longtitude
if ( abs(long1 - long2) <= 0 | abs(lat1 - lat2)<=0 | abs(lat1)>90 | abs(lat2)>90 | abs(long1)>180 | abs(long2)>180 ){
    #msgBox <- tkmessageBox(title = "Error message box",
                           #message = "Both the range of latitude and longtitude should greater than 1!", icon = "info", type = "ok")
  print("Error in latitudes or longtitudes")
}

# for given polygon get buildings
q0 <- opq(bbox = c( min(long1,long2), min(lat1,lat2), max(long1,long2),max(lat1,lat2))) # Chiswick Eyot in London, U.K.
q1 <- add_osm_feature(q0, key = 'building')
q2 <- add_osm_feature(q0, key = 'highway')
x <- osmdata_sp(q1)
y <- osmdata_sp(q2)
leaflet() %>% addTiles() %>% addPolygons(data=x$osm_polygons, color="blue")
#%>%addPolygons(data=y$osm_polygons, color="blue") 





# Generate features
building_poly <- x$osm_polygons
with_type <- which(!is.na(building_poly$building) ) 
# levels(building_poly$building)
# Subset only those with type
structures_poly_selected <- building_poly[with_type,]
#structure_points <- SpatialPointsDataFrame(rgeos::gCentroid(structures_poly_selected,byid=TRUE), structures_poly_selected@data, match.ID=FALSE)

# Q: how to deal with "yes?"
building_type <- ifelse(structures_poly_selected$building == "detached" | 
                          structures_poly_selected$building == "house" | 
                          structures_poly_selected$building == "hut" | 
                          structures_poly_selected$building == "residential" |
                          structures_poly_selected$building == "Teacher_Housing",0,1 )
areas <- sapply(slot(structures_poly_selected, "polygons"), function(x) sapply(slot(x, "Polygons"), slot, "area"))
poly_coords <- sapply(slot(structures_poly_selected, "polygons"), function(x) sapply(slot(x, "Polygons"), slot, "coords"))
building_corrds <- t(sapply(slot(structures_poly_selected, "polygons"), slot, "labpt"))
building_area <- sapply(slot(structures_poly_selected, "polygons"), slot, "area")
dist_to_nearest <- nn2(building_corrds,building_corrds, eps=0,searchtype="priority", k=5)$nn.dists[,2]
area_of_nearest <- building_area[nn2(building_corrds,building_corrds, eps=0,searchtype="priority", k=2)$nn.idx[,2]]
poly_area<-sapply(building_area, max) #Some buildings have multiple polygons so choose largest
n_poly <- sapply(areas, length)
poly_complexity <- sapply(poly_coords, function(x){length(unlist(x))})
poly_complexity_of_nearest <- poly_complexity[nn2(building_corrds,building_corrds, eps=0,searchtype="priority", k=2)$nn.idx[,2]]
n_poly_of_nearest <- n_poly[nn2(building_corrds,building_corrds, eps=0,searchtype="priority", k=2)$nn.idx[,2]]
#road_areas <- sapply(slot(road_poly, "polygons"), slot, "area")
#road_coords <- t(sapply(slot(road_poly, "polygons"), slot, "labpt"))

# Calculate distance to road
get_coords <- function(x){x@Lines} 
road_coords <- sapply(y$osm_lines@lines, get_coords)
road_coords <- sapply(road_coords, coordinates)
road_coords_all <- do.call("rbind", road_coords)
dist_to_nearest_road <- nn2(road_coords_all, building_corrds, eps=0,searchtype="priority", k=2)$nn.dists[,2]
lng <- building_corrds[,1]
lat <- building_corrds[,2]

# get country
points = data.frame(lon=c(long1, long1, long2, long2), lat=c(lat1, lat2, lat1, lat2))
coords2country = function(points)
{#ref: https://stackoverflow.com/questions/14334970/convert-latitude-and-longitude-coordinates-to-country-name-in-r
  countriesSP <- getMap(resolution='low')
  #countriesSP <- getMap(resolution='high') #you could use high res map from rworldxtra if you were concerned about detail
  
  # convert our list of points to a SpatialPoints object
  
  # pointsSP = SpatialPoints(points, proj4string=CRS(" +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
  
  #setting CRS directly to that from rworldmap
  pointsSP = SpatialPoints(points, proj4string=CRS(proj4string(countriesSP)))  
  
  # use 'over' to get indices of the Polygons object containing each point 
  indices = over(pointsSP, countriesSP)
  #indices$ADMIN #country
  indices$ISO3 
}
countries_loc <- as.character(coords2country(points))
print(countries_loc)
region_map <-  read.csv('africa_abbr.csv')
region_out <- sapply(list(countries_loc)[[1]],function(x) region_map$Region[which(region_map$Abbr== x)])
# Choose the region is the polygon major locates in
region <- rep(as.factor(levels(region_out))[which.max(table(region_out))],length(structures_poly_selected ))


# load model and mkae predictions
# setwd("C:\\Users\\Esther\\Desktop\\URAP")
model = load('africa_rf_mod2.RData')

# prection
building_sprayable <- factor(building_type,levels=c(0,1))
all_structure_points_test <- data.frame(region, poly_area, dist_to_nearest, area_of_nearest, dist_to_nearest_road,n_poly,poly_complexity,poly_complexity_of_nearest,n_poly_of_nearest,lng,lat)
test_data_predictions <- predict(rf_mod2, all_structure_points_test, predict.all=T)

# evaluation
test_data_predictions_prob <- apply(test_data_predictions$individual, 1, function(x){mean(as.numeric(x))})
#auc_test <- AUC(test_data_predictions_prob, building_sprayable)
pred_test <- prediction(test_data_predictions_prob, building_sprayable)
#perf_test <- performance(pred_test,"tpr","fpr")
#perf_test <- performance(pred_test,"auc")
accuracy <- sum(test_data_predictions[[1]] ==building_sprayable)/length(test_data_predictions[[1]])
expected_acc <- 

# output
point.df <- data.frame(
  lat = lat,
  lng = lng,
  size = rep(1,length(test_data_predictions[[1]])),
  color = ifelse(test_data_predictions[[1]]==1,'red','limegreen')
)

m <- leaflet(point.df) %>%
  addTiles() %>%
  #setView(lng = mean(long1,long2), lat = mean(lat1,lat2), zoom = 12)
  fitBounds(long1, lat1, long2, lat2)
m %>% addCircles(data = point.df, lng = ~lng, lat = ~lat,color=~color,radius = ~size,fill = TRUE)
  