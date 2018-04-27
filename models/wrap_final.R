# Overarching function that receives polygon and 
# gives back buildings with predictions
#install.packages("osmdata")
#install.packages('rworldmap')
library(osmdata)
library(leaflet)
require(tcltk)
require(sp)
require(RANN)
require(rworldmap)
library(devtools)
library(htmlwidgets)
library(webshot)
library(mapview)
library(randomForest)
library(stringr)
library(pryr)
#install_github("wch/webshot")
#webshot::install_phantomjs()
#devtools::install_github("r-spatial/mapview@develop")
setwd("C:\\Users\\Esther_Su\\Desktop\\URAP\\Data")
region_map <-  read.csv('africa_abbr.csv')
input_type = readline("Pleae enter the input type number: 1.Two longitudes/latitudes 2.Coordinates 3.Region name(e.g.RWA)?") 

if (input_type == "1"){ # input longtitudes and latitudes
  long1 <- readline("What is the first number of longitude?") 
  long2 <- readline("What is the second number of longitude?")
  lat1 <- readline("What is the first number of latitude?") 
  lat2 <- readline("What is the second number of latitude?")
  long1 = min(as.numeric(long1),as.numeric(long2))
  long2 = max(as.numeric(long1),as.numeric(long2))
  lat1 = min(as.numeric(lat1),as.numeric(lat2))
  lat2 = max(as.numeric(lat1),as.numeric(lat2))
  # test 
  # long1 =  30.5
  # long2 = 30.0
  # lat1 = -2.35
  # lat2 = -2.00
  
  # checking latitude and longtitude
  if ( abs(long1 - long2) <= 0 | abs(lat1 - lat2)<=0 | abs(lat1)>90 | abs(lat2)>90 | abs(long1)>180 | abs(long2)>180 ){
    #msgBox <- tkmessageBox(title = "Error message box",
    #message = "Both the range of latitude and longtitude should greater than 1!", icon = "info", type = "ok")
    print("Error in latitudes or longtitudes")
  } 
  
  # for given polygon get buildings
  #q0 <- opq(bbox = c( min(long1,long2), min(lat1,lat2), max(long1,long2),max(lat1,lat2)), timeout = 1000, memsize=8000000000) 
  q0 <- opq(bbox = c( min(long1,long2), min(lat1,lat2), max(long1,long2),max(lat1,lat2)))
} else if(input_type == "2"){ #coordinate
  flpath <- readline("Please enter the filname with filepath?")
  header <- readline("Does the file contain a header: 1.True, 2.False?")
  # test
  flpath <- "coordinates.csv"
  coordinates <- read.csv(flpath,header= ifelse(header=="1",TRUE,FALSE)) 
  q0 <- opq(bbox =as.matrix(coordinates))
}else if(input_type == "3"){ #country name
  cnames <-readline("Please enter the country name or ISO code?")
  #test
  #bb = getbb('RWA') # abbr
  #bb = getbb('Rwanda') #countries
  bb = getbb(cnames)
  q0 <- opq(bbox = bb)
}else{print("error")}
q1 <- add_osm_feature(q0, key = 'building')
q2 <- add_osm_feature(q0, key = 'highway')
x <- osmdata_sp(q1)
#object_size(x)
y <- osmdata_sp(q2)
#(m0 <- leaflet() %>% addTiles() %>% addPolygons(data=x$osm_polygons, color="blue"))
(m0 <- leaflet() %>% addProviderTiles(providers$OpenStreetMap) %>% addPolygons(data=x$osm_polygons, color="blue"))
#somePlace <- ggmap::geocode("Vienna")
#mapview(someplace)
# m0$width <- 874
# m0$height <- 700
# saveWidget(m0, 'm0.html', selfcontained = FALSE)
# mapshot(m0, file = paste0(getwd(),"/","original_graph.jpg"))
# saveWidget(m0, "original_graph.html")
# webshot("original_graph.html", file = paste0(getwd(),"/","original_graph.png"),
#         cliprect = "viewport")
#%>%addPolygons(data=y$osm_polygons, color="blue") 

# Generate features
building_poly <- x$osm_polygons
with_type <- which(!is.na(building_poly$building) & building_poly$building != "yes") 
orig_type <- ifelse(building_poly$building[with_type] == "detached" | 
                    building_poly$building[with_type] == "house" | 
                    building_poly$building[with_type] == "hut" | 
                    building_poly$building[with_type] == "residential" |
                    building_poly$building[with_type] == "Teacher_Housing",0,1)       
to_predict <- which(building_poly$building == "yes")
print(paste0("There are total ", length(with_type)+length(to_predict), " buildings within the polygon."))
print(paste0("There are ", length(with_type), " buildings with types."))
print(paste0("There are ",length(to_predict)," buildings require to predict their types!"))
#unique(structures_poly_selected$building)
# levels(building_poly$building)
# Subset only those with type
structures_poly_selected <- building_poly[to_predict,]

orig_building_coord <- t(sapply(slot(building_poly[with_type,], "polygons"), slot, "labpt"))
orig_long <- orig_building_coord[,1]
orig_lat <- orig_building_coord[,2]
orig_w_type <- cbind(orig_long,orig_lat,orig_type,rep(0,length(orig_type)))

# building_type <- ifelse(structures_poly_selected$building == "detached" | 
#                           structures_poly_selected$building == "house" | 
#                           structures_poly_selected$building == "hut" | 
#                           structures_poly_selected$building == "residential" |
#                           structures_poly_selected$building == "Teacher_Housing",0,1 )
# getdata
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
points= as.matrix(coordinates)
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
region_out <- sapply(list(countries_loc)[[1]],function(x) region_map$Region[which(region_map$Abbr== x)])
country_names <- sapply(list(countries_loc)[[1]],function(x) region_map$Sub_Region [which(region_map$Abbr== x)])
country_names <- unique(country_names)
print(paste0("The area to predict is located in ",country_names))
# Choose the region is the polygon major locates in
region <- rep(as.factor(levels(region_out))[which.max(table(region_out))],length(structures_poly_selected ))


# load model and mkae predictions
# setwd("C:\\Users\\Esther\\Desktop\\URAP")
#model = load('africa_rf_mod2.RData')
model = load('africa_rf_mod_201804.RData')

# prection
all_structure_points_test <- data.frame(region, poly_area, dist_to_nearest, area_of_nearest, dist_to_nearest_road,n_poly,poly_complexity,poly_complexity_of_nearest,n_poly_of_nearest,lng,lat)
test_data_predictions <- predict(rf_mod, all_structure_points_test, predict.all=T)
test_data_predictions_prob <- apply(test_data_predictions$individual, 1, function(x){mean(as.numeric(x))})
opt_cutoff = 0.194
test_data_predictions_output <- as.numeric(test_data_predictions_prob>opt_cutoff)
final_temp <-  cbind(all_structure_points_test$lng,all_structure_points_test$lat,test_data_predictions_output,rep(1,length(test_data_predictions_output)))
final_result <- data.frame(rbind(final_temp, orig_w_type ))
colnames(final_result) <- c("Longtitude","Latitude","Sprayable","Prediction")
rownames(final_result) <- c(rownames(all_structure_points_test),rownames(orig_w_type))
#building_sprayable <- factor(building_type,levels=c(0,1))

# evaluation
africa_stat <- read.csv(paste0(getwd(),"/","agg_output_201804.csv"))
loc = sapply(str_to_lower(country_names), function(x) which(x ==africa_stat$countries..i..))
out_stat = t(africa_stat[loc,2:length(africa_stat)])

#test_data_predictions_prob <- apply(test_data_predictions$individual, 1, function(x){mean(as.numeric(x))})
#auc_test <- AUC(test_data_predictions_prob, building_sprayable)
#pred_test <- prediction(test_data_predictions_prob, building_sprayable)
#perf_test <- performance(pred_test,"tpr","fpr")
#perf_test <- performance(pred_test,"auc")
# accuracy <- sum(test_data_predictions[[1]] ==building_sprayable)/length(test_data_predictions[[1]])


# output
point.df <- data.frame(
  lat = final_result$Latitude,
  lng = final_result$Longtitude,
  size = rep(1, nrow(final_result) ),
  predict_ind = final_result$Prediction,
  color = c(ifelse(final_result[final_result[,4] ==1,3]==1,'red','limegreen'),
            ifelse(final_result[final_result[,4] ==0,3]==1,'pink','lightgreen') )
 )

m <- leaflet(point.df) %>%
  addProviderTiles(providers$OpenStreetMap) %>%
  #setView(lng = mean(long1,long2), lat = mean(lat1,lat2), zoom = 12)
  fitBounds(long1, lat1, long2, lat2)
(m_final<- m %>% addCircles(data = point.df, lng = point.df$lng, lat = point.df$lat,color=point.df$color,radius = point.df$size,fill = TRUE) )

print(out_stat)
write.csv(out_stat,file= "prediction_statistcs.csv")  #return the expected performance of the prediction
write.csv(final_result,file= "prediction_output.csv") #return sprayable or not

## save result to png
saveWidget(m_final, "precition_graph.html", selfcontained = FALSE)
webshot("precition_graph.html", file = "precition_graph.png",
        cliprect = "viewport")
#mapshot(m_final, file = "precition_graph.jpg")

  