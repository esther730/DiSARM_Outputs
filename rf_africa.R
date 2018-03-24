# install.packages('caret')
# install.packages('raster')
# install.packages('glmnet')
# install.packages('fields')
# install.packages('rgdal')
# install.packages('rgeos')
# install.packages('leaflet')
# install.packages('maptools')
# install.packages('RANN')
# install.packages('gbm')
# install.packages('randomForest')
# install.packages('SuperLearner')
# install.packages('glmnet')
# install.packages('polspline')
# install.packages('mgcv')
# install.packages('hydroTSM')
# install.packages('earth')
# install.packages('party')
# install.packages('ggplot2')
# install.packages('e1071')
# install.packages('wesanderson')
# install.packages('raster')
# install.packages('cvAUC')
# install.packages('Hmisc')
# install.packages('gtools')
# install.packages('osmar')
# install.packages('stringr')
# install.packages('XML')
# install.packages('RCurl')
# install.packages('curl')
# install.packages('rlist')
# install.packages('httr')
# install.packages('foreach')
# install.packages('parallel')
# install.packages('doParallel')
# install.packages('readr')
# install.packages('ROCR')

library(caret)  # Classification and Regression Training
library(raster) # Geographic Data Analysis and Modeling
library(glmnet) # Lasso and Elastic-Net Regularized Generalized Linear Models
library(fields) # Tools for Spatial Data
library(rgdal) # Bindings for the 'Geospatial' Data Abstraction Library
library(rgeos) # Interface to Geometry Engine - Open Source ('GEOS')
library(leaflet) # Create Interactive Web Maps with the JavaScript 'Leaflet' Library
library(maptools) # Tools for Reading and Handling Spatial Objects
library(RANN) # Fast Nearest Neighbour Search (Wraps ANN Library) Using L2 Metric
library(gbm) # Generalized Boosted Regression Models
library(randomForest) # Breiman and Cutler's Random Forests for Classification and Regression
library(SuperLearner) # mplements the super learner prediction method and contains a library of prediction algorithms to be used in the super learner.
library(glmnet) # Lasso and Elastic-Net Regularized Generalized Linear Models
library(polspline) # Polynomial Spline Routines
library(mgcv) # Mixed GAM Computation Vehicle with Automatic Smoothness Estimation
library(hydroTSM) # Time Series Management, Analysis and Interpolation for Hydrological Modelling
library(earth) # Multivariate Adaptive Regression Splines
library(party) # A Laboratory for Recursive Partytioning
library(ggplot2)
library(e1071) # Misc Functions of the Department of Statistics, Probability Theory Group (Formerly: E1071), TU Wien
library(wesanderson) # A Wes Anderson Palette Generator
library(raster) # Reading, writing, manipulating, analyzing and modeling of gridded spatial data
library(cvAUC) # cross-validated area under the ROC curve (AUC) estimators
library(Hmisc) # Harrell Miscellaneous, data analysis, high-level graphics, utility operations, functions for computing sample size and power, importing and annotating datasets, imputing missing values, advanced table making, variable clustering, character string manipulation, conversion of R objects to LaTeX and html code, and recoding variables.
library(gtools) # Various R Programming Tools

require(osmdata)
require(osmar) # OpenStreetMap and R

require(stringr)
require(XML)
require(RCurl)
require(curl)
require(rlist)
require(httr)
require(foreach)
require(parallel) 
require(doParallel)
#require(readr)
require(ROCR)
library(sp)
library(lattice) # required for trellis.par.set():



#nCores <- as.integer(Sys.getenv("SLURM_CPUS_ON_NODE"))
#registerDoParallel(nCores)

website1 <- "http://download.geofabrik.de/africa.html" 
content <- readLines(website1)
Sys.sleep(sample(5, 1))
link <- getHTMLLinks(content)
pos <- sapply(link, function(x) gregexpr(".zip",x)[[1]][1]!=-1 )
links <- link[pos]
nafrica <- sum(pos)

#country = 'angola'

get_data <- function(structures_poly, training=TRUE,country){
  # Subset only those with type
  with_type <- which(!is.na(structures_poly$type))
  
  if(training==TRUE){
    structures_poly_selected <- structures_poly[with_type,]
  }else{
    structures_poly_selected <- structures_poly[-with_type,]  
  }
  
  structure_points <- SpatialPointsDataFrame(rgeos::gCentroid(structures_poly_selected,byid=TRUE), structures_poly_selected@data, match.ID=FALSE)
  
  # add country #Esther
  structure_points$country <- rep(country,length(with_type))
  
  # Calculate distance to road
  get_coords <- function(x){x@Lines} 
  road_coords <- sapply(roads@lines, get_coords)
  road_coords <- sapply(road_coords, coordinates)
  road_coords_all <- do.call("rbind", road_coords)
  
  structure_points$dist_to_nearest_road <- nn2(road_coords_all, structure_points@coords, eps=0,searchtype="priority", k=2)$nn.dists[,2]
 
  
  # Esther water, region, size
  # water_areas <- sapply(slot(water, "polygons"), slot, "area")
  # water_coords <- t(sapply(slot(water, "polygons"), slot, "labpt"))
  # 
  # structure_points$dist_to_nearest_water <- nn2(water_coords, structure_points@coords, eps=0,searchtype="standard", k=1)$nn.dists[,1]
  # structure_points$nearest_water_area <- water_areas[nn2(water_coords, structure_points@coords, eps=0,searchtype="standard", k=1)$nn.idx[,1] ]
  # structure_points$size <- rep(length(with_type),length(with_type))
  africa_info <- read.csv(paste0(getwd(),"/","africa_abbr.csv"))
  structure_points$region <- rep(africa_info$Region[i],length(with_type))
  
  # Calculate areas
  areas <- sapply(slot(structures_poly_selected, "polygons"), function(x) sapply(slot(x, "Polygons"), slot, "area"))
  structure_points$poly_area<-sapply(areas, max) #Some buildings have multiple polygons so choose largest
  structure_points$n_poly <- sapply(areas, length)

  # Calculate polygon 'complexity'
  poly_coords <- sapply(slot(structures_poly_selected, "polygons"), function(x) sapply(slot(x, "Polygons"), slot, "coords"))
  structure_points$poly_complexity <- sapply(poly_coords, function(x){length(unlist(x))})
  
  # Calculate distance to nearest structure
  structure_points$dist_to_nearest <- nn2(structure_points@coords,structure_points@coords, eps=0,searchtype="priority", k=5)$nn.dists[,2]
  
  # Calculate features of nearest structure
  structure_points$area_of_nearest <- structure_points$poly_area[nn2(structure_points@coords,structure_points@coords, eps=0,searchtype="priority", k=2)$nn.idx[,2]]
  structure_points$poly_complexity_of_nearest <- structure_points$poly_complexity[nn2(structure_points@coords,structure_points@coords, eps=0,searchtype="priority", k=2)$nn.idx[,2]]
  structure_points$n_poly_of_nearest <- structure_points$n_poly[nn2(structure_points@coords,structure_points@coords, eps=0,searchtype="priority", k=2)$nn.idx[,2]]
  structure_points$lng <- structure_points@coords[,1]
  structure_points$lat <- structure_points@coords[,2]
  
  
  return(structure_points)
}

optimal_cutoff <- function(predictions, observations){
  
  perc_corr_class_sprayable <- NULL
  perc_corr_class_not_sprayable <- NULL
  
  true_num_sprayable <- sum(observations)
  true_num_not_sprayable <- sum(observations==0)
  
  for(i in 1:1000){
    predicted_class_loop <- ifelse(predictions>=i/1000,1,0)
    perc_corr_class_sprayable <- c(perc_corr_class_sprayable, 
                                   round(sum(observations==1 & predicted_class_loop==1,na.rm=T)/true_num_sprayable,3))
    
    perc_corr_class_not_sprayable <- c(perc_corr_class_not_sprayable, 
                                       round(sum(observations==0 & predicted_class_loop==0,na.rm=T)/true_num_not_sprayable,3))
  }
  
  plot_data <- data.frame(threshold = rep(seq(0,1,length.out=1000),2),
                          perc_corr_class = c(perc_corr_class_sprayable, perc_corr_class_not_sprayable),
                          group = c(rep("Sprayable",1000), rep("Not sprayable",1000)))
  
  ggplot(data=plot_data, aes(x=threshold, y=perc_corr_class, group=group, color=group)) + 
    geom_line(size=1.5) + xlab("Threshold") + ylab("Proportion correctly classified")
  
  # Identify where the lines overlap
  #opt_threshold <- (which(perc_corr_class_sprayable == perc_corr_class_not_sprayable)/1000)[1]
  opt_threshold <- which.min((perc_corr_class_sprayable - perc_corr_class_not_sprayable)^2)[1]/1000
  
  
  # Table predicted class against observed using different thresholds
  predicted_class <- ifelse(predictions>=opt_threshold,1,0) # Threshold of 0.7 gives 95% sensitivity in both
  table(predicted_class, observations)
  return(list(plot_data = plot_data, 
              opt_cutoff = opt_threshold,
              confusion_matrix = table(predicted_class, observations)))
}

# data = links[41]
#setwd("C:\\Users\\Esther\\Desktop\\URAP")
#single <- function(data){
cnt = 0
spray_ratio_list = rep(NA,55)
#i = 39
# for (i in c(39:40)){
for (i in c(1:nafrica)){
  file <- basename(links[i])
  country <- gsub('-latest-free.shp.zip','',file)
  #download.file(paste0("http://download.geofabrik.de/",links[i]),file)
  # unzip(file, exdir = paste0(getwd(),"/osm_data/",country,"_OSM"))
  # if (file.info(file)$size) {}
  #file.remove(file) #due to memory reason
   # structures_poly <- rgdal::readOGR(dsn = paste0(getwd(),"/osm_data/",country,"_OSM"),
   #                           layer = "gis.osm_buildings_a_free_1")
  # 
  # save(structures_poly,file = paste0(getwd(),'/',country,"_structures_poly.RData"))
   load(paste0(getwd(),'/',country,"_structures_poly.RData"))
  #sum(is.na(structures_poly$type))/length(structures_poly$type)
  # slot(structures_poly, "data")
  # summary(structures_poly)
  # trellis.par.set(sp.theme())
 
  # if (length(structures_poly)<=3000){}
  # roads <- rgdal::readOGR(dsn = paste0(getwd(),'/osm_data/',country,"_OSM"),
  #                  layer = "gis.osm_roads_free_1")
  # save(roads,file = paste0(getwd(),'/',country,"_roads.RData"))
    load(paste0(getwd(),'/',country,"_roads.RData"))
  #water <- rgdal::readOGR(dsn = paste0(getwd(),'/osm_data/',country,"_OSM"), 
  #                        layer = "gis.osm_water_a_free_1")
  
  #summary(roads)
  structure_points <- get_data(structures_poly, training=TRUE,country=country)
  min_train_size <- 500
  if (length(structure_points) < min_train_size){
    add_structure_point_index= sample(c(1:length(structure_points)),min_train_size-length(structure_points),replace=TRUE)
    structure_points <- rbind(structure_points,structure_points[add_structure_point_index,])
  }
  max_train_size <- 5000
  if (length(structure_points) > max_train_size){
    reduce_structure_point_index= sample(c(1:length(structure_points)),max_train_size,replace=FALSE)
    structure_points <- structure_points[reduce_structure_point_index,]
  }
  
 
  
  
  structure_points$sprayable <- ifelse(structure_points$type == "detached" | 
                                         structure_points$type == "house" | 
                                         structure_points$type == "hut" | 
                                         structure_points$type == "residential" |
                                         structure_points$type == "Teacher_Housing",0,1)
  # Esther
  structure_points$size <- ifelse(nrow(structure_points)<500,0,ifelse(nrow(structure_points)<5000,1,2))
  #prediction_structure_points <- get_data(structures_poly,training=FALSE)
  save(structure_points,file = paste0(getwd(),'/',country,"_structure_points_201803.RData"))
  
  set.seed(1981)
  spray_ratio = sum(structure_points$sprayable==1)/sum(structure_points$sprayable==0)
  spray_ratio_list[i] = spray_ratio
  if( sum(structure_points$sprayable==1)== 0){
    random_non_sprayables <- sample(which(structure_points$sprayable==0),sum(structure_points$sprayable==0)*0.9)
    structure_points <- structure_points[c(random_non_sprayables),]
    structure_points_test  <- structure_points[-c(random_non_sprayables),]
  }else if(sum(structure_points$sprayable==0)== 0){
    random_sprayables <- sample(which(structure_points$sprayable==1),sum(structure_points$sprayable==1)*0.9)
    structure_points <- structure_points[c(random_sprayables),]
    structure_points_test  <- structure_points[-c(random_sprayables),]
  }else if ( spray_ratio<0.2 ){
    random_sprayables <- sample(which(structure_points$sprayable==1),sum(structure_points$sprayable==1)*0.9)
    random_non_sprayables <- sample(which(structure_points$sprayable==0),sum(structure_points$sprayable==0)*0.5)
    structure_points <- structure_points[c(random_sprayables, random_non_sprayables),]
    structure_points_test  <- structure_points[-c(random_sprayables, random_non_sprayables),]
  }else if(spray_ratio>5){
      random_sprayables <- sample(which(structure_points$sprayable==1),sum(structure_points$sprayable==1)*0.5)
      random_non_sprayables <- sample(which(structure_points$sprayable==0),sum(structure_points$sprayable==0)*0.9)
      structure_points <- structure_points[c(random_sprayables, random_non_sprayables),]
      structure_points_test  <- structure_points[-c(random_sprayables, random_non_sprayables),]
    }else{
      random_sprayables <- sample(which(structure_points$sprayable==1),sum(structure_points$sprayable==1)*0.9)
      random_non_sprayables <- sample(which(structure_points$sprayable==0),sum(structure_points$sprayable==0)*0.9)
      structure_points <- structure_points[c(random_sprayables,random_non_sprayables),]
      structure_points_test <- structure_points[-c(random_sprayables,random_non_sprayables),]
    }

  
 
  cnt <- cnt + 1
  if (cnt < 2){
    all_structure_points <-  structure_points
    all_structure_points_test <-  structure_points_test
  }else{
    all_structure_points <- rbind(all_structure_points,structure_points)
    
    all_structure_points_test <- rbind(all_structure_points_test,structure_points_test)
  }
  #do.call(file.remove, list(list.files(paste0(getwd(),"/osm_data/",country,"_OSM"), full.names = TRUE)))
}
#all_structure_points$country <- as.factor(all_structure_points$country)
#all_structure_points_test$country <- as.factor(all_structure_points_test$country)
#size_split <- pretty(all_structure_points$size,30)
#all_structure_points$size_group <- cut(all_structure_points$size,size_split)
#all_structure_points_test$size_group <- cut(all_structure_points_test$size,size_split)
save(all_structure_points,file = paste0(getwd(),'/',"africa","_all_structure_points_201803.RData"))
save(all_structure_points_test,file = paste0(getwd(),'/',"africa","_all_structure_points_test_201803.RData"))
# 50,000
# min max of each country 
# increase by country 
rf_mod <- randomForest(as.factor(sprayable) ~ region +
                         # size_group + 
                         #dist_to_nearest_water + nearest_water_area + 
                          poly_area + 
                           dist_to_nearest + 
                           area_of_nearest + 
                           dist_to_nearest_road + 
                           n_poly +
                           poly_complexity +
                           poly_complexity_of_nearest +
                           n_poly_of_nearest,
                         #nodesize = 5,
                         ntree = 1000,
                         data = all_structure_points,family="binomial")
#save(rf_mod,file = paste0(getwd(),'/',country,"_rf_mod.RData"))
save(rf_mod,file = paste0(getwd(),'/',"africa","_rf_mod_201803.RData"))
auc_train <- AUC(rf_mod$votes[,2], all_structure_points$sprayable)
impotance_data <- varImpPlot(rf_mod)
#write.csv(impotance_data, file = paste0(getwd(),'/',country,"_impotance_data_OSM.csv"))
write.csv(impotance_data, file = paste0(getwd(),'/',"africa","_impotance_data_OSM_201803.csv"))
opt_cutoff_rf_mod <- optimal_cutoff(rf_mod$votes[,2], all_structure_points$sprayable)
write.csv(opt_cutoff_rf_mod, file = paste0(getwd(),'/',"africa","_opt_cutoff_rf_mod_OSM_201803.csv"))
test_data_predictions <- predict(rf_mod, all_structure_points_test, predict.all=T)
test_data_predictions_prob <- apply(test_data_predictions$individual, 1, function(x){mean(as.numeric(x))})
auc_test <- AUC(test_data_predictions_prob, all_structure_points_test$sprayable)
pred_test <- prediction(test_data_predictions_prob, all_structure_points_test$sprayable)
perf_test <- performance(pred_test,"tpr","fpr")
#jpeg(paste0(getwd(),'/',country,"_test_roc_plot.jpg")) 
dev.set()
dev.new()
jpeg(paste0(getwd(),'/',"africa","_test_roc_plot_201803.jpg")) 
plot(perf_test,col="black",lty=2, lwd=3)
legend(0.3,0.3,paste(c("AUC = "),auc_test,sep=""),border="white",cex=1.3,box.col = "white")
dev.off()
  
confusion_table <- table(as.numeric(test_data_predictions_prob>opt_cutoff_rf_mod$opt_cutoff), all_structure_points_test$sprayable)
output <- cbind(spray_ratio_list,auc_train,auc_test,matrix(confusion_table))
#write.csv(output, file = paste0(getwd(),'/',country,"_output.csv"))
write.csv(output, file = paste0(getwd(),'/',"africa_output_201803.csv"))
  #do.call(file.remove, list(list.files(paste0(getwd(),"/osm_data/",country,"_OSM"), full.names = TRUE)))
  
  #return(output)


# i=40
# n <- nafrica
# #nsub <- 5
# result <- foreach(
#   i = 1:n,
#   .packages = c("stringr","ROCR"), 
#   .combine = rbind,              
#   .verbose = TRUE) %dopar% {  
#     outputs <- single(links[i])
#     
#     }
# 
# results <- data.frame( result )
# countries <-sapply(links, function(x) gsub('-latest-free.shp.zip','',basename(x)))
# colnames(results) <- countries
# write.csv(results, file = paste0(getwd(),'/',"africa_output.csv"))
proc.time()


#train = load(paste0(getwd(),'/',"africa","_all_structure_points.RData"))
#test = load(paste0(getwd(),'/',"africa","_all_structure_points_test.RData"))
countries <-sapply(links, function(x) gsub('-latest-free.shp.zip','',basename(x)))
#for (i in c(39:41)){
for (i in c(1:nafrica)){
  structure_points <- all_structure_points[all_structure_points$country==countries[[i]],]
  structure_points_test <- all_structure_points_test[all_structure_points_test$country==countries[[i]],]
  #auc_train <- AUC(rf_mod$votes[,2], all_structure_points$sprayable)
  #impotance_data <- varImpPlot(rf_mod)
  #write.csv(impotance_data, file = paste0(getwd(),'/',country,"_impotance_data_OSM.csv"))
  #opt_cutoff_rf_mod <- optimal_cutoff(rf_mod$votes[,2], all_structure_points$sprayable)
  #write.csv(opt_cutoff_rf_mod, file = paste0( getwd(),'/',countries[[i]],"_opt_cutoff_rf_mod_OSM.csv") ) 
  train_data_predictions <- predict(rf_mod, structure_points, predict.all=T)
  train_data_predictions_prob <- apply(train_data_predictions$individual, 1, function(x){mean(as.numeric(x))})
  auc_train <- AUC(train_data_predictions_prob, structure_points$sprayable)
  pred_train <- prediction(train_data_predictions_prob, structure_points$sprayable)
  perf_train <- performance(pred_train,"tpr","fpr")
  
  
  test_data_predictions <- predict(rf_mod, structure_points_test, predict.all=T)
  test_data_predictions_prob <- apply(test_data_predictions$individual, 1, function(x){mean(as.numeric(x))})
  auc_test <- AUC(test_data_predictions_prob, structure_points_test$sprayable)
  pred_test <- prediction(test_data_predictions_prob, structure_points_test$sprayable)
  perf_test <- performance(pred_test,"tpr","fpr")
  dev.cur()
  dev.set()
  dev.new()
  jpeg(paste0(getwd(),'/',countries[[i]],"_test_roc_plot.jpg")) 
  plot(perf_test,col="black",lty=2, lwd=3)
  legend(0.3,0.3,paste(c("AUC = "),auc_test,sep=""),border="white",cex=1.3,box.col = "white")
  dev.off()
  dev.set()
  dev.new()
  jpeg(paste0(getwd(),'/',countries[[i]],"_train_roc_plot.jpg")) 
  plot(perf_train,col="black",lty=2, lwd=3)
  legend(0.3,0.3,paste(c("AUC = "),auc_train,sep=""),border="white",cex=1.3,box.col = "white")
  dev.off()
  confusion_table <- table(as.numeric(test_data_predictions_prob>opt_cutoff_rf_mod$opt_cutoff), structure_points_test$sprayable)  
  output <- c(auc_train,auc_test,matrix(confusion_table))
  write.csv(output, file = paste0(getwd(),'/',countries[[i]],"_output.csv"))
}
proc.time()

# library(osmdata)
# library(leaflet)

# get_building_given_polygon <- funciton()
# q0 <- opq(bbox = c( -122.41, 37.77, -122.42,37.78)) # Chiswick Eyot in London, U.K.
# q1 <- add_osm_feature(q0, key = 'building')
# x <- osmdata_sp(q1)
# leaflet() %>% addTiles() %>% addPolygons(data=x$osm_polygons, color="red")

