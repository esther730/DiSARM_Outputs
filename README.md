# MASS Model 

## Walkthrough

  ### Querying and Preprocessing
  The input can be three formats: two longtitudes/latitudes, multiple longtitude/latitudes(matrix), a given region name(e.g.country name)
   ```sh
  # type 1
  long,lat
  30.5	-2.35
  30.0	-2.00
  
  # type 2
  coordinates <- read.csv("coordinates.csv",header=FALSE)
  
  # type 3
  bb <-  getbb('RWA')
   ```
  
  - Given a polygon, preprocess it into a rectangle
  ```sh
  # type 1
  q0 <- opq(bbox = c( min(long1,long2), min(lat1,lat2), max(long1,long2),max(lat1,lat2)))
  
  # type 2
  q0 <- opq(bbox =as.matrix(coordinates))
  
  # type 3
  q0 <- opq(bbox = bb)
  ```
  
  - Get osm labels(“building”, “highway”) through Overpass query
  ```sh
  q1 <- add_osm_feature(q0, key = 'building')
  q2 <- add_osm_feature(q0, key = 'highway')
  x <- osmdata_sp(q1)
  y <- osmdata_sp(q2)
  ```
  
   ### Feature Generation
  - **region**  
  Calssify a building into different regions of Africa(six categorical variables: North, Middle, Western, South, Eastern, Other)
  ```sh
  region <- rep(as.factor(levels(region_out))[which.max(table(region_out))],length(structures_poly_selected ))
  ```
  
  - **poly_area**  
  For each builidng, calculate its area (choose the largest one if a building has multiple polygons)  
  ```sh
  building_area <- sapply(slot(structures_poly_selected, "polygons"), slot, "area")
  poly_area<-sapply(building_area, max)
  ```
  
  - **dist_to_nearest**  
  For each building, calculate its distant to nearest k=5 neighboring buildings   
  ```sh
  building_area <- sapply(slot(structures_poly_selected, "polygons"), slot, "area")
  dist_to_nearest <- nn2(building_corrds,building_corrds, eps=0,searchtype="priority", k=5)$nn.dists[,2]
  ```
  
  - **area_of_nearest**
  For each builidng, calculate the avaerage areas of its k=2 neighbours
  ```sh
  area_of_nearest <- building_area[nn2(building_corrds,building_corrds, eps=0,searchtype="priority", k=2)$nn.idx[,2]]
  ```
  
  - **dist_to_nearest_road**  
  For each builidng, calculate distance between it and nearest k=2 roads
  ```sh
  get_coords <- function(x){x@Lines} 
  road_coords <- sapply(y$osm_lines@lines, get_coords)
  road_coords <- sapply(road_coords, coordinates)
  road_coords_all <- do.call("rbind", road_coords)
  dist_to_nearest_road <- nn2(road_coords_all, building_corrds, eps=0,searchtype="priority", k=2)$nn.dists[,2]
  ```
  
  - **n_poly**  
  Number of polygons a building has (usually 1, some buildings have multiple polygons)   
  ```sh
  n_poly <- sapply(areas, length)
  ```
  
  - **poly_complexity**  
  The Length of the corrdinate list of the building
  ```sh
  poly_complexity <- sapply(poly_coords, function(x){length(unlist(x))})
  ```
  
  - **poly_complexity_of_nearest**   
  The Length of the corrdinate list of the nearest k=2 neighboring buildings
  ```sh
  poly_complexity_of_nearest <- poly_complexity[nn2(building_corrds,building_corrds, eps=0,searchtype="priority", k=2)$nn.idx[,2]]
  ```
  
  - **n_poly_of_nearest**  
  Number of polygons the nearest k =2 buildings have (usually 1, some buildings have multiple polygons)
  ```sh
  n_poly_of_nearest <- n_poly[nn2(building_corrds,building_corrds, eps=0,searchtype="priority", k=2)$nn.idx[,2]]
  ```
  
  - **lng**   
  Building Longtitude
  ```sh
  lng <- building_corrds[,1]
  ```
  
  - **lat**   
  Building Latitude
  ```sh
  lat <- building_corrds[,2]
  ```
 
  
  ### Training
  - training size for each country was controled between 500 to 10,000 by resampling and random sampling without replacement.
  - 80% of data(319,237 records) will be used for training and the rest part(63,782 records) for testing. 
  
  - training uses RandomForest   
  ```sh
  rf_mod <- randomForest(as.factor(sprayable) ~ region +
                         poly_area + 
                         dist_to_nearest + 
                         area_of_nearest + 
                         dist_to_nearest_road + 
                         n_poly +
                         poly_complexity +
                         poly_complexity_of_nearest +
                         n_poly_of_nearest,
                         ntree = 1000,
                         data = all_structure_points,family="binomial")
  ```
  
  
  ### Prediction
  - Predict with opt_cutoff = 0.194 (calculated by optimal_cutoff function)
  ```sh
 est_data_predictions <- predict(rf_mod, all_structure_points_test, predict.all=T)
 test_data_predictions_prob <- apply(test_data_predictions$individual, 1, function(x){mean(as.numeric(x))})
 opt_cutoff = 0.194
 test_data_predictions_output <- as.numeric(test_data_predictions_prob>opt_cutoff)
 final_result = cbind(all_structure_points_test,test_data_predictions[[1]])
  ```
  
  - Evaluate  
    For a input polygon, the function will identify where the polygon locates, and return the statistics for the corrsponding       countries. Current metrics included: Confusion matrix measures(accuracy, precision, sensitivity, specificity, prevalence), PPV, FDR, NPV, FOR, NER, F1 Score, Cohen Kappa, MCC, AUC

  
  - Output  
  [in plan] Outputs a csv file of a matrix of columns(IDs, lng, lat, predictions) and a html/jpg map with building types on it. (Green points mean non-sprayable, and red ones mean spreayable.)


 
## Examples
  - For example, if the input polygon is located in Rwanda, the outputs will be the following three parts: 
  1. Prediciton:
  ```sh
  ID	        lng	          lat	      Sprayable Prediction
  106059259	  30.10411	   -2.027780	1          1
  530254380	  30.09169	   -2.143496	1          0
  ```
  
  2. Statisitics:
  accuracy: 0.9899875  
  precision_value: 0.9612591  
  sensitivity: 1  
  specificity: 0.9866778  
  prevalence:	0.2484355  
  PPV:	0.9612591  
  FDR:	0.03874092  
  NPV:	1  
  FOR:	0  
  NER:	0.2484355  
  F1 score:	0.9802469  
  cohen_kappa:	0.9735446  
  MCC:	0.9738855  
  AUC:	0.9933389 
  ```sh
  |              | 0              | 1             |
  | 0            | 1185           | 16            |
  | 1            | 0              | 397           |
  ```
  3. Visulization: (Optional)
  
  ![Screenshot](https://github.com/locational/urap-maas-models/blob/master/models/sample_output/precition_graph.png)

  
