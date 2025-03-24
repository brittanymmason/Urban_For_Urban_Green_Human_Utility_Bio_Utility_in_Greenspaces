###########################################################################
# Alternative Calculations for Biodiversity Index
###########################################################################

library(readr)
library(raster)
library(sf)
library(dplyr)
library(data.table)
library(randomForest)
library(missForest)
library(caret)
library(ggplot2)
library(htmltools)
library(scales)
library(mgcv)

# Read in and prepare data -------------------------------------------------------

# load in iNaturalist data by greenspace data created from the the "1_Summarize_iNaturalist_Data_by_Greenspace"
inat <- read_csv("Data/iNaturalist_park.csv")

# load park shapefiles
parks <- st_read("Data/shapefiles/urban_parks.shp", crs = 4326)

# combine files
parks_inat <- left_join(parks %>% select(-mncplty, -grns__2), 
                        inat %>% select(-geometry), 
                        by="poly_id")

# load in raster files exported from Google Earth Engine, which will be used as covariates in the models
treecover <- raster("Data/rasters/treecover.TIF")
water <- raster("Data/rasters/water.TIF")
impervious <- projectRaster(raster("Data/rasters/impervious.tif"), crs="+proj=longlat +datum=WGS84")
veg <- raster("Data/rasters/Percent_NonTree_Vegetation.TIF", band=2) # band 2 has the non-vegetated data

## Calculate average values per raster object -------------------------------------------------------
## Tree Cover ---------------------------------------------------------------------------------------

# extract raster values, calculate the mean for each value, then add the data to the park shapefile
tree.values <- raster::extract(treecover$Percent_Tree_Cover, parks_inat)

# calculate the mean for each value and add this to the parks files
parks_inat$treecover <- unlist(lapply(tree.values, FUN=mean, na.rm=TRUE))

# tree cover is the average percentage of pixel which is covered by tree canopy for each park

## non-tree vegetation ---------------------------------------------------------------------------------------

# extract raster values
veg.values <- raster::extract(veg, parks_inat)

# calculate the mean for each value
parks_inat$vegetation <- unlist(lapply(veg.values, FUN=mean, na.rm=TRUE))

# non-tree vegetation is the average percentage of pixel which is covered by non-tree vegetation for each park

## water ---------------------------------------------------------------------------------------

# extract raster values
water.values <- raster::extract(water, parks_inat)

# the water raster is coded 1=land, 2=water, 4=snow/ice, 200=cloud shadow, 201=cloud. We are 
# only interested in 2 here. I manually checked the raster and found that in our study area, cloud 
# shadows and clouds did not block water features, so we can ignore those pixels. 
# We need to calculate the percentage of pixels that contain water compared to total pixels in shapefile.
parks_inat$water <- unlist(lapply(water.values, function(x){length(which(x %in% c(2)))/length(x)}))

# water is the percentage of water coverage in each park

## impervious ---------------------------------------------------------------------------------------

# extract raster values
impervious.values <- raster::extract(impervious, parks_inat)

# calculate the mean for each value
parks_inat$impervious <- unlist(lapply(impervious.values, FUN=mean, na.rm=TRUE))

# impervious is the average percentage of pixel which is covered by developed impervious surface for each park

# Method 2: Random Forest with All Real Data -------------------------------------------------------

# pull relevant data
real_data <- parks_inat %>%
                as.data.frame() %>%
                select(poly_id, species_count_observations, total_observations, total_observers,
                       treecover, water, impervious, vegetation) %>%
                rename(richness=species_count_observations) %>%
                filter(complete.cases(richness))

# log transform the response and remove poly_id
rf_data <- real_data %>% 
                mutate(richness=log(richness)) %>%
                select(-poly_id)

# set testing and training data
set.seed(23)
samp <- sample(nrow(rf_data), 0.8*nrow(rf_data))
train <- rf_data[samp,]
test <- rf_data[-samp,]

control <- trainControl(method = "cv", number = 10, search ="grid")
rf_train <- train(richness~., data=train, method="rf", trControl=control)
print(train)

# find best mtry
tuneGrid <- expand.grid(.mtry=c(1:10))
rf_mtry <- train(richness~., data=train, method="rf", tuneGrid=tuneGrid,
                 trControl=control, importance=TRUE, nodesize=14, ntree=300)
print(rf_mtry)

best_mtry <- rf_mtry$bestTune$mtry

# find best maxnodes
store_maxnode <- list()
tuneGrid <- expand.grid(.mtry = best_mtry)
for (maxnodes in c(30:10)) {
  set.seed(23)
  rf_maxnode <- train(richness~.,
                      data = train,
                      method = "rf",
                      tuneGrid = tuneGrid,
                      trControl = control,
                      importance = TRUE,
                      nodesize = 14,
                      maxnodes = maxnodes,
                      ntree = 300)
  current_iteration <- toString(maxnodes)
  store_maxnode[[current_iteration]] <- rf_maxnode
}
results_mtry <- resamples(store_maxnode)
summary(results_mtry)

# now find the best ntrees
store_maxtrees <- list()
for (ntree in c(250, 300, 350, 400, 450, 500, 550, 600, 800, 1000, 2000)) {
  set.seed(23)
  rf_maxtrees <- train(richness~.,
                       data = train,
                       method = "rf",
                       metric = "RMSE",
                       tuneGrid = tuneGrid,
                       trControl = control,
                       importance = TRUE,
                       nodesize = 14,
                       maxnodes = 24,
                       ntree = ntree)
  key <- toString(ntree)
  store_maxtrees[[key]] <- rf_maxtrees
}
results_tree <- resamples(store_maxtrees)
summary(results_tree)

# run model
fit_rf <- train(richness~.,
                train,
                method = "rf",
                metric = "RMSE",
                tuneGrid = tuneGrid,
                trControl = control,
                importance = TRUE,
                nodesize = 14,
                ntree = 350,
                maxnodes = 26)

# evaluate model
prediction <- predict(fit_rf, test)

df <- data.frame(pred = prediction, real=test$richness)

varImp(fit_rf)

# let's check model performance
par(mar=c(5,5,2,2))
plot(prediction, test$richness, xlab="Observed Richness", ylab="Predicted Richness",
     cex=2, lwd=2, cex.lab=2, cex.axis=2)
abline(0,1, lwd=2)

# R2 value
# function to calculate R2
RSQUARE = function(y_actual,y_predict){
  cor(y_actual,y_predict)^2
}

RSQUARE(test$richness, as.numeric(prediction))

mean((prediction-test$richness)^2)

# now run the random forest on the entire data set
fit_rf <- randomForest(richness~ treecover + water + impervious + vegetation + 
                         total_observations + total_observers,
                       rf_data,
                       method = "rf",
                       metric = "RMSE")

varImpPlot(fit_rf)


# add poly_id to rf_data
rf_data$poly_id <- parks_inat[complete.cases(parks_inat$species_count_observations),]$poly_id

## Predict richness for set number of observations and observers -------------------------------------------------------

# scale the data, so observations are 1,000. We will leave the values for treecover, water, impervious, and 
# vegetation alone since those are percentage of area covered by each habitat type which will not change
scaling_value <-  1000 / rf_data$total_observations
scaled <- rf_data %>%
  mutate(total_observations = total_observations * scaling_value,
         total_observers = total_observers * scaling_value) %>%
  select(-richness)

# predict new richness values from set observations and observers
prediction <- predict(fit_rf, scaled)

# add to scaled data frame
scaled$predicted.richness <- prediction

# add poly_id back
scaled$poly_id <- real_data$poly_id

# to create output table, add real number of observations and number of observers back into the data
outputdata <- left_join(scaled %>% select(poly_id, predicted.richness), 
                        inat, 
                        by="poly_id")

# change NA values for total observations and total observers to 0 and clean up data
outputdata <- outputdata %>%
  as.data.frame() %>%
  mutate(
    total_observations = if_else(is.na(total_observations), 0, total_observations),
    total_observers = if_else(is.na(total_observers), 0, total_observers)
  ) %>%
  select(poly_id, predicted.richness, total_observations, total_observers) %>%
  mutate(predicted_richness_method2=exp(outputdata$predicted.richness))
  
# add observed richness
output_method2 <- inner_join(outputdata, 
                      real_data %>% select(poly_id, richness), 
                      by="poly_id")

# Method 3: GAM model of all real data ------------------------------------------

# log transform the response variable
real_data <- real_data %>%
                mutate(log_richness=log(richness))

# fit with GAM
richness_obs <- gam(log_richness ~ s(total_observations, bs="cr", k=30) +
                             s(total_observers, bs="cr") + s(treecover, bs="cr") + 
                             s(water, bs="cr", k=20) +
                             s(impervious, bs="cr") + s(vegetation, bs="cr"), 
                           method="REML",
                           data = real_data)
summary(richness_obs)
gam.check(richness_obs)

# scale the values based on 1,000 total observations
scaling_value <- 1000/real_data$total_observations
scaled <- real_data %>%
  dplyr::mutate(total_observations=total_observations*scaling_value,
                total_observers=total_observers*scaling_value) %>%
  select(poly_id, total_observations, total_observers, treecover, water, impervious, 
         vegetation)

# predict the richness based on the scaled data
scaled$predicted_richness_gam <- predict(richness_obs, newdata = scaled, type = "response")

# clean up the data
output_method3 <- scaled %>%
  dplyr::select(predicted_richness_gam) %>%
  dplyr::mutate(predicted_richness_method3=exp(predicted_richness_gam),
                poly_id=real_data$poly_id) %>%
  dplyr::select(-predicted_richness_gam)

# Method 4: Random Forest for parks with greater than 50 iNat observations -------------------------------------------------------

# pull relevant data
real_data_over50 <- parks_inat %>%
  as.data.frame() %>%
  select(poly_id, species_count_observations, total_observations, 
         total_observers,treecover, water, impervious, vegetation) %>%
  filter(total_observations > 50) %>%
  rename(richness=species_count_observations)

# prepare data for random forest model
rf_data50 <- real_data_over50 %>%
                mutate(richness=log(richness)) %>%
                select(-poly_id)

# test and train data
set.seed(23)
samp <- sample(nrow(rf_data50), 0.8*nrow(rf_data50))
train <- rf_data50[samp,]
test <- rf_data50[-samp,]

control <- trainControl(method = "cv", number = 10, search ="grid")
rf_train <- train(richness~., data=train, method="rf", trControl=control)
print(train)

# find best mtry
tuneGrid <- expand.grid(.mtry=c(1:10))
rf_mtry <- train(richness~., data=train, method="rf", tuneGrid=tuneGrid,
                 trControl=control, importance=TRUE, nodesize=14, ntree=300)
print(rf_mtry)

best_mtry <- rf_mtry$bestTune$mtry

# find best maxnodes
store_maxnode <- list()
tuneGrid <- expand.grid(.mtry = best_mtry)
for (maxnodes in c(30:10)) {
  set.seed(23)
  rf_maxnode <- train(richness~.,
                      data = train,
                      method = "rf",
                      tuneGrid = tuneGrid,
                      trControl = control,
                      importance = TRUE,
                      nodesize = 14,
                      maxnodes = maxnodes,
                      ntree = 300)
  current_iteration <- toString(maxnodes)
  store_maxnode[[current_iteration]] <- rf_maxnode
}
results_mtry <- resamples(store_maxnode)
summary(results_mtry)

# now find the best ntrees
store_maxtrees <- list()
for (ntree in c(250, 300, 350, 400, 450, 500, 550, 600, 800, 1000, 2000)) {
  set.seed(23)
  rf_maxtrees <- train(richness~.,
                       data = train,
                       method = "rf",
                       metric = "RMSE",
                       tuneGrid = tuneGrid,
                       trControl = control,
                       importance = TRUE,
                       nodesize = 14,
                       maxnodes = 24,
                       ntree = ntree)
  key <- toString(ntree)
  store_maxtrees[[key]] <- rf_maxtrees
}
results_tree <- resamples(store_maxtrees)
summary(results_tree)

# run model
fit_rf <- train(richness~.,
                train,
                method = "rf",
                metric = "RMSE",
                tuneGrid = tuneGrid,
                trControl = control,
                importance = TRUE,
                nodesize = 14,
                ntree = 550,
                maxnodes = 24)

# evaluate model
prediction <- predict(fit_rf, test)

df <- data.frame(pred = prediction, real=test$richness)

varImp(fit_rf)

# let's check model performance
par(mar=c(5,5,2,2))
plot(prediction, test$richness, xlab="Observed Richness", ylab="Predicted Richness",
     cex=2, lwd=2, cex.lab=2, cex.axis=2)
abline(0,1, lwd=2)

# R2 value
RSQUARE(test$richness, as.numeric(prediction))

mean((prediction-test$richness)^2)

# now run the random forest on the entire data set
fit_rf <- randomForest(richness~ treecover + water + impervious + vegetation + 
                         total_observations + total_observers,
                       rf_data50,
                       method = "rf",
                       metric = "RMSE")

varImpPlot(fit_rf)


# add poly_id to rf_data
rf_data50$poly_id <- real_data_over50$poly_id

## Predict richness for set number of observations and observers -------------------------------------------------------

# scale the data, so observations are 1,000. We will leave the values for treecover, water, impervious, and vegetation alone since those
# are percentage of area covered by each habitat type which will not change
scaling_value <- 1000/rf_data50$total_observations
scaled <- rf_data50 %>%
  dplyr::mutate(total_observations=total_observations*scaling_value,
                total_observers=total_observers*scaling_value)

# remove columns we do not want included in the calculation
pred_data <- scaled %>% dplyr::select(-richness,  -poly_id)

# predict new richness values from set observations and observers
prediction <- predict(fit_rf, pred_data)

# add to scaled data frame
scaled$predicted.richness <- prediction

# create output table
output_method4 <- scaled %>%
  select(poly_id, predicted.richness) %>%
  mutate(predicted_richness_method4=exp(scaled$predicted.richness),
         poly_id=real_data_over50$poly_id) %>%
  select(-predicted.richness)

# Method 5: Model the data to get richness ------------------------------------------

# log transform species richness
real_data_over50$log_richness <- log(real_data_over50$richness)

# GAM
richness_obs <- gam(log_richness ~ s(total_observations, bs="cr", k=10) +
                      s(total_observers, bs="cr") + s(treecover, bs="cr") + 
                      s(water, bs="cr", k=20) +
                      s(impervious, bs="cr") + s(vegetation, bs="cr"), 
                    method="REML",
                    data = real_data_over50)
summary(richness_obs)
gam.check(richness_obs)

scaling_value <- 1000/real_data_over50$total_observations
scaled <- real_data_over50 %>%
  dplyr::mutate(total_observations=total_observations*scaling_value,
                total__observers=total_observers*scaling_value) %>%
  select(poly_id, total_observations, total_observers, treecover, water, impervious, 
         vegetation)

scaled$predicted_richness_gam_ho <- predict(richness_obs, newdata = scaled[,2:length(scaled)], type = "response")

output_method5 <- scaled %>%
  dplyr::select(poly_id, predicted_richness_gam_ho) %>%
  dplyr::mutate(predicted_richness_method5=exp(predicted_richness_gam_ho)) %>%
  dplyr::select(-predicted_richness_gam_ho)

# Combine the data and prepare for export ---------------------------------

output_combined <- output_method2 %>%
  left_join(output_method3, by = "poly_id") %>%
  left_join(output_method4, by = "poly_id") %>%
  left_join(output_method5, by = "poly_id") %>%
  rename(predicted_richness_method1=predicted.richness,
         species_count=richness) %>%
  select(poly_id, total_observations, total_observers, species_count,
         predicted_richness_method1, predicted_richness_method2,
         predicted_richness_method3, predicted_richness_method4, 
         predicted_richness_method5) 

write.csv(output_combined, "Data/bio_use_predictors/alternative_bio_utility_measures.csv")
