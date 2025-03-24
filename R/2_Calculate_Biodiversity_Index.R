#############################################################################
# Calculate Biodiversity Use of Greenspaces from Citizen Science Data
#############################################################################

library(readr)
library(raster)
library(sf)
library(dplyr)
library(terra)
library(data.table)
library(randomForest)
library(missForest)
library(caret)
library(ggplot2)
library(leaflet)
library(htmltools)

# Read in and prepare data -------------------------------------------------------

# load in iNaturalist data by greenspace data created from the the "1_Summarize_iNaturalist_Data_by_Greenspace"
inat <- read_csv("Data/iNaturalist_park.csv")

# load the shapefile that contains Broward greenspaces
parks <- st_read("Data/shapefiles/urban_parks.shp", crs = 4326)

# load in raster files exported from Google Earth Engine, which will be used as covariates in the model
treecover <- raster("Data/rasters/treecover.TIF")
water <- raster("Data/rasters/water.TIF")
impervious <- projectRaster(raster("Data/rasters/impervious.tif"), crs="+proj=longlat +datum=WGS84")
veg <- raster("Data/rasters/Percent_NonTree_Vegetation.TIF", band=2) # band 2 has the non-vegetated data

# Calculate average values per raster object -------------------------------------------------------
## Tree Cover ---------------------------------------------------------------------------------------

# extract raster values, calculate the mean for each value, then add the data to the park shapefile
tree.values <- raster::extract(treecover$Percent_Tree_Cover, parks)

# calculate the mean for each value and add this to the parks files
parks$treecover <- unlist(lapply(tree.values, FUN=mean, na.rm=TRUE))

# tree cover is the average percentage of pixel which is covered by tree canopy for each park

## non-tree vegetation ---------------------------------------------------------------------------------------

# extract raster values
veg.values <- raster::extract(veg, parks)

# calculate the mean for each value
parks$vegetation <- unlist(lapply(veg.values, FUN=mean, na.rm=TRUE))

# non-tree vegetation is the average percentage of pixel which is covered by non-tree vegetation for each park

## water ---------------------------------------------------------------------------------------

# extract raster values
water.values <- raster::extract(water, parks)

# the water raster is coded 1=land, 2=water, 4=snow/ice, 200=cloud shadow, 201=cloud. We are 
# only interested in 2 here. I manually checked the raster and found that in our study area, cloud 
# shadows and clouds did not block water features, so we can ignore those pixels. 
# We need to calculate the percentage of pixels that contain water compared to total pixels in shapefile.
parks$water <- unlist(lapply(water.values, function(x){length(which(x %in% c(2)))/length(x)}))

# water is the percentage of water coverage in each park

## impervious ---------------------------------------------------------------------------------------

# extract raster values
impervious.values <- raster::extract(impervious, parks)

# calculate the mean for each value
parks$impervious <- unlist(lapply(impervious.values, FUN=mean, na.rm=TRUE))

# impervious is the average percentage of pixel which is covered by developed impervious surface for each park

# Predict bio-use using a random forest model -------------------------------------------------------

# combine data from parks shapefile with inat data frame
inat.pred <- right_join(inat, parks, by="poly_id")

# pull the data of interest and condense parks with more than one row of data
inat.pred <- inat.pred %>% 
                dplyr::select(poly_id, species_count_observations, park_size_ha, total_observations,
                              total_observers, treecover, water, impervious, vegetation) %>% 
                group_by(poly_id) %>%
                dplyr::summarise(richness=first(species_count_observations), park_size_ha=sum(park_size_ha),
                   total_observations=first(total_observations), total_number_of_observers=first(total_observers),
                   treecover=mean(treecover), water=sum(water), impervious=mean(impervious), 
                   vegetation=mean(vegetation))

# Random Forest  -----------------------------------------

# pull relevant data
rf_data <- inat.pred %>%
  dplyr::select(treecover, water, impervious, vegetation, 
         richness, total_observations, total_number_of_observers) %>%
  filter(complete.cases(richness))

# log transform the response to make it more normally distributed
rf_data$richness <- log(rf_data$richness)

# test and train data
set.seed(23) # set seed for reproducibility
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
                ntree = 300,
                maxnodes = 26)

# evaluate model
prediction <- predict(fit_rf, test)

# compare model predictions to test data
df <- data.frame(pred = prediction, real=test$richness)

varImp(fit_rf)

# let's check model performance
par(mar=c(5,5,2,2))
plot(prediction, test$richness, xlab="Observed Richness", ylab="Predicted Richness",
     cex=2, lwd=2, cex.lab=2, cex.axis=2)
abline(0,1, lwd=2)

# function to calculate R2
RSQUARE = function(y_actual,y_predict){
  cor(y_actual,y_predict)^2
}

# R2 value
RSQUARE(test$richness, as.numeric(prediction))

mean((prediction-test$richness)^2)

# now run the random forest on the entire data set without including park size
fit_rf_nps <- randomForest(richness~ treecover + water + impervious + vegetation + 
                             total_observations + total_number_of_observers,
                           rf_data,
                           method = "rf",
                           metric = "RMSE")

varImpPlot(fit_rf_nps)

# add poly_id to rf_data
rf_data$poly_id <- inat.pred[complete.cases(inat.pred$richness),]$poly_id

## Predict richness for set number of observations and observers -------------------------------------------------------

# scale the data, so observations are 1,000. We will leave the values for treecover, water, impervious, and 
# vegetation alone since those are percentage of area covered by each habitat type which will not change
scaling_value <-  1000 / rf_data$total_observations
scaled <- rf_data %>%
  mutate(total_observations = total_observations * scaling_value,
         total_number_of_observers = total_number_of_observers * scaling_value) %>%
  dplyr::select(-poly_id)

# predict new richness values from set observations and observers
scaled$predicted.richness <- predict(fit_rf_nps, scaled)

# add poly_id back
scaled$poly_id <- rf_data$poly_id

## Predict missing variables -------------------------------------------------------------------------------------------

# add in missing data back into the data frame. Multiply park size by 1000 since we added observation number as 1000
all <- inat.pred %>% 
           dplyr::filter(is.na(richness)) %>%
           dplyr::mutate(predicted.richness=NA, 
                         total_observations=1000, 
                         park_size_ha=log(exp(park_size_ha)*1000)) %>% 
           dplyr::select(-richness, -total_observations, -total_number_of_observers,
                          -park_size_ha)

# prepare the data for missForest
mf_data <- bind_rows(scaled, all) %>%
  dplyr::select(-poly_id) %>%
  as.matrix()

# run missForest
mf_calc <- missForest(mf_data, verbose=TRUE)

# create data frame and poly_id back in
mf_pred <- as.data.frame(mf_calc$ximp)
mf_pred$poly_id <- bind_rows(scaled, all)$poly_id

# add true/false if data was predicted during this step to visualize the data
mf_pred$predicted.data <- mf_pred$poly_id %in% unique(all$poly_id)

# back transform predicted.richness
mf_pred$predicted.richness <- exp(mf_pred$predicted.richness) 

# create data frame for parks with iNaturalist data (scaled) and parks without iNaturalist data (predicted)
scaled <- mf_pred[mf_pred$predicted.data==FALSE,]
predicted <- mf_pred[mf_pred$predicted.data==TRUE,]

# Save plots to a file
png("Figures/Supplemental_Figures/FigureA1.png", width=1200, height=1200, res=150)
par(mfrow=c(2,2), cex=1.3,mar=c(4, 4, 2.5, 2.5))
# let's make plots to see how well our data predictions match the real data predictions
# Plot 1: Tree Cover vs Predicted Richness
plot(predicted.richness ~ treecover, data=scaled, ylab="Species Richness", xlab="Tree Cover")
points(predicted.richness ~ treecover, col="red", data=predicted)
mtext("a", side=3, line=1, at=par("usr")[1], adj=0, cex=1.5) # Left align title

# Plot 2: Non-Tree Vegetation Cover vs Predicted Richness
plot(predicted.richness ~ vegetation, data=scaled, ylab="Species Richness", xlab="Non-Tree Vegetation Cover")
points(predicted.richness ~ vegetation, col="red", data=predicted)
mtext("b", side=3, line=1, at=par("usr")[1], adj=0, cex=1.5)

# Plot 3: Water vs Predicted Richness
plot(predicted.richness ~ water, data=scaled, ylab="Species Richness", xlab="Water")
points(predicted.richness ~ water, col="red", data=predicted)
mtext("c", side=3, line=1, at=par("usr")[1], adj=0, cex=1.5)

# Plot 4: Impervious Cover vs Predicted Richness
plot(predicted.richness ~ impervious, data=scaled, ylab="Species Richness", xlab="Impervious Cover")
points(predicted.richness ~ impervious, col="red", data=predicted)
mtext("d", side=3, line=1, at=par("usr")[1], adj=0, cex=1.5)

dev.off()

## Leave one out to test miss forest model -------------------------------

# To test the robustness of the miss forest model, we will run a leave-one-out cross
# validation test
test_data <- bind_rows(scaled, all) %>%
  dplyr::select(-poly_id, -predicted.data) %>%
  filter(complete.cases(predicted.richness))

loo <- data.frame(predicted.richness=test_data$predicted.richness,
                  missForest.richness="")
for (i in 1:nrow(loo)){
  test_data <- as.data.frame(test_data)
  loo.richness <- ifelse(test_data$predicted.richness==test_data$predicted.richness[i], NA, test_data$predicted.richness)
  test_data2 <- test_data
  test_data2$predicted.richness <- loo.richness
  test_data2 <- as.matrix(test_data2)
  mf_test <- missForest(test_data2, verbose=TRUE)
  loo$missForest.richness[i] <- as.numeric(mf_test$ximp[i,8])
}

# convert the miss Forest predicted richness to numeric
loo$missForest.richness <- as.numeric(loo$missForest.richness)

# plot the relationship
par(mar=c(5,5,2,2))
plot(loo$missForest.richness, loo$predicted.richness, xlab="Scaled Richness", ylab="Predicted Scaled Richness", 
     cex=2, lwd=2, cex.lab=2, cex.axis=2)
abline(0,1, lwd=2)

# calculate R2
RSQUARE(loo$predicted.richness, loo$missForest.richness)

# Create Output Table from the Results ---------------------------------------------------------

# create output table with predicted species richness
# add real number of observations and number of observers back into the data
outputdata <- left_join(mf_pred %>% dplyr::select(predicted.richness, poly_id), 
                        inat, 
                        by="poly_id")
outputdata$total_observations[is.na(outputdata$total_observations)] <- 0
outputdata$total_observers[is.na(outputdata$total_observers)] <- 0

output <- data.frame(poly_id = outputdata$poly_id, predicted_richness = outputdata$predicted.richness, 
                     iNat_observations = outputdata$total_observations, iNat_users = outputdata$total_observers,
                     park_size=outputdata$park_size_ha)

write.csv(output, "Data/bio_use_predictors/bio_use_predictions.csv")

# create a supplemental table, summarizing this data
output_s <- output %>%
  separate(poly_id, into=c("municipality", "park name"), sep="_") %>%
  mutate(
    predicted_richness = round(predicted_richness, 2),
    park_size = round(park_size, 2)
  ) %>%
  dplyr::select(municipality, `park name`, park_size, iNat_observations, iNat_users, predicted_richness)

colnames(output_s) <- c("Municipality", "Park Name", "Park Size (ha)", 
                        "Number of iNaturalist Observations", "Number of iNaturalist Users",
                        "Predicted Richness")

write_csv(output_s, "Data/bio_use_predictors/data.summary.csv")
