# script is the main analysis script

# packages
library(tidyverse)
library(scales)
library(sf)
library(ggpubr)
library(margins)
library(jtools)
library(corrgram)
library(lares)
library(corrplot)
library(cluster)
library(Rtsne)
library(ggeffects)

# read in human utility file
human_utility <- read_csv("Data/human_utility.csv") %>%
  mutate(total_utility=rowSums(.[3:10])) %>%
  unite(poly_id, Municipality, `Park name`)

# read in file produced in 2_Calculate_Biodiversity_Index that contains the predicted richness value
biodiversity_utility <- read_csv("Data/bio_use_predictors/alternative_bio_utility_measures.csv")

# read in urban park data
greenspaces <- st_read("Data/shapefiles/urban_parks.shp") %>%
  mutate(greenspace_area_m2=as.numeric(st_area(.)))

# now merge data into one dataset
analysis_data <- greenspaces %>%
  st_set_geometry(NULL) %>%
  left_join(., biodiversity_utility, by="poly_id") %>%
  left_join(., human_utility, by="poly_id") %>%
  mutate(scaled_human_utility=scales::rescale(total_utility)) %>%
  mutate(predicted_richness_method1=scales::rescale(predicted_richness_method1)) %>%
  mutate(total_observers=scales::rescale(total_observers)) %>%
  mutate(predicted_richness_method2=scales::rescale(predicted_richness_method2)) %>%
  mutate(predicted_richness_method3=scales::rescale(predicted_richness_method3)) %>%
  mutate(predicted_richness_method4=scales::rescale(predicted_richness_method4)) %>%
  mutate(predicted_richness_method5=scales::rescale(predicted_richness_method5))

# Biodiversity Predicted  -----------------------------------

## Method 1  -----------------------------------

# filter the data 
method1_data <- analysis_data %>%
                    filter(complete.cases(predicted_richness_method1))

# calculate predicted biodiversity and human utility trend
model <- lm(predicted_richness_method1 ~ scaled_human_utility + log(greenspace_area_m2),
             data=method1_data, family="gaussian")
summary(model)

# Calculate R-squared
rss <- sum(model$residuals^2)
tss <- sum((method1_data$predicted_richness_method1 - mean(method1_data$predicted_richness_method1, na.rm=TRUE))^2)
r_squared <- 1 - (rss/tss)

# Number of observations and predictors
n <- nrow(analysis_data)
p <- length(model$coefficients) - 1

# Calculate adjusted R-squared
adj_r_squared <- 1 - ((1 - r_squared) * (n - 1) / (n - p - 1))

# Print adjusted R-squared
print(adj_r_squared)

# get slope values
slopes_scaled_human_utility <- ggpredict(model, terms = "scaled_human_utility")
df_slopes_scaled_human_utility1 <- as.data.frame(slopes_scaled_human_utility)

# plot the data
ggplot() +
  geom_point(data = method1_data, aes(y = predicted_richness_method1, x = scaled_human_utility, color = log10(greenspace_area_m2)), size = 3) +
  geom_line(data = df_slopes_scaled_human_utility1, aes(x = x, y = predicted), color = "blue") +
  geom_ribbon(data = df_slopes_scaled_human_utility1, aes(x = x, ymin = conf.low, ymax = conf.high), fill = "blue", alpha = 0.2) +
  xlab("Human Utility") +
  ylab("Biodiversity") +
  ggtitle("Method 1") +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  scale_colour_gradientn(colours = c("#D55E00", "#E69F00", "#F0E442", "#009E74", "#56B4E9", "#0072B2")) +
  labs(colour = expression(Log(Greenspace ~ Area ~ (m^2)))) +
  theme(legend.position = "bottom", legend.key.width = unit(0.6, "inch"))

## Method 2  -----------------------------------

# filter the data 
method2_data <- analysis_data %>%
  filter(complete.cases(predicted_richness_method2))

# calculate predicted biodiversity and human utility trend
model <- lm(predicted_richness_method2 ~ scaled_human_utility + log(greenspace_area_m2),
            data=method2_data, family="gaussian")
summary(model)

# Calculate R-squared
rss <- sum(model$residuals^2)
tss <- sum((method2_data$predicted_richness_method2 - mean(method2_data$predicted_richness_method2, na.rm=TRUE))^2)
r_squared <- 1 - (rss/tss)

# Number of observations and predictors
n <- nrow(analysis_data)
p <- length(model$coefficients) - 1

# Calculate adjusted R-squared
adj_r_squared <- 1 - ((1 - r_squared) * (n - 1) / (n - p - 1))

# Print adjusted R-squared
print(adj_r_squared)

# get slope values
slopes_scaled_human_utility <- ggpredict(model, terms = "scaled_human_utility")
df_slopes_scaled_human_utility2 <- as.data.frame(slopes_scaled_human_utility)

# plot the data
ggplot() +
  geom_point(data = method2_data, aes(y = predicted_richness_method2, x = scaled_human_utility, color = log10(greenspace_area_m2)), size = 3) +
  geom_line(data = df_slopes_scaled_human_utility2, aes(x = x, y = predicted), color = "blue") +
  geom_ribbon(data = df_slopes_scaled_human_utility2, aes(x = x, ymin = conf.low, ymax = conf.high), fill = "blue", alpha = 0.2) +
  xlab("Human Utility") +
  ylab("Biodiversity") +
  ggtitle("Method 2") +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  scale_colour_gradientn(colours = c("#D55E00", "#E69F00", "#F0E442", "#009E74", "#56B4E9", "#0072B2")) +
  labs(colour = expression(Log(Greenspace ~ Area ~ (m^2)))) +
  theme(legend.position = "bottom", legend.key.width = unit(0.6, "inch"))

## Method 3  -----------------------------------

# filter the data 
method3_data <- analysis_data %>%
  filter(complete.cases(predicted_richness_method3))

# calculate predicted biodiversity and human utility trend
model <- lm(predicted_richness_method3 ~ scaled_human_utility + log(greenspace_area_m2),
            data=method3_data, family="gaussian")
summary(model)

# Calculate R-squared
rss <- sum(model$residuals^2)
tss <- sum((method3_data$predicted_richness_method3 - mean(method3_data$predicted_richness_method3, na.rm=TRUE))^2)
r_squared <- 1 - (rss/tss)

# Number of observations and predictors
n <- nrow(analysis_data)
p <- length(model$coefficients) - 1

# Calculate adjusted R-squared
adj_r_squared <- 1 - ((1 - r_squared) * (n - 1) / (n - p - 1))

# Print adjusted R-squared
print(adj_r_squared)

# get slope values
slopes_scaled_human_utility <- ggpredict(model, terms = "scaled_human_utility")
df_slopes_scaled_human_utility3 <- as.data.frame(slopes_scaled_human_utility)

# plot the data
ggplot() +
  geom_point(data = method3_data, aes(y = predicted_richness_method3, x = scaled_human_utility, color = log10(greenspace_area_m2)), size = 3) +
  geom_line(data = df_slopes_scaled_human_utility3, aes(x = x, y = predicted), color = "blue") +
  geom_ribbon(data = df_slopes_scaled_human_utility3, aes(x = x, ymin = conf.low, ymax = conf.high), fill = "blue", alpha = 0.2) +
  xlab("Human Utility") +
  ylab("Biodiversity") +
  ggtitle("Method 3") +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  scale_colour_gradientn(colours = c("#D55E00", "#E69F00", "#F0E442", "#009E74", "#56B4E9", "#0072B2")) +
  labs(colour = expression(Log(Greenspace ~ Area ~ (m^2)))) +
  theme(legend.position = "bottom", legend.key.width = unit(0.6, "inch"))

## Method 4  -----------------------------------

# filter the data 
method4_data <- analysis_data %>%
  filter(complete.cases(predicted_richness_method4))

# calculate predicted biodiversity and human utility trend
model <- lm(predicted_richness_method4 ~ scaled_human_utility + log(greenspace_area_m2),
            data=method4_data, family="gaussian")
summary(model)

# Calculate R-squared
rss <- sum(model$residuals^2)
tss <- sum((method4_data$predicted_richness_method4 - mean(method4_data$predicted_richness_method4, na.rm=TRUE))^2)
r_squared <- 1 - (rss/tss)

# Number of observations and predictors
n <- nrow(analysis_data)
p <- length(model$coefficients) - 1

# Calculate adjusted R-squared
adj_r_squared <- 1 - ((1 - r_squared) * (n - 1) / (n - p - 1))

# Print adjusted R-squared
print(adj_r_squared)

# get slope values
slopes_scaled_human_utility <- ggpredict(model, terms = "scaled_human_utility")
df_slopes_scaled_human_utility4 <- as.data.frame(slopes_scaled_human_utility)

# plot the data
ggplot() +
  geom_point(data = method4_data, aes(y = predicted_richness_method4, x = scaled_human_utility, color = log10(greenspace_area_m2)), size = 3) +
  geom_line(data = df_slopes_scaled_human_utility4, aes(x = x, y = predicted), color = "blue") +
  geom_ribbon(data = df_slopes_scaled_human_utility4, aes(x = x, ymin = conf.low, ymax = conf.high), fill = "blue", alpha = 0.2) +
  xlab("Human Utility") +
  ylab("Biodiversity") +
  ggtitle("Method 4") +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  scale_colour_gradientn(colours = c("#D55E00", "#E69F00", "#F0E442", "#009E74", "#56B4E9", "#0072B2")) +
  labs(colour = expression(Log(Greenspace ~ Area ~ (m^2)))) +
  theme(legend.position = "bottom", legend.key.width = unit(0.6, "inch"))

## Method 5  -----------------------------------

# filter the data 
method5_data <- analysis_data %>%
  filter(complete.cases(predicted_richness_method5))

# calculate predicted biodiversity and human utility trend
model <- lm(predicted_richness_method5 ~ scaled_human_utility + log(greenspace_area_m2),
            data=method5_data, family="gaussian")
summary(model)

# Calculate R-squared
rss <- sum(model$residuals^2)
tss <- sum((method5_data$predicted_richness_method5 - mean(method5_data$predicted_richness_method5, na.rm=TRUE))^2)
r_squared <- 1 - (rss/tss)

# Number of observations and predictors
n <- nrow(analysis_data)
p <- length(model$coefficients) - 1

# Calculate adjusted R-squared
adj_r_squared <- 1 - ((1 - r_squared) * (n - 1) / (n - p - 1))

# Print adjusted R-squared
print(adj_r_squared)

# get slope values
slopes_scaled_human_utility <- ggpredict(model, terms = "scaled_human_utility")
df_slopes_scaled_human_utility5 <- as.data.frame(slopes_scaled_human_utility)

# plot the data
ggplot() +
  geom_point(data = method5_data, aes(y = predicted_richness_method5, x = scaled_human_utility, color = log10(greenspace_area_m2)), size = 3) +
  geom_line(data = df_slopes_scaled_human_utility5, aes(x = x, y = predicted), color = "blue") +
  geom_ribbon(data = df_slopes_scaled_human_utility5, aes(x = x, ymin = conf.low, ymax = conf.high), fill = "blue", alpha = 0.2) +
  xlab("Human Utility") +
  ylab("Biodiversity") +
  ggtitle("Method 5") +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  scale_colour_gradientn(colours = c("#D55E00", "#E69F00", "#F0E442", "#009E74", "#56B4E9", "#0072B2")) +
  labs(colour = expression(Log(Greenspace ~ Area ~ (m^2)))) +
  theme(legend.position = "bottom", legend.key.width = unit(0.6, "inch"))

## Combine methods into one plot -----------------------------------------------------------------------------------------------

# combine all the analysis data and make it into long format
all_analysis <- analysis_data %>%
  dplyr::select(scaled_human_utility, 
         greenspace_area_m2,
         predicted_richness_method1,
         predicted_richness_method2,
         predicted_richness_method3,
         predicted_richness_method4,
         predicted_richness_method5) %>%
  rename(`Method 1`=predicted_richness_method1,
         `Method 2`=predicted_richness_method2,
         `Method 3`=predicted_richness_method3,
         `Method 4`=predicted_richness_method4,
         `Method 5`=predicted_richness_method5) %>%
  pivot_longer(cols = 3:7, 
               names_to = "analysis_type", 
               values_to = "biodiversity_utility") %>%
  as.data.frame()

# make the analysis type a factor
all_analysis$analysis_type <- factor(all_analysis$analysis_type, levels=unique(all_analysis$analysis_type))

# combine slope data
df_slopes_scaled_human_utility1 <- df_slopes_scaled_human_utility1 %>%
  mutate(analysis_type="Method 1")
df_slopes_scaled_human_utility2 <- df_slopes_scaled_human_utility2 %>%
  mutate(analysis_type="Method 2")
df_slopes_scaled_human_utility3 <- df_slopes_scaled_human_utility3 %>%
  mutate(analysis_type="Method 3")
df_slopes_scaled_human_utility4 <- df_slopes_scaled_human_utility4 %>%
  mutate(analysis_type="Method 4")
df_slopes_scaled_human_utility5 <- df_slopes_scaled_human_utility5 %>%
  mutate(analysis_type="Method 5")

df_slopes_scaled_human_utility_all <- rbind(df_slopes_scaled_human_utility1,
                                            df_slopes_scaled_human_utility2,
                                            df_slopes_scaled_human_utility3,
                                            df_slopes_scaled_human_utility4,
                                            df_slopes_scaled_human_utility5)

# plot the data
ggplot() +
  # Points colored by log10 of greenspace_area_m2 (continuous)
  geom_point(data = all_analysis, 
             aes(y = biodiversity_utility, 
                 x = scaled_human_utility, 
                 color = log10(greenspace_area_m2)),
             size = 3) +
  
  # Add the prediction line for each analysis type (use linetype or group)
  geom_line(data = df_slopes_scaled_human_utility_all, 
            aes(x = x, y = predicted, linetype = analysis_type, group = analysis_type), 
            color = "blue", size=1.3) +
  
  # Add the confidence interval ribbon for each analysis type (fill by analysis_type)
  geom_ribbon(data = df_slopes_scaled_human_utility_all, 
              aes(x = x, ymin = conf.low, ymax = conf.high, fill = analysis_type), 
              alpha = 0.2) +
  
  # Labels and title
  xlab("Human Utility") +
  ylab("Biodiversity") +
  theme_classic() +
  
  # Customize text and legend
  theme(text = element_text(size = 20)) +
  labs(colour = expression(Log(Greenspace ~ Area ~ (m^2)))) +
  
  # Customize legend position and size
  theme(legend.position = "bottom", 
        legend.key.width = unit(0.6, "inch")) +
  
  scale_fill_manual(values = c("Method 1" = "grey50", 
                               "Method 2" = "grey50", 
                               "Method 3" = "grey50", 
                               "Method 4" = "grey50",
                               "Method 5" = "grey50"), 
                    name = "Method") + 
  
  # Use a continuous color scale for greenspace area
  scale_colour_gradientn(colours = c("#D55E00", "#E69F00", "#F0E442", "#009E74", "#56B4E9", "#0072B2"))

ggsave("Figures/Supplemental_Figures/FigureB1.jpeg", height=8, width=8, dpi=300)

# Model the Data -----------------------------------------------------------------------------------------

## Method 1 -----------------------------------------------------------------------------------------------

# Let's see how each human utility features influences bio-use
binary_data <- analysis_data %>% 
  dplyr::select(`Pavilion/Picnic Area`, `Kids Playground`,
                `Body of Water`, `Jog/Walk Path`,
                `Athletic Facility`, `Nature Preserve`,
                `Dog Park`, `Indoor/Outdoor Fitness Center`)
binary_data <- binary_data %>% mutate(across(everything(), factor))
binary_data <- binary_data %>% mutate(across(everything(), ~ ifelse(. == 0, "No", "Yes")))
binary_data <- binary_data %>%
  mutate(park_size=log(analysis_data$greenspace_area_m2),
         predicted_richness_method1=analysis_data$predicted_richness_method1)

# model the binary attributes
model_binary1 <- lm(predicted_richness_method1 ~ `Pavilion/Picnic Area` + `Kids Playground` + `Body of Water` +
                     `Jog/Walk Path` + `Athletic Facility` + `Nature Preserve` + `Dog Park` + `Indoor/Outdoor Fitness Center` +
                     park_size, 
                   data=binary_data)
summary(model_binary1)
anova(model_binary1)
AIC(model_binary1)
plot(model_binary1)

a <- effect_plot(model_binary1, pred=`Pavilion/Picnic Area`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Pavilion/Picnic Area")
b <- effect_plot(model_binary1, pred=`Kids Playground`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Kids Playground*")
c <- effect_plot(model_binary1, pred=`Body of Water`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Body of Water*")
d <- effect_plot(model_binary1, pred=`Jog/Walk Path`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Jog/Walk Path")
e <- effect_plot(model_binary1, pred=`Athletic Facility`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Athletic Facility**")
f <- effect_plot(model_binary1, pred=`Nature Preserve`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Nature Preserve**")
g <- effect_plot(model_binary1, pred=`Dog Park`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Dog Park")
h <- effect_plot(model_binary1, pred=`Indoor/Outdoor Fitness Center`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Indoor/Outdoor Fitness Center")

figure <- ggarrange(a + rremove("ylab"), b + rremove("ylab"), c + rremove("ylab"), d + rremove("ylab"), 
                    e + rremove("ylab"), f + rremove("ylab"), g + rremove("ylab"), h + rremove("ylab"), ncol=2, nrow=4)
annotate_figure(figure, 
                top = textGrob("Method 1", gp = gpar(cex = 1.5, fontface = "bold")),
                left=textGrob("Biodiversity", rot=90, vjust=1, gp = gpar(cex = 1.2))) 

## Method 2 -----------------------------------------------------------------------------------------------

# Let's see how each human utility features influences bio-use
binary_data <- analysis_data %>% 
  dplyr::select(`Pavilion/Picnic Area`, `Kids Playground`,
                `Body of Water`, `Jog/Walk Path`,
                `Athletic Facility`, `Nature Preserve`,
                `Dog Park`, `Indoor/Outdoor Fitness Center`)
binary_data <- binary_data %>% mutate(across(everything(), factor))
binary_data <- binary_data %>% mutate(across(everything(), ~ ifelse(. == 0, "No", "Yes")))
binary_data <- binary_data %>%
  mutate(park_size=log(analysis_data$greenspace_area_m2),
         predicted_richness_method2=analysis_data$predicted_richness_method2)

# model the binary attributes
model_binary2 <- lm(predicted_richness_method2 ~ `Pavilion/Picnic Area` + `Kids Playground` + `Body of Water` +
                      `Jog/Walk Path` + `Athletic Facility` + `Nature Preserve` + `Dog Park` + `Indoor/Outdoor Fitness Center` +
                      park_size, 
                    data=binary_data)
summary(model_binary2)
anova(model_binary2)
AIC(model_binary2)
plot(model_binary2)

a <- effect_plot(model_binary2, pred=`Pavilion/Picnic Area`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Pavilion/Picnic Area")
b <- effect_plot(model_binary2, pred=`Kids Playground`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Kids Playground*")
c <- effect_plot(model_binary2, pred=`Body of Water`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Body of Water*")
d <- effect_plot(model_binary2, pred=`Jog/Walk Path`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Jog/Walk Path")
e <- effect_plot(model_binary2, pred=`Athletic Facility`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Athletic Facility**")
f <- effect_plot(model_binary2, pred=`Nature Preserve`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Nature Preserve**")
g <- effect_plot(model_binary2, pred=`Dog Park`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Dog Park")
h <- effect_plot(model_binary2, pred=`Indoor/Outdoor Fitness Center`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Indoor/Outdoor Fitness Center")

figure <- ggarrange(a + rremove("ylab"), b + rremove("ylab"), c + rremove("ylab"), d + rremove("ylab"), 
                    e + rremove("ylab"), f + rremove("ylab"), g + rremove("ylab"), h + rremove("ylab"), ncol=2, nrow=4)
annotate_figure(figure, 
                top = textGrob("Method 2", gp = gpar(cex = 1.5, fontface = "bold")),
                left=textGrob("Biodiversity", rot=90, vjust=1, gp = gpar(cex = 1.2))) 

## Method 3 -----------------------------------------------------------------------------------------------

# Let's see how each human utility features influences bio-use
binary_data <- analysis_data %>% 
  dplyr::select(`Pavilion/Picnic Area`, `Kids Playground`,
                `Body of Water`, `Jog/Walk Path`,
                `Athletic Facility`, `Nature Preserve`,
                `Dog Park`, `Indoor/Outdoor Fitness Center`)
binary_data <- binary_data %>% mutate(across(everything(), factor))
binary_data <- binary_data %>% mutate(across(everything(), ~ ifelse(. == 0, "No", "Yes")))
binary_data <- binary_data %>%
  mutate(park_size=log(analysis_data$greenspace_area_m2),
         predicted_richness_method3=analysis_data$predicted_richness_method3)

# model the binary attributes
model_binary3 <- lm(predicted_richness_method3 ~ `Pavilion/Picnic Area` + `Kids Playground` + `Body of Water` +
                      `Jog/Walk Path` + `Athletic Facility` + `Nature Preserve` + `Dog Park` + `Indoor/Outdoor Fitness Center` +
                      park_size, 
                    data=binary_data)
summary(model_binary3)
anova(model_binary3)
AIC(model_binary3)
plot(model_binary3)

a <- effect_plot(model_binary3, pred=`Pavilion/Picnic Area`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Pavilion/Picnic Area")
b <- effect_plot(model_binary3, pred=`Kids Playground`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Kids Playground*")
c <- effect_plot(model_binary3, pred=`Body of Water`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Body of Water*")
d <- effect_plot(model_binary3, pred=`Jog/Walk Path`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Jog/Walk Path")
e <- effect_plot(model_binary3, pred=`Athletic Facility`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Athletic Facility**")
f <- effect_plot(model_binary3, pred=`Nature Preserve`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Nature Preserve**")
g <- effect_plot(model_binary3, pred=`Dog Park`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Dog Park")
h <- effect_plot(model_binary3, pred=`Indoor/Outdoor Fitness Center`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Indoor/Outdoor Fitness Center")

figure <- ggarrange(a + rremove("ylab"), b + rremove("ylab"), c + rremove("ylab"), d + rremove("ylab"), 
                    e + rremove("ylab"), f + rremove("ylab"), g + rremove("ylab"), h + rremove("ylab"), ncol=2, nrow=4)
annotate_figure(figure, 
                top = textGrob("Method 3", gp = gpar(cex = 1.5, fontface = "bold")),
                left=textGrob("Biodiversity", rot=90, vjust=1, gp = gpar(cex = 1.2))) 

## Method 4 -----------------------------------------------------------------------------------------------

# Let's see how each human utility features influences bio-use
binary_data <- analysis_data %>% 
  dplyr::select(`Pavilion/Picnic Area`, `Kids Playground`,
                `Body of Water`, `Jog/Walk Path`,
                `Athletic Facility`, `Nature Preserve`,
                `Dog Park`, `Indoor/Outdoor Fitness Center`)
binary_data <- binary_data %>% mutate(across(everything(), factor))
binary_data <- binary_data %>% mutate(across(everything(), ~ ifelse(. == 0, "No", "Yes")))
binary_data <- binary_data %>%
  mutate(park_size=log(analysis_data$greenspace_area_m2),
         predicted_richness_method4=analysis_data$predicted_richness_method4)

# model the binary attributes
model_binary4 <- lm(predicted_richness_method4 ~ `Pavilion/Picnic Area` + `Kids Playground` + `Body of Water` +
                      `Jog/Walk Path` + `Athletic Facility` + `Nature Preserve` + `Dog Park` + `Indoor/Outdoor Fitness Center` +
                      park_size, 
                    data=binary_data)
summary(model_binary4)
anova(model_binary4)
AIC(model_binary4)
plot(model_binary4)

a <- effect_plot(model_binary4, pred=`Pavilion/Picnic Area`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Pavilion/Picnic Area")
b <- effect_plot(model_binary4, pred=`Kids Playground`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Kids Playground*")
c <- effect_plot(model_binary4, pred=`Body of Water`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Body of Water*")
d <- effect_plot(model_binary4, pred=`Jog/Walk Path`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Jog/Walk Path")
e <- effect_plot(model_binary4, pred=`Athletic Facility`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Athletic Facility**")
f <- effect_plot(model_binary4, pred=`Nature Preserve`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Nature Preserve**")
g <- effect_plot(model_binary4, pred=`Dog Park`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Dog Park")
h <- effect_plot(model_binary4, pred=`Indoor/Outdoor Fitness Center`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Indoor/Outdoor Fitness Center")

figure <- ggarrange(a + rremove("ylab"), b + rremove("ylab"), c + rremove("ylab"), d + rremove("ylab"), 
                    e + rremove("ylab"), f + rremove("ylab"), g + rremove("ylab"), h + rremove("ylab"), ncol=2, nrow=4)
annotate_figure(figure, 
                top = textGrob("Method 4", gp = gpar(cex = 1.5, fontface = "bold")),
                left=textGrob("Biodiversity", rot=90, vjust=1, gp = gpar(cex = 1.2))) 

## Method 5 -----------------------------------------------------------------------------------------------

# Let's see how each human utility features influences bio-use
binary_data <- analysis_data %>% 
  dplyr::select(`Pavilion/Picnic Area`, `Kids Playground`,
                `Body of Water`, `Jog/Walk Path`,
                `Athletic Facility`, `Nature Preserve`,
                `Dog Park`, `Indoor/Outdoor Fitness Center`)
binary_data <- binary_data %>% mutate(across(everything(), factor))
binary_data <- binary_data %>% mutate(across(everything(), ~ ifelse(. == 0, "No", "Yes")))
binary_data <- binary_data %>%
  mutate(park_size=log(analysis_data$greenspace_area_m2),
         predicted_richness_method5=analysis_data$predicted_richness_method5)

# model the binary attributes
model_binary5 <- lm(predicted_richness_method5 ~ `Pavilion/Picnic Area` + `Kids Playground` + `Body of Water` +
                      `Jog/Walk Path` + `Athletic Facility` + `Nature Preserve` + `Dog Park` + `Indoor/Outdoor Fitness Center` +
                      park_size, 
                    data=binary_data)
summary(model_binary5)
anova(model_binary5)
AIC(model_binary5)
plot(model_binary5)

a <- effect_plot(model_binary5, pred=`Pavilion/Picnic Area`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Pavilion/Picnic Area")
b <- effect_plot(model_binary5, pred=`Kids Playground`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Kids Playground*")
c <- effect_plot(model_binary5, pred=`Body of Water`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Body of Water*")
d <- effect_plot(model_binary5, pred=`Jog/Walk Path`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Jog/Walk Path")
e <- effect_plot(model_binary5, pred=`Athletic Facility`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Athletic Facility**")
f <- effect_plot(model_binary5, pred=`Nature Preserve`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Nature Preserve**")
g <- effect_plot(model_binary5, pred=`Dog Park`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Dog Park")
h <- effect_plot(model_binary5, pred=`Indoor/Outdoor Fitness Center`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Indoor/Outdoor Fitness Center")

figure <- ggarrange(a + rremove("ylab"), b + rremove("ylab"), c + rremove("ylab"), d + rremove("ylab"), 
                    e + rremove("ylab"), f + rremove("ylab"), g + rremove("ylab"), h + rremove("ylab"), ncol=2, nrow=4)
annotate_figure(figure, 
                top = textGrob("Method 5", gp = gpar(cex = 1.5, fontface = "bold")),
                left=textGrob("Biodiversity", rot=90, vjust=1, gp = gpar(cex = 1.2))) 

## Combine plots -----------------------------------------------

# Initialize a list to store models
models <- list()

# Define the response variables
response_vars <- c("predicted_richness_method1", 
                   "predicted_richness_method2", 
                   "predicted_richness_method3", 
                   "predicted_richness_method4", 
                   "predicted_richness_method5")

# Prepare the data
binary_data <- analysis_data %>% 
  mutate(across(c(`Pavilion/Picnic Area`, `Kids Playground`,
                  `Body of Water`, `Jog/Walk Path`,
                  `Athletic Facility`, `Nature Preserve`,
                  `Dog Park`, `Indoor/Outdoor Fitness Center`), 
                ~ factor(ifelse(. == 0, "No", "Yes")))) %>%
  mutate(park_size = log(greenspace_area_m2))

# Loop through each dataset
for (i in 1:5) {
  # Ensure the response variable exists in the dataset
  if (response_vars[i] %in% names(binary_data)) {
    # Fit the model
    model_formula <- as.formula(paste(response_vars[i], "~ `Pavilion/Picnic Area` + `Kids Playground` + `Body of Water` +
                                       `Jog/Walk Path` + `Athletic Facility` + `Nature Preserve` + `Dog Park` + 
                                       `Indoor/Outdoor Fitness Center` + park_size"))
    models[[i]] <- lm(model_formula, data = binary_data)
  } else {
    print(paste("Response variable", response_vars[i], "does not exist in dataset", i))
  }
}

# Get predictions
plot_data <- data.frame()

for (i in 1:5) {
  if (!is.null(models[[i]])) {
    # Ensure factors maintain the same levels for predictions
    new_data <- expand.grid(
      `Pavilion/Picnic Area` = factor(c("No", "Yes"), levels = levels(binary_data$`Pavilion/Picnic Area`)),
      `Kids Playground` = factor("Yes", levels = levels(binary_data$`Kids Playground`)),
      `Body of Water` = factor("Yes", levels = levels(binary_data$`Body of Water`)),
      `Jog/Walk Path` = factor("Yes", levels = levels(binary_data$`Jog/Walk Path`)),
      `Athletic Facility` = factor("Yes", levels = levels(binary_data$`Athletic Facility`)),
      `Nature Preserve` = factor("Yes", levels = levels(binary_data$`Nature Preserve`)),
      `Dog Park` = factor("Yes", levels = levels(binary_data$`Dog Park`)),
      `Indoor/Outdoor Fitness Center` = factor("Yes", levels = levels(binary_data$`Indoor/Outdoor Fitness Center`)),
      park_size = mean(binary_data$park_size, na.rm = TRUE)  # Use mean park size
    )
    
    # Get predictions with confidence intervals
    pred <- predict(models[[i]], newdata = new_data, interval = "confidence")
    
    # Add predictions to data frame
    new_data <- cbind(new_data, as.data.frame(pred))
    new_data$dataset <- paste("Method", i)
    
    # Combine with previous data
    plot_data <- rbind(plot_data, new_data)
  }
}

# Plot the predictions using ggplot2
(pp <- ggplot(plot_data, aes(x = `Pavilion/Picnic Area`, y = fit, color = dataset)) +
  geom_line(aes(group = dataset), size = 1) +  # Plot the predicted values, grouping by dataset
  # Add vertical lines for the lower and upper error bounds
  geom_segment(aes(x = `Pavilion/Picnic Area`, xend = `Pavilion/Picnic Area`, 
                   y = lwr, yend = upr), 
               alpha = 0.5, size = 1, 
               arrow = arrow(ends = "both", type = "closed", length = unit(0.05, "inches"))) +  # Adjust color and transparency
  labs(title = "Pavilion/Picnic Area",
       x = "",
       y = "Predicted Response",
       color = NULL) +  # Set color to NULL to remove legend title
  theme_minimal() +
    scale_color_manual(values = c("red", "gold2", "cyan3", "blue3", "purple")) +  # Use a color palette for better distinction
  theme(legend.position = "right"))

# Kids playground

# Create an empty data frame for combined predictions
plot_data <- data.frame()

# Assuming your original data has the appropriate levels for factors
original_levels <- levels(binary_data$`Kids Playground`)

for (i in 1:5) {
  if (!is.null(models[[i]])) {
    # Create a new dataframe for predictions with the same factor levels
    new_data <- expand.grid(
      `Pavilion/Picnic Area` = "Yes",  # Ensure this matches original factor levels
      `Kids Playground` = factor(original_levels),
      `Body of Water` = "Yes",
      `Jog/Walk Path` = "Yes",
      `Athletic Facility` = "Yes",
      `Nature Preserve` = "Yes",
      `Dog Park` = "Yes",
      `Indoor/Outdoor Fitness Center` = "Yes",
      park_size = mean(binary_data$park_size, na.rm = TRUE)  # Use the mean or a relevant value
    )
    
    # Get predictions with confidence intervals
    pred <- predict(models[[i]], newdata = new_data, interval = "confidence")
    
    # Create a data frame with predictions and intervals
    new_data <- cbind(new_data, as.data.frame(pred))
    new_data$dataset <- paste("Method", i)  # Add dataset identifier
    
    # Combine with previous data
    plot_data <- rbind(plot_data, new_data)
  }
}


# Plot the predictions using ggplot2
(kp <- ggplot(plot_data, aes(x = `Kids Playground`, y = fit, color = dataset)) +
    geom_line(aes(group = dataset), size = 1) +  # Plot the predicted values, grouping by dataset
    # Add vertical lines for the lower and upper error bounds
    geom_segment(aes(x = `Kids Playground`, xend = `Kids Playground`, 
                     y = lwr, yend = upr), 
                 alpha = 0.5, size = 1, 
                 arrow = arrow(ends = "both", type = "closed", length = unit(0.05, "inches"))) +  # Adjust color and transparency
    labs(title = "Kids Playground",
         x = "",
         y = "Predicted Response",
         color = NULL) +  # Set color to NULL to remove legend title
    theme_minimal() +
    scale_color_manual(values = c("red", "gold2", "cyan3", "blue3", "purple")) +  # Use a color palette for better distinction
    theme(legend.position = "right"))

# Body of Water

# Create an empty data frame for combined predictions
plot_data <- data.frame()

# Assuming your original data has the appropriate levels for factors
original_levels <- levels(binary_data$`Body of Water`)

for (i in 1:5) {
  if (!is.null(models[[i]])) {
    # Create a new dataframe for predictions with the same factor levels
    new_data <- expand.grid(
      `Pavilion/Picnic Area` = "Yes",  # Ensure this matches original factor levels
      `Kids Playground` = "Yes",
      `Body of Water` = factor(original_levels),
      `Jog/Walk Path` = "Yes",
      `Athletic Facility` = "Yes",
      `Nature Preserve` = "Yes",
      `Dog Park` = "Yes",
      `Indoor/Outdoor Fitness Center` = "Yes",
      park_size = mean(binary_data$park_size, na.rm = TRUE)  # Use the mean or a relevant value
    )
    
    # Get predictions with confidence intervals
    pred <- predict(models[[i]], newdata = new_data, interval = "confidence")
    
    # Create a data frame with predictions and intervals
    new_data <- cbind(new_data, as.data.frame(pred))
    new_data$dataset <- paste("Method", i)  # Add dataset identifier
    
    # Combine with previous data
    plot_data <- rbind(plot_data, new_data)
  }
}


# Plot the predictions using ggplot2
(bw <- ggplot(plot_data, aes(x = `Body of Water`, y = fit, color = dataset)) +
    geom_line(aes(group = dataset), size = 1) +  # Plot the predicted values, grouping by dataset
    # Add vertical lines for the lower and upper error bounds
    geom_segment(aes(x = `Body of Water`, xend = `Body of Water`, 
                     y = lwr, yend = upr), 
                 alpha = 0.5, size = 1, 
                 arrow = arrow(ends = "both", type = "closed", length = unit(0.05, "inches"))) +  # Adjust color and transparency
    labs(title = "Body of Water",
         x = "",
         y = "Predicted Response",
         color = NULL) +  # Set color to NULL to remove legend title
    theme_minimal() +
    scale_color_manual(values = c("red", "gold2", "cyan3", "blue3", "purple")) +  # Use a color palette for better distinction
    theme(legend.position = "right"))

# Job/Walk Path

# Create an empty data frame for combined predictions
plot_data <- data.frame()

# Assuming your original data has the appropriate levels for factors
original_levels <- levels(binary_data$`Jog/Walk Path`)

for (i in 1:5) {
  if (!is.null(models[[i]])) {
    # Create a new dataframe for predictions with the same factor levels
    new_data <- expand.grid(
      `Pavilion/Picnic Area` = "Yes",  # Ensure this matches original factor levels
      `Kids Playground` = "Yes",
      `Body of Water` = "Yes",
      `Jog/Walk Path` = factor(original_levels),
      `Athletic Facility` = "Yes",
      `Nature Preserve` = "Yes",
      `Dog Park` = "Yes",
      `Indoor/Outdoor Fitness Center` = "Yes",
      park_size = mean(binary_data$park_size, na.rm = TRUE)  # Use the mean or a relevant value
    )
    
    # Get predictions with confidence intervals
    pred <- predict(models[[i]], newdata = new_data, interval = "confidence")
    
    # Create a data frame with predictions and intervals
    new_data <- cbind(new_data, as.data.frame(pred))
    new_data$dataset <- paste("Method", i)  # Add dataset identifier
    
    # Combine with previous data
    plot_data <- rbind(plot_data, new_data)
  }
}


# Plot the predictions using ggplot2
(jw <- ggplot(plot_data, aes(x = `Jog/Walk Path`, y = fit, color = dataset)) +
    geom_line(aes(group = dataset), size = 1) +  # Plot the predicted values, grouping by dataset
    # Add vertical lines for the lower and upper error bounds
    geom_segment(aes(x = `Jog/Walk Path`, xend = `Jog/Walk Path`, 
                     y = lwr, yend = upr), 
                 alpha = 0.5, size = 1, 
                 arrow = arrow(ends = "both", type = "closed", length = unit(0.05, "inches"))) +  # Adjust color and transparency
    labs(title = "Jog/Walk Path",
         x = "",
         y = "Predicted Response",
         color = NULL) +  # Set color to NULL to remove legend title
    theme_minimal() +
    scale_color_manual(values = c("red", "gold2", "cyan3", "blue3", "purple")) + # Use a color palette for better distinction
    theme(legend.position = "right"))

# Athletic Facility

# Create an empty data frame for combined predictions
plot_data <- data.frame()

# Assuming your original data has the appropriate levels for factors
original_levels <- levels(binary_data$`Athletic Facility`)

for (i in 1:5) {
  if (!is.null(models[[i]])) {
    # Create a new dataframe for predictions with the same factor levels
    new_data <- expand.grid(
      `Pavilion/Picnic Area` = "Yes",  # Ensure this matches original factor levels
      `Kids Playground` = "Yes",
      `Body of Water` = "Yes",
      `Jog/Walk Path` = "Yes",
      `Athletic Facility` = factor(original_levels),
      `Nature Preserve` = "Yes",
      `Dog Park` = "Yes",
      `Indoor/Outdoor Fitness Center` = "Yes",
      park_size = mean(binary_data$park_size, na.rm = TRUE)  # Use the mean or a relevant value
    )
    
    # Get predictions with confidence intervals
    pred <- predict(models[[i]], newdata = new_data, interval = "confidence")
    
    # Create a data frame with predictions and intervals
    new_data <- cbind(new_data, as.data.frame(pred))
    new_data$dataset <- paste("Method", i)  # Add dataset identifier
    
    # Combine with previous data
    plot_data <- rbind(plot_data, new_data)
  }
}


# Plot the predictions using ggplot2
(af <- ggplot(plot_data, aes(x = `Athletic Facility`, y = fit, color = dataset)) +
    geom_line(aes(group = dataset), size = 1) +  # Plot the predicted values, grouping by dataset
    # Add vertical lines for the lower and upper error bounds
    geom_segment(aes(x = `Athletic Facility`, xend = `Athletic Facility`, 
                     y = lwr, yend = upr), 
                 alpha = 0.5, size = 1, 
                 arrow = arrow(ends = "both", type = "closed", length = unit(0.05, "inches"))) +  # Adjust color and transparency
    labs(title = "Athletic Facility",
         x = "",
         y = "Predicted Response",
         color = NULL) +  # Set color to NULL to remove legend title
    theme_minimal() +
    scale_color_manual(values = c("red", "gold2", "cyan3", "blue3", "purple")) + # Use a color palette for better distinction
    theme(legend.position = "right"))

# Nature Preserve

# Create an empty data frame for combined predictions
plot_data <- data.frame()

# Assuming your original data has the appropriate levels for factors
original_levels <- levels(binary_data$`Nature Preserve`)

for (i in 1:5) {
  if (!is.null(models[[i]])) {
    # Create a new dataframe for predictions with the same factor levels
    new_data <- expand.grid(
      `Pavilion/Picnic Area` = "Yes",  # Ensure this matches original factor levels
      `Kids Playground` = "Yes",
      `Body of Water` = "Yes",
      `Jog/Walk Path` = "Yes",
      `Athletic Facility` = "Yes",
      `Nature Preserve` = factor(original_levels),
      `Dog Park` = "Yes",
      `Indoor/Outdoor Fitness Center` = "Yes",
      park_size = mean(binary_data$park_size, na.rm = TRUE)  # Use the mean or a relevant value
    )
    
    # Get predictions with confidence intervals
    pred <- predict(models[[i]], newdata = new_data, interval = "confidence")
    
    # Create a data frame with predictions and intervals
    new_data <- cbind(new_data, as.data.frame(pred))
    new_data$dataset <- paste("Method", i)  # Add dataset identifier
    
    # Combine with previous data
    plot_data <- rbind(plot_data, new_data)
  }
}


# Plot the predictions using ggplot2
(np <- ggplot(plot_data, aes(x = `Nature Preserve`, y = fit, color = dataset)) +
    geom_line(aes(group = dataset), size = 1) +  # Plot the predicted values, grouping by dataset
    # Add vertical lines for the lower and upper error bounds
    geom_segment(aes(x = `Nature Preserve`, xend = `Nature Preserve`, 
                     y = lwr, yend = upr), 
                 alpha = 0.5, size = 1, 
                 arrow = arrow(ends = "both", type = "closed", length = unit(0.05, "inches"))) +  # Adjust color and transparency
    labs(title = "Nature Preserve",
         x = "",
         y = "Predicted Response",
         color = NULL) +  # Set color to NULL to remove legend title
    theme_minimal() +
    scale_color_manual(values = c("red", "gold2", "cyan3", "blue3", "purple")) +  # Use a color palette for better distinction
    theme(legend.position = "right"))

# Dog Park

# Create an empty data frame for combined predictions
plot_data <- data.frame()

# Assuming your original data has the appropriate levels for factors
original_levels <- levels(binary_data$`Dog Park`)

for (i in 1:5) {
  if (!is.null(models[[i]])) {
    # Create a new dataframe for predictions with the same factor levels
    new_data <- expand.grid(
      `Pavilion/Picnic Area` = "Yes",  # Ensure this matches original factor levels
      `Kids Playground` = "Yes",
      `Body of Water` = "Yes",
      `Jog/Walk Path` = "Yes",
      `Athletic Facility` = "Yes",
      `Nature Preserve` = "Yes",
      `Dog Park` = factor(original_levels),
      `Indoor/Outdoor Fitness Center` = "Yes",
      park_size = mean(binary_data$park_size, na.rm = TRUE)  # Use the mean or a relevant value
    )
    
    # Get predictions with confidence intervals
    pred <- predict(models[[i]], newdata = new_data, interval = "confidence")
    
    # Create a data frame with predictions and intervals
    new_data <- cbind(new_data, as.data.frame(pred))
    new_data$dataset <- paste("Method", i)  # Add dataset identifier
    
    # Combine with previous data
    plot_data <- rbind(plot_data, new_data)
  }
}


# Plot the predictions using ggplot2
(dp <- ggplot(plot_data, aes(x = `Dog Park`, y = fit, color = dataset)) +
    geom_line(aes(group = dataset), size = 1) +  # Plot the predicted values, grouping by dataset
    # Add vertical lines for the lower and upper error bounds
    geom_segment(aes(x = `Dog Park`, xend = `Dog Park`, 
                     y = lwr, yend = upr), 
                 alpha = 0.5, size = 1, 
                 arrow = arrow(ends = "both", type = "closed", length = unit(0.05, "inches"))) +  # Adjust color and transparency
    labs(title = "Dog Park",
         x = "",
         y = "Predicted Response",
         color = NULL) +  # Set color to NULL to remove legend title
    theme_minimal() +
    scale_color_manual(values = c("red", "gold2", "cyan3", "blue3", "purple")) +  # Use a color palette for better distinction
    theme(legend.position = "right"))

# Indoor/Outdoor Fitness Center

# Create an empty data frame for combined predictions
plot_data <- data.frame()

# Assuming your original data has the appropriate levels for factors
original_levels <- levels(binary_data$`Indoor/Outdoor Fitness Center`)

for (i in 1:5) {
  if (!is.null(models[[i]])) {
    # Create a new dataframe for predictions with the same factor levels
    new_data <- expand.grid(
      `Pavilion/Picnic Area` = "Yes",  # Ensure this matches original factor levels
      `Kids Playground` = "Yes",
      `Body of Water` = "Yes",
      `Jog/Walk Path` = "Yes",
      `Athletic Facility` = "Yes",
      `Nature Preserve` = "Yes",
      `Dog Park` = "Yes",
      `Indoor/Outdoor Fitness Center` = factor(original_levels),
      park_size = mean(binary_data$park_size, na.rm = TRUE)  # Use the mean or a relevant value
    )
    
    # Get predictions with confidence intervals
    pred <- predict(models[[i]], newdata = new_data, interval = "confidence")
    
    # Create a data frame with predictions and intervals
    new_data <- cbind(new_data, as.data.frame(pred))
    new_data$dataset <- paste("Method", i)  # Add dataset identifier
    
    # Combine with previous data
    plot_data <- rbind(plot_data, new_data)
  }
}


# Plot the predictions using ggplot2
(fc <- ggplot(plot_data, aes(x = `Indoor/Outdoor Fitness Center`, y = fit, color = dataset)) +
    geom_line(aes(group = dataset), size = 1) +  # Plot the predicted values, grouping by dataset
    # Add vertical lines for the lower and upper error bounds
    geom_segment(aes(x = `Indoor/Outdoor Fitness Center`, xend = `Indoor/Outdoor Fitness Center`, 
                     y = lwr, yend = upr), 
                 alpha = 0.5, size = 1, 
                 arrow = arrow(ends = "both", type = "closed", length = unit(0.05, "inches"))) +  # Adjust color and transparency
    labs(title = "Indoor/Outdoor Fitness Center",
         x = "",
         y = "Predicted Response",
         color = NULL) +  # Set color to NULL to remove legend title
    theme_minimal() +
    scale_color_manual(values = c("red", "gold2", "cyan3", "blue3", "purple")) +  # Use a color palette for better distinction
    theme(legend.position = "right"))

figure <- ggarrange(pp + rremove("ylab"), kp + rremove("ylab"), bw + rremove("ylab"), jw + rremove("ylab"), 
                    af + rremove("ylab"), np + rremove("ylab"), dp + rremove("ylab"), fc + rremove("ylab"), ncol=2, nrow=4)
annotate_figure(figure, 
                left=textGrob("Biodiversity", rot=90, vjust=1, gp = gpar(cex = 1.2))) 
ggsave("Figures/Supplemental_Figures/FigureB2.jpg", height=10, width=8, units="in", dpi=300)

