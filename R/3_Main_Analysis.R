# script is the main analysis script

# packages
library(ggplot2)
library(readr)
library(tidyverse)
library(scales)
library(sf)
library(ggpubr)
library(margins)
library(jtools)
library(grid)
library(corrgram)
library(lares)
library(corrplot)
library(cluster)
library(Rtsne)
library(ggeffects)
library(effects)
library(broom)
library(reshape2)
library(patchwork)
library(stringr)

# read in human utility file
human_utility <- read_csv("Data/human_utility.csv") %>%
  mutate(total_utility=rowSums(.[3:10])) %>%
  unite(poly_id, Municipality, `Park name`)

# read in file produced in 2_Calculate_Biodiversity_Index that contains the predicted richness value
biodiversity_utility <- read_csv("Data/bio_use_predictors/bio_use_predictions.csv")

# read in urban park data
greenspaces <- st_read("Data/shapefiles/urban_parks.shp") %>%
  mutate(greenspace_area_m2=as.numeric(st_area(.)))

# merge the three data sources into one data frame and rescale human utility and biodiversity
# measures to be between 0 and 1
analysis_data <- greenspaces %>%
  st_set_geometry(NULL) %>%
  left_join(., biodiversity_utility, by="poly_id") %>%
  left_join(., human_utility, by="poly_id") %>%
  mutate(scaled_human_utility=scales::rescale(total_utility)) %>%
  mutate(scaled_biodiversity_utility=scales::rescale(predicted_richness)) %>%
  mutate(scaled_iNat_users=scales::rescale(iNat_users)) 

# Comparison Plots of Human Utility, Biodiversity Utility, and Park Size --------

# compare park size and human utility
ggplot(analysis_data, aes(x=greenspace_area_m2, y=total_utility))+
  geom_point()+
  geom_smooth(method="lm")+
  scale_x_log10()+
  theme_minimal()

# histogram of greenspace size
ggplot(analysis_data, aes(x = (greenspace_area_m2/10000)))+ 
  geom_histogram(fill = "darkseagreen3", color = "black")+  
  xlab("Greenspace area (ha)") + 
  ylab("Number of greenspaces")+
  theme_classic() +
  scale_x_log10() + 
  theme(text=element_text(size=25)) 

## Visualize relationships in the data -----------------------------------

# Now visualize total utility and predicted richness vs park size
(ps_ns <- ggplot(analysis_data, aes(x=greenspace_area_m2, y=scaled_biodiversity_utility))+
  geom_point()+
  geom_smooth(method="lm")+
  xlab("Park Size (m2)")+
  ylab("Biodiversity")+
  ggtitle("Predicted Bio-use with Park Size") +
  scale_x_log10() + 
   theme_minimal())

(ps_ut <- ggplot(analysis_data, aes(x=greenspace_area_m2, y=total_utility))+
    geom_point()+
    geom_smooth(method="lm")+
    xlab("Park Size (m2)")+
    ylab("Human Utility")+
    ggtitle("Predicted Bio-use with Park Size") +
    scale_x_log10() + 
    theme_minimal())

ps_ut / ps_ns 

ggsave("Figures/Supplemental_Figures/Figure2A.jpeg", height=10, width=8, units="in")

# plot biodiversity utility and human utility
(ps_s <- ggplot(analysis_data, aes(y=scaled_biodiversity_utility, x=total_utility))+
  geom_point()+
  xlab("Human utility")+
  ylab("Biodiversity utility")+
  ggtitle("Predicted Bio-use with Park Size") +
  theme_minimal() +
  geom_smooth())

# plot human utility by biodiversity utility and color points by greenspace area
(wops_sps <- ggplot(analysis_data, aes(y=scaled_biodiversity_utility, x=scaled_human_utility, color=log10(greenspace_area_m2)))+
    geom_point(size=3)+
    xlab("Human Utility")+
    ylab("Biodiversity")+
    theme_classic()+
    theme(text=element_text(size=20)) +
    scale_colour_gradientn(colours = c("#D55E00", "#E69F00", "#F0E442", "#009E74", "#56B4E9", "#0072B2")) +
    labs(colour=expression(Log(Greenspace~Area~(m^2)))) + 
    theme(legend.position = "bottom",
          legend.key.width = unit(0.6, "inch")))

# add trend line from GLM
model <- glm(scaled_biodiversity_utility ~ scaled_human_utility + log(greenspace_area_m2),
             data=analysis_data, family="gaussian")
summary(model)

# get predictions
slopes_scaled_human_utility <- ggpredict(model, terms = "scaled_human_utility")
df_slopes_scaled_human_utility <- as.data.frame(slopes_scaled_human_utility)

# make the previous plot with the glm trend line
ggplot() +
  geom_point(data = analysis_data, aes(y = scaled_biodiversity_utility, x = scaled_human_utility, color = log10(greenspace_area_m2)), size = 3) +
  geom_line(data = df_slopes_scaled_human_utility, aes(x = x, y = predicted), color = "blue") +
  geom_ribbon(data = df_slopes_scaled_human_utility, aes(x = x, ymin = conf.low, ymax = conf.high), fill = "blue", alpha = 0.2) +
  xlab("Human Utility") +
  ylab("Biodiversity") +
  ggtitle("Method 1") +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  scale_colour_gradientn(colours = c("#D55E00", "#E69F00", "#F0E442", "#009E74", "#56B4E9", "#0072B2")) +
  labs(colour = expression(Log(Greenspace ~ Area ~ (m^2)))) +
  theme(legend.position = "bottom", legend.key.width = unit(0.6, "inch"))

ggsave("Figures/Figure3.jpeg", height=8, width=8, dpi=300)

## Biodiversity Predicted with Number of iNat Users -----------------------------------

# Now visualize total utility and predicted richness vs park size
(users_ns <- ggplot(analysis_data, aes(x=greenspace_area_m2, y=iNat_users))+
   geom_point()+
   geom_smooth(method="lm")+
   xlab("Park Size (m2)")+
   ylab("Number of iNaturalist Users")+
   ggtitle("iNat Users") +
   theme_minimal())

# plot biodiversity utility and human utility
(users_s <- ggplot(analysis_data, aes(y=scaled_iNat_users, x=total_utility))+
    geom_point()+
    xlab("Human utility")+
    ylab("Number of iNaturalist Users")+
    ggtitle("iNat Users") +
    theme_minimal() +
    geom_smooth())

# plot human utility by number of iNaturalist users and color code by greenspace area
(users_sps <- ggplot(analysis_data, aes(y=scaled_iNat_users, x=scaled_human_utility, color=log10(greenspace_area_m2)))+
    geom_point()+
    xlab("Human utility")+
    ylab("Number of iNaturalist Users")+
    ggtitle("iNat Users") +
    theme_minimal()+
    scale_colour_gradientn(colours = terrain.colors(10)))


# Physical Attribute Plots ------------------------------------------------

# read in the data
human_utility <- read_csv("Data/human_utility.csv") %>%
  mutate(total_utility=rowSums(.[3:10])) %>%
  unite(poly_id, Municipality, `Park name`)

human_utility <- human_utility %>%
  pivot_longer(cols = -poly_id, 
               names_to = "attribute", 
               values_to = "value") %>%
  mutate(value=ifelse(value==1, "Presence", "Absence"))

(appearence <- ggplot(human_utility, aes(x = str_wrap(attribute, width=15), fill = factor(value)))+
    geom_bar(position = "dodge", color = "black")+
    scale_fill_manual(values=c("salmon2", "darkseagreen3")) +
    labs(x = "Physical attribute", y = "Count", fill = "")+
    theme_classic()+
    ylab("Number of green spaces") + 
    theme(text=element_text(size=20),
          legend.position = "top",
          axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))) 

hist <- ggplot(analysis_data, aes(x = total_utility))+ 
  geom_histogram(binwidth = 1, fill = "darkseagreen3", color = "black")+  
  xlab("Number of physical attributes")+
  ylab("Number of green spaces")+
  theme_classic() +
  theme(text=element_text(size=20)) 

# plot together using patchwork
combined_plot <- hist + appearence 

# Customize layout (optional)
combined_plot <- combined_plot + plot_layout(ncol = 2) +
  plot_annotation(tag_levels = c("a", "b"))

# Display the combined plot
combined_plot

ggsave("Figures/Figure2.jpeg", height=8, width=14, dpi=300)

# Correlations between Human Utility Types --------

jpeg("Figures/Supplemental_Figures/FigureA3.jpeg", , width = 3000, height = 3000 * sqrt(2))

# create a correlogram to compare human utility attributes
analysis_data %>% 
  dplyr::select(`Pavilion/Picnic Area`, `Kids Playground`, `Body of Water`, `Jog/Walk Path`,
         `Athletic Facility`, `Nature Preserve`, `Dog Park`, `Indoor/Outdoor Fitness Center`) %>%
  cor() %>%
  corrplot(method="color", order='hclust', tl.col='black', tl.cex=5, cl.cex=5)

# let's check significance of the correlations
pval <- psych::corr.test(analysis_data %>% 
                           dplyr::select(`Pavilion/Picnic Area`, `Kids Playground`, `Body of Water`, `Jog/Walk Path`,
                                         `Athletic Facility`, `Nature Preserve`, `Dog Park`, `Indoor/Outdoor Fitness Center`), 
                         adjust="none")$p

# let's format the p-values 
pos <- expand.grid(1:ncol(pval), ncol(pval):1)
text(pos, format.pval(pval, digits=2, eps=0.001), cex=4)

dev.off()

# calculate the matrix cross product of the variables
mcp <- crossprod(as.matrix(analysis_data %>% 
                             dplyr::select(`Pavilion/Picnic Area`, `Kids Playground`, `Body of Water`, `Jog/Walk Path`,
                                    `Athletic Facility`, `Nature Preserve`, `Dog Park`, `Indoor/Outdoor Fitness Center`)))

mcp[upper.tri(mcp, diag=TRUE)] <- 0 
subset(as.data.frame.table(mcp), Freq > 0)

# correlation parameters 
corrs <- as.data.frame(psych::corr.test(analysis_data %>% 
                                          dplyr::select(`Pavilion/Picnic Area`, `Kids Playground`, `Body of Water`, `Jog/Walk Path`,
                                                        `Athletic Facility`, `Nature Preserve`, `Dog Park`, `Indoor/Outdoor Fitness Center`), 
                                        adjust="none")$r
)

corr.test_result <- psych::corr.test(analysis_data %>% 
                                 dplyr::select(`Pavilion/Picnic Area`, `Kids Playground`, `Body of Water`, `Jog/Walk Path`,
                                               `Athletic Facility`, `Nature Preserve`, `Dog Park`, `Indoor/Outdoor Fitness Center`), 
                               adjust="none")

corr.test_result$ci

# Model the Data --------

## Bio-use, Human Utility, and Park Size --------

# We want to model biouse using human utility and park size
# To start, let's look at the data trends
hist(analysis_data$scaled_biodiversity_utility)
hist(analysis_data$scaled_human_utility)
hist(analysis_data$greenspace_area_m2)

# use a glm to model biodiversity utility by human utility and greenspace area
model <- glm(scaled_biodiversity_utility ~ scaled_human_utility + log(greenspace_area_m2),
             data=analysis_data, family="gaussian")
summary(model)
anova(model)
AIC(model)
plot(model)

# Calculate R-squared
rss <- sum(model$residuals^2)
tss <- sum((analysis_data$scaled_biodiversity_utility - mean(analysis_data$scaled_biodiversity_utility))^2)
r_squared <- 1 - (rss/tss)

# Number of observations and predictors
n <- nrow(analysis_data)
p <- length(model$coefficients) - 1

# Calculate adjusted R-squared
adj_r_squared <- 1 - ((1 - r_squared) * (n - 1) / (n - p - 1))

# Print adjusted R-squared
print(adj_r_squared)

plot(effect("scaled_human_utility", model), main = "Slope Line for Scaled Human Utility")
plot(effect("log(greenspace_area_m2)", model), main = "Slope Line for Log Greenspace Area")

# let's try removing park size and running the model to see if that improves model fit
model2 <- glm(scaled_biodiversity_utility ~ scaled_human_utility,
             data=analysis_data, family="gaussian")
summary(model2)
anova(model2)
AIC(model2)
plot(model2)
# based on the AIC, this does not improve model fit

# Calculate R-squared
rss <- sum(model2$residuals^2)
tss <- sum((analysis_data$scaled_biodiversity_utility - mean(analysis_data$scaled_biodiversity_utility))^2)
r_squared <- 1 - (rss/tss)

# Number of observations and predictors
n <- nrow(analysis_data)
p <- length(model$coefficients) - 1

# Calculate adjusted R-squared
adj_r_squared <- 1 - ((1 - r_squared) * (n - 1) / (n - p - 1))

# Print adjusted R-squared
print(adj_r_squared)

# let's check the model without human utility
model3 <- glm(scaled_biodiversity_utility ~ log(greenspace_area_m2),
             data=analysis_data, family="gaussian")
summary(model3)
anova(model3)
AIC(model3)
plot(model3)

# Calculate R-squared
rss <- sum(model3$residuals^2)
tss <- sum((analysis_data$scaled_biodiversity_utility - mean(analysis_data$scaled_biodiversity_utility))^2)
r_squared <- 1 - (rss/tss)

# Number of observations and predictors
n <- nrow(analysis_data)
p <- length(model3$coefficients) - 1

# Calculate adjusted R-squared
adj_r_squared <- 1 - ((1 - r_squared) * (n - 1) / (n - p - 1))

# Print adjusted R-squared
print(adj_r_squared)

## Human Utility and Park Size --------

# We are also interested in the relationship between human utility and greenspaces, so let's create another
# model with human utility as the response
model_hu <- glm(scaled_human_utility ~ log(greenspace_area_m2), data=analysis_data, family="gaussian")
summary(model_hu)
anova(model_hu)
AIC(model_hu)
plot(model_hu)

# Calculate R-squared
rss <- sum(model_hu$residuals^2)
tss <- sum((analysis_data$scaled_human_utility - mean(analysis_data$scaled_human_utility))^2)
r_squared <- 1 - (rss/tss)

# Number of observations and predictors
n <- nrow(analysis_data)
p <- length(model_hu$coefficients) - 1

# Calculate adjusted R-squared
adj_r_squared <- 1 - ((1 - r_squared) * (n - 1) / (n - p - 1))

# Print adjusted R-squared
print(adj_r_squared)

## Bio-Use by Binary Human Utility Features --------

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
                           scaled_biodiversity_utility=analysis_data$scaled_biodiversity_utility)

model_binary <- lm(scaled_biodiversity_utility ~ `Pavilion/Picnic Area` + `Kids Playground` + `Body of Water` +
                    `Jog/Walk Path` + `Athletic Facility` + `Nature Preserve` + `Dog Park` + `Indoor/Outdoor Fitness Center` +
                     park_size, 
                   data=binary_data)
summary(model_binary)
anova(model_binary)
AIC(model_binary)
plot(model_binary)

adj.r2 <- summary(model_binary)
adj.r2$adj.r.squared

# Plot Bio-Use by Binary Human Utility Features --------

(a <- effect_plot(model_binary, pred=`Pavilion/Picnic Area`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Pavilion/Picnic Area"))
(b <- effect_plot(model_binary, pred=`Kids Playground`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Kids Playground*"))
(c <- effect_plot(model_binary, pred=`Body of Water`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Body of Water*"))
(d <- effect_plot(model_binary, pred=`Jog/Walk Path`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Jog/Walk Path"))
(e <- effect_plot(model_binary, pred=`Athletic Facility`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Athletic Facility**"))
(f <- effect_plot(model_binary, pred=`Nature Preserve`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Nature Preserve**"))
(g <- effect_plot(model_binary, pred=`Dog Park`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Dog Park"))
(h <- effect_plot(model_binary, pred=`Indoor/Outdoor Fitness Center`, ylim=c(0.17, 0.45), cat.geom="line",  y.label="Bio-use Model Prediction", x.label="", main.title="Indoor/Outdoor Fitness Center"))

figure <- ggarrange(a + rremove("ylab"), b + rremove("ylab"), c + rremove("ylab"), d + rremove("ylab"), 
                    e + rremove("ylab"), f + rremove("ylab"), g + rremove("ylab"), h + rremove("ylab"), ncol=2, nrow=4)
annotate_figure(figure, 
                left=textGrob("Biodiversity", rot=90, vjust=1, gp = gpar(cex = 1.2))) 
ggsave("Figures/Figure4.jpg", height=10, width=8, units="in", dpi=300)

# Summary of Data ---------------------------------------------------------

# summary of greenspace area
summary(analysis_data$greenspace_area_m2/10000)

# summary of total_utility
summary(analysis_data$total_utility)

# calculate the percentage occurence of each human utility attribute
percentage <- data.frame(picnic=sum(analysis_data$`Pavilion/Picnic Area`),
                         playground=sum(analysis_data$`Kids Playground`),
                         water=sum(analysis_data$`Body of Water`),
                         walk=sum(analysis_data$`Jog/Walk Path`),
                         athletic=sum(analysis_data$`Athletic Facility`),
                         nature=sum(analysis_data$`Nature Preserve`),
                         dog=sum(analysis_data$`Dog Park`),
                         fitness=sum(analysis_data$`Indoor/Outdoor Fitness Center`))

percentage %>% 
  mutate(across(everything(), ~ ./rowSums(percentage)*100))

# calculate greenspaces per municipality
summary(analysis_data %>% 
          group_by(mncplty) %>%
          summarise(count=n()))

# Make supplementry table -------------------------------------------------

# make a supplementry file summarizing biodiveristy utility and human utility for each park
hum_ut_table <- analysis_data %>%
                    dplyr::select(poly_id, scaled_human_utility)

sup_table_join <- left_join(biodiversity_utility, hum_ut_table, by=c("poly_id"))

sup_table_join <- sup_table_join %>%
  separate(poly_id, into = c("Municipality", "Park_Name"), sep = "_", extra = "merge") %>%
  mutate(
    scaled_human_utility = round(scaled_human_utility, 2),
  ) %>%
  rename(`Park Size (ha)`=park_size,
         `Number of iNaturalist Observations`=iNat_observations,
         `Number of iNaturalist Users`=iNat_users,
         `Predicted Richness`=predicted_richness,
         `Human Utility Attribute Index`=scaled_human_utility)  %>%
  mutate(Municipality=ifelse(Municipality=="Broward", "Undefined", Municipality))

write_csv(sup_table_join, "Data/bio_use_predictors/supplemental_table.csv")
