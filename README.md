This repository contains all the code and files used in the analyses reported in "Urban greenspaces benefit both human utility and biodiversity."

# R Folder

Within this folder are three scripts to (1) summarize iNaturalist data by greenspace, (2) calculate the biodiversity index value, and (3) conduct the main analyses presented in the paper. The numbers calculated in these scripts may not perfectly align with those presented in the paper due to differences in randomization. However, values and figures will produce the same results as presented.

Additionally contained in this folder is a sub folder titled "Alternative_bio_use_calculation_methods." Here there are two scripts that are used to (1) calculate alternative methods to obtain the biodiversity index and (2) run the main analyses individually for each method.

# Data folder

**Broward_polygons_metadata:** A CSV file containing the municipality and park name.

**Clean_GBIF:** A RDS file containing Research Grade iNaturalist observations sourced from the Global Biodiversity Information Facility (GBIF).

**Human_utility:** A CSV file containing data on the presence (1) or absence (0) of eight human utility attributes.

**iNaturalist_data:** A RDS file containing all iNaturalist data, sourced from www.iNaturalist.org.

**iNaturalist_park:** The CSV output file from 1_Summarize_iNaturalist_Data_by_Greenspace.R. This contains the summarized iNaturalist data for each park.

**bio_use_predictors:** A folder containing CSV files that have biodiversity utility prediction values.

**rasters:** A folder containing raster files of impervious cover, percent non-tree cover, tree cover, and water cover.

**shapefiles:** A folder containing shapefiles on Broward County municipalities and urban parks.

# Figures folder

This folder contains all figures presented in the paper and outputted by R scripts.
