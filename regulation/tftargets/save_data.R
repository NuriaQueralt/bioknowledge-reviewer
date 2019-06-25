# This script is to save each object data list into a csv file.
# imports
library(jsonlite)

# load data
rm(list = ls())
load("data/tftargets.rda")

# save each dataset into a different file
datasets <- ls()
for (data in datasets){
	write(toJSON(get(data), pretty=TRUE),paste0("data/",tolower(data),".json"))
}

