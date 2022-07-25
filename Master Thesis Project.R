
Script of R Coding for this study
'Time series analysis of the deciduous forest degradation in the Berlin-Brandenburg region based on Sentinel-2 data'

Prepared by: Kazi Sajjad Hossain
             M.Sc.Student, Environmental Planning 
             Technical University of Berlin 
################################################################################
                     Cloud Masking
Source: This code was adapted from https://blogs.fu-berlin.de/reseda/quality-band/
################################################################################

#install necessary packages and libraries
install.packages("raster")
install.packages("rgdal")

library(raster)
library(rgdal)

#Setting up the working directrory
setwd("D:\\Thesis\\Resampling\\L2A_T33UUU_20190512T102029") 

#Open and plot the image
sen2 <- stack("T33UUU_20190512T102029_B01_13_resampled.tif")
plot(sen2)
plotRGB(sen2, 4,3,2, stretch="lin")

#Seperate spectral bands and classification (band 1 to 12 is the mulispectral Senitnel- 2 scene, band 13 is the quality classification band)
sen2_bands <- sen2[[-13]]
sen2_mask <- sen2[[13]]
plot(sen2_mask)
sen2_mask_combi <- sen2_mask

#masking the pixels
sen2_mask_combi[sen2_mask == 3 |sen2_mask == 8 | sen2_mask == 9 |sen2_mask == 10 |sen2_mask == 11] <- NA
plot(sen2_mask_combi)
writeRaster(sen2_mask_combi, "sen2_mask.tif")

#Apply mask
sen2_bands_masked <- raster::mask(sen2_bands,sen2_mask_combi)
plotRGB(sen2_bands_masked, 4,3,2,stretch="lin")
writeRaster(sen2_bands_masked, "T33UUU_20190512T102029_B01_13_resampled_masked_10m.tif")

Note: The above code was adapted for other images as well


#################################################################################
Calculation of indices for all the images together using 'for' loop in RStudio
#################################################################################

#setting working directory
setwd("D:\\Thesis\\Resampling\\All_Images_2017_2020\\Only_tif")

#NDVI
File <- list.files() #make a list of all files path
print(File)

for(i in File){
  sen2 <- stack(i)
  str1 = gsub(".tif", "_NDVI",i) #changing the file name 
  str2 = '.tif'
  filename= paste(str1,str2) #adding two strings as one string
  VI <- function(img, a, b) {
    NIR <- img[[a]] 
    Red <- img[[b]] 
    vi <- (NIR-Red)/ (NIR+Red)
    return(vi)
  }
  
  ndvi <- VI(sen2,8,4)
  plot(ndvi)
  writeRaster(ndvi,filename)
}

#DWSI
File <- list.files() #make a list of all files path
print(File)

for(i in File){
  sen2 <- stack(i)
  str1 = gsub(".tif", "_DWSI",i) #changing the file name 
  str2 = '.tif'
  filename= paste(str1,str2) #adding two strings as one string
  
  VI <- function(img, a, b,c,d) {
    NIR <- img[[a]]
    Green <- img[[b]]
    SWIR <- img[[c]]
    Red <- img[[d]]
    vi <- (NIR+Green)/(SWIR+Red)
    return(vi)
  }
  
  dwsi <- VI(sen2,8,3,11,4)
  plot(dwsi)
  writeRaster(dwsi,filename)
}

#NDMI
File <- list.files() #make a list of all files path
print(File)


for(i in File){
  sen2 <- stack(i)
  str1 = gsub(".tif", "_NDMI",i) #changing the file name 
  str2 = '.tif'
  filename= paste(str1,str2) #adding two strings as one string
  VI <- function(img, a, b) {
    NIR <- img[[a]] 
    SWIR <- img[[b]] 
    vi <- (NIR-SWIR)/ (NIR+SWIR)
    return(vi)
  }
  
  ndmi <- VI(sen2,8,11)
  plot(ndmi)
  writeRaster(ndmi,filename)
}

#PSRI
File <- list.files() #make a list of all files path
print(File)

for(i in File){
  sen2 <- stack(i)
  str1 = gsub(".tif", "_PSRI",i) #changing the file name 
  str2 = '.tif'
  filename= paste(str1,str2) #adding two strings as one string
  VI <- function(img, a, b,c) {
    Blue <- img[[a]] 
    Red <- img[[b]] 
    Red_edge_1 <- img[[c]] 
    
    
    vi <- (Red-Blue)/ Red_edge_1
    return(vi)
  }
  
  psri <- VI(sen2,2,4,6)
  plot(psri)
  writeRaster(psri,filename)
}


################################################################################
Extracting pixel value as a dataframe with defined extent (shapefile) from an 
arbitrary numberof image layers of arbitrary resolutions.This code is only for 
one image with three layers.It is also adapted for rest of the images 

Source: This code was adapted from  https://www.youtube.com/watch?v=Z2gMTVumJZo
################################################################################

#setting up working directory
setwd('D:\\Thesis\\Resampling\\Indices\\DWSI') #path of the raster file

install.packages("sf")
install.packages("snow")
#activate necessary packages
library(raster)
library(sf) # shapefile packaage that should come with rgdal
library(rgdal)
library(snow) #for parallel processing though this may be outdated
#try this package that does not require snow (https://github.com/isciences/exactextractr)

##name of shapefile containing polygons to use to extract pixel values
shapefile_filename = 'D:\\Thesis\\Resampling\\AOI\\sample_areas.shp' 

#set the fieldnames with unique ID/name and the class field
uid_fieldname = 'main_bioto'   
class_fieldname = 'name' 



#image with any number of bands whose pixels will be extracted
im_filename = 'T33UUU_DWSI_June.tif' ##A composite image of the month June for 2018, 2019, 2020 

Note: This code also adapted for the composite image 'T33UUU_DWSI_July.tif' of july 2018,2019,2020 

# Loading shapefile
polys <- readOGR(shapefile_filename)


#displaying key information like how many features and crs
polys
#just show the attribute column names
names(polys@data)

#Loading images
im <- stack(im_filename)
num_layers <- nlayers(im)

# S4 method for Raster,SpatialPolygons
#small=TRUE allows for extraction of pixels that are larger than the polygon
#df=TRUE outputs as data frame
#normalizeWeights=FALSE and weights=TRUE returns a column saying what percent of that pixel is in the poly
#nl = number of layers in stack from which you want to extract
#beginCluster(num_nodes) do not use for now

values_df_all <- extract(im, polys, small=TRUE, 
                         df=TRUE, weights=TRUE, normalizeWeights=FALSE,
                         nl=num_layers)
#create dummy df that has the extract IDs (1:number of polygons) aligned with
#the polygon unique IDs/name
#first subset the polys shapefile by column, only finding the UniqueID/Unique name and class name cols
uid_and_class_data <- polys[,(names(polys)==uid_fieldname | names(polys)==class_fieldname)]
#this should have one column that is the polygon UniqueID/Unique name, and one that is class code
poly_uid_df <- data.frame(uid_and_class_data@data)
#add the 1:n column to the the linking ID
num_rows <- data.frame("ID"=1:nrow(poly_uid_df))
poly_uid_df <- cbind(num_rows, poly_uid_df)

#this outer join gets the species and crown UID information attached
#in a one-to-many fashion to the pixel level results of values_df_all
all_crowns_df<-merge(x=values_df_all,y=poly_uid_df,by="ID",all=TRUE)

#exporting the df as csv file in working directory
write.csv(all_crowns_df, file = 'dwsi_sample_3dates_all_poly_june.csv')







################################################################################ 
                     #Calculating boxplot statitics 
################################################################################


install.packages("ggpubr")

library(ggpubr)
library(tidyr)

setwd("D:\\Thesis\\R coding")

#NDVI
data<- read.csv('sample_15_statistics_for_single_dates - jun_july_NDVI.csv')
print(data)


ggboxplot(data,'Date','pixel_value')+
  stat_compare_means(label.x = 1.5,method = 'wilcox.test',
                     comparisons =  list(c("06.06.2018", "03.06.2019"), 
                                         c("06.06.2018", "02.06.2020"),
                                         c("03.06.2019", "02.06.2020"),
                                         c("31.07.2018", "26.07.2019"),
                                         c("31.07.2018", "30.07.2020"),
                                         c("26.07.2019", "30.07.2020"),
                                         c("06.06.2018", "31.07.2018"),
                                         c("03.06.2019", "26.07.2019"),
                                         c("02.06.2020", "30.07.2020")))

##only_significant(NDVI) p???0.01
ggboxplot(data,'Date','pixel_value')+
  stat_compare_means(label.x = 1.5,method = 'wilcox.test',
                     comparisons =  list(c("03.06.2019", "02.06.2020")))

#DWSI
data<- read.csv('sample_15_statistics_for_single_dates - jun_july_DWSI.csv')
print(data)


ggboxplot(data,'Date','pixel_value')+
  stat_compare_means(label.x = 1.5,method = 'wilcox.test',
                     comparisons =  list(c("06.06.2018", "03.06.2019"), 
                                         c("06.06.2018", "02.06.2020"),
                                         c("03.06.2019", "02.06.2020"),
                                         c("31.07.2018", "26.07.2019"),
                                         c("31.07.2018", "30.07.2020"),
                                         c("26.07.2019", "30.07.2020"),
                                         c("06.06.2018", "31.07.2018"),
                                         c("03.06.2019", "26.07.2019"),
                                         c("02.06.2020", "30.07.2020")))

#only_significant(DWSI) p???0.01

ggboxplot(data,'Date','pixel_value')+
  stat_compare_means(label.x = 1.5,method = 'wilcox.test',
                     comparisons =  list(c("06.06.2018", "02.06.2020"),
                                         c("31.07.2018", "26.07.2019"),
                                         c("02.06.2020", "30.07.2020")))




#NDMI
data<- read.csv('sample_15_statistics_for_single_dates - jun_july_NDMI.csv')
print(data)


ggboxplot(data,'Date','pixel_value')+
  stat_compare_means(label.x = 1.5,method = 'wilcox.test',
                     comparisons =  list(c("06.06.2018", "03.06.2019"), 
                                         c("06.06.2018", "02.06.2020"),
                                         c("03.06.2019", "02.06.2020"),
                                         c("31.07.2018", "26.07.2019"),
                                         c("31.07.2018", "30.07.2020"),
                                         c("26.07.2019", "30.07.2020"),
                                         c("06.06.2018", "31.07.2018"),                                         
                                         c("03.06.2019", "26.07.2019"),
                                         c("02.06.2020", "30.07.2020")))
#only_significant(NDMI)  p???0.01
ggboxplot(data,'Date','pixel_value')+
  stat_compare_means(label.x = 1.5,method = 'wilcox.test',
                     comparisons =  list(c("31.07.2018", "26.07.2019"),
                                         c("03.06.2019", "26.07.2019")))

#PSRI
data<- read.csv('sample_15_statistics_for_single_dates - jun_july_PSRI.csv')
print(data)


ggboxplot(data,'Date','pixel_value')+
  stat_compare_means(label.x = 1.5,method = 'wilcox.test',
                     comparisons =  list(c("06.06.2018", "03.06.2019"), 
                                         c("06.06.2018", "02.06.2020"),
                                         c("03.06.2019", "02.06.2020"),
                                         c("31.07.2018", "26.07.2019"),
                                         c("31.07.2018", "30.07.2020"),
                                         c("26.07.2019", "30.07.2020"),
                                         c("06.06.2018", "31.07.2018"),                                         
                                         c("03.06.2019", "26.07.2019"),
                                         c("02.06.2020", "30.07.2020")))

#only_significant(PSRI)  p???0.01
ggboxplot(data,'Date','pixel_value')+
  stat_compare_means(label.x = 1.5,method = 'wilcox.test',
                     comparisons =  list(c("31.07.2018", "26.07.2019"),
                                         c("03.06.2019", "26.07.2019")))
