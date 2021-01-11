# pkgTest is a helper function to load packages and install packages only when they are not installed yet.
pkgTest <- function(x)
{
  if (x %in% rownames(installed.packages()) == FALSE) {
    install.packages(x, dependencies= TRUE)
  }
  library(x, character.only = TRUE)
}
neededPackages <- c("strucchange","zoo", "bfast", "raster", "leaflet", "MODISTools")
for (package in neededPackages){pkgTest(package)}
# Function to create time series object
# val_array: data array for one single pixel (length is number of time steps)
# time_array: array with dates at which raster data is recorded (same length as val_array)
timeser <- function(val_array, time_array) {
  z <- zoo(val_array, time_array) # create zoo object
  yr <- as.numeric(format(time(z), "%Y")) # extract the year numbers
  jul <- as.numeric(format(time(z), "%j")) # extract the day numbers (1-365)
  delta <- min(unlist(tapply(jul, yr, diff))) # calculate minimum time difference (days) between observations
  zz <- aggregate(z, yr + (jul - 1) / delta / 23) # aggregate into decimal year timestamps
  (tso <- as.ts(zz)) # convert into timeseries object
  return(tso)
}
# Downloading the NDVI data, starting from 2000-01-01
VI <- mt_subset(product = "MOD11A2 ",
                site_id = "nl_gelderland_loobos",
                band = "250m_16_days_NDVI",
                start = "2000-01-01",
                km_lr = 2,
                km_ab = 2,
                site_name = "testsite",
                internal = TRUE,
                progress = FALSE)

# Downloading the pixel reliability data, starting from 2000-01-01
QA <- mt_subset(product = "MOD13Q1",
                site_id = "nl_gelderland_loobos",
                band = "250m_16_days_pixel_reliability",
                start = "2000-01-01",
                km_lr = 2,
                km_ab = 2,
                site_name = "testsite",
                internal = TRUE,
                progress = FALSE)

# convert df to raster
VI_r <- mt_to_raster(df = VI)
QA_r <- mt_to_raster(df = QA)
## polish the data
# set all values to be set at <0 or >1 NA after creating mask on pixel reliability flag
m <- QA_r
m[(QA_r < 0 | QA_r > 1)] <- NA # continue working with QA 0 (good data), and 1 (marginal data)

# apply the mask to the NDVI raster
VI_m <- mask(VI_r, m, maskvalue=NA, updatevalue=NA)

# plot the first image
plot(m,1) # plot mask
plot(VI_m,1) # plot  NDVI raster
# extract data for selected pixels from the cleaned raster.
click(VI_m, , xy=TRUE, id=TRUE, n=1, cell=TRUE,)

library(leaflet)
r <- raster(VI_m,1)
pal <- colorNumeric(c("#ffffff", "#4dff88", "#004d1a"), values(r),
                    na.color = "transparent")

map <- leaflet() %>% addTiles() %>%
  addRasterImage(r, opacity = 0.8, colors = pal,) %>%
  addLegend(pal = pal, values = values(r),
            title = "NDVI")
map
## at a specific pixel for example 1 row, complete left-hand side check VI data. 
## raster dimensions are: 33x33

px <- 78 # pixel number so adjust this number to select the center pixel
tspx <- timeser(as.vector(VI_m[px]),as.Date(names(VI_m), "X%Y.%m.%d")) # convert pixel 1 to a time series
plot(tspx, main = 'NDVI') # NDVI time series cleaned using the "reliability information"

bfm1 <- bfastmonitor(tspx, response ~ trend + harmon, order = 3, start = c(2019,1)) # Note: the first observation in 2019 marks the transition from 'history' to 'monitoring'
plot(bfm1)
dates <- as.Date(names(VI_m), "X%Y.%m.%d")

# here we define the function that we will apply across the brick using the calc function:
bfmRaster = function(pixels)
{
  tspx <- timeser(pixels, dates) # create a timeseries of all pixels
  bfm <- bfastmonitor(tspx, response ~ trend + harmon, order = 3, start = c(2019,1)) # run bfast on all pixels
  return(c(bfm$breakpoint, bfm$magnitude)) 
}

# calc function 
bfmR <- calc(VI_m, bfmRaster)
names(bfmR) <- c('time of break', 'magnitude of change')
plot(bfmR) # resulting time and magnitude of change
plot(bfmR,1)
click(VI_m, id=FALSE, xy=FALSE, cell=TRUE, n=1)
px <- 460 # pixel number so adjust this number to select the center pixel
tspx <- timeser(as.vector(VI_m[px]),as.Date(names(VI_m), "X%Y.%m.%d")) # convert pixel 1 to a time series
plot(tspx, main = 'NDVI') # NDVI time series cleaned using the "reliability information"

tspx[tspx < 0] <- NA
bfm <- bfastmonitor(tspx, response ~ trend + harmon, order = 3, start = c(2019,1))
plot(bfm)

pkgTest("remotes")
remotes::install_github("bfast2/bfast1")
detache("package:bfast", unload=TRUE)
library(bfast) # always restart your R session in case of any errors experienced then rerun the above code.
breaks <- bfast0n(tspx, response ~ trend + harmon, order = 3)
breaks
## 
##   Segment partition: Optimal 1
## 
## Call:
## breakpoints. Formula (formula = formula, data = data_pp)
## 
## Observation number at break points:
## NA 
## 
## Corresponding points to breakdates:
## NA
px <- 92
tspx <- timeser(as.vector(VI_m[px]),as.Date(names(VI_m), "X%Y.%m.%d"))
breaks <- bfast0n(tspx, response ~ trend + harmon, order = 3)
breaks
## 
##   Segment partition- Optimal 2: 
## 
## Call:
## breakpoints.formula(formula = formula, data = data_pp)
## 
## Observation number breaakpoints:
## 0327 
## 
## Corresponding points to breakdates:
## 0.80147
dates.no.na <- as.numeric(time(tspx))
dates.no.na[is.na(tspx)] <- NA
dates.no.na <- na.omit(dates.no.na)
dates.no.na[breaks$breakpoints[1]]
# Plot tspx it
abline(v=dates.no.na[breaks$breakpoints[1]], col="red")
plot(tspx)
library(bfast)
## a demo of ndvi series of time:
ndvi <- ts(rowSums(simts$time.series))
tsp(ndvi) <- tsp(simts$time.series)
## input variable for the sinus and cosinus functions
f <- 23
w <- 1/f
tl <- 1:length(ndvi)
## 3th order harmonic model
co <- cos(2 * pi * tl * w)
si <- sin(2 * pi * tl * w)
co2 <- cos(2 * pi * tl * w * 2)
si2 <- sin(2 * pi * tl * w * 2)
co3 <- cos(2 * pi * tl * w * 3)
si3 <- sin(2 * pi * tl * w * 3)
# fit the seasonal model using linear regression
fitm<- lm(ndvi~co+si+co2+si2+co3+si3) 
predm <- fitted(fitm) ## predict based on the modelfit
plot(co, type = "l", ylab = "cos and sin")
lines(si, type = "l", lty = 2)
predm <- ts(as.numeric(predm), start=c(2000,4), frequency=23) 
plot(ndvi, lwd = 3, col = "grey", ylab = "NDVI")
lines(predm, type = "l", col = "red") # fitted

