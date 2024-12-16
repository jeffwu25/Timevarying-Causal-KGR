--------------------------------------------------------------------------------
County
--------------------------------------------------------------------------------
county/tl_2019_us_county/:
This is a folder that contains the shapefile for CONUS counties in 2019. Files were downloaded from the US Census Bureau TIGER/Line Shapefiles website (https://www.census.gov/cgi-bin/geo/shapefiles/index.php). R users may also use the `tigris` package.

--------------------------------------------------------------------------------
county/smokePM2pt5_predictions_daily_county_20060101-20201231.rds:
This is a file that contains a data frame with the final set of daily smoke PM2.5 predictions on smoke days at the county level from January 1, 2006 to December 31, 2020 for the contiguous US. County-level smoke PM2.5 predictions are aggregated from smoke PM2.5 predictions at the 10 km resolution using population and area of intersection-weighted averaging. The 'GEOID' column in this file corresponds to the 'GEOID' column in the county shapefile.

All rows in this file are predictions on smoke days. Predictions on non-smoke days are by construction 0 ug/m^3 and not included in this file. A smoke PM2.5 prediction of 0 in this file means that the county-day did have a smoke day but did not have elevated PM2.5. The full set of smoke PM2.5 predictions on both smoke days and non-smoke days can be obtained by setting the smoke PM2.5 prediction to 0 on county-days in the counties and in the January 1, 2006-December 31, 2020 date range that are not in this file. For example, the R code below returns the full set of smoke PM2.5 predictions:

library(lubridate)
library(sf)
library(dplyr)
library(tidyr)

# Load smokePM predictions on smoke days
preds = readRDS("./final/county/smokePM2pt5_predictions_daily_county_20060101-20201231.rds")

# Load counties
counties = read_sf("./final/county/tl_2019_us_county")

# Load full set of dates
dates = seq.Date(ymd("20060101"), ymd("20201231"), by = "day")

# Get full combination of county-days
# Warning: this may require a large amount of memory
out = expand.grid(GEOID = counties$GEOID, date = dates)

# Match smokePM predictions on smoke days to county-days
out = left_join(out, preds, by = c("GEOID", "date"))

# Predict 0 for remaining county-days, which are non-smoke days
out = mutate(out, smokePM_pred = replace_na(smokePM_pred, 0))

--------------------------------------------------------------------------------
county/smokePM2pt5_predictions_daily_county_20060101-20201231.csv:
This is the same as smokePM2pt5_predictions_daily_county_20060101-20201231.rds, except it is saved as a CSV file.
