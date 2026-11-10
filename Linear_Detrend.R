library(ggplot2)
setwd('/home/preston/Projects/Daycent/Crop_Optimization/Crop_Yield')

files=list.files()
files=files[-c(5)]

yields=lapply(files, read.csv)

files=substr(files, 1, nchar(files)-11)
names(yields)=files

yield_avg <- lapply(yields, function(x) {
  aggregate(x[,3:5], by = list(Year = x$Year), FUN = mean)
})


detrend_to_2019 <- function(x) {
  for (col in c("grain_c", "straw_c", "roots_c")) {
    fit <- lm(x[[col]] ~ x$Year)
    slope <- coef(fit)[2]
    x[[col]] <- x[[col]] + slope * (2019 - x$Year)
  }
  return(x)
}

yields_2019 <- lapply(yields, detrend_to_2020)

yield_avg <- lapply(yields_2019, function(x) {
  aggregate(x[,3:5], by = list(Year = x$Year), FUN = mean)
})

#check slopes
for (crop in names(yields)) {
  cat("\n---", crop, "---\n")
  for (col in c("grain_c", "straw_c", "roots_c")) {
    slope_before <- coef(lm(yields[[crop]][[col]] ~ yields[[crop]]$Year))[2]
    slope_after  <- coef(lm(yields_2019[[crop]][[col]] ~ yields_2019[[crop]]$Year))[2]
    cat(col, "| slope before:", round(slope_before, 4),
        "| slope after:", round(slope_after, 4), "\n")
  }
}  


#export detrended
setwd('/home/preston/Projects/Daycent/Crop_Optimization/Crop_Yield')

files=list.files()
files=files[-c(3,6)]

setwd('/home/preston/Projects/Daycent/Crop_Optimization/Crop_Yield/detrend')

mapply(function(data, name) {
  write.csv(data, file = name, row.names = FALSE)
}, yields_2019, files)