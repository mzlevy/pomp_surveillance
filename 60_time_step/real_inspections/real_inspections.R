length(which(all_active_data$INSP_POSITIVA == 1))
all_active_data$YEAR <- sapply(all_active_data$FECHA, function(x) strsplit(x, "-")[[1]][1])
unique(all_active_data$YEAR)
all_active_data$MONTH <- sapply(all_active_data$FECHA, function(x) strsplit(x, "-")[[1]][2])
unique(all_active_data$MONTH)

length(which(all_staticInsp$INSP_POSITIVA == 1))
all_staticInsp$YEAR <- sapply(all_staticInsp$FECHA, function(x) strsplit(x, "/")[[1]][3])
unique(all_staticInsp$YEAR)
all_staticInsp$MONTH <- sapply(all_staticInsp$FECHA, function(x) strsplit(x, "/")[[1]][1])
unique(all_staticInsp$MONTH)

years <- c("2017", "2018", "2019")
months <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
data <- list()
data_dex <- 1
for (year in years) {
  for (month in months) {
    static_year <- NA
    if (year == "2017") {
      static_year <- "17"
    } else if (year == "2018") {
      static_year <- "18"
    } else if (year == "2019") {
      static_year <- "19"
    }
    static_month <- NA
    if (month == "01") {static_month <- "1"}
    else if (month == "02") {static_month <- "2"}
    else if (month == "03") {static_month <- "3"}
    else if (month == "04") {static_month <- "4"}
    else if (month == "05") {static_month <- "5"}
    else if (month == "06") {static_month <- "6"}
    else if (month == "07") {static_month <- "7"}
    else if (month == "08") {static_month <- "8"}
    else if (month == "09") {static_month <- "9"}
    else if (month == "10") {static_month <- "10"}
    else if (month == "11") {static_month <- "11"}
    else if (month == "12") {static_month <- "12"}
    
    if (!(month == "07" & year == "2017")) {
      count_I <- 0
      temp_active <- all_active_data[which(all_active_data$YEAR == year & all_active_data$MONTH == month),]
      temp_static <- all_staticInsp[which(all_staticInsp$YEAR == static_year & all_staticInsp$MONTH == static_month),]
      count_I <- length(which(temp_active$INSP_POSITIVA == 1)) + length(which(temp_static$INSP_POSITIVA == 1))
    }
    
    data[[data_dex]] <- count_I
    data_dex <- data_dex + 1
  }
}
data <- do.call(rbind, data)
data <- data[1:34]
data <- cbind(1:34, data)
data <- cbind(data, data[,2])
colnames(data) <- c("day", "true", "reports")
data <- as.data.frame(data)
write.csv(data, "~/PETM-shiny/pomp/real_infestation_data.csv")

# Argument to run in pomp
data <- read.csv("~/PETM-shiny/pomp/real_infestation_data.csv", row.names = FALSE)
chagas <- data
