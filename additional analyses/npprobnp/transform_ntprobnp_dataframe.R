names(datED)[14] <- "ntprobnp"
write.csv(datED,"ntprobnp_results.csv")
data <- read.csv("ntprobnp_results.csv")
data$ntprobnp <- as.character(data$ntprobnp) %>% as.numeric()
require(tidyverse)

data$yob <- data$b_patient_code %>% gsub(pattern = "m-(....)-.*",
                             replacement = "\\1")
data$yob <- paste0(data$yob,"-01-01")
data$yob <- as.Date(data$yob)

data$Datum.der.BE <- as.Date(data$Datum.der.BE, origin = "1900-01-01")
time_length(difftime(data$Datum.der.BE, data$yob), "years")


data %>%
  ggplot(aes(x = q_114_hypotension_collapse_v5, y = ntprobnp))+
  geom_boxplot()

data %>%
    ggplot(aes(x = nach.Müller, y = ntprobnp))+
  geom_boxplot()

