#library(plyr)
require(tidyverse)
#require(dplyr)

## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  #library(plyr)

  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }

  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )

  # Rename the "mean" column
  datac <- rename(datac, c("mean" = measurevar))

  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult

  return(datac)
}


d <- readxl::read_xlsx("analysis/data/raw_data/NORA Konferenz Auswertung Fragebogen.xlsx")
d %>% str()
names(d)
clean <- d %>%
  select(1:44) %>%
  tidyr::gather(key = "symptom",value = "VAS",6:44) %>%
  mutate(VAS = as.numeric(VAS),
         experience = as.numeric(`Erfahrung in jahren`),
         symptom = factor(symptom))
clean$VAS[667] <- 55.6
clean[which(clean$symptom == "death"&clean$VAS <70),]
clean <- clean %>%
  filter(ID != 23)

clean %>%
  ggplot(aes(reorder(symptom,VAS),VAS))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 70,hjust = 1))

MEAN <- clean %>%
  filter(complete.cases(.)) %>%
  #split(.$symptom)
  group_by(symptom) %>%
  #summarize()
  summarise(mean_VAS = mean(VAS,na.rm = T),
            weighted_mean_VAS = weighted.mean(x = VAS, w=experience,na.rm = T),
            sd_VAS = sd(VAS, na.rm=T))
MEAN %>%
  ggplot(aes(reorder(symptom,weighted_mean_VAS),weighted_mean_VAS))+
           geom_point()+
  theme(axis.text.x = element_text(angle = 70,hjust = 1))+
  geom_errorbar(aes(ymin = weighted_mean_VAS-sd_VAS, ymax =weighted_mean_VAS+ sd_VAS))+
  labs(x = "Symptoms", y = "severity according to weighted mean VAS")


MEAN %>%
  ggplot(aes(reorder(symptom,sd_VAS),sd_VAS))+
  geom_point()+
  theme(axis.text.x = element_text(angle = 70,hjust = 1))

clean %>%
  filter(symptom =="angioedema") %>%
  ggplot(aes(Land,experience))+
  geom_jitter()+
  theme(axis.text.x = element_text(angle = 70,hjust = 1))




