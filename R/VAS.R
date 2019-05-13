require(tidyverse)
#require(dplyr)
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
  summarize(mean_VAS = mean(VAS,na.rm = T),
            weighted_mean_VAS = weighted.mean(x = VAS, w=experience,na.rm = T),
            sd_VAS = sd(VAS))
MEAN %>%
  ggplot(aes(reorder(symptom,weighted_mean_VAS),weighted_mean_VAS))+
           geom_point()+
  theme(axis.text.x = element_text(angle = 70,hjust = 1))

MEAN %>%
  ggplot(aes(reorder(symptom,sd_VAS),sd_VAS))+
  geom_point()+
  theme(axis.text.x = element_text(angle = 70,hjust = 1))

clean %>%
  filter(symptom =="angioedema") %>%
  ggplot(aes(Land,experience))+
  geom_jitter()+
  theme(axis.text.x = element_text(angle = 70,hjust = 1))




