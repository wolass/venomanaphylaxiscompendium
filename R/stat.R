#####Setup####
require(dplyr)
require(magrittr)
require(ggplot2)
require(forestplot)
require(summarytools)
require(vcd)
require(MatchIt)
require(tidyr)
library(ggradar)
library(tibble)
library(scales)
require(ggpubr)
require(rlang)
library(forcats)
require(fmsb)
require(purrr)


#### GET THE DATA ########
# Note the path that we need to use to access our data files when rendering this document
# data5 <-foreign::read.spss("analysis/data/raw_data/anaphylaxis_registry_15_mar_2019_mergedSD.sav",
#                            use.value.labels = TRUE,
#                            to.data.frame = T,
#                            trim.factor.names = T,
#                            trim_values = T)

#load('../RefractoryAnaOrg/analysis/data/raw_data/data4.R')
load("data.R") # here the anascore is done using the severity_analysis R script

# update the anascore according to nora consensus
#source("R/anascore_vas.R")

data4 <-data

manual_greens <- rev(c("#E5F5E0",
                       "#C7E9C0",
                       "#A1D99B",
                       "#74C476",
                       "#41AB5D",
                       "#238B45",
                       "#005A32"))

####### General Functions #############
source("R/functions.R")

#### Data cleaning ####
data4 %<>% correctLabels1()

data4$reaction_type_brown <- AssessAnaphylaxisDefinition(data4)

data4$severity_brown <- severity_brown(data4)

other_insects <- c("horse fly",
                   "bumble-bee",
                   "mosquito",
                   "insects")

data4$d_insect_gr4 <- as.character(data4$q_340_insects)
data4$d_insect_gr4[which(data4$q_340_insects%in%other_insects)] <- "other"
data4$d_insect_gr4 %<>% factor
data4%<>%
  mutate(d_AAI_prescribed = ifelse(q_540_why_autoinj_v5 %in% c("prescribed, available, but not used",
                                                               "prescribed, but not available"),
                                   "yes",
                                   ifelse(q_540_why_autoinj_v5 =="not prescribed",
                                          "no",
                                          NA)
  ))

data4$d_111_sum_wf <- data$d_111_sum %>% char_to_num_levels()
data4$d_112_sum_wf <- data$d_112_sum %>% char_to_num_levels()
data4$d_113_sum_wf <- data$d_113_sum %>% char_to_num_levels()
data4$d_114_sum_wf <- data$d_114_sum %>% char_to_num_levels()

data4$q_200_height_v7 %<>% as.character %>% as.numeric
data4$q_200_weight_v7 %<>% as.character %>% as.numeric

data4$d_severity_rm %<>% as.numeric()
data4$d_organ_sum %<>% car::recode("'none' = '0';
                                   'one' = '1';
                                   'two' = '2';
                                   'three' = '3';
                                   'four' = '4';
                                   'unknown' = NA") %>% as.numeric()
data4$q_116_VAS_v7 %<>% as.numeric()
data4$q_120_time_between_v4%<>% as.numeric()
data4$q_131_biphasic_time_v4 %<>% car::recode("'b_version<4' = NA;
                                              'no biphasic reaction' = NA;
                                              'unknown' = NA;
                                              '4 – 12 hours after initial symptoms' = '1';
                                              '12 – 24 hours' = '2';
                                              'more than 24 hours' = '3';
                                              '2 -4 hours (V4)' = '0'") %>% as.numeric()
levels(data4$q_141_fatal_time_v5) %<>%  {c(NA,.[-1])}
data4$q_141_fatal_time_v5 %<>% as.numeric()
levels(data4$q_561_hospital_admission_v6) %<>% {c(NA,.[2:3],NA)}
levels(data4$q_562_intensive_care_v6)%<>% {c(NA,.[2:3],NA)}
levels(data4$q_632_autoinj_number_v7) %<>% {c(.[1:3],NA)}
levels(data4$q_422_stress) <- c("no","yes")
levels(data4$q_340_insects) <-c(levels(data4$q_340_insects)[1:7],"other")

### REstricted ring and messmer (conformat with the RUEFF work)
data4$d_severity_rmr <- ifelse(data4$d_severity_rm%in%c(3,4),
                               "severe",
                               "mild")
### ATOPIC VARIABLE (posible many NAs)
data4$atopy <- ifelse(data4$q_410_rhinitis_cur=="yes"|
                        data4$q_410_rhinitis_prev_v5=="yes"|
                        data4$q_410_asthma_cur=="yes"|
                        data4$q_410_asthma_prev_v5=="yes"|
                        data4$q_410_ad_cur=="yes"|
                        data4$q_410_ad_prev_v5=="yes",
                      "yes",
                      "no") %>% factor()

#### Cardiologic grouping var
data4$cc_cardio <- ifelse(data4$q_410_cardio_cur=="yes"|
                            data4$q_410_cardio_prev_v5=="yes"
                          ,
                          "yes",
                          "no") %>% factor()

### Add tryptase coategory based on RUEFF (cut off 8)
data4 %<>%
  mutate(tryp_cat=ifelse(q_212_tryptase_value_v5<8,
                         "low",
                         "high"))
###  var sc. route of administration ####
data4$route_sc <- rep(FALSE,length(data4$b_submitdate))
data4$route_sc[data4$d_330_drug_group %in%
                 c("biologics",
                   #"xray_cm",
                   "local_anaesthetics",
                   #"chemotherapeutics",
                   "immunisation"#,
                   #"volume replacement",
                   #"muscle relaxant"
                 )] <- TRUE
data4$route_sc[data4$d_elicitor_gr5 == "insects"] <- T


##### FINAL DATABASE RDB ######

rdb <- data4[data4$reaction_type_brown=="anaphylaxis",]
countries <-rdb %>%
  filter(d_elicitor_gr5=="insects") %>%
  select(d_centres_country) %>%
  group_by(d_centres_country) %>%
  summarize(n=n())%>% arrange(desc(n))

#### Grouping variable insects vs no insects #####
grouping <- ifelse(rdb$d_elicitor_gr5=="insects",
                   "insects",
                   ifelse(is.na(rdb$d_elicitor_gr5)|rdb$d_elicitor_gr5=="unkown",
                          NA,
                          "other")
                   ) %>%
  factor
rdb$grouping <- grouping
rdbp <-rdb
rdbp$grouping<- grouping



##### DIAGRAM ##############
source("R/make_flow.R")
make_flow()


###Statistic###########

### Binomial trigger either insects or other
testInsectsbinomial <- makeTests(groups = "grouping",rdb=rdb) %>% arrange(pval)

testInsectsbinomialTab <- testInsectsbinomial %>% #filter(pval<1e-30) %>%
  select(variableName,counts_1,counts_2,fraq_1,fraq_2,
         pval,
         section) %>%
         {split(.,.$section)}

age_sex_matched <- rdb %>%
  match_patients(grouping_var = "grouping",
                 matching_vars = c("b_sex","d_age"),
                 df=T)

tests_matched <- makeTests(groups ="grouping",
                           rdb =age_sex_matched)

tests_matchedTab <- tests_matched %>% filter(pval<1e-30) %>%
  select(variableName,counts_1,counts_2,fraq_1,fraq_2,
         pval,
         section) %>%
         {split(.,.$section)}

testInsectsbinomialTab2 <- testInsectsbinomialTab
testInsectsbinomialTab <- tests_matchedTab

rdb$groupYJ <- ifelse(rdb$d_insect_gr4 == "yellow jacket",
       "y-j",
       "other")

testsYJ <- makeTests("groupYJ",rdb = rdb) %>%
  arrange(pval) %>%
  select(variableName,counts_1,counts_2,fraq_1,fraq_2,
         pval,
         section) %>%
         {split(.,.$section)}



###### FOREST PLOT ####

source("R/make_all_forest_plots.R")
#make_all_forest_plots() # important chached figures!!!

##### convert to function ######



# cramerFun(data =rdbp,
#           grouping = "grouping",
#           vars = c("q_114",
#                    "q_112_incontinence",
#                    "q_114_dizziness",
#                    "q_112_abdominal_pain",
#                    "q_114_loss_of_consciousness",
#                    "q_114_hypotension_collapse_v5",
#                    "q_114_reductions_of_alertness",
#                    "q_111_angioedema",
#                    "q_113_wheezing_expiratory_distress_v5",
#                    "q_111_pruritus"))


temp <- cramerFun(data =rdbp,
               vars =c( "q_114",
               "q_112_incontinence",
               "q_114_dizziness",
               "q_112_abdominal_pain",
               "q_114_loss_of_consciousness",
               "q_114_hypotension_collapse_v5",
               "q_114_reductions_of_alertness",
               "q_111_angioedema",
               "q_113_wheezing_expiratory_distress_v5"),
               grouping = "grouping")


#### Demographics####

demoTab <- cbind(n = rdb$b_sex[rdb$d_elicitor_gr5=="insects"] %>% summary(),
                 Age = rdb$d_age[rdb$d_elicitor_gr5=="insects"] %>%
                   split(rdb$b_sex[rdb$d_elicitor_gr5=="insects"]) %>%
                   lapply(.,function(x){mean(x) %>% signif(3)}),
                 Cardiologic = f4(rdb$q_410_cardio_cur),
                 DM = f4(rdb$q_410_diab_cur_v6),
                 `Food allergy` = f4(rdb$q_410_foodallergy_cur_v6),
                 Mastocytosis = f4(rdb$q_410_masto_cur),
                 Malignancy = f4(rdb$q_410_malig_cur),
                 `Atopic dermatitis` = f4(rdb$q_410_ad_cur),
                 `Rhinitis` = f4(rdb$q_410_rhinitis_cur),
                 Asthma = f4(rdb$q_410_asthma_cur),
                 `tryptase [median]` = rdb$q_212_tryptase_value_v5 [rdb$d_elicitor_gr5=="insects"]%>%
                   split(rdb$b_sex[rdb$d_elicitor_gr5=="insects"]) %>%
                   lapply(median,na.rm=T),
                 `Insects as elicitors` = rdb$grouping %>%
                   split(rdb$b_sex) %>%
                   lapply(function(x){(length(which(x=="insects"))/length(x)*100) %>% roP})
)

demo_age_tab <- cbind(n = rdb$d_age_gr2[rdb$d_elicitor_gr5=="insects"] %>% summary(),
                 Age = rdb$d_age[rdb$d_elicitor_gr5=="insects"] %>%
                   split(rdb$d_age_gr2[rdb$d_elicitor_gr5=="insects"]) %>%
                   lapply(.,function(x){mean(x) %>% signif(3)}),
                 Cardiologic = f5(rdb$q_410_cardio_cur),
                 DM = f5(rdb$q_410_diab_cur_v6),
                 `Food allergy` = f5(rdb$q_410_foodallergy_cur_v6),
                 Mastocytosis = f5(rdb$q_410_masto_cur),
                 Malignancy = f5(rdb$q_410_malig_cur),
                 `Atopic dermatitis` = f5(rdb$q_410_ad_cur),
                 `Rhinitis` = f5(rdb$q_410_rhinitis_cur),
                 Asthma = f5(rdb$q_410_asthma_cur),
                 `tryptase [median]` = rdb$q_212_tryptase_value_v5 [rdb$d_elicitor_gr5=="insects"]%>%
                   split(rdb$d_age_gr2[rdb$d_elicitor_gr5=="insects"]) %>%
                   lapply(median,na.rm=T),
                 `Insects as elicitors` = rdb$grouping %>%
                   split(rdb$d_age_gr2) %>%
                   lapply(function(x){(length(which(x=="insects"))/length(x)*100) %>% roP})
)


#### plot MOR#####
#this is a grid extra figure and not required
#plot_MOR <- plot_mor()

#### Countries ####

t1 <- table(rdb$d_centres_country,rdb$q_340_insects)
t2 <- t1 %>% prop.table(1) %>% {round(.*100,1)}


#### Previous ANA ####
rdb$grouping <- grouping
# only visualize the p vaues - not necessary
#rdb$q_160_ever_react %>% table(rdb$grouping) %>% {.[1:2,]} %>% assocstats()

#rdb$q_1622_ever_mild_v4 %>% table(grouping) %>% {.[1:2,]}%>% assocstats()

#rdb$q_1621_ever_severe_v4 %>% table(grouping) %>% {.[1:2,]} %>% assocstats()

#### Previous Reacktions ########
# Check for previous
repeatedPat<- rdb %>%
  #filter(reaction_type_brown=="anaphylaxis",
  #       d_elicitor_gr5 =="insects") %>%
  group_by(b_patient_code) %>% summarise(n =n()) %>%
  arrange(desc(n)) %>% filter(n>1) %>% select(b_patient_code) %>%
  pull() %>% as.character()


tempDF<-filter(rdb,b_patient_code %in% repeatedPat) %>% #rowwise() %>%
  mutate(Rdate =  ifelse(substr(as.character(b_reactiondate),1,2)=="00",
           paste0("15",substr(as.character(b_reactiondate),3,10)),
           as.character(b_reactiondate))
         ) %>% mutate(Rdate = ifelse(substr(Rdate,4,5)=="00",
                                paste0(substr(Rdate,1,3),
                                      "06",
                                      substr(Rdate,6,10)),
                                Rdate)%>% as.Date(format = "%d.%m.%Y")) %>%
  # select(b_reactiondate,Rdate)
  select(#variableSelectionTab(rdb) %>%
          # filter(section=="symptoms") %>%
           #select(variableName) %>% pull() %>% as.character(),
         b_patient_code,d_severity_rm,Rdate,q_340_insects,d_elicitor_gr5
  ) %>% arrange(b_patient_code,Rdate) %>%
  data.frame()

iCode <- tempDF %>% split(tempDF$b_patient_code,drop = T) %>% lapply(function(x){
  any(x$d_elicitor_gr5=="insects")
}) %>% unlist %>% {names(.)[which(.==T)]}

tempDF %>% filter(b_patient_code%in%iCode)
repeated_reaction_in_Reg <- split(tempDF,tempDF$b_patient_code,drop = T) %>% lapply(function(x){
    ifelse(x$d_severity_rm[1]<x$d_severity_rm[1],
           "smaler",ifelse(x$d_severity_rm[1]==x$d_severity_rm[2],
                           "equal",
                           "greater"))
  }) %>% unlist() %>% factor %>% summary()


iCodeDF <- tempDF %>% filter(b_patient_code%in%iCode) %>%
  group_by(b_patient_code) %>% summarize(first = q_340_insects[1],
                                         first_sev = d_severity_rm[1],
                                         second = q_340_insects[2],
                                         second_sev = d_severity_rm[2]) %>%
  mutate(difference = ifelse(first_sev < second_sev, "greater",
                             ifelse(first_sev == second_sev, "same","lower")))
 greater_same <- iCodeDF %>%
   group_by(difference) %>% summarize(n = n())

iCodeDF %>% filter(difference=="greater")

### PLOT Symptom differences ####

# plotSympt <- ggplot(testInsectsbinomialTab$symptoms %>%
#                       tidyr::gather(key = "group",value = "Fraction", 4:5) %>%
#                       arrange(pval),
#                     aes(reorder(variableName,pval),Fraction, fill = group))+
#   geom_bar(stat = "identity", position = "dodge",na.rm = T)+
#   theme_classic()+
#   theme(axis.text.x = element_text(angle = 50,hjust =1))
#
# png(width = 400*2,height = 430*2,res = 300, filename = "symptomsG.png",pointsize = 7)
# tempplot <- testInsectsbinomialTab$symptoms %>%
#   tidyr::gather(key = "group",value = "Fraction", 4:5) %>%
#   arrange(pval) %>%
#   {.[c(3,4,11:12,7:8,9:10,13,14),]} %>%
#   mutate(variableName = car::recode(
#     variableName,recodes = "'q_114_dizziness'='dizziness';
#     'q_114_hypotension_collapse_v5'='collapse';
#     'q_114_loss_of_consciousness'='loss of consciousness';
#     'q_112_abdominal_pain'='abdominal pain';
#     'q_113_wheezing_expiratory_distress_v5'='expiratory distress'"),
#     group = car::recode(group,"'fraq_1'='IVA';
#                         'fraq_2'='non-IVA'")) %>%
#                         {mutate(.,variableName = factor(.$variableName,
#                                                         levels=unique(.$variableName)))} %>%
#
# ggplot(
#                     aes(variableName,Fraction, fill = group))+
#   geom_bar(stat = "identity", position = "dodge",na.rm = T)+
#   theme_classic()+
#   theme(axis.text.x = element_text(angle = 30,hjust =1),
#         legend.position = c(1,1),
#         legend.justification = c(1,1))+
#   labs(x = "Symptoms",y = "Proportion of cases",fill="")+
#   scale_fill_manual(values = rev(c("#E5F5E0", "#74C476","#005A32")))
# dev.off()
#

##### RADAR CHART SPIDER ####

temp1 <- lapply(list(quo(q_111),
            quo(q_112),
            quo(q_113),
            quo(q_114),
            quo(atopy)),
       funtempspider) %>%
  do.call(what = rbind) %>%
  as.data.frame()

temp <- temp1 %>%
  filter(grouping == "insects",
         q_111 == "yes")

temp2 <- temp1 %>%
  filter(grouping == "other",
         q_111 == "yes")

temp <- temp %>% t()

colnames(temp) <- temp[5,]

tempdf <- temp %>% data.frame(stringsAsFactors = F) %>%
{rbind(rep(1,length(.)),rep(0,length(.)),.[2,],temp2 %>% t() %>% {.[2,]})}
tempdf<- tempdf %>%
  mutate_if(is.character,as.numeric)
tempdf <- age_sex_matched %>%
  group_by(grouping) %>%
  summarize(tryptase = mean(q_212_tryptase_value_v5,na.rm = T)) %>%
    {rbind(11,0,.[,2])} %>%
    cbind(tempdf,.)
 radarchart(tempdf,
             #axistype=1,
           #seg=5, %>%
           #plty=1,
           #title="(axis=1, 5 segments, with specified vlabels)",
           #vlcex=0.5
           vlabels = c("skin",
                      "gastrologic",
                      "respiratory",
                      "cardiac",
                      "atopic",
                      "tryptase"))


# Split the plot into kids and grown ups AGE division - 12
# Get the proper table

symptTabKids <- makeTests(rdb=rdb[rdb$d_age<13,],
          groups = "grouping") %>%
          {split(.,.$section)} %>%
  .$symptoms %>% mutate(subset = "children")

symptTabAdults <- makeTests(rdb=rdb[rdb$d_age>13,],
                                  groups = "grouping") %>%
                                  {split(.,.$section)} %>%
                          .$symptoms %>% mutate(subset = "adults")
ex <- full_join(symptTabKids,symptTabAdults)

# ggplot(ex %>% filter(pval<1e-10) %>%
#          tidyr::gather(key = "group",value = "Fraction", c("fraq_1","fraq_2")) %>%
#          arrange(pval),
#        aes(reorder(variableName,pval),Fraction, fill = group))+
#   geom_bar(stat = "identity", position = "dodge",na.rm = T)+
#   facet_grid(.~subset)+
#   theme_classic()+
#   theme(axis.text.x = element_text(angle = 50,hjust =1))


hypotensionPlot <-
  ggplot(ex %>%
           filter(variableName =="q_114_hypotension_collapse_v5") %>%
           tidyr::gather(key = "group",
                         value = "Fraction",
                         c("fraq_1",
                           "fraq_2")) %>%
           mutate(group = car::recode(group,
                  "'fraq_1'='IVA';
                   'fraq_2'='non-IVA'")) %>%
         arrange(pval),
       aes(#reorder(variableName,pval),
          subset,
         Fraction,
           fill = group))+
  geom_bar(stat = "identity", position = "dodge",na.rm = T)+
  #facet_grid(.~subset)+
  scale_fill_manual(values = rev(c("#E5F5E0", "#74C476", "#005A32")))+
  theme_classic()+

  #theme(axis.text.x = element_blank(),
  #      axis.ticks.x = element_blank())+
  labs(y = "proportion [%]", x = "hypotension")+
   theme(axis.text.x = element_text(angle = 50,hjust =1),
         legend.position = "top")



#### Mastocytose patienten  as a subroup #####



#### MEasures of Association ############
vcd::assocstats(table(rdb$q_114_hypotension_collapse_v5,
                      grouping)[1:2,])
vcd::assocstats(table(rdb$q_114,
                      grouping)[1:2,])

crammerElicitorBinvsSympt <-
  lapply(variableSelectionTab(rdb) %>% filter(section=="symptoms",level==3) %>%
           pull(variableName) %>% as.character(),function(x){
             assocstats(rdb[,x] %>%
                          table(grouping) %>% .[1:2,1:2])[[4]]
           }) %>%
  unlist %>%
  {data.frame(CV=.,
              var=variableSelectionTab(rdb) %>%
                filter(section=="symptoms",level==3) %>%
                pull(variableName))} %>%
  arrange(desc(CV))


### REmove the unknown group ####

### Check associated variables
checkVarTab(data = rdb,
              var = as.character(crammerElicitorBinvsSympt[3,2]),
              groupby = "grouping")

#### The most associated differences
AssociatedVars <- testInsectsbinomial %>%
  arrange(desc(Cramer)) %>%
  select(-c(2:5,11:16)) %>%
  filter(Cramer>0.25)

#### Crammer Plot ####
##### Plot the most Crammers V in ALL patients, Choldren, Adults and
# testInsectsbinomialTab$symptoms
# symptTabKids
# symptTabAdults
# rdb$grouping <- grouping
# all = makeTests(rdb=rdb,
#                 grouping = rdb$grouping) %>%
#                 {split(.,.$section)} %>%
#   .$symptoms %>% mutate(subset = "all")
# ex <- full_join(all,
#   symptTabKids) %>%
#   full_join(symptTabAdults) %>% arrange(desc(Cramer)) %>%# select(variableName) %>% pull
#   filter(variableName%in%c("q_114",
#                            "q_114_dizziness",
#                            "q_114_reductions_of_alertness",
#                            "q_114_loss_of_consciousness",
#                            "q_114_hypotension_collapse_v5",
#                            "q_112_abdominal_pain",
#                            "q_112_vomiting"))
# ggplot(ex %>%
#          tidyr::gather(key = "group",value = "Fraction", c("fraq_1","fraq_2")) %>%
#          arrange(pval),
#        aes(reorder(variableName,pval),Fraction, fill = group))+
#   geom_bar(stat = "identity", position = "dodge",na.rm = T)+
#   facet_grid(.~subset)+
#   theme_classic()+
#   theme(axis.text.x = element_text(angle = 50,hjust =1))





cramerPlot <-
  ex %>%
  filter(!is.na(Cramer),
         variableName%in%ex$variableName[which(ex$Cramer>0.1)]
         ) %>%
  ggplot( aes(variableName,Cramer,fill=subset))+
  geom_bar(stat = "identity", position = "dodge")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 50,hjust =1))

### Tryptase plot#####


# plot.proportions(data = rdb,
#                  varx = "d_centres_country",
#                  vary = "q_340_insects",
#                  minN = 20)


##### reaction severity after SIT ######
### rationale
# Do patients with previous SIT had less severe reactions than these without previous SIT?
### data and variables
dt <- rdb %>%
  filter(d_elicitor_gr5=="insects",
         q_610_sit_prior_v5 %in% c("no","yes")) %>%
  group_by(q_610_sit_prior_v5,d_severity_rm) %>%
  summarize(n = n()) %>%
  mutate(d_severity_rm = factor(d_severity_rm)) %>%
  ggbarplot(x = "q_610_sit_prior_v5",
            y = "n",
            fill = "d_severity_rm",
            position = position_fill())

dt3 <- rdb %>%
  filter(d_elicitor_gr5=="insects",
         q_610_sit_prior_v5 %in% c("no","yes")) %>%
  mutate(severity_brown = factor(severity_brown)) %>%
  group_by(q_610_sit_prior_v5,severity_brown) %>%
  summarize(n = n()
            )


fisher.test(matrix(dt3$n,ncol=2,byrow=T))

### Plot
### stat test
fisher.test(matrix(dt$data$n,ncol=4,byrow=T)[,2:4])
### coonclusion
# There were no statisitcally significant differences in these groups
### discussion
# Probably there is a potential difference in the tests when we look into
# a previous reaction severity and next reaction severity

require(purrr)
##### F1 repeate reaction severity after SIT ######
F1 <- list()
### rationale
F1$rationale<- "Do patients with previous SIT had less severe reactions than these without previous SIT? we compare these reactions in the same patients."
### libraries
F1$libs <- "purrr"
### data and variables
# patients with at leat two documented reactions in our database
# split-apply-combine using dplyr,tidyr,purrr
F1$data[["pat_ids"]]<-
  data4 %>%
  filter(d_elicitor_gr5=="insects") %>%
  select(b_patient_code) %>%
  group_by(b_patient_code) %>%
  summarise(reps =n()) %>%
  filter(reps>1) %>%
  select(b_patient_code) %>% pull() %>% as.character()
F1$data[["dt"]]  <-
  data4 %>%
  filter(d_elicitor_gr5=="insects",
         #q_610_sit_prior_v5 %in% c("no","yes"),
         b_patient_code %in% F1$data$pat_ids) %>%
  select(b_patient_code,d_severity_rm,b_reactiondate) %>%
  mutate(Rdate =  ifelse(substr(as.character(b_reactiondate),1,2)=="00",
                         paste0("15",substr(as.character(b_reactiondate),3,10)),
                         as.character(b_reactiondate))
  ) %>%
  mutate(Rdate = ifelse(substr(Rdate,4,5)=="00",
                            paste0(substr(Rdate,1,3),
                                   "06",
                                   substr(Rdate,6,10)),
                            Rdate)%>%
           as.Date(format = "%d.%m.%Y")
         ) %>%
  group_by(b_patient_code) %>%
  #do(ifelse(.$Rdate==max(.$Rdate),"last","prev"))
  tidyr::nest() %>%
  mutate(reaction = purrr::map(data,function(df){
    ifelse(df$Rdate==max(df$Rdate),
           "last",
           "prev")
  })) %>%
  tidyr::unnest()

F1$data[["dt2"]] <- F1$data$dt %>%
  group_by(b_patient_code) %>%
  tidyr::nest() %>%
  mutate(compare_r = map(data,function(df){
  ifelse(df$d_severity_rm[df$reaction=="last"]<df$d_severity_rm[df$reaction=="prev"],
         "reduction",
         "no reduction")
})) %>%
  tidyr::unnest(compare_r) %>%
  group_by(compare_r) %>%
  summarize(n=n())


### Plot
F1$vis <- F1$data$dt%>% ggplot(aes(reaction,y=factor(d_severity_rm)))+
  geom_point()+
  geom_line(aes(group=b_patient_code),
            position = position_jitter(height = 0.05))
### stat test
F1$tests <- NULL#fisher.test(matrix(dt$n,ncol=4,byrow=T)[,2:4])
### coonclusion
F1$conclusion <-  "More patients showed no reduction in severity after in repeated reaction after SIT"
### discussion
F1$discussion <- "Low patient numbers could have prevented an adequate analysis. This is very plausible as many of these reactions which happen in the meantime could be really mild and therefore not cosidered anaphylaxis and therefore not reported in the registry"

source("R/findingOrientedResearch.R")

### Adrenaline managment in insect cases ###
# F2 Adrenaline use Corellation with severity? #####
F2<- list(rationale = "We saw that patients who were treated for anaphylaxis after insect sting less often recived adrenaline than patients who were treated for anaphylaxis due to other triggers. We want to evaluate the cause of this observation",
             libs = "")

F2$plot[["RM"]] <- rdb %>%
  filter(!is.na(d_severity_rm),
         !is.na(d_522_adren_agg),
         !is.na(grouping),
         d_522_adren_agg!="unknown") %>%
  group_by(d_severity_rm,grouping) %>%
  summarize(n =mean(d_522_adren_agg=="yes")) %>%
  ggplot(aes(d_severity_rm,n,fill = grouping))+
  geom_bar(stat="identity",position="dodge")


F2$test <- rdb %>%
  filter(!is.na(d_severity_rm),
         !is.na(d_522_adren_agg),
         !is.na(grouping),
         d_522_adren_agg!="unknown") %>%
  {glm(d_522_adren_agg~grouping*d_severity_rm,family = "binomial",data =.)} %>%
  summary()


F2$plot[["brown"]] <- rdb %>%
  filter(!is.na(severity_brown),
         !is.na(d_522_adren_agg),
         !is.na(grouping),
         d_522_adren_agg!="unknown") %>%
  group_by(severity_brown,grouping) %>%
  summarize(n =mean(d_522_adren_agg=="yes")) %>%
  ggplot(aes(severity_brown,n,fill = grouping))+
  geom_bar(stat="identity",position="dodge")


F2$plot[["insectvsother"]] <- rdb %>%
  filter(!is.na(d_severity_rm),
         !is.na(d_522_adren_agg),
         !is.na(grouping),
         d_522_adren_agg!="unknown") %>%
   # select(#d_elicitor_gr5,
  #    d_522_adren_agg,
   #   d_severity_rm,
    #  grouping)
  #pull() %>%
  ggplot(aes(d_severity_rm,fill=d_522_adren_agg))+
  geom_bar(position = "fill")+
  facet_grid(.~grouping)
  #group_by(grouping,severity) %>%
  #summarize()
  # c(data.tab = F2$data.d %>%
  #     group_by(grouping,d_severity_rm,d_522_adren_agg)
  #   ) %>%
  # c(vis = list(plot="1")
  #   ) %>%
  # c(conclusion = "This is dumb test") %>%
  # c(discussion = "test could be ok?"
  #   )

F3 <- list()
F3$rationale <- "The treatment groups need to be adjusted for the availibility of autoinjectors so that the actual use can be calcualted"
F3$plot
#rdb$
#ggplot(aes())

#### Propensity Score MAtching for Adrenaline use #####
#' The functition to match samples based on the minimal set of predictors
#' provide data, grouping variables and predictor variables as well as others..
prop_fun <- function(data, grouping_var, predictor_vars, other_vars=NULL,
                     method = "optimal",
                     samplesize = 200) {
  if(is.logical(grouping_var)){
    temp <- data %>%
      select(!!!predictor_vars,!!!other_vars) %>%
      mutate(grouping = grouping_var) %>%
      filter(complete.cases(.))
  } else {
  temp <- data %>%
    select(!!!grouping_var,!!!predictor_vars,!!!other_vars) %>%
    filter(complete.cases(.)) %>%
    mutate(grouping = ifelse(get(grouping_var) == "insects", T,
                           F))
  }
  match_stats <- matchit(formula = paste0(grouping_var,
                                          "~",
                                          paste0(predictor_vars,
                                                 collapse = "+"))%>%
                           as.formula(),
                         data=data.frame(temp[sample(1:length(temp[,1]),
                                                     samplesize),]),
                         method = method,
                         ratio = 1,
                         na.rm=T)
  return(list(pre_data = temp,
              post_data = match.data(match_stats),
              stats = match_stats))#%>%
}#str()

temp <- prop_fun(data = rdb,
         grouping_var = "grouping",
         predictor_vars = c("b_sex","d_age"),
         other_vars = c("q_116_VAS_v7"),
         samplesize = 1200
         )
temp$post_data %>%
ggplot(aes(grouping,y=q_116_VAS_v7))+
  geom_violin()

temp$post_data %>%
{kruskal.test(.$q_116_VAS_v7,.$grouping)}
rdb %>%
{kruskal.test(.$q_116_VAS_v7,.$grouping)}

ggplot(temp$post_data, aes(grouping, q_116_VAS_v7))+
  geom_boxplot()


####3 Here is the adrenalin data analysis
#temp <-
rdb %>%
  #filter(d_AAI_prescribed!="no") %>%
  group_by(grouping, d_AAI_prescribed,d_522_adren_agg) %>%
  summarize(n = n()) %>%
  ggplot(aes(y=n, x = d_AAI_prescribed,fill=d_522_adren_agg))+
  geom_bar(stat = "identity",position ="fill")+
  facet_grid(.~grouping)



ANAscore_matched <-
  match_patients(rdb %>%
                   filter(d_522_adren_agg %in% c("yes", "no")),
                 "grouping",
                 c("ANAscore","d_age"),
                 grouping_var_level = "insects",
                 T)
testANAscoreMatched <- makeTests(groups = "grouping",
                                 rdb=ANAscore_matched) %>%
  #arrange(pval) %>%
  #filter(pval<1e-30) %>%
  select(variableName,counts_1,counts_2,fraq_1,fraq_2,
         pval,
         section) %>%
         {split(.,.$section)}


countYesManagment <- function(var){
  tbl <- ANAscore_matched %>%
    group_by(grouping,!! sym(var))  %>%
    summarize(n = n()) %>%
    filter(!!sym(var) == "yes")
  data.frame(variable = names(tbl)[2],
             IVA = tbl[1,3],
             nonIVA = tbl[2,3])
}

plotManagement <- testANAscoreMatched$management$variableName %>%
  map(countYesManagment) %>%
  do.call(what = rbind)  %>%
  tidyr::gather(key = "grouping", value = "positive",2:3) %>%
  filter(variable %in% c(
    "d_522_adren_agg",
    "q_522_antih_iv",
    "q_522_cortico_iv",
    "q_522_o2",
    "q_522_volume",
    "q_561_hospital_admission_v6",
    "q_562_intensive_care_v6",
    "q_522_beta2_inhal"
  )) %>%
  mutate(
    variable = car::recode(
      variable,
      recodes = "'q_522_cortico_iv'='corticoids iv.';
            'q_522_antih_iv'='antihistamines iv.';
            'd_522_adren_agg'='adrenaline iv./im.';
            'q_522_adren_iv'='adrenaline iv.';
            'q_522_o2'='100% oxygen';
            'q_522_beta2_inhal'='beta-2 mimetics inh.';
      'q_522_volume'='volume iv.';
      'q_561_hospital_admission_v6'='hospital admission';
      'q_562_intensive_care_v6'='intensive care'"),
    grouping = ifelse(grouping =="n", "IVA","non-IVA")
    ) %>%
  ggplot(aes(reorder(variable, -positive),positive/1976,fill=grouping))+
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    legend.position = c(1, 1),
    legend.justification = c(1, 1)
  ) +
  labs(x = "Therapy", y = "proportion", fill = "") +
  scale_fill_manual(values = rev(c("#E5F5E0", "#74C476", "#005A32"))) +
  theme(axis.title.x = element_blank())+ scale_y_continuous(labels = scales::percent_format(accuracy = 1))


backup <- testANAscoreMatched$management[, -7] %>%
  tidyr::gather(value = "Proportion",
                key = "Group",
                4:5) %>%
  filter(
    #counts_1>100|pval<0.05
    variableName %in% c(
      "q_521_autoinj_v5",
      "q_521_beta2_v5",
      "q_521_antih_v5",
      "q_521_cortic_v5",
      "d_522_adren_agg",
      "q_522_adren_im",
      "q_522_adren_iv")
    #   "d_520_adren1",
    #   "q_562_intensive_care_v6",
    #   "q_522_adren_im",
    #   "q_522_adren_iv",
    #   "q_552_beta2_inhal",
    #   "q_522_beta2_inhal",
    #   "q_561_hospital_admission_v6",
    #   "q_550_2nd_v5",
    #   "q_552_volume_v5",
    #   "q_552_cortico_iv_v5",
    #   "q_552_antih_oral_v5"
    # )
  ) %>%
  mutate(
    variableName = car::recode(
      variableName,
      recodes = "'q_552_cortico_iv'='corticoids iv.';
      'q_522_antih_iv'='antihistamines iv.';
      'q_522_adren_im'='adrenaline im.';
      'q_522_adren_iv'='adrenaline iv.';
      'q_522_o2'='100% oxygen';
      'q_522_beta2_inhal'='beta-2 mimetics inh.'"
          ),
          Group = car::recode(Group, "'fraq_1'='IVA';
                              'fraq_2'='non-IVA'")
          ) %>%
        ggplot(aes(reorder(variableName, -Proportion), Proportion, fill =
                     Group)) +
        geom_bar(stat = "identity", position = "dodge") +
        theme_classic() +
        theme(
          axis.text.x = element_text(angle = 20, hjust = 1),
          legend.position = c(1, 1),
          legend.justification = c(1, 1)
        ) +
        labs(x = "Therapy", y = "proportion", fill = "") +
        scale_fill_manual(values = rev(c("#E5F5E0", "#74C476", "#005A32"))) +
  theme(axis.title.x = element_blank())


upper_panel <-ggpubr::ggarrange(
      # ANAscore_matched %>%
      #   mutate(grouping = recode(grouping,
      #                            'insects'='IVA',
      #                            'other'='non-IVA'),
      #          d_522_adren_agg = recode(d_522_adren_agg,
      #                                  'no'='no adrenaline',
      #                                  'yes'='adrenaline given') ) %>%
      #                    group_by(grouping,
      #                             #d_AAI_prescribed,
      #                             d_522_adren_agg) %>%
      #                    summarize(n = n()) %>%
      #                    ggplot(aes(y = n, fill = grouping, x =
      #                                 d_522_adren_agg)) +
      #                    geom_bar(stat = "identity",position = "dodge")+
      #                    theme_classic()+
      #                    theme(legend.position = "right",
      #                          axis.text.x = element_text(angle = 30,hjust = 1))+
      #                    labs(y ="number of cases",
      #                         x = "Treatment with adrenaline",
      #                         fill = "Elicitor")+
      #                    scale_fill_manual(values = rev(c("#E5F5E0", "#74C476","#005A32")))+
      #   theme(axis.title.x = element_blank()),

      plotManagement,

ANAscore_matched %>%
  filter(!is.na(q_540_why_autoinj_v5)) %>%
  group_by(grouping,
           q_540_why_autoinj_v5) %>%
  summarize(n = n()) %>%
  group_by(grouping) %>%
  nest() %>%
  mutate(summed =
  map(data,.f=function(tab){
    rep(sum(tab$n),4)
  })
  ) %>% unnest() %>%
  mutate(proportion = n/summed,
         grouping = recode(grouping,
                           'insects'='IVA',
                            'other'='non-IVA'))  %>%
  ggplot(aes(y = proportion, fill = grouping, x =
               q_540_why_autoinj_v5,
             label = n)) +
  geom_bar(stat = "identity",position="dodge")+
  geom_text(aes(label=n),
           position=position_dodge(width = 1),
           size=4,
           vjust=0)+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle=20, hjust = 1))+
  labs(y ="proportion",
       x = "Reason for\nnot administering adrenaline",
       fill = "Elicitor")+
  scale_fill_manual(values = rev(c("#E5F5E0", "#74C476","#005A32")))+
  theme(axis.title.x = element_blank())+
  ylim(0,.99),
ncol = 2,
common.legend = F,
widths = c(1,0.7),
align = "h",
labels = c("A","B"))

lower_panel <- ggpubr::ggarrange(
  ANAscore_matched  %>%
    ggplot(aes(x = d_age, fill = grouping)) +
    geom_density(alpha = 0.5)+
    labs(x = "Age [years]",fill = "elicitor")+
    theme_classic()+
    theme(legend.position = "none")+
    scale_fill_manual(values = rev(c("#c4c4c4", "#000000")))+
    labs(y ="density",x = "Age [years]",fill = "")+
    geom_segment(aes(x = 22, xend =22, y = 0, yend = 0.02),
                 linetype = 2),

  ANAscore_matched %>%
    group_by(b_sex,grouping) %>%
    ggplot(aes(fill=grouping,x=b_sex))+
    geom_bar(position = "fill")+
    theme_classic()+
    theme(legend.position = "none")+
    scale_fill_manual(values = rev(c("#c4c4c4", "#454545")))+
    labs(y ="proportion",x = "Sex",fill = "elicitor"),

  ANAscore_matched  %>%
    #select(grouping,d_severity_rm) %>%
    group_by(grouping,d_severity_rm) %>%
    summarize(n = n()) %>%
    ggplot(aes(x = as.factor(d_severity_rm),y=n, fill = grouping)) +
    geom_bar(stat="identity",position="dodge")+
    theme_classic()+
    labs(x = "Severity grade [R&M]",fill = "R&M",y="cases [n]")+
    theme(legend.position = "none")+
    scale_fill_manual(values = rev(c("#c4c4c4", "#454545"))),

common.legend = T,
legend = "right",
align = "h",
ncol = 3,
labels = c("B","C","D"))


# png(res=300,filename = "adrenaline.png",pointsize = 6,width = 3500,height=1900)
# grid.draw(fig_adrenuse1)
# dev.off()

test_adrenuse1<- ANAscore_matched %>%
  group_by(grouping,d_522_adren_agg) %>%
  summarize(n = n()) %>% {matrix(.$n,ncol = 2)} %>%
  chisq.test()

###### Test if me see less skin symptoms in mastocytosis patients #####
df <- data4 %>%
  filter(d_elicitor_gr5=="insects") %>%
  group_by(q_111,q_410_masto_cur) %>%
  summarize(n()) %>%
  spread(key = q_410_masto_cur,value = `n()`) %>%
  as.data.frame()

o <- prop.table(as.matrix(df[1:2,2:3]),2)*100
o <- as.matrix(df[1:2,2:3])
chisq.test(o)
rownames(o) <- c("no skin symptoms", "skin symptoms present")
colnames(o) <- c("no masto","masto present")




### route of administration severity ####
rdb %>%
  filter(route_sc==T | d_elicitor_gr5=="insects") %>%
  #{table(.$d_severity_rm,.$d_elicitor_gr5)}
  group_by(d_severity_rm,d_elicitor_gr5) %>%
  summarise(n = n()) %>%
  ggplot(aes(fill=d_elicitor_gr5,x = d_severity_rm,y=n))+
  geom_bar(stat = "identity",position = "fill")


## matchit according to route of administration )
# options("optmatch_max_problem_size" = Inf)
temp <- prop_fun(data = rdb,
                 grouping_var = "grouping",
                 predictor_vars = c("b_sex","d_age"),
                 other_vars = c("d_severity_rm","route_sc")
)


temp$post_data %>% group_by(grouping,route_sc) %>%
  summarize(n())


#Severity analysis####

#### Symptoms severity #####
symptomsdf <- data.frame( symptoms=c("q_111_angioedema",                        "q_111_erythema_flush_v5",
                                     "q_111_pruritus",                          "q_111_urticaria",
                                     "q_111_conjunctivitis_v5",                 "q_112_abdominal_pain",
                                     "q_112_abdominal_distention_v5",           "q_112_diarrhoea",
                                     "q_112_dysphagia_v5",                      "q_112_vomiting",
                                     "q_112_incontinence",                      "q_112_nausea",
                                     "q_113_respiratory_arrest",                "q_113_dyspnea",
                                     "q_113_chest_tightness_v5",                "q_113_cough_v5",
                                     "q_113_change_in_voice_v5",                "q_113_throat_tightness_v5",
                                     "q_113_wheezing_expiratory_distress_v5",   "q_113_rhinitis_v5",
                                     "q_113_stridor_inspiratory"   ,            "q_114_loss_of_consciousness",
                                     "q_114_hypotension_collapse_v5",           "q_114_chest_pain_angina_v5",
                                     "q_114_palpitations_cardiac_arrythmia_v5", "q_114_cardiac_arrest",
                                     "q_114_dizziness",                         "q_114_tachycardia",
                                     "q_114_reductions_of_alertness",           "q_115_dysarthria_v6",
                                     "q_115_dysphonia_v6",                      "q_115_hot_sweat_tremble_v6",
                                     "q_115_tingle_hands_feet_paresthesia_v6",  "q_115_sight_disorder_v6",
                                     "q_115_agony_v6",                          "q_115_cyanosis_pallor_v6",
                                     "q_140_fatal"                            ),
                          weights=c(3,2,1,2,2,3,1,4,3,4,12,3,8,4,1,20,4,2,8,8,1,20,18,20,10,6,18,12,10,18,2,8,3,15,20,6,50),
                          grades=c(1,1,1,1,1,1,1,2,1,2,3,1,2,2,1,4,2,1,2,2,1,4,3,4,2,2,3,3,2,3,1,2,1,3,4,2,5),
                          myVas = c(449,74,56,114,132,170,54,173,319,628,
                                    740,133,968,244,292,88,285,516,479,30,
                                    401,887,846,516,580,983,552,214,770,298,
                                    272,303,505,653,708,829,1000)
)

#data <- data5
### ANASCORE######
#MAke the anascore variable as the VAS points (0 - 1000) of the symptoms
# Use only the maximum severe symptom.
#
# data$ANAscore <- (1:length(data[,1])) %>%
#   map(function(x){
#     data[x,] %>%
#       select(as.character(symptomsdf[,1])) %>%
#       tidyr::gather(key="symptom",value = "val") %>%
#       cbind(point = symptomsdf$myVas) %>%
#       filter(val=="yes") %>%
#       summarize(score = max(point))
#   }) %>% unlist()
#
# data$ANAscore[data$q_130_biphasic_v4 =="yes"] <-
#   data$ANAscore[data$q_130_biphasic_v4 =="yes"]+150
# load("data.R")

# add_anascore_points <-function(x,points) {
#   data$ANAscore[which(x =="yes")] <<-
#     data$ANAscore[which(x =="yes")]+points
#
# }
#
#
# add_anascore_points(data$d_552_adren_agg_v5,110)
# add_anascore_points(data$d_522_adren_agg,100)
# add_anascore_points(data$d_560_adren2_v5,130)
# add_anascore_points(data$q_562_intensive_care_v6,200)
# add_anascore_points(data$q_522_volume,50)
# add_anascore_points(data$q_552_volume_v5,65)
# add_anascore_points(data$q_561_hospital_admission_v6,100)
# add_anascore_points(data$q_522_dopamine,75)
# add_anascore_points(data$q_552_dopamine_v5_v5,85)

# save(data, file="data.R")
load("data.R")

data %>%
  ggplot(aes(d_severity_rm,as.numeric(q_116_VAS_v7)))+
  geom_violin()

data %>%
  filter(!is.na(d_severity_rm)) %>%
  ggplot(aes(d_severity_rm,as.numeric(q_116_VAS_v7)))+
  geom_boxplot() +
  labs(x="Ring and Messmer", y="VAS")

# rdb %>%
#   ggplot(aes(severity_brown,as.numeric(q_116_VAS_v7)))+
#   geom_violin()
#
#
# rdb %>%
# ggplot(aes(x = ANAscore,y = as.numeric(q_116_VAS_v7))) +
#   geom_point()
#
#
# rdb %>%
#   filter(!is.na(d_severity_rm)) %>%
#   ggplot(aes(x = ANAscore,y = as.numeric(q_116_VAS_v7))) +
#   geom_jitter()+
#   facet_grid(.~d_severity_rm)
#
#
# rdb %>%
#   filter(!is.na(d_severity_rm)) %>%
#   ggplot(aes(x = ANAscore,y = q_116_VAS_v7)) +
#   geom_jitter()+
#   facet_grid(severity_brown~d_severity_rm)
#
# rdb %>%
#   filter(!is.na(d_severity_rm),
#          !is.na(q_116_VAS_v7)) %>%
#   ggplot(aes(x = ANAscore,y = factor(q_116_VAS_v7))) +
#   geom_jitter()+
#   facet_grid(severity_brown~d_severity_rm)
#
#
# rdb %>%
#   filter(!is.na(d_severity_rm),
#          !is.na(q_116_VAS_v7)) %>%
#   ggplot(aes(x = ANAscore,y = factor(q_116_VAS_v7), color= q_140_fatal)) +
#   geom_jitter()+
#   facet_grid(relevel(severity_brown,"severe")~d_severity_rm)
#
# rdb %>%
#   filter(!is.na(d_severity_rm),
#          !is.na(q_116_VAS_v7)) %>%
#   ggplot(aes(x = ANAscore,y = factor(q_116_VAS_v7), color= q_113_respiratory_arrest)) +
#   geom_jitter()+
#   facet_grid(relevel(severity_brown,"severe")~d_severity_rm)
#

# purrr::map(symptomsdf[,1] %>% as.character,
#            function(x){
#   rdb %>%
#   select(q_116_VAS_v7,!!!x) %>%
#     group_by(q_116_VAS_v7) %>%
#     summarize(x = sum(get(x)=="yes")/n()) %>%
#                {.$q_116_VAS_v7[which(.$x == max(.$x,na.rm=T))]}
#            })


# data_lm=as.data.frame(rdb)
# # Fit a logistic regression model
# fit_glm=glm(paste("q_116_VAS_v7~",paste0(symptomsdf[,1],collapse = "+")) %>% as.formula()
#             ,data=data_lm)
# fit_glm_f <- MASS::stepAIC(fit_glm)
# # generate summary
# summary(fit_glm)
# # Using varImp() function
# library(caret)
# varImp(fit_glm_f) %>%
# {data.frame(var = rownames(.)[order(.,decreasing = T)],
#             imp = .[order(.,decreasing = T),])}
#
#
# #Import the random forest library and fit a model
# library(randomForest)
#
# fit_glm_f$formula %>%
#
#
# # data_lm %>%
# #   select(fit_glm_f$data %>% names(),q_116_VAS_v7) %>%
# #   tidyr::drop_na() %>%
# #           {randomForest(fit_glm_f$formula, data=.)}
#
# library(rpart)
# fit_rpart <- rpart::rpart(formula = fit_glm_f$formula,data = data_lm)
# # Create an importance based on mean decreasing gini
# # importance(fit_rf)
# # importance(fit_rpart)
# # compare the feature importance with varImp() function
# # varImp(fit_rf)
#
# # Create a plot of importance scores by random forest
# # varImpPlot(fit_rf)
#
#
# rdb$q_116_VAS_v7 %>% factor() %>% summary()
# data$q_116_VAS_v7 %>% factor() %>% summary()
# # gbm::gbm(fit_glm_f$formula,data = data)
#
#
# #Check if a symptom var is positive in a given case
# # data[1,] %>% select(as.character(symptomsdf[1,1])) %>%
# #   {ifelse(is.na(.),
# #           0,
# #           ifelse(.="yes"),
# #   )}
#
#
#
#
# #### Correlate ANASCORE wit VAS ####
# rdb %>%
#   filter(!is.na(d_severity_rm),
#          !is.na(q_116_VAS_v7)) %>%
#   ggplot(aes(x = ANAscore,y = factor(q_116_VAS_v7), color= q_113_respiratory_arrest)) +
#   geom_jitter()+
#   facet_grid(relevel(severity_brown,"severe")~d_severity_rm)

data %>%
  filter(!is.na(d_severity_rm),
         !is.na(q_116_VAS_v7)) %>%
  ggplot(aes(x = ANAscore,y = factor(q_116_VAS_v7), color= q_140_fatal)) +
  geom_jitter()+
  facet_grid(.~d_severity_rm)

data %>%
  filter(!is.na(d_severity_rm),
         !is.na(q_116_VAS_v7)) %>%
  ggplot(aes(x = ANAscore,y = factor(q_116_VAS_v7), color= q_140_fatal)) +
  geom_jitter()#+
#  facet_grid(.~d_severity_rm)

data %>%
  filter(#!is.na(d_severity_rm),
    !is.na(q_116_VAS_v7)) %>%
  ggplot(aes(x = d_age,y = as.numeric(q_116_VAS_v7)))+#, color= q_140_fatal)) +
  geom_jitter()+
  geom_smooth(aes(mean(as.numeric(q_116_VAS_v7))))

data %>%
  group_by(d_age) %>%
  summarize(vas = mean(as.numeric(q_116_VAS_v7),na.rm=T)) %>%
  ggplot(aes(d_age, vas))+
  geom_point()+
  geom_smooth()+
  labs(y = "mean VAS score", x = "Age [years]")


data %>%
  group_by(d_age) %>%
  summarize(vas = mean(as.numeric(ANAscore),na.rm=T)) %>%
  ggplot(aes(d_age, vas))+
  geom_point()+
  geom_smooth()+
  labs(y = "mean ANAscore", x = "Age [years]")
cowplot::plot_grid(
  data %>%
    #group_by(d_age) %>%
    #summarize(vas = mean(as.numeric(ANAscore),na.rm=T)) %>%
    ggplot(aes(d_elicitor_gr5, ANAscore,fill=b_sex))+
    geom_boxplot()+
    #geom_smooth()+
    labs(y = "mean ANAscore", x = "Age [years]")+
    theme(legend.position = "none"),

  data %>%
    #group_by(d_age) %>%
    #summarize(vas = mean(as.numeric(ANAscore),na.rm=T)) %>%
    ggplot(aes(d_elicitor_gr5, as.numeric(q_116_VAS_v7),fill=b_sex))+
    geom_boxplot()+
    #geom_smooth()+
    labs(y = "mean VAS", x = "Age [years]")
)

yesonly <- function(x){
  ifelse(is.na(x),
         "no",
         ifelse(x=="yes",
                "yes",
                "no"))
}

data %>%
  group_by(factor(ANAscore)) %>%
  summarise(mean(as.numeric(q_116_VAS_v7),na.rm=T))

data$q_116_VAS_v7 %<>% as.numeric(as.character())

data %>%
  select()


fit <- glm(q_116_VAS_v7~ANAscore+
             q_130_biphasic_v4+
             q_562_intensive_care_v6+
             q_522_volume+
             #q_530_adren2_time_in_min_v5+
             q_522_dopamine,
           data = data)
summary(fit)
#### Assign VAS to other variables ####



#### Classification problem ####

# require(class)
# cls <-data$q_116_VAS_v7
# data %>%
#   filter(!is.na(q_116_VAS_v7))
# knn(train = data,)

###### Select the features for severity models ######

testInsectsbinomial %>%
  filter(section=="cofactors") %>%
  arrange(desc(Cramer)) %>%
  select(1,8,9,10,pval) %>%{.[c(14,16:18,22,23,24,26,35,43),]} %>%
  pull(variableName)

# FEATURE SELECTION########
# The features should contain:
# 1. Atopic diseases groupped together
# 2. Used medis grouped according to their background and effect on the anaphylaxis
# 3. Age
# 4. Sex
# 5. elicitor
# 6. Therapy? Not feasable
# 7. Tryptase value

# Build the models for RM, ANASCORE, BROWN and VAS
# Severity as a response and clinical features as predictors (include tryptase)

# Model for RM should take into account the fact that we have multiple level response

# Check the collinearity of different assessment scales intra


#### Tryptase and cardio symptoms ####


# function
tryp_assoc <- function(symptom){
  quo(symptom)
  tryp_ROA <- age_sex_matched %>%
    filter(!is.na(tryp_cat),
           !!symptom != "unknown",
           !is.na(grouping)) %>%
    mutate(tryp_cat = tryp_cat %>% factor(levels=c("low","high"),labels = c("< 8","> 8"))) %>%
    group_by(tryp_cat,!!symptom,
             grouping) %>%
    summarize(n=n())

  tab <- tryp_ROA
  plot <- tryp_ROA %>%
    ggbarplot(x=as_name(symptom),
              fill = "tryp_cat",
              y = "n",
              facet.by = "grouping",
              position = position_fill(),
              palette = c( "#848484","#1f1f1f"))+
    labs(y = "proportion",
         fill = "tryptase")

  test <- tryp_ROA %>% spread(key=!!symptom,value = n) %>%
    group_by(grouping) %>%
    nest() %>%
    {map(.$data,function(x){
      x[,2:3] %>% chisq.test
    })}
  return(list(tab,plot,test))
}

cardiac_typtase_effect <- names(data) %>% broom::tidy() %>%
  filter(grepl(x,pattern = "q_114")) %>% pull() %>%
  map(function(x){
    tryp_assoc(as.name(x))
  })

output_plot_tryp <- ggarrange(
  cardiac_typtase_effect[[6]][[2]]+
    theme_classic()+
    labs(x = "cardiac arrest")+
    theme(legend.position = "none"#,
          #axis.title.x = element_text(angle= 10,hjust = 0.5,vjust = 1)
    ),
  cardiac_typtase_effect[[2]][[2]]+
    theme_classic()+
    labs(y = "proportion",
         x = "loss of consciousness",
         fill = "tryptase")+
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank()#,
          #axis.title.x = element_text(angle= 5,hjust = 1)
    ),
  widths = c(1,1.2),
  align = "hv",
  common.legend = T,
  legend = "top"
)

fit_try_cardiac <- glm(q_114_cardiac_arrest ~ tryp_cat+grouping, data = age_sex_matched %>%
                         filter(!is.na(tryp_cat),
                                q_114_cardiac_arrest != "unknown",
                                !is.na(grouping)) %>%
                         dplyr::select(q_114_cardiac_arrest,tryp_cat,grouping),
                       family = binomial)
summary(fit_try_cardiac)

#### Tryptase+cardio+insect ####

# function
tryp_assoc_insect <- function(symptom){
  tryp_ROA <- age_sex_matched %>%
    filter(!is.na(tryp_cat),
           !!symptom != "unknown",
           !is.na(d_insect_gr4),
           d_insect_gr4!="other") %>%
    mutate(tryp_cat = tryp_cat %>% factor(levels=c("low","high"))) %>%
    group_by(tryp_cat,!!symptom,
             d_insect_gr4) %>%
    summarize(n=n())

  tab <- tryp_ROA
  plot <- tryp_ROA %>%
    ggbarplot(x=as_name(symptom),
              fill = "tryp_cat",
              y = "n",
              facet.by = "d_insect_gr4",
              position = position_fill(),
              palette = "lancet")+
    labs(y = "proportion",
         fill = "tryptase")

  test <- tryp_ROA %>% spread(key=!!symptom,value = n) %>%
    group_by(d_insect_gr4) %>%
    nest()# %>%
    #{map(.$data,function(x){
    #  x[,2:3] %>% chisq.test
    #})}
  return(list(tab,plot,test))
}

cardiac_typtase_insect_effect <- names(data) %>% broom::tidy() %>%
  filter(grepl(x,pattern = "q_114")) %>% pull() %>%
  map(function(x){
    tryp_assoc_insect(as.name(x))
  })

cardiac_typtase_insect_effect[[6]][[2]]
# this shows thatthere was no effect between different insects
# especially between wasps and yellow jackets

#### ACE Inhibitors and insects-tryptase-axis ####


cardiac_ace <- function(s){
  age_sex_matched %>%
    filter(!is.na(q_423_ace),
           q_423!="unknown",
           !!s != "unknown",
           !is.na(grouping)) %>%
    group_by(q_423_ace,
             !!s,
             grouping) %>%
    summarize(n=n()) %>%
    ggbarplot(x = as_name(s),
              fill = "q_423_ace",
              y = "n",
              facet.by = "grouping",
              position = position_fill(),
              palette = "lancet")+
    labs(y = "proportion",
         fill = "ACE-Inhibitors")
}

cardiac_ace_plots <- names(data) %>% broom::tidy() %>%
  filter(grepl(x,pattern = "q_114")) %>% pull() %>%
  map(function(x){
    cardiac_ace(as.name(x))
  })

plot_ace_cardiacs <-
  ggarrange(cardiac_ace_plots[[1]]+
              labs(x = "cardiologic symptoms"),
            cardiac_ace_plots[[6]]+
              labs(x = "cardiac arrest")+
              theme(axis.title.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.text.y = element_blank()),
            widths = c(1, 0.8),
            common.legend = T)


#### Very important!!!! above!
# This has to do with the cascade of

cardiac_beta <- function(s){
  age_sex_matched %>%
    filter(!is.na(q_423_beta),
           q_423_beta!="unknown",
           !!s != "unknown",
           !is.na(grouping)) %>%
    group_by(q_423_beta,
             !!s,
             grouping) %>%
    summarize(n=n()) %>%
    ggbarplot(x = as_name(s),
              fill = "q_423_beta",
              y = "n",
              facet.by = "grouping",
              position = position_fill(),
              palette = "lancet")+
    labs(y = "proportion",
         fill = "Beta blockers")
}

cardiac_beta_plots <- names(data) %>% broom::tidy() %>%
  filter(grepl(x,pattern = "q_114")) %>% pull() %>%
  map(function(x){
    cardiac_beta(as.name(x))
  })

plot_beta_cardiacs <-
  ggarrange(cardiac_beta_plots[[6]]+
              labs(x = "cardiac arrest"),
            cardiac_beta_plots[[5]]+
              labs(x = "arrythmia")+
              theme(axis.title.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.text.y = element_blank()),
            cardiac_beta_plots[[4]]+
              labs(x = "chest pain / angina")+
              theme(axis.title.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.text.y = element_blank()),
            widths = c(1, 0.8,0.8),
            nrow = 1,
            ncol = 3,
            common.legend = T)


cardiac_masto <- function(s){
  age_sex_matched %>%
    filter(!is.na(q_410_masto_cur),
           q_410_masto_cur!="unknown",
           !!s != "unknown",
           !is.na(grouping)) %>%
    group_by(q_410_masto_cur,
             !!s,
             grouping) %>%
    summarize(n=n()) %>%
    ggbarplot(x = as_name(s),
              fill = "q_410_masto_cur",
              y = "n",
              facet.by = "grouping",
              position = position_fill(),
              palette = "lancet")+
    labs(y = "proportion",
         fill = "MAstocytosis")
}

cardiac_masto_plots <- names(data) %>% broom::tidy() %>%
  filter(grepl(x,pattern = "q_114")) %>% pull() %>%
  map(function(x){
    cardiac_masto(as.name(x))
  })

plot_beta_cardiacs <-
  ggarrange(cardiac_beta_plots[[6]]+
              labs(x = "cardiac arrest"),
            cardiac_beta_plots[[5]]+
              labs(x = "arrythmia")+
              theme(axis.title.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.text.y = element_blank()),
            cardiac_beta_plots[[4]]+
              labs(x = "chest pain / angina")+
              theme(axis.title.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.text.y = element_blank()),
            widths = c(1, 0.8,0.8),
            nrow = 1,
            ncol = 3,
            common.legend = T)




#### Very important!!!! above!
# This has to do with the cascade of


##### Figure MOR2  eli_green####
eli_green <-
  rdbp %>%
  select(b_reactiondate,
         grouping,
         q_340_insects,
         d_centres_country) %>%
  mutate(MOR = substr(b_reactiondate,4,5),
         q_340_insects = ifelse(
           is.na(q_340_insects),
           "non-IVA",
           as.character(q_340_insects)
           ) %>%
           factor(
             levels=c(
               "yellow jacket",
               "bee","hornet",
               "bumble-bee",
               "horsefly",
               "mosquito",
               "other",
               "non-IVA"))) %>%
  filter(
    !is.na(q_340_insects),
    MOR!="00",
    q_340_insects != "non-IVA") %>%
  mutate(MOR = car::recode(MOR,
                           "'01' = 'Jan';
                           '02' = 'Feb';
                           '03' = 'Mar';
                           '04' = 'Apr';
                           '05' = 'May';
                           '06' = 'Jun';
                           '07' = 'Jul';
                           '08' = 'Aug';
                           '09' = 'Sep';
                           '10' = 'Oct';
                           '11' = 'Nov';
                           '12' = 'Dec'")) %>%
  mutate(MOR = factor(MOR, levels = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))) %>%
  group_by(MOR,q_340_insects) %>%
  summarize(count = n()) %>%
  ggpubr::ggbarplot(x = "MOR",
    y = "count",
    fill = "q_340_insects")+
  labs(fill="Insect",x = "",y = "Cases [n]")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_blank(),
    #axis.ticks.x = element_blank(),
    legend.position = c(0.01,.98),
    legend.justification = c(0,1),
    #panel.background = element_rect(fill = "#c9ccc5"),
    legend.background = element_blank())+
  scale_fill_brewer(palette = "Greys",direction = -1)
  # scale_fill_manual(values = rev(c("#E5F5E0",
  #                                  "#C7E9C0",
  #                                  "#A1D99B",
  #                                  "#74C476",
  #                                  "#41AB5D",
  #                                  "#238B45",
  #                                  "#005A32")))

age_dens <- ggpubr::ggdensity(rdbp %>%
                    filter(!is.na(grouping)) %>%
                    mutate(grouping = ifelse(grouping == "insects", "IVA","non-IVA")) %>%
                    mutate(grouping = relevel(factor(grouping),"non-IVA")),
                  x = "d_age",
                  fill = "grouping",
                  position = "fill",
                  rug = T,
                  color = "grouping",
                  palette = "grey")+#manual_greens[c(4,1)])+
  theme(legend.position = "right")+
  labs(x = "Age [years]",fill = "", color = "")

#png(width = 400*2,height = 430*2,res = 300, filename = "elicitors_green.png",pointsize = 7)
#eli_green
#dev.off()


#png(width =700*2,height = 630*2,res = 300, filename = "IVAonly.png",pointsize = 10)
IVAonly <- rdb %>%
  mutate(q_340_insects = as.character(q_340_insects) %>%
           factor(levels=c("yellow jacket",
                           "bee",
                           "hornet",
                           "other")),
         d_centres_country = car::recode(d_centres_country,
                                        "'not EU' = 'Brazil'")) %>%
  plot.proportions("d_centres_country","q_340_insects",5)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title.x = element_blank())+
  labs(x = "Country",
       y="Proportion",
       fill = "Insect")+
    scale_fill_brewer(palette = "Greys",direction = -1)


right_panel <- ggpubr::ggarrange(age_dens,
                  IVAonly,
                  nrow = 2,
                  ncol = 1,
                  heights = c(1,2),
                  labels = c("B","C"))
fig_basic <- ggpubr::ggarrange(
  eli_green,
  right_panel,
  labels = c("A",""),
  widths = c(1,1.2)
)
# library(RColorBrewer)
# brewer.pal(8, "Greens")

#dev.off()

png(width =500*2,height = 630*2,res = 300, filename = "IVAcomplete.png",pointsize = 10)
gridExtra::grid.arrange(
  #cowplot::plot_grid(
  rdbp %>%
    select(b_reactiondate,grouping,q_340_insects,d_centres_country) %>%
    mutate(MOR = substr(b_reactiondate,4,5)) %>%
    filter(!is.na(q_340_insects),MOR!="00") %>%
    ggplot(aes(MOR,fill=q_340_insects))+
    geom_bar()+
    theme_classic()+
    labs(fill="Insect",x = "Month of the year")+
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = c(0.01,.98),
          legend.justification = c(0,1))+
    scale_fill_manual(values = rev(c("#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45","#005A32"))),
  rdbp %>%
    select(b_reactiondate,grouping,d_elicitor_gr5,d_centres_country) %>%
    mutate(MOR = substr(b_reactiondate,4,5),
           d_elicitor_gr5 = relevel(d_elicitor_gr5, "insects")) %>%
    filter(!is.na(grouping),MOR!="00") %>%
    group_by(MOR,grouping) %>%
    summarize(n = n()) %>%
    group_by(MOR) %>%
    summarise(prop = n[1]/sum(n)) %>%
    ggplot(aes(MOR,prop))+
    geom_bar(stat="identity")+
    theme_classic()+
    labs(x = "Month of the year",y = "Fraction of IVA")+
    scale_fill_manual(values = rev(c("#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45","#005A32"))),
  #theme(legend.position = c(0.01,.98),
  #      legend.justification = c(0,1)),
  heights = c(1,0.5),
  ncol = 1
)
dev.off()

require(ggpubr)

#### Figure Symptoms #####
severity_joined <- full_join(
  rdb %>%
    mutate(d_severity_rm = d_severity_rm %>% factor(),
           grouping = car::recode(grouping,
                                  "'insects'='IVA';
                                  'other'='non-IVA'")
           ) %>%
    filter(!is.na(d_severity_rm),
           !is.na(grouping)) %>%
    group_by(d_severity_rm,
             grouping) %>%
    summarize(n =n()),
  rdb %>%
    mutate(d_severity_rm = d_severity_rm %>% factor()) %>%
    mutate(grouping = ifelse(d_insect_gr4 == "yellow jacket",
                                 "yellow jacket",
                                 "other")) %>%
    filter(!is.na(d_severity_rm),
           !is.na(d_insect_gr4)) %>%
    group_by(d_severity_rm,
             grouping) %>%
    summarize(n = n()),

  by = c("d_severity_rm","grouping","n")
) %>% full_join(
  rdb %>%
    mutate(d_severity_rm = d_severity_rm %>% factor(),
           grouping = d_age_gr2) %>%
    filter(!is.na(d_severity_rm)
    ) %>%
    group_by(d_severity_rm,
             grouping) %>%
    summarize(n = n()),
  by = c("d_severity_rm","grouping","n")
) %>%
  data.frame(subset = c(rep("Elicitor",8),
                    rep("Species",8),
                    rep("Age group",8)))

severity_joined_brown <- full_join(
  rdb %>%
    mutate(d_severity_rm = severity_brown %>% factor(),
           grouping = car::recode(grouping,
                                  "'insects'='IVA';
                                  'other'='non-IVA'")
           ) %>%
    filter(!is.na(d_severity_rm),
           !is.na(grouping)) %>%
    group_by(d_severity_rm,
             grouping) %>%
    summarize(n =n()),
  rdb %>%
    mutate(d_severity_rm = severity_brown %>% factor()) %>%
    mutate(grouping = ifelse(d_insect_gr4 == "yellow jacket",
                             "yellow jacket",
                             "other")) %>%
    filter(!is.na(d_severity_rm),
           !is.na(d_insect_gr4)) %>%
    group_by(d_severity_rm,
             grouping) %>%
    summarize(n = n()),

  by = c("d_severity_rm","grouping","n")
    ) %>% full_join(
      rdb %>%
        mutate(d_severity_rm = severity_brown %>% factor(),
               grouping = d_age_gr2) %>%
        filter(!is.na(d_severity_rm)
        ) %>%
        group_by(d_severity_rm,
                 grouping) %>%
        summarize(n = n()),
      by = c("d_severity_rm","grouping","n")
    ) %>%
  data.frame(subset = c(rep("Elicitor",4),
                        rep("Species",4),
                        rep("Age group",4)))



check <- rbind(severity_joined %>%
  filter(grouping %in% c("IVA","non-IVA")) %>%
  spread(grouping,n) %>%
  {matrix(data = cbind(.$IVA,.$`non-IVA`),
          nrow = 4,ncol=2, byrow = F)} %>%
  t() %>%
  chisq.test() %>%
  broom::tidy(),


severity_joined %>%
  filter(subset == "Species") %>%
  spread(grouping,n) %>%
  {matrix(data = cbind(.[,3],.[,4]),
          nrow = 4,ncol=2, byrow = F)} %>%
  t() %>%
  chisq.test()%>%
  broom::tidy(),

severity_joined %>%
  filter(subset == "Age group") %>%
  spread(grouping,n) %>%
  {matrix(data = cbind(.[,3],.[,4]),
          nrow = 4,ncol=2, byrow = F)} %>%
  t() %>%
  chisq.test()%>%
  broom::tidy()
) %>%
  mutate(subset = c("Elicitor","Species","Age group"))


colnames(tempdf) <- c("skin",
                      "gastrologic",
                      "respiratory",
                      "cardiac",
                      "atopic",
                      "tryptase")
rownames(tempdf) <- c("max","min","IVA","non-IVA")

names(tempdf)[2] <- "gastrointestinal"
########## FIGURE GOES HERE #######
fig_symptoms <- ggarrange(
  # Symtoms
  ggplot(tests_matched %>%
           tidyr::gather(key = "group",
                         value = "Fraction",
                         c("fraq_1","fraq_2")) %>%
           filter(!is.na(Fraction),
                  section =="symptoms",
                  counts_1<2000,
                  pval < 1e-9
           ) %>%
           mutate(variableName = car::recode(variableName,
                                             "'q_112_nausea' = 'nausea';
                                             'q_114_loss_of_consciousness' = 'unconsciousness';
                                             'q_113_throat_tightness_v5' = 'throat tightness';
                                             'q_114_reductions_of_alertness' = 'reduced alertness';
                                             'q_112_incontinence' = 'incontinence';
                                             'q_112_abdominal_pain' = 'abdominal pain';
                                             'q_112_diarrhoea' = 'diarrhoea';
                                             'q_113_wheezing_expiratory_distress_v5' = 'expiratory distress';
                                             'q_113_rhinitis_v5' = 'rhinitis';
                                             'q_114_dizziness' = 'dizziness';
                                             'q_114_hypotension_collapse_v5' = 'hypotension'"
           ),
           group = car::recode(
             group,
             "'fraq_1' = 'IVA';
             'fraq_2' = 'non-IVA'"
           )),
         aes(reorder(variableName,desc(counts_1)),Fraction, fill = group))+
    geom_bar(stat = "identity", position = "dodge",na.rm = T)+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 30,hjust =1, size = 15),
          legend.position = c(0.99,0.99),
          legend.justification = c(1,1),
          axis.title.x.bottom = element_blank())+
    labs(x = "symptom",y = "proportion [%]", fill = "elicitor")+
    scale_fill_manual(values = c("#1f1f1f", "#848484"))+
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)),

  tempdf  %>%
  {.[,c(2,1,5,3,4,6)]} %>%
    rownames_to_column( var = "group" ) %>%
    mutate_at(vars(-group),funs(rescale)) %>%
    filter(group%in%c("IVA","non-IVA")) %>%
    ggradar(axis.label.offset = c(1.1),
            axis.label.size = 5,
            group.point.size = 3,
            legend.text.size = 7,
            grid.label.size = 0)+
    theme(legend.position = "none")+
    scale_color_manual(values = c("#1f1f1f", "#848484")),

  common.legend = F,
  labels = c("A","B"),
  nrow = 1, ncol =2,
  widths = c(1,1.1))

#### Cramer Plot just for reference ####
plot.Cramer(cramerFun(data =rdbp,
                      vars=c("q_114",
                             "q_112_incontinence",
                             "q_114_dizziness",
                             "q_112_abdominal_pain",
                             "q_114_loss_of_consciousness",
                             "q_114_hypotension_collapse_v5",
                             "q_114_reductions_of_alertness",
                             "q_111_angioedema",
                             "q_113_wheezing_expiratory_distress_v5"),
                      grouping = "grouping"))

plot.Cramer(supraCramerFun(data =rdbp,
               vars =c( "q_114",
               "q_112_incontinence",
               "q_114_dizziness",
               "q_112_abdominal_pain",
               "q_114_loss_of_consciousness",
               "q_114_hypotension_collapse_v5",
               "q_114_reductions_of_alertness",
               "q_111_angioedema",
               "q_113_wheezing_expiratory_distress_v5"),
               grouping = "grouping",
               subset = "d_age_gr2"))

##### COFACTORS FIGURE!!!! ####

require(ggpubr)
cof_fig<- ggarrange(
  ggplot()+
    background_image(png::readPNG("analysis/figures/figForestfinalrmr.png")),

  ggarrange(
    age_sex_matched %>%
    filter(!is.na(d_severity_rmr),
           !is.na(tryp_cat)) %>%
    mutate(d_severity_rmr = d_severity_rmr %>% factor()) %>%
    group_by(grouping,d_severity_rmr,tryp_cat) %>%
    summarize(n = n()) %>%
    ggbarplot(x= "tryp_cat",
              fill = "d_severity_rmr",
              y = "n",
              facet.by = "grouping",
              position = position_fill(reverse = F),
              palette = c("#848484","#1f1f1f"))+
    labs(y = "proportion",
         x = "tryptase levels",
         fill = "severity"),

    ggarrange(
      cardiac_typtase_effect[[6]][[2]]$data %>%
      ggbarplot(x = "q_114_cardiac_arrest",
                y = "n",
                fill = "tryp_cat",
                facet.by = "grouping",
                position = position_fill(),
                palette = c("#848484","#1f1f1f"))+
        labs(x = "cardiac arrest",
             fill = "tryptase [ng/ml]")+
      theme(axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.position = "top"#,
            #axis.title.x = element_text(angle= 5,hjust = 1)
      ),#+
      #theme(
      #  legend.text = element_blank(),
      #  legend.title = element_blank(),
      #  legend.key = element_rect(fill = "white")
      #), #+
      #scale_fill_discrete(
      #  guide = guide_legend(override.aes = list(color = "white"))
      #),#+
        #theme(legend.position = "none"#,
              #axis.title.x = element_text(angle= 10,hjust = 0.5,vjust = 1)
        #),
      cardiac_typtase_effect[[2]][[2]]$data %>%
      ggbarplot(x = "q_114_loss_of_consciousness",
                y = "n",
                fill = "tryp_cat",
                facet.by = "grouping",
                position = position_fill(),
                palette = c("#848484","#1f1f1f"))+
        labs(y = "proportion",
             x = "loss of consciousness",
             fill = "tryptase")+
        theme(axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              axis.line.y = element_blank(),
              axis.ticks.y = element_blank(),
              legend.position = c(.1,.99),
              legend.direction = "horizontal"#,
              #axis.title.x = element_text(angle= 5,hjust = 1)
        )+
      labs(fill = "tryptase [ng/ml]",
           x = "loss of consciousness"),
      widths = c(1,1),
      common.legend = T
    )+
      theme(plot.margin = unit(c(0.8,0,0,0),"lines")),

    nrow = 1,
    ncol = 2,
    widths = c(1,1.6)
    ),
  nrow = 2,
  ncol=1,
  heights = c(1,0.7),
  labels = c("A","B")
)

#### Heatmap Symptom+Therapy ####
#### For all cases
require(gplots)
# x <- rdb %>%
#   select(starts_with("q_1"),
#          starts_with("q_5")
#          ) %>%
#   map(function(x){
#     ifelse(is.na(x), 0,ifelse(x=="yes",1,0))
#   }) %>%
#   do.call(what = cbind) %>% as.matrix()
# x <- x[,which(apply(x,2,sum)>0)]

# dendr <- x %>% dist(method = "binary") %>%
#   hclust(method = "ward.D2") %>%
#   as.dendrogram #%>%
#   #ladderize %>%
#   #color_branches(k=4)
#
# heatmap.2(x,
#           density.info = "none",
#             trace = "none",
#           Rowv = dendr,
#           RowSideColors = rdb %>%
#             mutate(group = ifelse(grouping=="insects",
#                                   "yellow",
#                                   ifelse(grouping =="other",
#                                          "brown",
#                                          "green"))) %>%
#             select(group) %>% pull()
#             )
#
#
# heatmap.2(x,
#           density.info = "none",
#           trace = "none",
#           Rowv = dendr,
#           RowSideColors = rdb %>%
#             mutate(group = car::recode(d_elicitor_gr5,
#                                        "'insects'='yellow';
#                                        'food' = 'brown';
#                                        'drugs' = 'green';
#                                        'other' = 'blue';
#                                        'unknown' = 'grey';
#                                        'unkown' = 'grey'")) %>%
#             select(group) %>% pull() %>% as.character()
# )
#
#
x <- age_sex_matched %>%
  select(#starts_with("q_1"),
         #-q_112,
         #-q_113,
         #-q_114,
         #-starts_with("q_16"),
        starts_with("q_5"),
         -"q_522_unknown"
  ) %>%
  map(function(x){
    ifelse(is.na(x), 0,ifelse(x=="yes",1,0))
  }) %>%
  do.call(what = cbind) %>% as.matrix()

library(randomForest)
fit_rf=randomForest(grouping~., data=data.frame(x,grouping = age_sex_matched$grouping ))
rf_plot_treatment <- varImpPlot(fit_rf)
rf_ggplot_treatment <- rf_plot_treatment %>% data.frame(treatment = rownames(.),.) %>%
  filter(MeanDecreaseGini >10) %>%
  mutate(treatment = treatment %>% recode_variables()) %>%
  ggdotchart(x = "treatment",
             y = "MeanDecreaseGini",
             #color = "organ",
             palette = "lancet",
             sorting = "descending",
             add = "segments",
             #add.params = list(color = "organ"),
             rotate = T,
             #group = "organ",
             dot.size = 4,
             size = 1.5,
             font.label = list(size = 8, color = "black"),
             #label = "Cramer"
             ggtheme = theme_pubr())

vars <- importance(fit_rf) %>%
  broom::tidy() %>%
  arrange(desc(MeanDecreaseGini)) %>%
  filter(MeanDecreaseGini > 10)

x <- age_sex_matched %>%
  select(vars$.rownames) %>%
  map(function(x){
    ifelse(is.na(x), 0,ifelse(x=="yes",1,0))
  }) %>%
  do.call(what = cbind) %>% as.matrix()

#x <- x[,which(apply(x,2,sum)>5)]




#
# heatmap.2(t(x),
#           density.info = "none",
#           trace = "none",
#           #Rowv = dendr,
#           ColSideColors = age_sex_matched %>%
#             mutate(group = car::recode(d_elicitor_gr5,
#                                        "'insects'='yellow';
#                                        'food' = 'brown';
#                                        'drugs' = 'green';
#                                        'other' = 'blue';
#                                        'unknown' = 'grey';
#                                        'unkown' = 'grey'")) %>%
#             select(group) %>% pull() %>% as.character(),
#           distfun= function(x) dist(x, method="binary"),
#           hclustfun=function(x) hclust(x, method="ward.D2")
# )


heatmap.2(t(x),
          density.info = "none",
          trace = "none",
          #Rowv = dendr,
          ColSideColors = age_sex_matched %>%
            mutate(group = ifelse(grouping=="insects",
                                  "green",
                                  ifelse(grouping =="other",
                                         "white",
                                         "gray"))) %>%
            select(group) %>% pull() %>% as.character(),
          distfun= function(x) dist(x, method="binary"),
          hclustfun=function(x) hclust(x, method="ward.D")
)

# require(heatmap.plus)
# heatmap.plus(t(x),
#           density.info = "none",
#           trace = "none",
#           ColSideColors = age_sex_matched %>%
#             mutate(group = ifelse(grouping=="insects",
#                                  "green",
#                                  ifelse(grouping =="other",
#                                         "white",
#                                         "gray")),
#                    severity = car::recode(
#                      d_severity_rm %>% as.factor(),
#                      "'1' = 'white';
#                      '2' = 'yellow';
#                      '3' = 'orange';
#                      '4' = 'red'"
#                    )) %>%
#             select(group,severity)  %>% as.matrix(),
#           distfun= function(x) dist(x, method="canberra"),
#           hclustfun=function(x) hclust(x, method="ward.D2")
# )



age_sex_matched %>%
  select(starts_with("q_1")
    # #-q_112,
    # #-q_113,
    # #-q_114,
    # #-starts_with("q_16"),
    # starts_with("q_5"),
    # -"q_522_unknown"
  ) %>%
  map(function(x){
    ifelse(is.na(x), 0,ifelse(x=="yes",1,0))
  }) %>%
  do.call(what = cbind) %>% as.matrix()


pca1 <- age_sex_matched %>%
  select(#starts_with("q_1"),
         starts_with("q_5")) %>%
  map(function(x){
    ifelse(is.na(x), 0,ifelse(x=="yes",1,0))
  }) %>%
  do.call(what = cbind) %>%
  as.data.frame() %>%

  {prcomp(.,center =T, scale. =F)}

# devtools::install_github("vqv/ggbiplot")
# require(ggbiplot)
# pca1 %>%
#   ggbiplot(pca1,
#            choices = 1:2)

# install.packages("ggfortify")
require(ggfortify)
autoplot(pca1, data=age_sex_matched, colour = "grouping" )

autoplot(prcomp(x, center = F,scale. = T),
         data=age_sex_matched,
         colour = "grouping")

#### Bradykinin effects ####
##We may suspect the effect of bradykinin due
#to the action on both cardiologic and gastrointestinal organs

# measures of associacian
# vcd::assocstats(table(age_sex_matched$q_,
#                       age_sex_matched$q_112_incontinence)[1:2,1:2])

#### Atopic and insect stings #####

o <-  rdb %>%
  filter(grouping=="insects") %>%
  select(b_case_id,atopy,d_age,b_sex) %>%
  {filter(.,complete.cases(.))} %>%
  mutate(grouping = ifelse(atopy=="yes",1,0)) %>%
  matchit(formula = as.formula("grouping~d_age+b_sex"),
  method = "nearest",
  ratio = 1) %>%
  match.data() %>%
  select(b_case_id) %>%
  pull() %>%
  {rdb[rdb$b_case_id %in% .,]}

o %>%
  group_by(atopy) %>%
  summarize(unconsciousness = sum(q_114_loss_of_consciousness=="yes",na.rm=T)) %>%
  ggbarplot(x = "atopy",
            #fill = "tryp_cat",
            y = "unconsciousness")

table(o$atopy,o$q_114_loss_of_consciousness)[,1:2] %>% chisq.test()

o %>%
  group_by(atopy) %>%
  summarize(cardiologic = sum(q_114=="yes",na.rm=T)) %>%
  ggbarplot(x = "atopy",
            #fill = "tryp_cat",
            y = "cardiologic")



##### Heatmap of therapy and symptoms ####
require(rlang)
phi_f <- function(data,
         grouping_var,
         grouping_var_val,
         var1,
         var2){
         data %>% #age_sex_matched %>%
  filter(!!sym(grouping_var)==grouping_var_val) %>%  #filter(grouping =="insects") %>%
  select(!!sym(var1),!!sym(var2)) %>%  #select(d_520_adren1,q_111) %>%
  #mutate(!!sym(var1) =
  #          ifelse(is.na(d_520_adren1),NA,
  #                 ifelse(d_520_adren1=="yes","yes",
  #                   ifelse(d_520_adren1=="no","no",NA)))) %>%
  {.[complete.cases(.),]} %>%
  {table(.[,1],.[,2])[1:2,1:2]} %>%
  assocstats() %>%
  .$phi
}

# phi_f(data = age_sex_matched,
#       grouping_var = "grouping",
#       grouping_var_val = "insects",
#       var1 = "d_520_adren1",
#       var2 = "q_111")
multiphi_f <- function(data,
                       grouping_var,
                       grouping_var_val,
                       vars1,
                       vars2){
  o <- map(vars2,function(treat){
    map(vars1,function(sympt){
      phi_f(data,
            grouping_var,
            grouping_var_val,
            var1 = sympt,
            var2 = treat)
    }) %>%
      do.call(what = rbind)
  }) %>% do.call(what = cbind)
  colnames(o) <- vars2
  rownames(o) <- vars1
  return(o)
}
# multiphi_f(
#   data = age_sex_matched,
#   grouping_var = "grouping",
#   grouping_var_val = "insects",
#   vars1 = c("q_111","q_112"),
#   vars2 = c("d_520_adren1")
# )
#
# matrix_treat_sympt<- multiphi_f(
#   data = age_sex_matched,
#   grouping_var = "grouping",
#   grouping_var_val = "insects",
#   vars1 = age_sex_matched %>%
# {names(.)[grepl(pattern = "q_11",x = names(.))]} %>%
#   .[c(1:8,10:40)],
# vars2 = c(age_sex_matched %>%
# {names(.)[grepl(pattern = "q_52",x = names(.))]} %>%
#   .[c(1:5,7:16,18:20,22:23,25:26)],
# age_sex_matched %>%
# {names(.)[grepl(pattern = "q_55",x = names(.))]} %>%
#   .[c(4:9,10,12,14:16,18,21)]
# ))
#
# h_insects <- heatmap.2(matrix_treat_sympt,
#           density.info = "none",
#           trace = "none",
#           distfun= function(x) dist(x, method="euclidean"),
#           hclustfun=function(x) hclust(x, method="ward.D"))
#
# h_insects$carpet %>% apply(MARGIN = 1,sum,na.rm = T) %>% {which(.>1.4)} %>% names()

matrix_treat_sympt<- multiphi_f(
  data = age_sex_matched,
  grouping_var = "grouping",
  grouping_var_val = "insects",
  vars1 = age_sex_matched %>%
  {names(.)[grepl(pattern = "q_11",x = names(.))]} %>%
    .[c(1:8,10:40)],
  vars2 = c("q_521_autoinj_v5",
            "q_522_adren_iv",
            "q_552_adren_iv_v5",
            "q_521_other_v5",
            "q_521_antih_v5"  ,
            "q_522_adren_im",
            "q_522_antih_oral",
            "q_552_adren_im_v5",
            "q_552_beta2_inhal_v5",
            "q_521_cortic_v5",
            "q_552_antih_oral_v5",
            "q_552_cortico_oral_v5",
            "q_552_volume_v5",
            "q_552_adren_inhal_v5",
            "q_552_beta2_iv_v5",
            "q_552_o2_v5",
            "q_522_o2",
            "q_522_beta2_inhal",
            "q_521_beta2_v5",
            "q_552_cortico_iv_v5",
            "q_552_antih_iv_v5" )
  )

matrix_treat_sympt_other<- multiphi_f(
  data = age_sex_matched,
  grouping_var = "grouping",
  grouping_var_val = "other",
  vars1 = age_sex_matched %>%
  {names(.)[grepl(pattern = "q_11",x = names(.))]} %>%
    .[c(1:8,10:40)],
  vars2 = c("q_521_autoinj_v5",
            "q_522_adren_iv",
            "q_552_adren_iv_v5",
            "q_521_other_v5",
            "q_521_antih_v5"  ,
            "q_522_adren_im",
            "q_522_antih_oral",
            "q_552_adren_im_v5",
            "q_552_beta2_inhal_v5",
            "q_521_cortic_v5",
            "q_552_antih_oral_v5",
            "q_552_cortico_oral_v5",
            "q_552_volume_v5",
            "q_552_adren_inhal_v5",
            "q_552_beta2_iv_v5",
            "q_552_o2_v5",
            "q_522_o2",
            "q_522_beta2_inhal",
            "q_521_beta2_v5",
            "q_552_cortico_iv_v5",
            "q_552_antih_iv_v5" )
)

# require(grid)
# require(gridG)
# grab.grob <- function(){
#   grid.echo()
#   grid.grab()
# }
#
#
# gl <- heatmap.2(matrix_treat_sympt,
#           density.info = "none",
#           trace = "none",
#           distfun= function(x) dist(x, method="euclidean"),
#           hclustfun=function(x) hclust(x, method="ward.D"),
#           dendrogram = "row",
#           key = F,keysize = 0.5)
# grab.grob()
#
# gridExtra::grid.arrange(
#   ncol = 2, nrow = 1,
#   heatmap.2(matrix_treat_sympt,
#             density.info = "none",
#             trace = "none",
#             distfun= function(x) dist(x, method="euclidean"),
#             hclustfun=function(x) hclust(x, method="ward.D"),
#             dendrogram = "row",
#             key = F,keysize = 0.5),
# heatmap.2(matrix_treat_sympt_other,
#           density.info = "none",
#           trace = "none",
#           distfun= function(x) dist(x, method="euclidean"),
#           hclustfun=function(x) hclust(x, method="ward.D"),
#           dendrogram = "row",
#           margins = c(9,12),key = F,keysize = 0.5
#           )
# )
#
# drawGridHeatmap  <- function(hm) {
#   plot(hm)
#   grab.grob()
# }
#
# gl <- lapply(list(heatmap.2(matrix_treat_sympt,
#                             density.info = "none",
#                             trace = "none",
#                             distfun= function(x) dist(x, method="euclidean"),
#                             hclustfun=function(x) hclust(x, method="ward.D"),
#                             dendrogram = "row",
#                             key = F,keysize = 0.5),
#                   heatmap.2(matrix_treat_sympt_other,
#                             density.info = "none",
#                             trace = "none",
#                             distfun= function(x) dist(x, method="euclidean"),
#                             hclustfun=function(x) hclust(x, method="ward.D"),
#                             dendrogram = "row",
#                             margins = c(9,12),key = F,keysize = 0.5
#                   )), drawGridHeatmap)
#
#
# library(gridGraphics)
# grab_grob <- function(){
#   grid.echo()
#   grid.grab()
# }
# heatmap.2(matrix_treat_sympt,
#            density.info = "none",
#            trace = "none",
#            distfun= function(x) dist(x, method="euclidean"),
#            hclustfun=function(x) hclust(x, method="ward.D"),
#            dendrogram = "row",
#            key = F,keysize = 0.5)
# g <- grab_grob()
# grid.newpage()
# heatmap.2(matrix_treat_sympt_other,
#           density.info = "none",
#           trace = "none",
#           distfun= function(x) dist(x, method="euclidean"),
#           hclustfun=function(x) hclust(x, method="ward.D"),
#           dendrogram = "row",
#           margins = c(9,12),key = F,keysize = 0.5
# )
# g2 <- grab_grob()
# grid.newpage()
#
# lay <- grid.layout(nrow = 1, ncol=2)
# pushViewport(viewport(layout = lay))
# grid.draw(editGrob(g, vp=viewport(layout.pos.row = 1,
#                                   layout.pos.col = 1, clip=F)))
# grid.draw(editGrob(g2, vp=viewport(layout.pos.row = 1,
#                                   layout.pos.col = 2, clip=F)))
# upViewport(1)

# ord <- hclust( dist(matrix_treat_sympt, method = "euclidean"), method = "ward.D" )$order
# ord
#
# matrix_treat_sympt %>%
#   as_tibble(rownames = "symptoms") %>%
#   gather(key="treatment",value = "phi",2:23) %>%
# ggplot(aes(treatment,symptoms,fill=phi))+
#   geom_tile()



rownames(matrix_treat_sympt) %<>% recode_variables()
colnames(matrix_treat_sympt) %<>% recode_variables()

rownames(matrix_treat_sympt_other) %<>% recode_variables()
colnames(matrix_treat_sympt_other) %<>% recode_variables()


require(heatmaply)
matrix_venom <- matrix_treat_sympt %>%
  scale() %>%
  heatmaply(hclust_method = "ward.D",
            dist_method = "euclidean",
            #Colv = "none",
            k_row = 2,
            k_col =2,
            colors = c("#FFFFFF","#000066"),
            plot_metod = "ggplot",
            return_ppxpy=TRUE
            )
matrix_other <- matrix_treat_sympt_other %>%
  scale() %>%
  heatmaply(hclust_method = "ward.D",
            dist_method = "euclidean",
            #Colv = "none",
            k_row = 2,
            k_col = 2,
            colors = c("#FFFFFF","#000066"),
            plot_metod = "ggplot"
            #return_ppxpy=TRUE
  )

### refractory cases

yj_bee <- rdb %>%
  mutate(d_severity_rm = d_severity_rm %>% factor()) %>%
  mutate(grouping = ifelse(d_insect_gr4 == "yellow jacket",
                           "yellow jacket",
                           ifelse(d_insect_gr4 =="bee", "bee","other"))) %>%
  filter(!is.na(d_severity_rm),
         !is.na(d_insect_gr4)) %>%
  group_by(d_severity_rm,
           grouping) %>%
  summarize(n = n())

