rcalc <- function(var,level="yes"){
  tab <- table(data3[,var],data3$rANA)
  if(level %in% (data3[,var] %>% levels)){
    pn <- tab[level,"no"]
    pp <- tab[level,"yes"]
    l <- level
  } else {
    pn <- tab[2,"no"]
    pp <- tab[2,"yes"]
    l <- levels(data3[,var])[2]
  }
  nn <- data3$rANA %>% factor %>% summary %>% {.[1]}-pn
  np <- data3$rANA %>% factor %>% summary %>% {.[2]}-pp
  eval <- c(pn,pp,nn,np) %>% matrix(byrow = T,ncol = 2)
  prt <- (prop.table(eval,2)*100) %>% {signif(.,3)}
  p <- eval %>% fisher.test() %>% .$p.value %>% roP
  return(c(var,prt[1,1],prt[1,2],p,l))
}


roP<- function(x){
  ifelse(x > 0.001,
         round(x,3),ifelse(x==0,0,0.0001)
  )
}

pval<- function(x){
  ifelse(x > 0.001,
         paste0("p = ",round(x,3)),
         "p < 0.001")
}

f1 <- function(x){
  names(x) <- NULL
  return(do.call(t.test,x)$p.value %>% roP)
}
f2 <- function(x){
  x %>% table(data3$rANA) %>% {.[1:2,]} %>% fisher.test() %>% {.$p.value} %>% roP
}
f3 <- function(x){
  x %>%
    lapply(function(x){(length(which(x=="yes"))/length(x)*100) %>% roP})
}




#' Match patients
#'
#' Function match_patients takes a data frame from ANAreg and otuputs matched
#' case_id according to a given matching variable. The various grouping variables
#' have to still be implemented. Right now only the "grouping" variable is recognized
#'
#' @param grouping_var has to be a binomial var  of 1s and 0s
#' @param data is a data frame derived from ANAreg
#' @param matching_vars variable name to which we shouldmatch the two groups
#' @return vector of case ids that can be used for further analysis
#'
#'
match_patients <- function(data,
                           grouping_var,
                           matching_vars,
                           df=F,
                           grouping_var_level = "insects"
){
  o <- data %>%
    select(b_case_id,!!!grouping_var,!!!matching_vars) %>%
    {filter(.,complete.cases(.))} %>%
    mutate(grouping = ifelse(grouping==grouping_var_level,1,0)) %>%
    matchit(formula = as.formula(paste0(grouping_var,
                                        "~",
                                        paste0(matching_vars,
                                               collapse = "+"))
    ),
    method = "nearest",
    ratio = 1) %>%
    match.data() %>%
    select(b_case_id) %>% pull()

  if(df == T){
    o <- data[data$b_case_id %in% o,]
  }
  return(o)
}


unknown_to_na <- function(vect){
  ifelse(vect == "unknown",
         NA,
         ifelse(grepl(vect,pattern = "no"),
                0,
                vect)
  )
}


char_to_num_levels <- function(variable){
  levels(variable) <-  levels(variable) %>%
    car::recode("'one' = 1;
                'two' = 2;
                'three' = 3;
                'four' = 4;
                'five' = 5;
                'six' = 6;
                'seven' = 7;
                'eight' = 8;
                'nine' = 9") %>%
    unknown_to_na()
  return(variable %>%
           as.character() %>%
           as.numeric())
}

recode_variables <- function(x){
  car::recode(x,
              recodes =
                "'q_111'='cutaneous symptoms';
              'q_111_angioedema'= 'angioedema';
              'q_111_erythema_flush_v5' = 'flush';
              'q_111_pruritus' = 'pruritus';
              'q_111_urticaria' = 'urticaria';
              'q_111_conjunctivitis_v5' = 'cojunctivitis';
              'q_112' = 'gastrointestinal symptoms';
              'q_114_dizziness'='dizziness';
              'q_114_loss_of_consciousness'='loss of consciousness';
              'q_112_abdominal_pain'='abdominal pain';
              'q_113_wheezing_expiratory_distress_v5'='expiratory distress';
              'q_522_volume'='volume iv.';
              'q_561_hospital_admission_v6'='hospital admission';
              'q_562_intensive_care_v6'='intensive care';
              'q_112_nausea' = 'nausea';
              'q_113_throat_tightness_v5' = 'throat tightness';
              'q_114_reductions_of_alertness' = 'reduced allertness';
              'q_112_incontinence' = 'incontinence';
              'q_112_abdominal_pain' = 'abdominal pain';
              'q_112_diarrhoea' = 'diarrhoea';
              'q_113_wheezing_expiratory_distress_v5' = 'expiratory distress';
              'q_113_rhinitis_v5' = 'rhinitis';
              'q_114_dizziness' = 'dizziness';
              'q_114_hypotension_collapse_v5' = 'hypotension';
              'q_112_dysphagia_v5' = 'dysphagia';
              'q_112_vomiting' = 'vomiting';
              'q_113' = 'respiratory symptoms';
              'q_113_respiratory_arrest' = 'respiratory arrest';
              'q_113_dyspnea' = 'dyspnea';
              'q_113_chest_tightness_v5' = 'chest tightness';
              'q_113_cough_v5' = 'cough';
              'q_113_change_in_voice_v5' = 'change in voice';
              'q_113_stridor_inspiratory' = 'inspiratory stridor';
              'q_114' = 'cardiologic symptoms';
              'q_114_chest_pain_angina_v5' = 'chest pain';
              'q_114_palpitations_cardiac_arrythmia_v5' = 'cardiac arrythmia';
              'q_114_cardiac_arrest' = 'cardiac arrest';
              'q_114_tachycardia' = 'tachycardia';
              'q_115_dysarthria_v6' = 'dysarthria';
              'q_115_dysphonia_v6' = 'dysphonia';
              'q_115_hot_sweat_tremble_v6' = 'hot sweat';
              'q_115_tingle_hands_feet_paresthesia_v6' = 'paresthesia';
              'q_115_sight_disorder_v6' = 'sight diorder';
              'q_115_agony_v6' = 'agonizing pain';
              'q_115_cyanosis_pallor_v6' = 'cyanosis';
              'q_521_autoinj_v5' = 'adrenaline autoinjector';
              'q_522_adren_iv'   = 'adrenaline i.v. 1st';
              'q_552_adren_iv_v5'   = 'adrenaline i.v. 2nd';
              'q_521_other_v5'   = 'other drugs 1st';
              'q_521_antih_v5'   = 'antihistmines 1st';
              'q_522_adren_im'   = 'adrenaline i.m. 1st';
              'q_522_antih_oral'   = 'antihistamines p.o. 1st';
              'q_552_adren_im_v5'   = 'adrenaline i.m. 2nd';
              'q_552_beta2_inhal_v5'   = 'β2-agonists inh. 2nd';
              'q_521_cortic_v5'   = 'corticosteroids 1st';
              'q_552_antih_oral_v5'   = 'antihistamines p.o. 2nd';
              'q_552_cortico_oral_v5'   = 'corticosteroids p.o. 2nd';
              'q_552_volume_v5'   = 'volume i.v. 2nd';
              'q_552_adren_inhal_v5'   = 'adrenaline inh. 2nd';
              'q_552_beta2_iv_v5'   = 'β2-agonists i.v. 2nd';
              'q_552_o2_v5'   = '100% oxygen 2nd';
              'q_522_o2'   = '100% oxygen 1st';
              'q_522_beta2_inhal'   = 'β2-agonists inh. 1st';
              'q_521_beta2_v5'   = 'β2-agonists 1st';
              'q_552_cortico_iv_v5'   = 'corticosteroids i.v. 2nd';
              'q_552_antih_iv_v5'   = 'antihistamines i.v. 2nd';
              'q_522_antih_iv' = 'antihistamines i.v. 1st';
              'q_522_cortico_iv' = 'corticosteroids i.v. 1st';
              'q_522_other_v5' = 'other meds 1st';
              'q_522_cortico_oral' = 'corticosteroids p.o. 1st';
              'q_550_2nd_v5' = '2nd adrenaline dose';
              "
  )
}



#### cases selection ####
# Define a new variable - reaction type to exclude non ANA cases

correctLabels1 <- function(data) {
  data$q_120_time_between_v4 %<>% factor(levels = c("00 – 10 Minutes",
                                                    "11 - 30 Minutes",
                                                    "31 - 60 Minutes",
                                                    "61 - 120 Minutes",
                                                    "121 – 240 Minutes (2 – 4 hours)",
                                                    "more than 240 Minutes" ,
                                                    "more than 120 min (V4)"),
                                         labels=c("10","30",
                                                  "60","120","4h","more","4h-5"))
  return(data)
}



AssessAnaphylaxisDefinition <- function(data){
  sentinel_skin <- rep(F, nrow(data))
  sentinel_skin[which(data$q_111_angioedema=="yes"|
                        data$q_111_urticaria=="yes"|
                        data$q_111_erythema_flush_v5=="yes"|
                        data$q_111_pruritus=="yes")] <- T

  sentinel_respiratory <- rep(F, nrow(data))
  sentinel_respiratory[which(data$q_115_cyanosis_pallor_v6=="yes"|
                               data$q_113_stridor_inspiratory=="yes"|
                               data$q_113_wheezing_expiratory_distress_v5=="yes"|
                               data$q_113_dyspnea=="yes"|
                               data$q_113_respiratory_arrest=="yes")] <- T

  sentinel_cardio <- rep(F, nrow(data))
  sentinel_cardio[which(
    data$q_112_incontinence=="yes"|
      data$q_114_reductions_of_alertness=="yes"|
      data$q_114_loss_of_consciousness=="yes"|
      data$q_114_cardiac_arrest=="yes"|
      data$q_114_hypotension_collapse_v5=="yes")] <- T

  sentinel_gastro <- rep(F, nrow(data))
  sentinel_gastro[which(data$q_112_abdominal_pain=="yes"|
                          data$q_112_vomiting=="yes")] <- T

  reaction_type_brown <- rep(NA, length(data[,1]))
  # Correct the sum of organs involved.
  organ_sum <- rep(0,nrow(data))
  organ_sum <- organ_sum+ ifelse(!is.na(data$q_111),
                                 ifelse(data$q_111=="yes",1,0),
                                 0)
  organ_sum <- organ_sum+ ifelse(!is.na(data$q_112),
                                 ifelse(data$q_112=="yes",1,0),
                                 0)
  organ_sum <- organ_sum+ ifelse(!is.na(data$q_113),
                                 ifelse(data$q_114=="yes",1,0),
                                 0)
  organ_sum <- organ_sum+ ifelse(!is.na(data$q_114),
                                 ifelse(data$q_114=="yes",1,0),
                                 0)

  organ_sum %<>% factor()
  reaction_type_brown[which(organ_sum=="1"&data$d_111_sum%in%
                              c("one","two","three","four","five"))] <- "skin only"

  reaction_type_brown[which(sentinel_skin==T&(sentinel_cardio==T|sentinel_respiratory==T))] <- "anaphylaxis"

  levels(data$q_120_time_between_v4) <- levels(data$q_120_time_between_v4)[c(1:8,7)]

  # BELOW: I am not including the information about the time from exposure to reaction as it was largly missing. The cases are nevertheless supposed allergic reactions therefore the likelyhood of an allergen and occurence of these symptoms should be strict enough to be sufficient as a definition of anaphylaxis.

  reaction_type_brown[which(#data$q_120_time_between_v4%in%levels(data$q_120_time_between_v4)[3:7]&
    ((sentinel_skin==T&sentinel_cardio==T)|
       (sentinel_skin==T&sentinel_gastro==T)|
       (sentinel_skin==T&sentinel_respiratory==T)|
       (sentinel_cardio==T&sentinel_gastro==T)|
       (sentinel_cardio==T&sentinel_respiratory==T)|
       (sentinel_respiratory==T&sentinel_gastro==T)))] <- "anaphylaxis"
  reaction_type_brown[is.na(reaction_type_brown)] <-"ANA-def not met"

  reaction_type_brown <- as.factor(reaction_type_brown)
  return(reaction_type_brown)
}


# Add severity of anaphylaxis according to brown
severity_brown <- function(data){
  severity_brown <- rep(NA,length(data[,1]))
  severity_brown[which(data$reaction_type_brown=="anaphylaxis")]<-"mild"
  severity_brown[which(data$reaction_type_brown=="anaphylaxis"&
                         (data$q_115_cyanosis_pallor_v6=="yes"|
                            data$q_114_hypotension_collapse_v5=="yes"|
                            data$q_114_loss_of_consciousness=="yes"|
                            data$q_114_reductions_of_alertness=="yes"|
                            data$q_112_incontinence=="yes"|
                            data$q_140_fatal=="yes"|
                            data$q_114_cardiac_arrest=="yes"|
                            data$q_113_respiratory_arrest=="yes"))] <- "severe"
  severity_brown <- as.factor(severity_brown)
  return(severity_brown)
}
