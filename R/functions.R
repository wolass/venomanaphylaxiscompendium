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
    dplyr::select(b_case_id,!!!grouping_var,!!!matching_vars) %>%
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
    dplyr::select(b_case_id) %>% pull()

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

compare <- function(x,percent =T,rounding =1){
  if(percent == T){
    paste0(round(x[1]*100,rounding),
           "% vs. ",
           round(x[2]*100,rounding),"%")
  } else {
    paste0(round(x[1],rounding),
           " vs. ",
           round(x[2],rounding))
  }
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
              'q_521_antih_v5'   = 'antihistamines 1st';
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

#'This function makes a chisq test on selected variables given a grouping variable.
#'It outputs a p.value
#'
#'
#'
funChi <-function(var,grouping_var){
  var %>%
    data.frame(gr = grouping_var) %>%
    {chisq.test(.[,1],.[,2])$p.value}
}

funN1 <- function(var, grouping_var){
  var %>%
    data.frame(grouping_var) %>% table %>%
    {.[2,1]}
  #group_by(grouping_var) %>%
  #summarise(count=n())
}
funN2 <- function(var, grouping_var){
  var %>%
    data.frame(grouping_var) %>% table %>%
    {.[2,2]}
  #group_by(grouping_var) %>%
  #summarise(count=n())
}
funF1 <- function(var, grouping_var){
  var %>%
    data.frame(grouping_var) %>%
    table %>% {.[1:2,]} %>%
    prop.table(2) %>%
    {.[2,1]}
  #group_by(grouping_var) %>%
  #summarise(count=n())
}

funF2 <- function(var, grouping_var){
  var %>%
    data.frame(grouping_var) %>%
    table %>% {.[1:2,]} %>%
    prop.table(2) %>%
    {.[2,2]}
  #group_by(grouping_var) %>%
  #summarise(count=n())
}
funC <- function(var, grouping_var){
  var %>%
    data.frame(grouping_var) %>%
    table %>% {.[1:2,]} %>%
    assocstats() %>% {.$cramer}
}

#' This function takes grouping vector and a DATA FRame ONLY OF  variables to test
#'
ChiF <- function(variablesDF,grouping){
  # C <- lapply(variablesDF,funC,grouping)
  data.frame(variables = names(variablesDF),
             pval = lapply(variablesDF,funChi,grouping) %>% unlist,
             counts_1 = lapply(variablesDF,funN1,grouping) %>%unlist,
             counts_2 = lapply(variablesDF,funN2,grouping) %>%unlist,
             fraq_1 = lapply(variablesDF,funF1,grouping) %>%unlist,
             fraq_2 = lapply(variablesDF,funF2,grouping) %>%unlist,
             Cramer = lapply(variablesDF,funC,grouping) %>% unlist
  )
}

variableSelectionTab <- function(data){
  data.frame(variableName=data %>% names(),
             type = data %>%lapply(class) %>% unlist,
             level = data %>% lapply(function(x){levels(x) %>% length}) %>% unlist(),
             section = c(rep("general",17),
                         rep("symptoms",72-17),
                         rep("location",77-72),
                         rep("previous reactions",86-77),
                         rep("diagnostic",114-86),
                         rep("triggers",184-114),
                         rep("cofactors",234-184),
                         rep("management",307-234),
                         rep("prophylaxis",369-307),
                         rep("other",length(data[1,])-369))
  )
}

#' Do a bunch of statistial tests on a data frame
#' This function automatically chooses test and performs them on the anaphylaxis registry database
#' @param grouping Factor with two levels to group the cases into case control sets
#' @param rdb database on which to perform the test (you can restrict it)
#' @return Tbl
#' @export
makeTests <- function(groups,rdb){
  grouping <- rdb[,groups]
  variableSelectionTab <- variableSelectionTab(rdb)
  variableSelectionTab %<>% # Make division for tests
    rowwise() %>%
    mutate(fun = if(level==3|level==2){
      "chi"
    } else if(variableName%in%c("q_212_tryptase_value_v5","d_age","ANAscore")) {
      "t.test"
    } else if (variableName%in%c("d_severity_rm","d_organ_sum",
                                 "q_116_VAS_v7", "q_120_time_between_v4",
                                 "q_131_biphasic_time_v4","q_141_fatal_time_v5",
                                 "q_632_autoinj_number_v7")){
      "Kruskal-Wallis"
    } else {
      "not tested"
    })

  # Perform chi2 test
  varDF <- rdb[,variableSelectionTab %>% filter(fun=="chi") %>%
                 dplyr::select(variableName) %>% pull %>% as.character()]
  ts <- lapply(varDF, # Check which tables are givin an error with chisq test
               function(x){ # You may implement also fisher tests here later
                 table(x,grouping)[1:2,] # get rid of the unknown row
               }) %>%
    lapply(function(x){
      length(which(x==0))
    }) %>% lapply(function(x){
      ifelse(x==0,T,F) # Gives True output if there are no 0 counts in any table cells
    }) %>% unlist()
 # varDF[,ts]
  #ChiF(varDF[ts],grouping) # put this into the dataframe
  variableSelectionTab %<>% full_join(ChiF(varDF[,ts],grouping),
                                      by = c("variableName"="variables")
  )
  #perform t.tests
  varDF <- rdb[,variableSelectionTab %>% filter(fun=="t.test")
               %>% dplyr::select(variableName) %>% pull %>% as.character()]
  #chek if t.test possible
  varDF %>%  lapply(function(column){
    split(column,grouping) %>% lapply(function(x){
      is.na(x) %>% `!` %>% sum() %>% {ifelse(.<2,F,T)}
    }) %>% rbind %>% all()
  }) %>% cbind()
  grouping %<>% factor()
  ttsts <- lapply(varDF,function(x){
    obj <- t.test(x[grouping==levels(grouping)[1]],
                  x[grouping==levels(grouping)[2]],
                  na.rm=T)
    return(c(pval=obj$p.value,
             mean1=obj$estimate[1],
             mean2 = obj$estimate[2]))
  }) %>% do.call(what=rbind) %>% as.data.frame() %>%
  {data.frame(variableName=rownames(.),.,stringsAsFactors = F)}
  # Join with the maind DF wit h substituting the missing values
  variableSelectionTab %<>%
    full_join( ttsts)
    #mutate(pval = coalesce(pval.y, pval.x)) %>%
    #dplyr::select(-c(pval.x,pval.y))
  # Perform Kruskall
  varDF <- rdb[,variableSelectionTab %>% filter(fun=="Kruskal-Wallis")
               %>% dplyr::select(variableName) %>% pull %>% as.character()]
  # Chek if kruskal possible
  varDF %<>%  lapply(function(column){
    split(column,grouping) %>% lapply(function(x){
      is.na(x) %>% `!` %>% sum() %>% {ifelse(.<2,F,T)}
    }) %>% rbind %>% all()
  }) %>% unlist() %>% cbind() %>%
  {varDF[,.]}
  ktest <-lapply(varDF,function(x){
    p <- kruskal.test(x,grouping)$p.value
    m <- split(x,grouping) %>% lapply(function(x){
      #if(length(x)   ###checkifgood
      median(x,na.rm=T)
    }) %>% do.call(what=cbind)
    IQR <- split(x,grouping) %>% lapply(function(x){
      #if(length(x)   ###checkifgood
      IQR(x,na.rm=T)
    }) %>% do.call(what=cbind)
    return(c(pval=p,
             median1 = m[1],
             median2 = m[2],
             IQR1 = IQR[1],
             IQR2 = IQR[2]))
  }) %>% do.call(what = rbind) %>%
    data.frame() %>%
    {data.frame(variableName=rownames(.),
                .,
                stringsAsFactors = F)}
  variableSelectionTab %<>% full_join(ktest) #%>%
    #mutate(pval = coalesce(pval.y, pval.x)) %>%
    #dplyr::select(-c(pval.x,pval.y))
  return(variableSelectionTab)
}



## Creat a formatted data frame

makeOdds <- function(data,var1,var2){ # var1 to predyktor, var2 - outcome
  if(var1=="b_sex"){
    or <- questionr::odds.ratio(data[,var1],data[,var2])
  } else {
    or <- questionr::odds.ratio(data[,var1] %>% factor(levels=c("no","yes")),data[,var2])
    data.frame(V = var1,
               p = format(or$p,scientific = T,digits = 1),
               mean = or[[1]],
               lower = or[[2]],
               upper = or[[3]],
               stringsAsFactors = F)
  }
}
makeDF <- function(data,variables,grouping,outcome){
  o <-lapply(variables,function(x){
    grouping <- data[,grouping]
    grA <- which(grouping==levels(grouping)[1])
    grB <- which(grouping==levels(grouping)[2])
    rbind(
      #c(makeOdds(data =data, var1 = x, var2 = outcome),S="all"),
      c(makeOdds(data =data[grA,], var1 = x, var2 = outcome),S="VIA"),
      c(makeOdds(data =data[grB,], var1 = x, var2 = outcome),S="non-VIA"),
      NA
    )
  }) %>% do.call(what=rbind) %>% data.frame()
  return(rbind(rep(NA,length(o[1,])),o) %>% data.frame())
}


#makeDF(rdb,"q_114",grouping,"severity_brown")
makeTableText <- function(makeDF){
  rbind(c("Cofactor","Elicitor","p"),
        makeDF[-c(1),c(1,6,2)])
}

makeForestPlot <- function(data,variables,grouping,outcome,cLow=0.2,cHigh = 4){
  df <- makeDF(data,variables,grouping,outcome)
  forestplot(labeltext = makeTableText(df),
             mean = df[,"mean"] %>% unlist,
             lower = df[,"lower"] %>% unlist,
             upper = df[,"upper"] %>% unlist,
             new_page = TRUE,
             is.summary=c(T,rep(F,length(df[,1])-1)),
             clip=c(cLow,cHigh),
             xlog=TRUE,
             col=fpColors(box="royalblue",line="darkblue", summary="royalblue"))
}

or_plot <- function(my_result){
  my_result %>%
    mutate_all(unlist) %>%
    mutate(variable = forcats::fct_inorder(V)) %>% #make variable into factor (otherwise ggplot will make alphabetically)
    ##start ggplot()
    ggplot(aes(x = forcats::fct_rev(V),
               y = mean))+
    geom_hline(yintercept = 1,
               colour = 'grey')+ #1.0 reference line
    geom_errorbar(aes(ymin = lower, ymax = upper),
                  width  = 0.2,
                  colour = '#41ae76',
                  size   = 0.9) + #OR lower-upper errorbars
    ##geom_linerange(aes(ymin=or_lower, ymax=or_upper), width=0.2, colour='#41ae76', size=0.9) + #alternative to errorbars
    geom_point(shape  = 18,
               size   = 2.5,
               colour = '#41ae76') +
    ##explanatory labels
    geom_text(aes(label = S),
              vjust = -0.5) +
    ##reference level labels
    #geom_text(data = my_result %>%
    #            filter(is.na(or)) %>%
    #            mutate(or_reference = 1.0),
    #          aes(label = value, x = variable, y = or_reference), angle = 90, vjust = -0.2, colour = '#737373') +
    coord_flip() +
    scale_y_log10() + #log scale breaks=c(0.2, 0.5, 1, 1.5, 2.0, 3.0), breaks=c(0.01, 0.1, 0.5, 1, 2.0)
    theme_bw() +                          #theme
    theme(axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title   = element_blank(),
          #panel.border = element_blank(),
          axis.text.x  = element_text(size=12, colour='black')
          #,panel.grid.major.y = element_blank() #optinally only remove horizontal grid lines (keeping odds ratio ones)
          ,panel.grid.major = element_blank() #these remove all grid lines
          ,panel.grid.minor = element_blank()
    )

  #or_plot

}

cramerFun <- function(data,grouping,vars){

  lapply(data[,vars],function(x){
    DescTools::CramerV(table(x,
                             data[,grouping]),
                       conf.level = 0.80,
                       method = "ncchisqadj")
  }) %>% do.call(what=rbind) %>%  {data.frame(var=rownames(.),.)} #%>%
  #data.frame(var=vars)
}

supraCramerFun <- function(data,vars,grouping,subset = NULL){
  if(is.null(subset)){
    o <- cramerFun(data,grouping,vars)
  } else{
    o <- split(data,data[,subset]) %>%  lapply(function(x){
      cramerFun(x,grouping,vars) %>% data.frame()
    })
  }
  out <- list(NULL)
  for(x in names(o)){
    out[[x]] <- mutate(o[[x]],subset = x)
  }
  out %>% do.call(what=rbind) %>% data.frame()
}

plot.Cramer <- function(x){
  if(is.null(x$subset)){
    p <- ggplot(x, aes(var,Cramer.V))+
      geom_bar(stat = "identity")+
      geom_errorbar(aes(ymin= lwr.ci,
                        ymax = upr.ci),
                    width=0.2)
  } else {
    p <- ggplot(x, aes(var,Cramer.V,fill=subset))+
      geom_bar(stat = "identity",position = "dodge")+
      geom_errorbar(aes(ymin= lwr.ci,
                        ymax = upr.ci),
                    position = position_dodge(0.9),
                    width=0.2)
  }
  p + geom_segment(mapping = aes(x=0.5,xend = length(levels(var))+0.5,
                                 y = 0.1, yend= 0.1),
                   linetype = 2, color = "gray")+
    geom_segment(mapping = aes(x=0.5,xend = length(levels(var))+0.5,
                               y = 0.3, yend= 0.3),
                 linetype = 2, color = "gray")+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 45,hjust=1))
}


f4 <- function(x){
  x[rdb$d_elicitor_gr5=="insects"]%>%
    split(rdb$b_sex[rdb$d_elicitor_gr5=="insects"]) %>%
    f3
}
f5 <- function(x){
  x[rdb$d_elicitor_gr5=="insects"]%>%
    split(rdb$d_age_gr2[rdb$d_elicitor_gr5=="insects"]) %>%
    f3
}

funtempspider <- function(var){
  age_sex_matched %>%
    group_by(grouping,!!var) %>%
    summarize(n= n()) %>%
    mutate(variable = deparse(var)) %>%
    group_by(grouping) %>%
    nest() %>%
    mutate(prop = map(data,function(x){
      x$n/sum(x$n)
    })
    ) %>% unnest() %>%
    as.matrix()
}

plot_mor <-function(){
gridExtra::grid.arrange(
  #cowplot::plot_grid(
  rdbp %>%
    dplyr::select(b_reactiondate,grouping,q_340_insects,d_centres_country) %>%
    mutate(MOR = substr(b_reactiondate,4,5)) %>%
    filter(!is.na(q_340_insects),MOR!="00") %>%
    ggplot(aes(MOR,fill=q_340_insects))+
    geom_bar(position = "fill")+
    theme_classic()+
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none"
    )+
    ylab("Proportion"),
  rdbp %>%
    dplyr::select(b_reactiondate,grouping,q_340_insects,d_centres_country) %>%
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
          legend.justification = c(0,1)),
  rdbp %>%
    dplyr::select(b_reactiondate,grouping,d_elicitor_gr5,d_centres_country) %>%
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
    labs(x = "Month of the year",y = "Fraction of insect elicited ANA"),
  #theme(legend.position = c(0.01,.98),
  #      legend.justification = c(0,1)),
  heights = c(0.4,1,0.5),
  ncol = 1
)
}

funtempspider2 <- function(var){
  age_sex_matched %>%
    group_by(grouping) %>%
    summarize(n= mean(!!var, na.rm = T)) %>%
    mutate(variable = deparse(var)) %>%
    group_by(grouping) %>%
    nest() %>%
    mutate(prop = map(data,function(x){
      x$n/sum(x$n)
    })
    ) %>% unnest() %>%
    as.matrix()
}

#' Check the variable in a dataset.
#' @param data a data frame
#' @param var  variable should be a binomial fctor
#' @param groupby a binomial factor variable in the data frame to divide into two groups
#'
#' @export
checkVarTab <- function(data, var, groupby){
  data[!is.na(data[,var])&!is.na(data[,groupby]),] %>%
    #filter(!is.na(get(groupby)&!is.na(get(var)))) %>%
    group_by(get(groupby),get(var)) %>%
    summarize(n = n()) %>%
    mutate(perc = n/sum(n)*100)
}


tryptase_plot <- function(data3){
  ggplot(data3[!is.na(data3$q_340_insects),],aes(q_212_tryptase_value_v5,color=q_340_insects))+
    geom_density()+
    xlim(0,30)+
    geom_density(mapping = aes(q_212_tryptase_value_v5),data3[data3$d_elicitor_gr5!="insects",],fill="black",alpha =0.2)
}



# Plot countries proportions
plot.proportions <- function(data,varx,vary,minN){
  ns <- data %>% filter(grouping == "insects") %>% group_by(get(varx)) %>%
    summarize(n=n()) %>% filter(n>minN)
  l=length(ns$n)
  ggplot(data[!is.na(data[,vary])&
                data[,varx]%in%ns$`get(varx)`,],
         aes(get(varx),
             fill = get(vary)))+
    geom_bar(position = "fill")+
    theme(axis.text.x = element_text(angle = 45,hjust=1))+
    annotate(geom = "text",
             x = 1:l,
             y = rep(0.3,l),
             label = paste("n =",ns%>% {as.character(.$n)}),
             angle = 90)
}
