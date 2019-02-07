#Setup
require(dplyr)
require(magrittr)
require(ggplot2)
require(forestplot)
require(summarytools)
st_options('bootstrap.css', FALSE)
#### GET THE DATA ########
# Note the path that we need to use to access our data files when rendering this document
load('../RefractoryAnaOrg/analysis/data/raw_data/data4.R')


####### Functions #############
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

data4 %<>% correctLabels1()


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

data4$reaction_type_brown <- AssessAnaphylaxisDefinition(data4)

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
data4$severity_brown <- severity_brown(data4)

other_insects <- c("horse fly",
                   "bumble-bee",
                   "mosquito",
                   "insects")
data4$d_insect_gr4 <- as.character(data4$q_340_insects)
data4$d_insect_gr4[which(data4$q_340_insects%in%other_insects)] <- "other"
data4$d_insect_gr4 %<>% factor


##### DIAGRAM ##############
library(DiagrammeRsvg)
require(DiagrammeR)

# Create a node data frame
nodes <-
  create_node_df(
    n = 7,
    type = "a",
    label= c(paste0("All cases in the\nEuropean Anaphylaxis Registry\n",length(data4$b_submitdate)),
             paste0("Formal anaphylaxis definition met\n",
                    length(which(data4$reaction_type_brown=="anaphylaxis"))),
             paste0("Insects as elicitors\n",
                    length(which(data4$d_elicitor_gr5== "insects" &
                                   data4$reaction_type_brown=="anaphylaxis"))),
             paste0("Other elicitors\n",
                    length(which(data4$d_elicitor_gr5!= "insects"&
                                   data4$reaction_type_brown=="anaphylaxis"))),
             paste(paste0("Yellow-jackets = ",
                    length(which(data4$q_340_insects== "yellow jacket"&
                                   data4$reaction_type_brown=="anaphylaxis"))),
                paste0("Bees = ",
                    length(which(data4$q_340_insects== "bee"&
                                   data4$reaction_type_brown=="anaphylaxis"))),
                paste0("Hornets = ",
                    length(which(data4$q_340_insects== "hornet"&
                                   data4$reaction_type_brown=="anaphylaxis"))),
                paste0("Other insects = ",
                    length(which(data4$d_insect_gr4=="other" &
                                   data4$reaction_type_brown=="anaphylaxis"))),
                sep="\n"
             ),
             paste0("Elicitor unknown\n",
                    length(which(data4$d_elicitor_gr5=="unkown" &
                                   data4$reaction_type_brown=="anaphylaxis"))),
             paste0("Other known\nelicitors\n",
                    length(which(data4$d_elicitor_gr5!="insects"&
                                   data4$d_elicitor_gr5!="unkown" &
                                   data4$reaction_type_brown=="anaphylaxis")))
             ),

    #color = c("red", "green",
    #          "grey", "blue"),
    #value = c(3.5, 2.6, 9.4, 2.7),
    shape= "rectangle",
    width = c(2.5,3,1.5,1.5,1.8,1.2,1.2),
    x = c(2, 2, 0.9, 3.1,1.1,4.7,3.1),
    y = c(0,-1,-2,-2,-3,-2,-3),
    height = c(rep(0.6,4),1,0.6,0.6))

edges <- create_edge_df(from = c(1,2,2,3,4,4),
                        to =   c(2,3,4,5,6,7))
# Add the node data frame to the
# graph object to create a graph
# with nodes
create_graph() %>%
  add_node_df(node_df = nodes) %>%
  add_edge_df(edge_df = edges) %>%
  add_global_graph_attrs(
    attr = "splines",
    value = 1,
    attr_type = "graph") %>%
  add_global_graph_attrs(value="black",
                         attr = "fontcolor",
                         attr_type = "node") %>% #render_graph()
  export_graph(file_name = "analysis/figures/flow.png", file_type = "png", title = NULL,
               width = NULL, height = NULL)


###Outcomes############

# variableSelctionTab <- data.frame(variableName=rdb %>% names(),
#                                   type = rdb %>%lapply(class) %>% unlist,
#                                   level = rdb %>% lapply(function(x){levels(x) %>% length}) %>% unlist(),
#                                   section = c(rep("general",17),
#                                               rep("symptoms",72-17),
#                                               rep("location",77-72),
#                                               rep("previous reactions",86-77),
#                                               rep("diagnostic",114-86),
#                                               rep("triggers",184-114),
#                                               rep("cofactors",234-184),
#                                               rep("management",307-234),
#                                               rep("prophylaxis",369-307),
#                                               rep("other",379-369))
# )

#### Symptoms ####

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
  temp <- var %>%
    data.frame(grouping_var) %>% table %>%
    {.[2,1]}
    #group_by(grouping_var) %>%
    #summarise(count=n())
}
funN2 <- function(var, grouping_var){
  temp <- var %>%
    data.frame(grouping_var) %>% table %>%
    {.[2,2]}
  #group_by(grouping_var) %>%
  #summarise(count=n())
}
funF1 <- function(var, grouping_var){
  temp <- var %>%
    data.frame(grouping_var) %>%
    table %>% {.[1:2,]} %>%
    prop.table(2) %>%
    {.[2,1]}
  #group_by(grouping_var) %>%
  #summarise(count=n())
}

funF2 <- function(var, grouping_var){
  temp <- var %>%
    data.frame(grouping_var) %>%
    table %>% {.[1:2,]} %>%
    prop.table(2) %>%
    {.[2,2]}
  #group_by(grouping_var) %>%
  #summarise(count=n())
}

#' This function takes grouping vector and a DATA FRame ONLY OF  variables to test
#'
ChiF <- function(variablesDF,grouping){
  data.frame(variables = names(variablesDF),
             pval = lapply(variablesDF,funChi,grouping) %>% unlist,
             counts_1 = lapply(variablesDF,funN1,grouping) %>%unlist,
             counts_2 = lapply(variablesDF,funN2,grouping) %>%unlist,
             fraq_1 = lapply(variablesDF,funF1,grouping) %>%unlist,
             fraq_2 = lapply(variablesDF,funF2,grouping) %>%unlist
  )
}


### Data Cleanin #####

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

#FINAL DATABASE RDB ######
rdb <- data4[data4$reaction_type_brown=="anaphylaxis",]
countries <-rdb %>%
  filter(d_elicitor_gr5=="insects") %>%
  select(d_centres_country) %>%
  group_by(d_centres_country) %>%
  summarize(n=n())%>% arrange(desc(n))

grouping <- ifelse(rdb$d_elicitor_gr5=="insects","insects",
                   ifelse(is.na(rdb$d_elicitor_gr5)|rdb$d_elicitor_gr5=="unkown",NA,"other")) %>%
  factor

group_triggerBinomial <- ifelse(rdb$d_elicitor_gr5=="insects","insects",
                                ifelse(rdb$d_elicitor_gr5=="unkown",NA,"other"))


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


###Statistic###########

#' Do a bunch of statistial tests on a data frame
#' This function automatically chooses test and performs them on the anaphylaxis registry database
#' @param grouping Factor with two levels to group the cases into case control sets
#' @param rdb database on which to perform the test (you can restrict it)
#' @return Tbl
#' @export
makeTests <- function(grouping,rdb){
  variableSelectionTab <- variableSelectionTab(rdb)
  variableSelectionTab %<>% # Make division for tests
    rowwise() %>%
    mutate(fun = if(level==3|level==2){
      "chi"
    } else if(variableName%in%c("q_212_tryptase_value_v5","d_age")) {
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
                 select(variableName) %>% pull %>% as.character()]
  ts <- lapply(varDF, # Check which tables are givin an error with chisq test
               function(x){ # You may implement also fisher tests here later
                 table(x,grouping)
               }) %>%
    lapply(function(x){
      length(which(x==0))
    }) %>% lapply(function(x){
      ifelse(x==0,T,F)
    }) %>% unlist()
  varDF[,ts]
  #ChiF(varDF[ts],grouping) # put this into the dataframe
  variableSelectionTab %<>% full_join(ChiF(varDF[,ts],grouping),
                                      by = c("variableName"="variables")
  )
  #perform t.tests
  varDF <- rdb[,variableSelectionTab %>% filter(fun=="t.test")
               %>% select(variableName) %>% pull %>% as.character()]
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
  variableSelectionTab %<>% full_join( ttsts, by = "variableName") %>%
    mutate(pval = coalesce(pval.y, pval.x)) %>%
    select(-c(pval.x,pval.y))
  # Perform Kruskall
  varDF <- rdb[,variableSelectionTab %>% filter(fun=="Kruskal-Wallis")
               %>% select(variableName) %>% pull %>% as.character()]
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
  variableSelectionTab %<>% full_join(ktest, by = "variableName") %>%
    mutate(pval = coalesce(pval.y, pval.x)) %>%
    select(-c(pval.x,pval.y))
  return(variableSelectionTab)
}

### Binomial trigger either insects or other
testInsectsbinomial <- makeTests(grouping,rdb=rdb) %>% arrange(pval)
# ggplot(testInsectsbinomial,aes(variableName,pval))+
#   geom_point()+
#   scale_y_log10()

testInsectsbinomialTab <- testInsectsbinomial %>% filter(pval<1e-30) %>%
  select(variableName,counts_1,counts_2,fraq_1,fraq_2,
         pval,
         section) %>%
         {split(.,.$section)}


testYJ <- ifelse(rdb$d_insect_gr4 == "yellow jacket",
       "y-j",
       "other")

testsYJ <- makeTests(testYJ,rdb = rdb) %>%
  arrange(pval) %>%
  select(variableName,counts_1,counts_2,fraq_1,fraq_2,
         pval,
         section) %>%
         {split(.,.$section)}



###### FOREST PLOT ####
# data <- gdata::read.xls(xls = "../../Downloads/forestplotdata.xlsx", stringsAsFactors=FALSE)
# length(data[,1])
# ## Labels defining subgroups are a little indented!
# subgps <- c(4,5,8,9,12,13,16,17,20,21,24,25,28,29,32,33)
# data$Variable[subgps] <- paste("  ",data$Variable[subgps])
#
# ## Combine the count and percent column
# np <- ifelse(!is.na(data$Count), paste(data$Count," (",data$Percent,")",sep=""), NA)
#
# ## The rest of the columns in the table.
# tabletext <- cbind(c("Subgroup","\n",data$Variable),
#                    c("No. of Patients (%)","\n",np),
#                    c("4-Yr Cum. Event Rate\n PCI","\n",data$PCI.Group),
#                    c("4-Yr Cum. Event Rate\n Medical Therapy","\n",data$Medical.Therapy.Group),
#                    c("P Value","\n",data$P.Value))
#
#
# forestplot(labeltext=tabletext,
#           graph.pos=3,
#           mean=c(NA,NA,data$Point.Estimate),
#           lower=c(NA,NA,data$Low),
#           upper=c(NA,NA,data$High),
#           title="Hazard Ratio",
#           xlab="     <---PCI Better---    ---Medical Therapy Better--->",
#           hrzl_lines=list("3" = gpar(lwd=1, col="#99999922"),
#                           "7" = gpar(lwd=60, lineend="butt", columns=c(2:6), col="#99999922"),
#                           "15" = gpar(lwd=60, lineend="butt", columns=c(2:6), col="#99999922"),
#                           "23" = gpar(lwd=60, lineend="butt", columns=c(2:6), col="#99999922"),
#                           "28" = gpar(lwd=60, lineend="butt", columns=c(2:6), col="#99999922")),
#           txt_gp=fpTxtGp(label=gpar(cex=1.25),
#                         ticks=gpar(cex=1.1),
#                         xlab=gpar(cex = 1.2),
#                         title=gpar(cex = 1.2)),
#           col=fpColors(box="black", lines="black", zero = "gray50"),
#           zero=1,
#           cex=0.9,
#           lineheight = "auto",
#           boxsize=0.5,
#           colgap=unit(6,"mm"),
#           lwd.ci=2,
#           ci.vertices=TRUE,
#           ci.vertices.height = 0.4)

# # An example of how the exponential works
# test_data <- data.frame(coef=c(2.45, 0.43),
#                         low=c(1.5, 0.25),
#                         high=c(4, 0.75),
#                         boxsize=c(0.5, 0.5))
# row_names <- cbind(c("Name", "Variable A", "Variable B"),
#                    c("HR", test_data$coef))
# test_data <- rbind(rep(NA, 3), test_data)
#
# forestplot(labeltext = row_names,
#            test_data[,c("coef", "low", "high")],
#            is.summary=c(TRUE, FALSE, FALSE),
#            boxsize   = test_data$boxsize,
#            zero      = 1,
#            xlog      = TRUE,
#            col = fpColors(lines="red", box="darkred"))
#
library(forestplot)
# Cochrane data from the 'rmeta'-package
# cochrane_from_rmeta <-
#   structure(list(
#     mean  = c(NA, NA, 0.578, 0.165, 0.246, 0.700, 0.348, 0.139, 1.017, NA, 0.531),
#     lower = c(NA, NA, 0.372, 0.018, 0.072, 0.333, 0.083, 0.016, 0.365, NA, 0.386),
#     upper = c(NA, NA, 0.898, 1.517, 0.833, 1.474, 1.455, 1.209, 2.831, NA, 0.731)),
#     .Names = c("mean", "lower", "upper"),
#     row.names = c(NA, -11L),
#     class = "data.frame")
#
# tabletext<-cbind(
#   c("", "Study", "Auckland", "Block",
#     "Doran", "Gamsu", "Morrison", "Papageorgiou",
#     "Tauesch", NA, "Summary"),
#   c("Deaths", "(steroid)", "36", "1",
#     "4", "14", "3", "1",
#     "8", NA, NA),
#   c("Deaths", "(placebo)", "60", "5",
#     "11", "20", "7", "7",
#     "10", NA, NA),
#   c("", "OR", "0.58", "0.16",
#     "0.25", "0.70", "0.35", "0.14",
#     "1.02", NA, "0.53"))
#
# forestplot(tabletext,
#            cochrane_from_rmeta,new_page = TRUE,
#            is.summary=c(TRUE,TRUE,rep(FALSE,8),TRUE),
#            clip=c(0.1,2.5),
#            xlog=TRUE,
#            col=fpColors(box="royalblue",line="darkblue", summary="royalblue"))
#

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
    grA <- which(grouping==levels(grouping)[1])
    grB <- which(grouping==levels(grouping)[2])
    rbind(
      c(makeOdds(data =data, var1 = x, var2 = outcome),S="all"),
      c(makeOdds(data =data[grA,], var1 = x, var2 = outcome),S="Lev1"),
      c(makeOdds(data =data[grB,], var1 = x, var2 = outcome),S="Lev2"),
      NA
    )
  }) %>% do.call(what=rbind) %>% data.frame()
  return(rbind(rep(NA,length(o[1,])),o))
}


#makeDF(rdb,"q_114",grouping,"severity_brown")
makeTableText <- function(makeDF){
  rbind(c("Symptoms","Elicitors","p"),
        makeDF[-c(1),c(1,6,2)])
}

makeForestPlot <- function(data,variables,grouping,outcome){
  df <- makeDF(data,variables,grouping,outcome)
  forestplot(labeltext = makeTableText(df),
             mean = df[,"mean"] %>% unlist,
             lower = df[,"lower"] %>% unlist,
             upper = df[,"upper"] %>% unlist,
             new_page = TRUE,
             is.summary=c(T,rep(F,length(df[,1])-1)),
             #clip=c(0.1,2.5),
             xlog=TRUE,
             col=fpColors(box="royalblue",line="darkblue", summary="royalblue"))
}


FigForestCofactors <- makeForestPlot(rdb,
               testInsectsbinomial %>%
                 filter(section=="cofactors") %>%
                 arrange(pval) %>%
                 select(variableName) %>% pull() %>% {c(.[c(1,4,5,8)], "q_423_beta","q_422_stress","q_410_masto_cur" )},
               grouping,
               "severity_brown")



#### MEasures of Association ############
rdb$q_114_hypotension_collapse_v5
vcd::assocstats(table(rdb$q_114_hypotension_collapse_v5,
                      grouping))
require(vcd)

# Tu skonczyles - zrób assocstats lapply na wielu zmiennych z groupingiem
crammerElicitorBinvsSympt <-
lapply(variableSelectionTab(rdb) %>% filter(section=="symptoms",level==3) %>%
         pull(variableName),function(x){
           assocstats(rdb[,x] %>%
                      table(grouping))[[4]]
         }) %>%
  unlist %>%
  {data.frame(CV=.,
              var=variableSelectionTab(rdb) %>%
                filter(section=="symptoms",level==3) %>%
                pull(variableName))}

crammerElicitorBinvsSympt2 <-
  lapply(variableSelectionTab(rdb) %>% filter(section=="symptoms",level==3) %>%
           pull(variableName),function(x){
             assocstats(rdb[,x] %>%
                          table(testYJ))[[4]]
           }) %>%
  unlist %>%
  {data.frame(CV=.,
              var=variableSelectionTab(rdb) %>%
                filter(section=="symptoms",level==3) %>%
                pull(variableName))}


rdbp <-rdb
rdbp$grouping<- grouping

printx <- function(x){
  ggplot(rdbp,aes(get(x),fill=grouping))+
    geom_bar()
}

lapply(crammerElicitorBinvsSympt$var %>% as.character() %>% .[1:5],
       printx)
crammerElicitorBinvsSympt[4,]

#### Demographics####

demoTab <- cbind(n = rdb$b_sex[rdb$d_elicitor_gr5=="insects"] %>% summary(),
                 Age = rdb$d_age[rdb$d_elicitor_gr5=="insects"] %>%
                   split(rdb$b_sex[rdb$d_elicitor_gr5=="insects"]) %>%
                   lapply(.,function(x){mean(x) %>% signif(3)}),
                 Cardiologic = rdb$q_410_cardio_cur[rdb$d_elicitor_gr5=="insects"]%>%
                   split(rdb$b_sex[rdb$d_elicitor_gr5=="insects"]) %>%
                   f3,
                 DM = rdb$q_410_diab_cur_v6[rdb$d_elicitor_gr5=="insects"]%>%
                   split(rdb$b_sex[rdb$d_elicitor_gr5=="insects"]) %>%
                   f3,
                 `Food allergy` = rdb$q_410_foodallergy_cur_v6[rdb$d_elicitor_gr5=="insects"]%>%
                   split(rdb$b_sex[rdb$d_elicitor_gr5=="insects"]) %>%
                   f3,
                 Mastocytosis = rdb$q_410_masto_cur[rdb$d_elicitor_gr5=="insects"] %>%
                   split(rdb$b_sex[rdb$d_elicitor_gr5=="insects"]) %>% f3,
                 Malignancy = rdb$q_410_malig_cur[rdb$d_elicitor_gr5=="insects"] %>%
                   split(rdb$b_sex[rdb$d_elicitor_gr5=="insects"])  %>% f3,
                 `Atopic dermatitis` = rdb$q_410_ad_cur[rdb$d_elicitor_gr5=="insects"] %>%
                   split(rdb$b_sex[rdb$d_elicitor_gr5=="insects"])  %>% f3,
                 `tryptase [median]` = rdb$q_212_tryptase_value_v5 [rdb$d_elicitor_gr5=="insects"]%>%
                   split(rdb$b_sex[rdb$d_elicitor_gr5=="insects"]) %>%
                   lapply(median,na.rm=T)
)


#### Countries ####

t1 <- table(rdb$d_centres_country,rdb$q_340_insects)
t2 <- t1 %>% prop.table(1) %>% {round(.*100,1)}

rdb %>% group_by(d_centres_country,q_340_insects)%>%
  select(d_centres_country,q_340_insects) %>%
  summarize(n()) %>% tidyr::spread(key = d_centres_country,
                            value = `n()`)


#### Previous ANA ####
rdb$grouping <- grouping
rdb$q_160_ever_react %>% table(rdb$grouping) %>% assocstats()

rdb$q_1622_ever_mild_v4 %>% table(grouping) %>% assocstats()

rdb$q_1621_ever_severe_v4 %>% table(grouping) %>% assocstats()

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
  mutate(difference = ifelse(first_sev <second_sev, "greater",
                             ifelse(first_sev == second_sev, "same","greater")))
greter_same <- iCodeDF %>%
  group_by(difference) %>% summarize(n = n()) %>% pull() %>% as.numeric()

### Check what was the elicitor in the repeated cases


### PLOT Szmptom differences ####



plotSympt <- ggplot(testInsectsbinomialTab$symptoms %>%
                      tidyr::gather(key = "group",value = "Fraction", 4:5) %>%
                      arrange(pval),
                    aes(reorder(variableName,pval),Fraction, fill = group))+
  geom_bar(stat = "identity", position = "dodge",na.rm = T)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 50,hjust =1))

plotSymptKids <- ggplot(makeTests(rdb) %>%
                      tidyr::gather(key = "group",value = "Fraction", 4:5) %>%
                      arrange(pval),
                    aes(reorder(variableName,pval),Fraction, fill = group))+
  geom_bar(stat = "identity", position = "dodge",na.rm = T)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 50,hjust =1))

