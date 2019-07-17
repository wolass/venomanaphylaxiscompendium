### make all forest plots
make_all_forest_plots <- function(){


test1 <- makeDF(age_sex_matched,
                c( "q_410_cur",
                   "q_410_asthma_cur",
                   "q_410_masto_cur",
                   "q_410_cardio_cur",
                   "q_423_beta",
                   "q_423_ace",
                   "q_4211_exercise",
                   "q_422_stress",
                   "q_423_asa"),
                "grouping",
                "severity_brown")

test1[1] <- c(NA,"Concomitant disease",
              NA,NA,"Asthma",
              NA,NA,"Mastocytosis",
              NA,NA,"Cardiologic disease",
              NA,NA,"β-blockers",
              NA,NA,"ACE-I",
              NA,NA,"Exercise",
              NA,NA,"Stress",
              NA,NA,"ASA",NA,NA)
test1 <- test1[-c(2:4,20:25),]


test2 <- makeDF(age_sex_matched,
                c( "q_410_cur",
                   "q_410_masto_cur",
                   "q_410_asthma_cur",
                   "q_410_cardio_cur",
                   "q_423_beta",
                   "q_423_ace",
                   "q_4211_exercise",
                   "q_422_stress",
                   "q_423_asa"),
                "grouping",
                "d_severity_rmr")

test2[1] <- c(NA,"Concomitant disease",
              NA,NA,"Mastocytosis",
              NA,NA,"Asthma",
              NA,NA,"Cardiovascular dis.",
              NA,NA,"β-blockers",
              NA,NA,"ACE-I",
              NA,NA,"Exercise",
              NA,NA,"Stress",
              NA,NA,"ASA",NA,NA)
test2 <- test2[-c(2:4,20:25),]


#or_plot(test1)


png("analysis/figures/figForestfinalb.png",
    height = 800,
    width = 1000,
    #res = 300,
    pointsize = 40,
    units="px")
forestplot(labeltext = makeTableText(test1),
           mean = test1[,"mean"] %>% unlist,
           lower = test1[,"lower"] %>% unlist,
           upper = test1[,"upper"] %>% unlist,
           new_page = TRUE,
           is.summary=c(T,rep(F,length(test1[,1])-1)),
           clip=c(0.2,4),
           xlog=TRUE,
           col=fpColors(box="#1f1f1f",line="gray", summary="royalblue"),
           lwd.xaxis = 5,
           lwd.ci = 4,
           lwd.zero = 4)
dev.off()

png("analysis/figures/figForestfinalrmr.png",
    height = 640,
    width = 860,
    #res = 300,
    pointsize = 30
    #units="px"
)
forestplot(labeltext = makeTableText(test2)[,1:2],
           mean = test2[,"mean"] %>% unlist,
           lower = test2[,"lower"] %>% unlist,
           upper = test2[,"upper"] %>% unlist,
           new_page = TRUE,
           is.summary=c(T,rep(F,length(test2[,1])-1)),
           clip=c(0.2,4),
           xlog=TRUE,
           col=fpColors(box="#1f1f1f",line="gray", summary="royalblue"),
           lwd.xaxis = 5,
           lwd.ci = 4,
           lwd.zero = 4,
           hrzl_lines = T,
           grid = T,
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.8)))
dev.off()




png("analysis/figures/figForest.png")
makeForestPlot(age_sex_matched,
               c( "q_410_masto_cur",
                  "q_410_asthma_cur",
                  "q_410_cardio_cur",
                  "q_423_beta",
                  "q_423_ace",
                  "q_423_asa"),
               "grouping",
               "d_severity_rmr")
dev.off()

png("analysis/figures/kidsForest.png")
makeForestPlot(rdbp[rdbp$d_age<18,],
               testInsectsbinomial %>%
                 filter(section=="cofactors") %>%
                 arrange(pval) %>%
                 select(variableName) %>% pull() %>% {c(.[c(1,2,6,7,15,18
                 )])},
               "grouping",
               "severity_brown",
               cLow= 0.3,
               cHigh = 4)
dev.off()
png("analysis/figures/adultsForest.png")
makeForestPlot(rdbp[rdbp$d_age>18,],
               testInsectsbinomial %>%
                 filter(section=="cofactors") %>%
                 arrange(pval) %>%
                 select(variableName) %>% pull() %>% {c(.[c(1,6,16,41,18,27
                 )])},
               "grouping",
               "severity_brown",
               cLow = 0.5,
               cHigh = 2)
dev.off()
}
