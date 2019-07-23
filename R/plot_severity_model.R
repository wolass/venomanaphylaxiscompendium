plot_severity <- function(){
  output_severity <- lapply(possible_outcomes,
                            function(x){
                              ana_outcome_matched(data = data4,
                                                  outcome = x,
                                                  matching_var = "q_340_insects",
                                                  confounders_cat = c("tryp_cat", "cc_cardio"),
                                                  confounders_int = c("d_age"),
                                                  selected_predictors = c(1,2,4))
                            }
  )


  do_cof_tab <- function(x){
    output_severity[[x]]$fit %>% summary() %>% {.$coefficients}%>% tidy() %>% as.matrix()
  }

  temp <- lapply(1:5,do_cof_tab) %>% do.call(what = rbind) %>%
    as_tibble() %>%
    {mutate(.,
            Estimate = as.numeric(Estimate),
            Pr...t.. = as.numeric(Pr...t..)
    )} %>%
    {cbind(.,outcome = lapply(possible_outcomes,
                              rep,
                              nrow(.)/5) %>% unlist)
    } %>%
    filter(.rownames!="(Intercept)")
  gridExtra::grid.arrange(output_severity[[2]]$plot,
                          output_severity[[3]]$plot,
                          output_severity[[1]]$plot,
                          output_severity[[4]]$plot,
                          output_severity[[5]]$plot,
                          temp %>%
                            ggplot(aes(.rownames,Pr...t.., color=outcome))+
                            geom_segment(x=0.5,xend=4.5,y=log10(0.05),yend=log10(0.05))+
                            scale_y_log10()+
                            geom_jitter()+
                            theme(axis.text.x = element_text(angle=40, hjust=1))
  )
}
