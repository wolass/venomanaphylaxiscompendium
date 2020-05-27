# age_sex_matched %>%
#   mutate(grouping = car::recode(grouping, "'insects' = 'VIA';'other'='non-VIA'")) %>%
#   group_by(grouping,q_114,d_522_adren_agg) %>%
#   count() %>%
#   filter(!is.na(q_114),
#          !is.na(d_522_adren_agg),
#          q_114!="unknown",
#          d_522_adren_agg!="unknown")%>%
# ggplot(aes(x = q_114,y = n, fill = d_522_adren_agg))+
#   facet_grid(~grouping)+
#   geom_bar(stat = "identity") +
#   labs(x = "Cardiologic symptoms",
#        y = "Number of cases",
#        fill = "Epinephrine Tx")
