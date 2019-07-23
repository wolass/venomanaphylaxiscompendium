plot_severity_tryptase <- function(){
  ggarrange(
    rdb %>%
      filter(d_elicitor_gr5=="insects") %>%
      mutate(typ_cat=ifelse(q_212_tryptase_value_v5<8,
                            "low",
                            "high")) %>%
      filter(!is.na(typ_cat)) %>%
      ggplot(aes(typ_cat,q_116_VAS_v7))+
      geom_boxplot(),
    rdb %>%
      filter(d_elicitor_gr5=="insects") %>%
      mutate(typ_cat=ifelse(q_212_tryptase_value_v5<8,
                            "low",
                            "high")) %>%
      filter(!is.na(typ_cat)) %>%
      ggplot(aes(typ_cat,ANAscore))+
      geom_boxplot(),
    rdb %>%
      filter(d_elicitor_gr5=="insects",
             !is.na(d_severity_rm)) %>%
      mutate(typ_cat=ifelse(q_212_tryptase_value_v5<8,
                            "low",
                            "high"),
             d_severity_rm = factor(d_severity_rm)) %>%
      group_by(d_severity_rm,typ_cat) %>%
      filter(!is.na(typ_cat)) %>%
      summarize(n = n()) %>%
      ggplot(aes(typ_cat,n,fill=d_severity_rm))+
      geom_bar(stat="identity",position= "fill"),
    rdb %>%
      filter(d_elicitor_gr5!="insects") %>%
      mutate(typ_cat=ifelse(q_212_tryptase_value_v5<8,
                            "low",
                            "high")) %>%
      filter(!is.na(typ_cat)) %>%
      ggplot(aes(typ_cat,q_116_VAS_v7))+
      geom_boxplot(),
    rdb %>%
      filter(d_elicitor_gr5!="insects") %>%
      mutate(typ_cat=ifelse(q_212_tryptase_value_v5<8,
                            "low",
                            "high")) %>%
      filter(!is.na(typ_cat)) %>%
      ggplot(aes(typ_cat,ANAscore))+
      geom_boxplot(),
    rdb %>%
      filter(d_elicitor_gr5!="insects",
             !is.na(d_severity_rm)) %>%
      mutate(typ_cat=ifelse(q_212_tryptase_value_v5<8,
                            "low",
                            "high"),
             d_severity_rm = factor(d_severity_rm)) %>%
      group_by(d_severity_rm,typ_cat) %>%
      filter(!is.na(typ_cat)) %>%
      summarize(n = n()) %>%
      ggplot(aes(typ_cat,n,fill=d_severity_rm))+
      geom_bar(stat="identity",position= "fill"),
    ncol =3 , nrow = 2,
    widths = c(1,1.1,2))
}
