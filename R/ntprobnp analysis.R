require(tidyverse)
df <-read_csv("ntprobnp_results.csv")
glimpse(df)
df$ntprobnp[28] <- NA
df$ntprobnp <- as.numeric(df$ntprobnp)
ggpubr::ggscatter(df,
       x = "d_severity_rm",
       y = "ntprobnp")
ggpubr::ggscatter(df,
                  x = "q_114_hypotension_collapse_v5",
                  y = "ntprobnp")
ggpubr::ggscatter(df,
                  x = "q_212_tryptase_value_v5",
                  y = "ntprobnp")
load("data.R")
sum(data$b_patient_code %in% df$b_patient_code)

df$b_patient_code
sum(data$b_patient_code %>% as.character() %>% gsub(pattern=" ",replacement="") %in% df$b_patient_code)

fdf <- data %>%
  mutate(b_patient_code = data$b_patient_code %>% as.character() %>% gsub(pattern=" ",replacement="")) %>%
  filter(b_patient_code %in% df$b_patient_code) %>%
  select(b_patient_code,d_age,q_410_ad_cur,q_410_cardio_cur,q_114_loss_of_consciousness,q_111_urticaria,ANAscore,) %>%
  full_join(df,by="b_patient_code")

ggpubr::ggscatter(fdf,
                  x = "d_age",
                  y = "ntprobnp")

ggpubr::ggscatter(fdf,
                  x = "d_age",
                  y = "ntprobnp",
                  color = "q_410_cardio_cur")

ggpubr::ggscatter(fdf,
                  x = "q_212_tryptase_value_v5",
                  y = "ntprobnp",
                  color = "q_410_cardio_cur")


ggpubr::ggscatter(fdf,
                  x = "d_age",
                  y = "ntprobnp",
                  color = "q_114_hypotension_collapse_v5",
                  facet.by = "d_severity_rm")


ggpubr::ggscatter(fdf,
                  x = "d_age",
                  y = "ntprobnp",
                  color = "q_114_loss_of_consciousness",
                  facet.by = "d_severity_rm")

ggpubr::ggscatter(fdf,
                  x = "q_212_tryptase_value_v5",
                  y = "ntprobnp",
                  color = "q_111_urticaria",
                  facet.by = "d_severity_rm")

ggpubr::ggscatter(fdf,
                  x = "ANAscore",
                  y = "ntprobnp",
                  color = "q_410_cardio_cur")

