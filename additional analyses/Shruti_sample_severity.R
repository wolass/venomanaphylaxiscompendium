#Shruti serum severity anlysis
d1 <- readxl::read_excel("../../Downloads/Kopie von Seren_ANA_Insekten ab 172_für Shruti.xls")

res <- data %>%
  filter(data$b_patient_code %in% c("w-1969-l-09-b-05-26 "))

data %>%
  select(b_patient_code) %>%
  filter(data$b_patient_code %>% grepl(pattern = "m-1986-t"))
data$b_patient_code[1:6] %>% as.character()
data$b_patient_code %<>% gsub(pattern = "[ ]+",replacement = "")

res <- data %>%
  filter(b_patient_code %in% (d1 %>%
                                     filter(!is.na(`AR-ID`)) %>%
                                     select(`AR-ID`) %>% pull())) %>%
  select(b_patient_code,d_severity_rm,ANAscore)
res %>%
  write.csv(file = "file_172_severity.csv")



check<-  (d1 %>%
         filter(!is.na(`AR-ID`)) %>%
         select(`AR-ID`) %>% pull())
check[which (!(check %in% data$b_patient_code))]


##### Second file
d2 <- readxl::read_excel("../../Downloads/Seren_ANA_Insekten_03.01.2019 für Shruti.xlsx",skip = 1)

res <- data %>%
  filter(b_patient_code %in% (d2 %>%
                                filter(!is.na(`AR-ID`),
                                       `ANA-Code`>98) %>%
                                select(`AR-ID`) %>% pull())) %>%
  select(b_patient_code,d_severity_rm,ANAscore)
res %>%
  write.csv(file = "file_SEREN_ANA.csv")


