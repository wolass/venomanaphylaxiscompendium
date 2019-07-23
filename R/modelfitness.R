modelfiteness <- function(){
  ggplot(temp,aes(.rownames,Pr...t.., color=outcome))+
  geom_segment(x=0.5,xend=4.5,y=log10(0.05),yend=log10(0.05))+
  scale_y_log10()+
  geom_jitter(width = 0.2)+
  theme(axis.text.x = element_text(angle=40, hjust=1))+
  annotate("text",x=4,y = c(0.0001,0.00001,.000001,.0000001,.00000001)*10,
           label = c("n=1776","n=1780","n=1206","n=234","n=1780"),size = 3)+
  labs(x = "",y="p value")
}
