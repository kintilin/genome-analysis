time <- read.table("ltr_insertime_out.txt")
library(ggplot2)
den <- density(time$V2)
denx=den$x
deny=den$y
#p <- 
ggplot()+
   geom_density(mapping=aes(x=time$V2),fill="dodgerblue2",alpha=0.5)+
   geom_vline(xintercept = denx[which.max(deny)])+
   geom_segment(mapping=aes(xend=denx[which.max(deny)],yend=0.145,x=5,y=0.16),arrow = arrow(length=unit(0.2, "cm")))+
   annotate("text", x = 10, y = 0.16,label = paste(round(denx[which.max(deny)],4),"MYA"),size=4)+
   labs(x="LTR insertion time (Million Year Ago)",y="Density")+
   #xlim(c(0,15))+
   theme_classic()

ggsave("LTR_insert_time.pdf",p,width = 210*0.4 , height = 210*0.4,units="mm")
