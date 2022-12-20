#### Importing my data into R

library(readxl)
Alldata<-read_excel("FISHMicroscopyResults.xlsx")

library(ggplot2)

library(tidyverse)
library(ggrepel)

#### making a graph for analysis of spatial distribution of bacteria in the gut of the various conditions
###Figure 2B

Alldata%>%pivot_longer(c(BulbBiovolume,PBacBiovolume_um3,DBacBiovolume_um3),
                       names_to = "Location", values_to = "Biovolume")%>%
  mutate(LocationClean=recode(Location,"BulbBiovolume"="Bulb",
                              "DBacBiovolume_um3"="Distal","PBacBiovolume_um3"="Proximal"))%>%
  mutate(LocationClean=factor(LocationClean,levels=c("Bulb","Proximal","Distal")))%>%
ggplot(aes(x=LocationClean, y=Biovolume, fill=Condition)) + 
  theme(legend.position ="none")+
  guides(colour=FALSE)+ 
  geom_jitter(width=0.1, size=2, shape=21)+
  facet_wrap(~Condition, nrow=1)+
  geom_boxplot(alpha=0.6, width=0.5, show.legend = FALSE, 
               position=position_dodge(width=0.5, preserve = "total"))+
  scale_fill_manual(values = c("grey60","#a8ddb5","#0868ac"))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     breaks = c(0,10,100,1000,10000,100000,1000000,10000000,1000000000),
                     expand=c(0,0))+
  theme_classic()+
  theme(legend.position ="none",axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=15),
        strip.text = element_text(size=18),
        axis.text.y = element_text(hjust=1., size=15))+
  labs(x=NULL, y="Biovolume(um^3)", size=10, face= "bold")
ggsave("Figures/Figure2B_BacteriaBiogeaography.png",width = 10,height=5,dpi=400)
ggsave("Figures/Figure2B_BacteriaBiogeaography.pdf",width = 10,height=5)

#### making a graph for the mucus distributon in proximal and distal gut
### Normalized mucus intensity for intestine (proximal and Distal) for various conditions
###Figure 5E
PlotIntestineMucus<-Alldata%>%mutate(IntestineMucus=NormalizedPMucus+NormalizedDMucus) %>% 
ggplot(aes(x=Condition, y=IntestineMucus, label=ID, fill=Condition)) + 
  geom_point(size=3, shape=21, position=position_dodge(width=0.5))+
  geom_boxplot(alpha=0.6, width=0.5, show.legend = FALSE, 
               position=position_dodge(width=0.5, preserve = "total"))+
  scale_fill_manual(values = c("grey60","#a8ddb5","#0868ac"))+
  stat_compare_means(comparisons=list( c("Axenic", "Conv"), c("Conv", "Mix9"), c("Axenic", "Mix9")),
                     label="p.format")+
 # geom_text_repel(aes(group=ID), size=4)+
  theme_classic()+
  theme(legend.position ="none",axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(hjust=1., size=15))+
  labs(x=NULL, y="Normalized Mucus Intensity (A.U.)")
ggsave("Figures/Figure5E_IntestineMucus.png",width = 4,height=5,dpi=400)
ggsave("Figures/Figure5E_IntestineMucus.pdf",width = 4,height=5)

#### Total bacteria biovolume
###Figure1A

ggplot(Alldata, aes(x=Condition, y=TotalBacBiovolume, fill=Condition)) + 
  theme(legend.position ="none")+
  geom_jitter(width=0.1, size=3, shape=21)+
  geom_boxplot(alpha=0.6, width=0.5, show.legend = FALSE, 
               position=position_dodge(width=0.5, preserve = "total"))+
  scale_fill_manual(values = c("grey60","#a8ddb5","#0868ac"))+
  stat_compare_means(comparisons=list(c("Conv", "Mix9")),
                     label="p.format")+
  theme_classic()+
  theme(legend.position ="none",axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(hjust=1., size=15))+
  labs(x=NULL, y="Biovolume (um^3)", size=10, face= "bold")
  
ggsave("Figures/Figure1A_TotalBiovolume.png",width = 4,height=5,dpi=400)
ggsave("Figures/Figure1A_TotalBiovolume.pdf",width = 4,height=5)

#### TotalBacteria and TotalGutMucus association
###Figure5C
PlotFigure5C<-Alldata %>% filter(Condition!="Axenic") %>% 
ggplot(aes(x=Condition, y=BacteriavsMucus, fill=Condition)) + 
  theme(legend.position ="none")+
  guides(colour=FALSE)+ 
  geom_point(size=3, shape=21, position=position_dodge(width=0.5))+
  geom_boxplot(alpha=0.6, width=0.5, show.legend = FALSE, 
               position=position_dodge(width=0.5, preserve = "total"))+
  scale_fill_manual(values = c("#a8ddb5","#0868ac"))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     breaks = c(0,10,100,1000,10000,100000,1000000,10000000,1000000000),
                     expand=c(0,0))+
  theme_classic()+
  theme(legend.position ="none",axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(hjust=1., size=15))+
  labs(x=NULL, y="Total Bacteria/Normalized Mucus", size=9, face= "bold")

ggsave("Figures/Figure5C_NormBacterMucus.png",width = 3,height=5,dpi=400)
ggsave("Figures/Figure5C_NormBacterMucus.pdf",width = 3,height=5)

#### Normalized entire gut mucus analysis for conventional and mix 9
### Figure 5D

PlotFigur5D<-Alldata %>% filter(Condition!="Axenic") %>% 
  ggplot(aes(x=Condition, y=NormalizedGutMucus, fill=Condition)) + 
  theme(legend.position ="none")+
  guides(colour=FALSE)+ 
  geom_jitter(size=3, shape=21, width=0.1)+
  geom_boxplot(alpha=0.6, width=0.5, show.legend = FALSE, outlier.shape = NA)+
  scale_fill_manual(values = c("#a8ddb5","#0868ac"))+
  theme_classic()+
  theme(legend.position ="none",axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(hjust=1., size=15))+
  labs(x=NULL, y="Normalized Mucus Intesity (A.U)", size=9, face= "bold")

ggsave("Figures/Figure5D_GutNormMucus.png",width = 3,height=5,dpi=400)
ggsave("Figures/Figure5D_GutNormMucus.pdf",width = 3,height=5)

### Per Fish Bacteria and Mucus Analysis
### Figure 5A&B

PlotFigure5AB<-Alldata %>% filter(Condition!="Axenic") %>% 
  ggplot(aes(x=TotalBacBiovolume, y=NormalizedGutMucus, fill=Condition, shape=Condition)) + 
  geom_point(size=3)+
  scale_shape_manual(values=c(22,23))+
  scale_fill_manual(values = c("#a8ddb5","#0868ac"))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10),
                     breaks = c(0,10,100,1000,10000,100000,1000000,10000000,1000000000))+
  theme_classic()+
  theme(axis.title = element_text(size=18),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(hjust=1., size=15))+
  labs(x="Bacteria Biovolume (um^3)", y="Normalized Mucus Intesity (A.U)", size=9, face= "bold")

ggsave("Figures/Figure5AB_GutNormMucusBacPerFish.png",width = 5,height=5,dpi=400)
ggsave("Figures/Figure5AB_GutNormMucusBacPerFish.pdf",width = 5,height=5)



library(patchwork)
PlotFigure5AB+PlotFigure5C+PlotFigur5D+PlotIntestineMucus+plot_annotation(tag_levels = "A")
ggsave("Figures/Figure5.png",width = 10,height=10,dpi=400)
