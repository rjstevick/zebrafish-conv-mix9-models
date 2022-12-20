#### Importing my data into R

library(readxl)
Manuscript_Bacteria_Data_F <- read_excel("Manuscript_Bacteria_Data_F.xlsx")
View(Manuscript_Bacteria_Data_F)                                                        
Manuscript_Mucus_Data_F <- read_excel("Manuscript_Mucus_Data_F.xlsx")
View(Manuscript_Mucus_Data_F) 
Manuscript_Mucus_Sum_F <- read_excel("Manuscript_Mucus_Sum_F.xlsx")
View(Manuscript_Mucus_Sum_F) 
Manuscript_Bacteria_Biovolume_Sum_F <- read_excel("Manuscript_Bacteria_Biovolume_Sum_F.xlsx")
View(Manuscript_Bacteria_Biovolume_Sum_F)   

library(ggplot2)

library(tidyverse)
library(ggrepel)

#### making a graph for analysis of bacteria in the gut of the various condiitions

ggplot(Manuscript_Bacteria_Data_F, aes(x=Location, y=Biovolume, fill=Condition)) + 
  theme(legend.position ="none")+
  guides(colour=FALSE)+ 
  geom_jitter(width=0.1, size=2, shape=21)+
  facet_wrap(~Condition, nrow=1)+
  geom_boxplot(alpha=0.6, width=0.5, show.legend = FALSE, 
               position=position_dodge(width=0.5, preserve = "total"))+
  scale_x_discrete(limits = c("Bulb", "Proximal", "Distal"))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     breaks = c(0,10,100,1000,10000,100000,1000000,10000000,1000000000),
                     expand=c(0,0))+
  theme_classic()+
  theme(strip.text.x = element_text(size=9, face="bold"))+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle=50, hjust=1., size=10, face="bold"))+
  theme(axis.text.y = element_text(hjust=1., size=10, face="bold"))+
  labs(x=NULL, y="Biovolume(um^3)", size=10, face= "bold")

#### making a graph for the mucus distributon in proximal and distal gut

ggplot(Manuscript_Mucus_Data_F, aes(x=Location, y=Intensity_SurfaceArea, label=ID, fill=Condition)) + 
  theme(legend.position ="none")+
  guides(colour=FALSE)+ 
  geom_point(size=3, shape=21, position=position_dodge(width=0.5))+
  facet_wrap(~Condition, nrow=1)+
  geom_boxplot(alpha=0.6, width=0.5, show.legend = FALSE, 
               position=position_dodge(width=0.5, preserve = "total"))+
  geom_text_repel(aes(group=ID), size=4)+
  scale_x_discrete(limits = c("Proximal", "Distal"))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     breaks = c(0,10,100,1000,10000,100000,1000000,10000000,1000000000),
                     expand=c(0,0))+
  theme_classic()+
  theme(strip.text.x = element_text(size=9, face="bold"))+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle=50, hjust=1., size=10, face="bold"))+
  theme(axis.text.y = element_text(hjust=1., size=10, face="bold"))+
  labs(x=NULL, y="Normalized Mucus Intensity (A.U.)", size=10, face= "bold")

#### Graph for mucus sum

ggplot(Manuscript_Mucus_Sum_F, aes(x=Condition, y=Intensity_SurfaceArea, label=ID, fill=Condition)) + 
  geom_point(size=3, shape=21, position=position_dodge(width=0.5))+
  geom_boxplot(alpha=0.6, width=0.5, show.legend = FALSE, 
               position=position_dodge(width=0.5, preserve = "total"))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     breaks = c(0,10,100,1000,10000,100000,1000000,10000000,1000000000),
                     expand=c(0,0))+
  geom_text_repel(aes(group=ID), size=4)+
  theme_classic()+
  theme(legend.position ="none")+
  theme(axis.text.x = element_text(angle=50, hjust=1., size=10, face="bold"))+
  theme(axis.text.y = element_text(hjust=1., size=10, face="bold"))+
  labs(x=NULL, y="Normalized Mucus Intensity (A.U.)", size=10, face= "bold")

#### Graph for bacteria sum

ggplot(Totalbacteria_Volume, aes(x=Condition, y=TotalBacBiovolume, fill=Condition)) + 
  theme(legend.position ="none")+
  geom_jitter(width=0.1, size=3, shape=21)+
  geom_boxplot(alpha=0.6, width=0.5, show.legend = FALSE, 
               position=position_dodge(width=0.5, preserve = "total"))+
  theme_classic()+
  theme(legend.position ="none")+
  theme(axis.text.x = element_text(angle=50, hjust=1., size=10, face="bold"))+
  theme(axis.text.y = element_text(hjust=1., size=10, face="bold"))+
  labs(x=NULL, y="Biovolume (um^3)", size=10, face= "bold")
  
Totalbacteria_Volume <- read_excel("Totalbacteria Volume.xlsx")
#### Bacteria and Mucus association

ggplot(Manuscript_Mucus_Bacteria_Sum_F, aes(x=Type, y=Intensity_SurfaceArea, label=ID, fill=Condition)) + 
  theme(legend.position ="none")+
  guides(colour=FALSE)+ 
  geom_point(size=3, shape=21, position=position_dodge(width=0.5))+
  facet_wrap(~Condition, nrow=1)+
  geom_boxplot(alpha=0.6, width=0.5, show.legend = FALSE, 
               position=position_dodge(width=0.5, preserve = "total"))+
  geom_text_repel(aes(group=ID), size=2)+
  geom_line(aes(group=ID))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     breaks = c(0,10,100,1000,10000,100000,1000000,10000000,1000000000),
                     expand=c(0,0))+
  theme_classic()+
  theme(strip.text.x = element_text(size=9, face="bold"))+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle=50, hjust=1., size=8, face="bold"))+
  theme(axis.text.y = element_text(hjust=1., size=8, face="bold"))+
  labs(x=NULL, y=NULL, size=9, face= "bold")

Manuscript_Mucus_Bacteria_Sum_F <- read_excel("Manuscript_Mucus_Bacteria_Sum_F.xlsx")


#### Mucus and Bacteria - Proximal and Distal

ggplot(Manuscript_Mucus_Bacteria_Data_New_F, aes(x=Location, y=Intensity_SurfaceArea, label=ID, fill=Condition)) + 
  theme(legend.position ="none")+
  guides(colour=FALSE)+ 
  geom_point(size=3, shape=21, position=position_dodge(width=0.5))+
  facet_wrap(~Condition, nrow=1)+
  geom_boxplot(alpha=0.6, width=0.5, show.legend = FALSE, 
               position=position_dodge(width=0.5, preserve = "total"))+
  scale_x_discrete(limits = c("PB", "PM","DB","DM"))+
  geom_text_repel(aes(group=ID), size=4)+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     breaks = c(0,10,100,1000,10000,100000,1000000,10000000,1000000000),
                     expand=c(0,0))+
  theme_classic()+
  theme(strip.text.x = element_text(size=9, face="bold"))+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle=50, hjust=1., size=10, face="bold"))+
  theme(axis.text.y = element_text(hjust=1., size=10, face="bold"))+
  labs(x=NULL, y=NULL, size=10, face= "bold")

####Bacteria and Mucus Sum (Only proximal and Distal)

ggplot(Manuscript_Mucus_Bacteria_Data_New_today, aes(x=Location, y=Intensity, label=ID, fill=Condition)) + 
  theme(legend.position ="none")+
  guides(colour=FALSE)+ 
  geom_point(size=3, shape=21, position=position_dodge(width=0.5))+
  facet_wrap(~Condition, nrow=1)+
  geom_boxplot(alpha=0.6, width=0.5, show.legend = FALSE, 
               position=position_dodge(width=0.5, preserve = "total"))+
  geom_text_repel(aes(group=ID), size=4)+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     breaks = c(0,10,100,1000,10000,100000,1000000,10000000,1000000000),
                     expand=c(0,0))+
  theme_classic()+
  theme(strip.text.x = element_text(size=10, face="bold"))+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle=50, hjust=1., size=10, face="bold"))+
  theme(axis.text.y = element_text(hjust=1., size=10, face="bold"))+
  labs(x=NULL, y=NULL, size=9, face= "bold")

#### Full gut analysis for conventional and mix 9

ggplot(Manuscript_Mucus_Bacteria_Sum_F_AV,aes(x=ID, y = Intensity_SurfaceArea, fill=Type), group =Type) + 
  geom_col(position="dodge") + 
  facet_wrap(~Condition, ncol=1)+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     breaks = c(0,10,100,1000,10000,100000,1000000,10000000,1000000000),
                     expand=c(0,0), sec.axis = sec_axis(~ . *1, name = "Biovolume(um^3)"))+
  labs(y="Normalized Mucus Intensity (A.U.")

ggplot(Manuscript_Mucus_Bacteria_Sum_F_AV,aes(x=ID, y = Intensity_SurfaceArea, fill=Condition) + 
  geom_col(position="dodge") + 
  facet_wrap(~Type, ncol=1)+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     breaks = c(0,10,100,1000,10000,100000,1000000,10000000,1000000000),
                     expand=c(0,0), sec.axis = sec_axis(~ . *1, name = "Biovolume(um^3)"))+
  labs(y="Normalized Mucus Intensity (A.U.")

  ggplot(Bacteria_Mucus, aes(x=Condition, y=Intensity, fill=Condition, fill=condition)) + 
    theme(legend.position ="none")+
    geom_point(size=3, shape=21, position=position_dodge(width=0.5), show.legend = "False")+
    geom_boxplot(alpha=0.8, width=0.3, fill="grey20", position=position_dodge(width=0.9)+
    scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                       breaks = c(0,10,100,1000,10000,100000,1000000,10000000,1000000000,1000000000),
                       expand=c(0,0))+
    theme_classic()+
    theme(legend.position ="none")+
    theme(axis.text.x = element_text(angle=50, hjust=1., size=12, face="bold"))+
    theme(axis.text.y = element_text(hjust=1., size=12, face="bold"))+
    labs(x=NULL, y="Total bacteria biovolume/Total mucus", size=17, face= "bold")
  
  ggplot(Conventional_Biovol_Mucus, aes(x = Biovolume, y = NormalizedMucus)) +
    geom_point() 
  
  ggplot(Mix9_Biovol_Mucus, aes(x = Biovolume, y = MucusIntensity)) +
    geom_point() 
  
  ggplot(Bacteria_Mucus, aes(x=Condition, y=Intensity)) + 
    theme(legend.position ="none")+
    geom_point(size=3, shape=21, position=position_dodge(width=0.5), show.legend = "False")+
    geom_boxplot(alpha=0.8, width=0.3, fill="grey20", position=position_dodge(width=0.9)+
                   scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                                      breaks = c(0,10,100,1000,10000,100000,1000000,10000000,1000000000,1000000000),
                                      expand=c(0,0))+
                   theme_classic()+
                   theme(legend.position ="none")+
                   theme(axis.text.x = element_text(angle=50, hjust=1., size=12, face="bold"))+
                   theme(axis.text.y = element_text(hjust=1., size=12, face="bold"))+
                   labs(x=NULL, y="Total bacteria biovolume/Total mucus", size=17, face= "bold")

                 ggplot(scale_remaned, aes(x=conditional, y=log_scale))+
                   theme(legend.position ="none")+
                   geom_point(size=3, shape=21, position=position_dodge(width=0.5), show.legend = "False")+
                   geom_boxplot(alpha=0.8, width=0.3, fill="black", position=position_dodge(width=0.9)+                                theme_classic()+
                                  theme(legend.position ="none")+
                                  theme(axis.text.x = element_text(angle=50, hjust=1., size=12, face="bold"))+
                                  theme(axis.text.y = element_text(hjust=1., size=12, face="bold"))+
                                  labs(x=NULL, y="Normalized", size=17, face= "bold")
                                
                                
                                ggplot(scale_remaned, aes(x=conditional, y=log_scale)) + 
                                  theme(legend.position ="none")+
                                  geom_jitter(width=0.1, size=3)+
                                  geom_boxplot(alpha=0.6, width=0.5, show.legend = FALSE, 
                                               position=position_dodge(width=0.5, preserve = "total"))+
                                  theme_classic()+
                                  theme(legend.position ="none")+
                                  theme(axis.text.x = element_text(angle=50, hjust=1., size=10, face="bold"))+
                                  theme(axis.text.y = element_text(hjust=1., size=10, face="bold"))+
                                  labs(x=NULL, y="Biovolume (um^3)", size=10, face= "bold")
                                
                                ggplot(TotalIntestineMucus, aes(x=Condition, y=Mucus, fill=Condition)) + 
                                  geom_jitter(width=0.1, size=3, shape=21)+
                                  geom_boxplot(alpha=0.6, width=0.5, show.legend = FALSE, 
                                               position=position_dodge(width=0.5, preserve = "total"))+
                                  theme_classic()+
                                  theme(legend.position ="none")+
                                  theme(axis.text.x = element_text(angle=50, hjust=1., size=10, face="bold"))+
                                  theme(axis.text.y = element_text(hjust=1., size=10, face="bold"))+
                                  labs(x=NULL, y="Normalized Mucus Intensity (A.U.)", size=10, face= "bold")
                                  
                                  
        t.test(Bacteria_Mucus, mu = 0, alternative = "two.sided")
        
        ggplot(Mix9_variabilty, aes(y=Biovolume))+
          geom_histogram()+
          geom_jitter(aes(x=0.5), width=0.4, height=0, shape=21, color="black")+
          scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                             breaks = c(0,10,100,1000,10000,100000,1000000))+
          theme_pubr()+
          scale_color_manual(values=c("lightblue", "red3"))+
          scale_fill_manual(values=c("lightblue", "red3"))+
          theme(legend.position="top", text=element_text(color="white"),
                plot.background = element_blank(),panel.background = element_blank(),
                axis.text=element_text(color="white"),
                axis.title = element_text(face="bold"),
                legend.background = element_blank(),
                axis.ticks = element_line(color="white"), axis.line = element_line(color="white"),
                panel.grid.major.y = element_line(inherit.blank = FALSE, color="white"),
                panel.grid.minor.y = element_line(inherit.blank = FALSE, size = 0.1, color="white"))+
          labs(x="Number of fish", y="Bacteria biovolume (Âµm^3)", fill=NULL, color=NULL) 
        
        ggplot(Mix9_variabilty, aes(y=Biovolume_um3))+
          geom_jitter(color="lightblue", size=3)+
          scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                             breaks = c(0,10,100,1000,10000,100000,1000000))
          theme_pubr()+
          theme(legend.position="top", text=element_text(color="white"),
                plot.background = element_blank(),panel.background = element_blank(),
                axis.text=element_text(color="white"),
                axis.title = element_text(face="bold"),
                legend.background = element_blank(),
                axis.ticks = element_line(color="white"), axis.line = element_line(color="white"),
                axis.text.x = element_blank(),
                panel.grid.major.y = element_line(inherit.blank = FALSE, color="white"),
                panel.grid.minor.y = element_line(inherit.blank =FALSE, color="white")
        
                