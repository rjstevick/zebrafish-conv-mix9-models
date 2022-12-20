
# load libraries
library(readxl)
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(patchwork)

# set global theme
theme_set(theme_classic()+
             theme(legend.position ="none",
                   axis.title = element_text(size=18),
                   legend.title = element_text(size=18, face="bold"),
                   axis.text = element_text(size=15), 
                   legend.text = element_text(size=18),
                   strip.text = element_text(size=18), 
                   plot.title = element_text(hjust=0.5, size=22),
                   strip.background = element_rect(fill="grey80", color="transparent")))

# set seed so jitters will all be the same
set.seed(2022)

# import data
Alldata<-read_excel("FISHMicroscopyResults.xlsx")

# -----------------------------------------

### Total bacteria biovolume
### Figure 1A

ggplot(Alldata, aes(x=Condition, y=TotalBacBiovolume, shape=Condition, fill=Condition)) +
   geom_jitter(width=0.1, size=3)+
   geom_boxplot(alpha=0.6, width=0.5, outlier.shape=NA)+
   scale_fill_manual(values = c("grey60","#a8ddb5","#0868ac"))+
   scale_shape_manual(values=c(21,22,23))+
   stat_compare_means(comparisons=list(c("Conv", "Mix9")), label="p.format")+
   labs(x=NULL, y="Biovolume (um^3)")

ggsave("Figures/Figure1A_TotalBiovolume.png",width = 4,height=5,dpi=400)
ggsave("Figures/Figure1A_TotalBiovolume.pdf",width = 4,height=5)

# -----------------------------------------


### analysis of spatial distribution of bacteria in the gut of the various conditions
### Figure 2B

Alldata %>%
   # convert to long form
   pivot_longer(c(BulbBiovolume,PBacBiovolume_um3,DBacBiovolume_um3),
                names_to = "Location", values_to = "Biovolume")%>%
   # rename and reorder
   mutate(LocationClean=recode(Location,"BulbBiovolume"="Bulb",
                               "DBacBiovolume_um3"="Distal","PBacBiovolume_um3"="Proximal"),
          LocationClean=factor(LocationClean,levels=c("Bulb","Proximal","Distal")))%>%
   # start plotting
   ggplot(aes(x=LocationClean, y=Biovolume, fill=Condition, shape=Condition)) +
   geom_jitter(size=3, width=0.1)+
   geom_boxplot(alpha=0.6, width=0.5, outlier.shape = NA)+
   facet_wrap(~Condition, nrow=1)+
   scale_shape_manual(values=c(21,22,23))+
   scale_fill_manual(values = c("grey60","#a8ddb5","#0868ac"))+
   scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                      breaks = c(0,10,100,1000,10000,100000,1000000,10000000,1000000000))+
   labs(x=NULL, y="Biovolume(um^3)")

ggsave("Figures/Figure2B_BacteriaBiogeography.png",width = 10,height=5,dpi=400)
ggsave("Figures/Figure2B_BacteriaBiogeography.pdf",width = 10,height=5)

# -----------------------------------------


### Per Fish Bacteria and Mucus Analysis
### Figure 5A&B

PlotFigure5AB<-Alldata %>% filter(Condition!="Axenic") %>%
   ggplot(aes(x=TotalBacBiovolume, y=NormalizedGutMucus, fill=Condition, shape=Condition)) +
   geom_point(size=5)+
   scale_shape_manual(values=c(22,23))+
   scale_fill_manual(values = c("#a8ddb5","#0868ac"))+
   scale_x_continuous(trans=scales::pseudo_log_trans(base = 10),
                      breaks = c(0,10,100,1000,10000,100000,1000000,10000000,1000000000))+
   theme(legend.position="right")+
   labs(x="Bacteria Biovolume (um^3)", y="Normalized Mucus Intesity (A.U)", fill="Fish Model", shape="Fish Model")
PlotFigure5AB
ggsave("Figures/Figure5AB_GutNormMucusBacPerFish.png",width = 5,height=5,dpi=400)
ggsave("Figures/Figure5AB_GutNormMucusBacPerFish.pdf",width = 5,height=5)

# -----------------------------------------


### TotalBacteria and TotalGutMucus association
### Figure 5C
PlotFigure5C<-Alldata %>% filter(Condition!="Axenic") %>%
   ggplot(aes(x=Condition, y=BacteriavsMucus, shape=Condition, fill=Condition)) +
   geom_jitter(size=3, width=0.1)+
   geom_boxplot(alpha=0.6, width=0.5, outlier.shape = NA)+
   scale_shape_manual(values=c(22,23))+
   scale_fill_manual(values = c("#a8ddb5","#0868ac"))+
   scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                      breaks = c(0,10,100,1000,10000,100000,1000000,10000000,1000000000))+
   labs(x=NULL, y="Total Bacteria/Normalized Mucus", title="Entire gut")
PlotFigure5C

ggsave("Figures/Figure5C_NormBacterMucus.png",width = 3,height=5,dpi=400)
ggsave("Figures/Figure5C_NormBacterMucus.pdf",width = 3,height=5)

# -----------------------------------------


### Normalized entire gut mucus analysis for conventional and mix 9
### Figure 5D

PlotFigure5D <- Alldata %>%
   filter(Condition!="Axenic") %>%
   ggplot(aes(x=Condition, y=NormalizedGutMucus, shape=Condition, fill=Condition)) +
   geom_jitter(size=3, width=0.1)+
   geom_boxplot(alpha=0.6, width=0.5, outlier.shape = NA)+
   scale_fill_manual(values = c("#a8ddb5","#0868ac"))+
   scale_shape_manual(values=c(22,23))+
   stat_compare_means(comparisons=list(c("Conv", "Mix9")),
                      label="p.format")+
   labs(x=NULL, y="Normalized Mucus Intesity (A.U)", title="Entire gut")
PlotFigure5D

ggsave("Figures/Figure5D_GutNormMucus.png",width = 3,height=5,dpi=400)
ggsave("Figures/Figure5D_GutNormMucus.pdf",width = 3,height=5)

# -----------------------------------------


### mucus distribution in proximal and distal gut
### Normalized mucus intensity for intestine (proximal and distal) for various conditions
### Figure 5E

PlotFigure5E<-Alldata %>% 
   mutate(IntestineMucus=NormalizedPMucus+NormalizedDMucus) %>%
   ggplot(aes(x=Condition, y=IntestineMucus, shape=Condition, fill=Condition)) +
   geom_jitter(size=3,  width=0.1)+
   geom_boxplot(alpha=0.6, width=0.5, outlier.shape = NA)+
   scale_shape_manual(values=c(21,22,23))+
   scale_fill_manual(values = c("grey60","#a8ddb5","#0868ac"))+
   stat_compare_means(comparisons=list( c("Axenic", "Conv"), c("Conv", "Mix9"), c("Axenic", "Mix9")),
                      label="p.format")+
   labs(x=NULL, y="Normalized Mucus Intensity (A.U)", title="Intestine (excluding bulb)")
PlotFigure5E

ggsave("Figures/Figure5E_IntestineMucus.png",width = 4,height=5,dpi=400)
ggsave("Figures/Figure5E_IntestineMucus.pdf",width = 4,height=5)

# -----------------------------------------


### plot figure 5 all together

((PlotFigure5AB+guide_area()+ PlotFigure5C+ plot_layout(guides="collect"))/(PlotFigure5D+PlotFigure5E))+
   plot_annotation(tag_levels = "A")  & 
   theme(plot.tag = element_text(size = 22))
ggsave("Figures/Figure5_1.png",width = 12,height=10,dpi=400)

(PlotFigure5AB+guide_area()+plot_layout(guides="collect"))/
   (PlotFigure5C+PlotFigure5D+PlotFigure5E)+plot_annotation(tag_levels = "A")& 
   theme(plot.tag = element_text(size = 22))
ggsave("Figures/Figure5_2.png",width = 13,height=10,dpi=400)

