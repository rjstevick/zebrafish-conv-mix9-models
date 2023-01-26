
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
                   panel.background = element_blank(), plot.background = element_blank(),
                   axis.title = element_text(size=18),
                   legend.title = element_text(size=18, face="bold"),
                   axis.text = element_text(size=15),
                   legend.text = element_text(size=18),
                   strip.text = element_text(size=22),
                   plot.title = element_text(hjust=0.5, size=22),
                   strip.background = element_rect(fill="grey80", color="transparent")))

# set seed so jitters will all be the same
set.seed(2022)

# import data
Alldata<-read_excel("FISHMicroscopyResults.xlsx")

# -----------------------------------------

### Total bacteria biovolume
### Figure 1E

Plot1E <- ggplot(Alldata, aes(x=Condition, y=TotalBacBiovolume, shape=Condition, fill=Condition)) +
   geom_jitter(width=0.1, size=3)+
   geom_boxplot(alpha=0.6, outlier.shape=NA)+
   scale_fill_manual(values = c("grey60","#a8ddb5","#0868ac"))+
   scale_shape_manual(values=c(21,22,23))+
   scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                      breaks = c(0,100,10000,1000000),
                      labels = c(0,expression(10^2),expression(10^4),expression(10^6)),
                      limits=c(NA,4000000))+
   stat_compare_means(comparisons=list(c("Conv", "Mix9")), label="p.format")+
   labs(x=NULL, y=expression(Bacteria~biovolume~'('*µm^3*')'))
Plot1E
ggsave("Figures/Figure1E_TotalBiovolume.png",width = 4,height=3,dpi=400)
ggsave("Figures/Figure1E_TotalBiovolume.pdf",width = 4,height=3)

# -----------------------------------------


### analysis of spatial distribution of bacteria in the gut of the various conditions
### Figure 2B

Alldata %>%
   # convert to long form
   pivot_longer(c(BulbBiovolume,PBacBiovolume_um3,DBacBiovolume_um3),
                names_to = "Location", values_to = "Biovolume") %>%
   # rename and reorder
   mutate(LocationClean=recode(Location,
                               "BulbBiovolume"="Bulb",
                               "DBacBiovolume_um3"="Distal","PBacBiovolume_um3"="Proximal"),
          LocationClean=factor(LocationClean,levels=c("Bulb","Proximal","Distal"))) %>%
   # start plotting
   ggplot(aes(x=LocationClean, y=Biovolume, fill=Condition, shape=Condition)) +
   geom_jitter(size=3, width=0.1)+
   geom_boxplot(alpha=0.6, outlier.shape = NA)+
   facet_wrap(~Condition, nrow=1)+
   scale_shape_manual(values=c(21,22,23))+
   scale_fill_manual(values = c("grey60","#a8ddb5","#0868ac"))+
   scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                      breaks = c(0,100,10000,1000000),
                      labels = c(0,expression(10^2),expression(10^4),expression(10^6)))+
   labs(x=NULL, y=expression(Bacteria~biovolume~'('*µm^3*')'))

ggsave("Figures/Figure2B_BacteriaBiogeography.png",width = 8.5, height=4,dpi=400)
ggsave("Figures/Figure2B_BacteriaBiogeography.pdf",width = 8.5, height=4)

# -----------------------------------------


### Per Fish Bacteria and Mucus Analysis
### Figure 5E

PlotFigure5E<-Alldata %>% filter(Condition!="Axenic") %>%
   ggplot(aes(x=TotalBacBiovolume, y=NormalizedGutMucus, fill=Condition, shape=Condition)) +
   geom_point(size=6, alpha=0.8)+
   scale_shape_manual(values=c(22,23))+
   scale_fill_manual(values = c("#a8ddb5","#0868ac"))+
   scale_x_continuous(trans=scales::pseudo_log_trans(base = 10),
                      breaks = c(1000,10000,100000,1000000), limits=c(NA,1e6),
                      labels = c(expression(10^3),expression(10^4),expression(10^5),expression(10^6)))+
   theme(legend.position="top", panel.background = element_rect(fill="grey90", color="transparent"),
         panel.grid.major = element_line(color="white"))+
   labs(x=expression(Bacteria~biovolume~'('*µm^3*')'), y="Normalized Mucus Intensity (A.U)", fill="Fish Model", shape="Fish Model")
PlotFigure5E
ggsave("Figures/Figure5E_GutNormMucusBacPerFish.png",width = 5,height=5,dpi=400)
ggsave("Figures/Figure5E_GutNormMucusBacPerFish.pdf",width = 5,height=5)

# -----------------------------------------


### TotalBacteria and TotalGutMucus association
### Figure 5F
PlotFigure5F<-Alldata %>% filter(Condition!="Axenic") %>%
   ggplot(aes(x=Condition, y=BacteriavsMucus, shape=Condition, fill=Condition)) +
   geom_jitter(size=3, width=0.1)+
   geom_boxplot(alpha=0.6, outlier.shape = NA)+
   scale_shape_manual(values=c(22,23))+
   scale_fill_manual(values = c("#a8ddb5","#0868ac"))+
   scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                      breaks = c(0,10,100,1000,10000), limits=c(0,NA),
                      labels = c(0,expression(10^1),expression(10^2),
                                 expression(10^3),expression(10^4)))+
   labs(x=NULL, y="Total Bacteria/Normalized Mucus")
PlotFigure5F

ggsave("Figures/Figure5F_NormBacterMucus.png",width = 3,height=5,dpi=400)
ggsave("Figures/Figure5F_NormBacterMucus.pdf",width = 3,height=5)

# -----------------------------------------


### Normalized entire gut mucus analysis for conventional and mix 9
### Figure 5G

PlotFigure5G <- Alldata %>%
   filter(Condition!="Axenic") %>%
   ggplot(aes(x=Condition, y=NormalizedGutMucus, shape=Condition, fill=Condition)) +
   geom_jitter(size=3, width=0.1)+
   geom_boxplot(alpha=0.6, outlier.shape = NA)+
   scale_fill_manual(values = c("#a8ddb5","#0868ac"))+
   scale_shape_manual(values=c(22,23))+
   stat_compare_means(comparisons=list(c("Conv", "Mix9")),
                      label="p.format")+
   labs(x=NULL, y="Normalized Mucus Intensity (A.U)", title="Entire gut")
PlotFigure5G

ggsave("Figures/Figure5G_GutNormMucus.png",width = 3,height=5,dpi=400)
ggsave("Figures/Figure5G_GutNormMucus.pdf",width = 3,height=5)

# -----------------------------------------


### mucus distribution in proximal and distal gut
### Normalized mucus intensity for intestine (proximal and distal) for various conditions
### Figure 5GH

stat.test <- compare_means(value ~ Condition, group.by = "Location", (Alldata  %>%
                              mutate(IntestineMucus=NormalizedPMucus+NormalizedDMucus) %>%
                              pivot_longer(c(IntestineMucus, NormalizedGutMucus)) %>%
                              mutate(Location=recode(name,"IntestineMucus"="Intestine (excluding bulb)",
                                                     "NormalizedGutMucus"="Entire gut"),
                                     Location=factor(Location,levels=c("Entire gut","Intestine (excluding bulb)"))) %>%
                              filter(Condition!="Axenic" |  Location!="Entire gut")))


PlotFigureGH <- Alldata  %>%
   mutate(IntestineMucus=NormalizedPMucus+NormalizedDMucus) %>%
   pivot_longer(c(IntestineMucus, NormalizedGutMucus)) %>%
   mutate(Location=recode(name,"IntestineMucus"="Intestine (excluding bulb)",
                               "NormalizedGutMucus"="Entire gut"),
          Location=factor(Location,levels=c("Entire gut","Intestine (excluding bulb)"))) %>%
   filter(Condition!="Axenic" |  Location!="Entire gut") %>%
   ggplot(aes(x=Condition, y=value)) +
   facet_grid(.~Location, scales = "free_x", space = "free_x") +
   geom_jitter(aes(shape=Condition, fill=Condition), size=3, width=0.1)+
   geom_boxplot(aes(fill=Condition), alpha=0.6, outlier.shape = NA)+
   scale_shape_manual(values=c(21,22,23))+
   scale_fill_manual(values = c("grey60","#a8ddb5","#0868ac"))+
   scale_y_continuous(limits=c(NA, 1350))+
   geom_bracket(aes(xmin = group1, xmax = group2, label = signif(p, 2)),
      data = stat.test,
      y.position = c(1100, 1200, 1300, 500))+
   theme(panel.spacing = unit(0.8, "lines"))+
   labs(x=NULL, y="Normalized Mucus Intensity (A.U)")
PlotFigureGH

# -----------------------------------------


### plot figure 5 all together

((PlotFigure5E)/(PlotFigure5F+PlotFigureGH+plot_layout(widths=c(1,3))))+
   plot_annotation(tag_levels = list(c("E","F","G","H")))  &
   theme(plot.tag = element_text(size = 22))

ggsave("Figures/Figure5EFGH.png",width = 10.5,height=9,dpi=400)
ggsave("Figures/Figure5EFGH.pdf",width = 10.5,height=9)




# -----------------------------------------

### Figure S1

rawdata <- gsheet::gsheet2tbl("https://docs.google.com/spreadsheets/d/1bmXLR8tvQ_LFOVfsBdXfs2pJGfwwEGW8ohENEwANJOc/edit?pli=1#gid=761509211")

data <- rawdata %>% drop_na(Biovolume_um3) %>%
   arrange(Biovolume_um3)

biovolumeplot<-data %>%
   ggplot(aes(y=Biovolume_um3))+
   geom_histogram(color="white", fill="#0868ac")+
   scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                      breaks = c(0,10,100,1000,10000,100000,1000000,10000000),
                      labels = c(0,expression(10^1),expression(10^2),
                                 expression(10^3),expression(10^4),
                                 expression(10^5),expression(10^6),expression(10^7)))+
   scale_x_continuous(expand=c(0,0))+
   labs(x="Number of fish (n)", y=expression(Bacteria~biovolume~'('*µm^3*')'), fill=NULL, color=NULL)
biovolumeplot

ggsave("Figures/FigureS1_Mix9_biovolumehist.png", width=3, height=4.5, bg="transparent")
ggsave("Figures/FigureS1_Mix9_biovolumehist.pdf", width=3, height=4.5)
