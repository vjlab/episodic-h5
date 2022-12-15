library(readxl)
library(ggplot2)
library(tidyr)
library(dplyr)
library(lubridate)
library(forcats)
library(ggh4x)

setwd("/Users/ruopengxie/vjlab Dropbox/Ruopeng Xie/HPAIH5-resugence/03.analysis/outbreaks_dating")

df.tmrca <- read_excel("TMRCA_summary.xlsx")

#panzootic
df.panzootic <- df.tmrca %>% 
  filter(Lineage == "Panzootic") %>% 
  group_by(Segments) %>% 
  mutate(index = row_number(),
         monophyletic_clade = paste(Segments,index,sep="_"),
         TMRCA = as.Date(date_decimal(TMRCA,tz = "EST")),
         TMRCA_lower = as.Date(date_decimal(TMRCA_lower,tz = "EST")),
         TMRCA_upper = as.Date(date_decimal(TMRCA_upper,tz = "EST")),
         subtype = case_when(grepl("H5N1",subtype) ~ "H5N1",
                             grepl("H5N8",subtype) ~ "H5N8",
                             TRUE ~ "Others"))

# add tmrca and duration
unknown_duration <- rbind(df.panzootic[,c("TMRCA","num_samples","monophyletic_clade")] %>% dplyr::rename(,x_date = TMRCA),
                          df.panzootic[,c("first_sample_date","num_samples","monophyletic_clade")] %>% dplyr::rename(,x_date = first_sample_date)) %>% 
  mutate(type = "unknown") %>% 
  filter(num_samples > 100)

tmrca_CI_duration <- rbind(df.panzootic[,c("TMRCA_lower","num_samples","monophyletic_clade")] %>% dplyr::rename(,x_date = TMRCA_lower),
                           df.panzootic[,c("TMRCA_upper","num_samples","monophyletic_clade")] %>% dplyr::rename(,x_date = TMRCA_upper)) %>% 
  mutate(type = "95% HPD TMRCA") %>% 
  filter(num_samples > 100)

detected_duration <- rbind(df.panzootic[,c("last_sample_date","num_samples","monophyletic_clade")] %>% dplyr::rename(,x_date = last_sample_date),
                           df.panzootic[,c("first_sample_date","num_samples","monophyletic_clade")] %>% dplyr::rename(,x_date = first_sample_date)) %>% 
  mutate(type = "detected")

data.panzootic.input.line <- rbind(unknown_duration,detected_duration,tmrca_CI_duration) %>% 
  mutate(segments = factor(gsub("_.*","",monophyletic_clade),levels = c("NS","MP","NP","PA","PB1","PB2","N1","N8","HA"))) %>% 
  arrange(.,segments,num_samples) %>% 
  mutate(monophyletic_clade= factor(monophyletic_clade,levels= unique(monophyletic_clade)),
         num_samples = case_when(type == "95% HPD TMRCA" ~ num_samples*4,
                                 TRUE ~ num_samples),
         segments = as.character(segments),
         segments = case_when(type == "95% HPD TMRCA" ~ "Others",
                              TRUE ~ segments))

data.panzootic.input.point <- df.panzootic %>% 
  subset(,c("last_sample_date","monophyletic_clade","subtype","Segments"))

(ggplot()+
  geom_line(data=data.panzootic.input.line, aes(y = monophyletic_clade, x = x_date,group = interaction(monophyletic_clade,type),color = segments, alpha = type, size = num_samples))+
  geom_point(data=data.panzootic.input.point, aes(y = monophyletic_clade, x=as.Date(last_sample_date)+30,color=subtype,alpha = subtype),shape=8,size=1,stroke=1)+
  scale_color_manual(breaks = c("HA","N8","N1","PB2","PB1","PA","NP","MP","NS","H5N1","H5N8","Others"),
                     values=c("#c03927","#1f524a","#883c29","#cdad85","#ea6e46","#4586bc","#e493b7","#98d0ea","#51a39d","#883c29","#1f524a","#BABABA"))+
  scale_alpha_manual(values=c(1,1,0,1,0.6,0.4),
                     breaks = c("H5N1","H5N8","Others","detected","unknown","95% HPD TMRCA"))+
  scale_size_continuous(range = c(0.3,3),
                        breaks = c(100,500,1000,2000),
                        limits = c(min(data.panzootic.input.line$num_samples),max(data.panzootic.input.line$num_samples)))+
  scale_x_date("Date", date_breaks =  "2 years", date_labels = "%Y",guide = "axis_minor") +
  theme_classic()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y = element_blank(),
        legend.position = c(0.3,2)))

ggsave("tmrca_panzootic.pdf",width=20,height=15,units="cm")


#JKE
df.jke <- df.tmrca %>% 
  filter(Lineage == "JKE") %>% 
  group_by(Segments) %>% 
  mutate(index = row_number(),
         monophyletic_clade = paste(Segments,index,sep="_"),
         TMRCA = as.Date(date_decimal(TMRCA,tz = "EST")),
         TMRCA_lower = as.Date(date_decimal(TMRCA_lower,tz = "EST")),
         TMRCA_upper = as.Date(date_decimal(TMRCA_upper,tz = "EST")),
         subtype = case_when(grepl("H5N1",subtype) ~ "H5N1",
                             grepl("H5N8",subtype) ~ "H5N8",
                             TRUE ~ "Others"))

# add tmrca and duration
unknown_duration <- rbind(df.jke[,c("TMRCA","num_samples","monophyletic_clade")] %>% dplyr::rename(,x_date = TMRCA),
                          df.jke[,c("first_sample_date","num_samples","monophyletic_clade")] %>% dplyr::rename(,x_date = first_sample_date)) %>% 
  mutate(type = "unknown") %>% 
  filter(num_samples > 100)

tmrca_CI_duration <- rbind(df.jke[,c("TMRCA_lower","num_samples","monophyletic_clade")] %>% dplyr::rename(,x_date = TMRCA_lower),
                           df.jke[,c("TMRCA_upper","num_samples","monophyletic_clade")] %>% dplyr::rename(,x_date = TMRCA_upper)) %>% 
  mutate(type = "95% HPD TMRCA") %>% 
  filter(num_samples > 100)

detected_duration <- rbind(df.jke[,c("last_sample_date","num_samples","monophyletic_clade")] %>% dplyr::rename(,x_date = last_sample_date),
                           df.jke[,c("first_sample_date","num_samples","monophyletic_clade")] %>% dplyr::rename(,x_date = first_sample_date)) %>% 
  mutate(type = "detected")

data.jke.input.line <- rbind(unknown_duration,detected_duration,tmrca_CI_duration) %>% 
  mutate(segments = factor(gsub("_.*","",monophyletic_clade),levels = c("NS","MP","NP","PA","PB1","PB2","N1","N8","HA"))) %>% 
  arrange(.,segments,num_samples) %>% 
  mutate(monophyletic_clade= factor(monophyletic_clade,levels= unique(monophyletic_clade)),
         num_samples = case_when(type == "95% HPD TMRCA" ~ num_samples*4,
                                 TRUE ~ num_samples),
         segments = as.character(segments),
         segments = case_when(type == "95% HPD TMRCA" ~ "Others",
                              TRUE ~ segments))

data.jke.input.point <- df.jke %>% 
  subset(,c("last_sample_date","monophyletic_clade","subtype","Segments"))

(ggplot()+
    geom_line(data=data.jke.input.line, aes(y = monophyletic_clade, x = x_date,group = interaction(monophyletic_clade,type),color = segments, alpha = type, size = num_samples))+
    geom_point(data=data.jke.input.point, aes(y = monophyletic_clade, x=as.Date(last_sample_date)+30,color=subtype,alpha = subtype),shape=8,size=1,stroke=1)+
    scale_color_manual(breaks = c("HA","N8","N1","PB2","PB1","PA","NP","MP","NS","H5N1","H5N8","Others"),
                       values=c("#c03927","#1f524a","#883c29","#cdad85","#ea6e46","#4586bc","#e493b7","#98d0ea","#51a39d","#883c29","#1f524a","#BABABA"))+
    scale_alpha_manual(values=c(1,1,0,1,0.6,0.4),
                       breaks = c("H5N1","H5N8","Others","detected","unknown","95% HPD TMRCA"))+
    scale_size_continuous(range = c(0.3,3),
                          breaks = c(100,500,1000,2000),
                          limits = c(min(data.panzootic.input.line$num_samples),max(data.panzootic.input.line$num_samples)))+
    scale_x_date("Date", date_breaks =  "2 years", date_labels = "%Y",guide = "axis_minor",
                 limits = c(min(data.panzootic.input.line$x_date),max(data.panzootic.input.line$x_date))) +
    theme_classic()+
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line.y = element_blank(),
          legend.position = c(0.3,2)))

ggsave("tmrca_jke.pdf",width=20,height=5,units="cm")
