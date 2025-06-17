library(nimble)
library(sp)
library(rgdal)
library(rgeos)
library(sf)
library(spdep)
library(maptools)
library(spatstat)
library(DRHotNet)
library(ggplot2)
library(ggspatial)
library(RColorBrewer)
library(mgcv)
library(splines)
library(fields)
library(npreg)

# Set the working directory
# setwd("...")

# Create window for simulation
study_area=square(r = 1500)
study_area=as(study_area,"SpatialPolygons")
# Disaggregate square into three squares
study_area=as(st_make_grid(as(study_area,"sf"), cellsize = 500),"Spatial")
aux_grid=as(st_make_grid(as(study_area,"sf"), cellsize = 75),"Spatial")
grid_square=SpatialPolygonsDataFrame(aux_grid,data=data.frame(IdS=1:length(aux_grid)),match.ID = F)

# Plot kappa functions (Figure 1 of the paper)
aux_grid2=as(st_make_grid(as(study_area,"sf"), cellsize = 25),"Spatial")
aux_grid2_c=gCentroid(aux_grid2,byid = T)
for (kappa_choice in c("A","B","C")){
  
  kappa_value=c()
  if (kappa_choice=="A"){kappa_function <- function(x,y,t) (1/3000)*(x+y)*exp(-t/100);t_max=0;x_max=1500;y_max=1500}
  if (kappa_choice=="B"){kappa_function <- function(x,y,t) 0.75*as.numeric(y>=1000)*as.numeric(t>=150);t_max=150;x_max=0;y_max=1000}
  if (kappa_choice=="C"){kappa_function <- function(x,y,t) 0.75*(exp(-(sqrt((x-375)^2+(y-1125)^2))/500)+exp(-(sqrt((x-1125)^2+(y-375)^2))/500))*as.numeric((t>=50 & t<=100) | t>=150);t_max=200;x_max=1125;y_max=375}

  for (i in 1:length(aux_grid2)){
    kappa_value=c(kappa_value,kappa_function(x=aux_grid2_c@coords[i,1],y=aux_grid2_c@coords[i,2],t=t_max))
  }
  
  kappa_value_T=c()
  for (i in seq(0,365,0.01)){
    kappa_value_T=c(kappa_value_T,kappa_function(x=x_max,y=y_max,t=i))
  }
  
  aux=SpatialPolygonsDataFrame(aux_grid2,data.frame(kappa_value),match.ID = F)
  proj4string(aux)="+proj=utm +zone=30 ellps=WGS84"
  aux_gg=as(aux,"sf")
  proj4string(study_area)="+proj=utm +zone=30 ellps=WGS84"
  study_area_gg=as(study_area,"sf")
  
  if (kappa_choice=="A"){
    postscript(paste0("Figures/kappa_",kappa_choice,".eps"),family = "Helvetica")
    print(ggplot() +
            annotation_spatial(aux_gg) +
            layer_spatial(aux_gg,aes(fill=kappa_value),col="transparent")+ #,col="transparent"
            scale_fill_continuous(type = "viridis",name=expression(kappa[1](italic(bold(x)),italic(t)[italic(max)])))+
            theme_bw()+
            ggtitle("Offspring parameter function over space")+
            scale_x_continuous(expand=c(0,0),breaks=c(0,500,1000,1500))+
            scale_y_continuous(expand=c(0,0),breaks=c(0,500,1000,1500))+
            coord_sf(xlim = c(0, 1500), ylim = c(0, 1500), 
                     datum=st_crs(aux_gg))+
            labs(subtitle = paste0("Scenario ",ifelse(kappa_choice=="A",1,ifelse(kappa_choice=="B",2,3))))+
            theme(text = element_text(size = 23),
                  axis.title = element_text(face = "bold"),
                  title = element_text(colour = "black",face = "bold"),
                  legend.position = "right"))
    dev.off()
    
    df=data.frame(Time=seq(0,365,0.01),kappa_value_T)
    postscript(paste0("Figures/kappa_T_",kappa_choice,".eps"),family = "Helvetica")
    print(ggplot(data=df) +
            geom_line(aes(x=Time,y=kappa_value_T))+
            xlab("Time")+ylab(expression(kappa[1](italic(bold(x)[italic(max)]),italic(t))))+
            theme_bw()+
            scale_size_manual(values = c(rep(1,10),3))+
            scale_color_manual(values = c(rep("gray80",10),"gray40"))+
            ggtitle("Offspring parameter function over time")+
            labs(subtitle=paste0("Scenario ",ifelse(kappa_choice=="A",1,ifelse(kappa_choice=="B",2,3))))+
            theme(text = element_text(size = 23),
                  axis.title = element_text(face = "bold"),
                  axis.text = element_text(colour = "black"),
                  title = element_text(colour = "black",face = "bold"),
                  legend.position = "none"))
    dev.off()
  } else if (kappa_choice=="B"){
    postscript(paste0("Figures/kappa_",kappa_choice,".eps"),family = "Helvetica")
    print(ggplot() +
            annotation_spatial(aux_gg) +
            layer_spatial(aux_gg,aes(fill=kappa_value),col="transparent")+ #,col="transparent"
            scale_fill_continuous(type = "viridis",name=expression(kappa[2](italic(bold(x)),italic(t)[italic(max)])))+
            theme_bw()+
            ggtitle("Offspring parameter function over space")+
            scale_x_continuous(expand=c(0,0),breaks=c(0,500,1000,1500))+
            scale_y_continuous(expand=c(0,0),breaks=c(0,500,1000,1500))+
            coord_sf(xlim = c(0, 1500), ylim = c(0, 1500), 
                     datum=st_crs(aux_gg))+
            labs(subtitle = paste0("Scenario ",ifelse(kappa_choice=="A",1,ifelse(kappa_choice=="B",2,3))))+
            theme(text = element_text(size = 23),
                  axis.title = element_text(face = "bold"),
                  title = element_text(colour = "black",face = "bold"),
                  legend.position = "right"))
    dev.off()
    
    df=data.frame(Time=seq(0,365,0.01),kappa_value_T)
    postscript(paste0("Figures/kappa_T_",kappa_choice,".eps"),family = "Helvetica")
    print(ggplot(data=df) +
            geom_line(aes(x=Time,y=kappa_value_T))+
            xlab("Time")+ylab(expression(kappa[2](italic(bold(x)[italic(max)]),italic(t))))+
            theme_bw()+
            scale_size_manual(values = c(rep(1,10),3))+
            scale_color_manual(values = c(rep("gray80",10),"gray40"))+
            ggtitle("Offspring parameter function over time")+
            labs(subtitle=paste0("Scenario ",ifelse(kappa_choice=="A",1,ifelse(kappa_choice=="B",2,3))))+
            theme(text = element_text(size = 23),
                  axis.title = element_text(face = "bold"),
                  axis.text = element_text(colour = "black"),
                  title = element_text(colour = "black",face = "bold"),
                  legend.position = "none"))
    dev.off()
  } else if (kappa_choice=="C"){
    postscript(paste0("Figures/kappa_",kappa_choice,".eps"),family = "Helvetica")
    print(ggplot() +
            annotation_spatial(aux_gg) +
            layer_spatial(aux_gg,aes(fill=kappa_value),col="transparent")+ #,col="transparent"
            scale_fill_continuous(type = "viridis",name=expression(kappa[3](italic(bold(x)),italic(t)[italic(max)])))+
            theme_bw()+
            ggtitle("Offspring parameter function over space")+
            scale_x_continuous(expand=c(0,0),breaks=c(0,500,1000,1500))+
            scale_y_continuous(expand=c(0,0),breaks=c(0,500,1000,1500))+
            coord_sf(xlim = c(0, 1500), ylim = c(0, 1500), 
                     datum=st_crs(aux_gg))+
            labs(subtitle = paste0("Scenario ",ifelse(kappa_choice=="A",1,ifelse(kappa_choice=="B",2,3))))+
            theme(text = element_text(size = 23),
                  axis.title = element_text(face = "bold"),
                  title = element_text(colour = "black",face = "bold"),
                  legend.position = "right"))
    dev.off()
    
    df=data.frame(Time=seq(0,365,0.01),kappa_value_T)
    postscript(paste0("Figures/kappa_T_",kappa_choice,".eps"),family = "Helvetica")
    print(ggplot(data=df) +
            geom_line(aes(x=Time,y=kappa_value_T))+
            xlab("Time")+ylab(expression(kappa[3](italic(bold(x)[italic(max)]),italic(t))))+
            theme_bw()+
            scale_size_manual(values = c(rep(1,10),3))+
            scale_color_manual(values = c(rep("gray80",10),"gray40"))+
            ggtitle("Offspring parameter function over time")+
            labs(subtitle=paste0("Scenario ",ifelse(kappa_choice=="A",1,ifelse(kappa_choice=="B",2,3))))+
            theme(text = element_text(size = 23),
                  axis.title = element_text(face = "bold"),
                  axis.text = element_text(colour = "black"),
                  title = element_text(colour = "black",face = "bold"),
                  legend.position = "none"))
    dev.off()
  }
  
}

# Summarize results of the simulation study and plot average estimates of b(x) and b(t) (Figure 2 of the paper)
df_pars=c()
model="smoothprod"
df_WAIC=c()
for (kappa_choice in c("A","B","C")){
  
  df=c()
  df_S=c()
  for (seed in c(1:100)){
    
    print(seed)
    
    load(paste0("Models/self_sim_revised_kappa_",kappa_choice,"_","standard","_",seed,".rda"))
    WAIC_standard=mcmc.output$WAIC$WAIC
    load(paste0("Models/self_sim_revised_kappa_",kappa_choice,"_",model,"_",seed,".rda"))
    df_WAIC=rbind(df_WAIC,data.frame(kappa_choice,model="smoothprod",Replicate=seed,WAIC_standard=WAIC_standard,WAIC_prop=mcmc.output$WAIC$WAIC))
    
    # Plot estimates
    select=grep("b_S",colnames(mcmc.output$samples))[1:400]
    estimates_bS=as.numeric(mcmc.output$summary[select,1])
    df_S=rbind(df_S,data.frame(Replicate=seed,
                               Cell=1:length(estimates_bS),
                               Mean=estimates_bS))
    
    select=grep("b_T",colnames(mcmc.output$samples))
    estimates_bT=as.numeric(mcmc.output$summary[select,1])
    estimates_bT_L=as.numeric(mcmc.output$summary[select,4])
    estimates_bT_U=as.numeric(mcmc.output$summary[select,5])
    df=rbind(df,data.frame(Replicate=seed,
                           Time=0:(length(estimates_bT)-1),
                           Mean=estimates_bT,
                           Lo=estimates_bT_L,
                           Up=estimates_bT_U))
    
  }
  
  # Mean spatial correction
  mean=c()
  for (cell in unique(df_S$Cell)){
    mean=c(mean,median(df_S$Mean[df_S$Cell==cell]))
  }
  df_S=rbind(df_S,data.frame(Replicate=999,
                             Cell=1:length(estimates_bS),
                             Mean=mean))
  df_S$Replicate=as.factor(df_S$Replicate)
  
  aux=SpatialPolygonsDataFrame(grid_square,data.frame(estimates=df_S$Mean[df_S$Replicate=="999"]),match.ID = F)
  proj4string(aux)="+proj=utm +zone=30 ellps=WGS84"
  aux_gg=as(aux,"sf")
  proj4string(study_area)="+proj=utm +zone=30 ellps=WGS84"
  study_area_gg=as(study_area,"sf")
  
  postscript(paste0("Figures/spatial_correction_sim_",kappa_choice,".eps"),family = "Helvetica")
  print(ggplot() +
          annotation_spatial(aux_gg) +
          layer_spatial(aux_gg,aes(fill=estimates),col="transparent")+ #,col="transparent"
          scale_fill_continuous(type = "viridis",name=expression(italic(b)(italic(bold(x)))))+
          theme_bw()+
          ggtitle("Spatial correction of the productivity parameter")+
          scale_x_continuous(expand=c(0,0),breaks=c(0,500,1000,1500))+
          scale_y_continuous(expand=c(0,0),breaks=c(0,500,1000,1500))+
          coord_sf(xlim = c(0, 1500), ylim = c(0, 1500), 
                   datum=st_crs(aux_gg))+
          labs(subtitle = paste0("Scenario ",ifelse(kappa_choice=="A",1,ifelse(kappa_choice=="B",2,3))))+
          theme(text = element_text(size = 23),
                axis.title = element_text(face = "bold"),
                # axis.text = element_blank(),
                # axis.ticks = element_blank(),
                title = element_text(colour = "black",face = "bold"),
                legend.position = "right"))
  dev.off()
  
  # Mean temporal correction
  mean=c()
  lo=c()
  up=c()
  for (time in unique(df$Time)){
    mean=c(mean,median(df$Mean[df$Time==time]))
    lo=c(lo,quantile(df$Mean[df$Time==time],0.025))
    up=c(up,quantile(df$Mean[df$Time==time],0.975))
  }
  df=rbind(df,data.frame(Replicate=999,
                         Time=0:(length(estimates_bT)-1),
                         Mean=mean,
                         Lo=lo,
                         Up=up))
  df$Replicate=as.factor(df$Replicate)
  df_aux=df[df$Replicate==999,]
  
  postscript(paste0("Figures/temporal_correction_sim_",kappa_choice,".eps"),family = "Helvetica")
  print(ggplot(data=df) +
          geom_line(aes(x=Time,y=Mean,group=Replicate,col=Replicate,size=Replicate))+
          xlab("Time")+ylab(expression(italic(b)(italic(t))))+
          theme_bw()+
          scale_size_manual(values = c(rep(1,length(unique(df$Replicate))-1),3))+
          scale_color_manual(values = c(rep("gray80",length(unique(df$Replicate))-1),"gray40"))+
          ggtitle("Temporal correction of the productivity parameter")+
          labs(subtitle=paste0("Scenario ",ifelse(kappa_choice=="A",1,ifelse(kappa_choice=="B",2,3))))+
          theme(text = element_text(size = 23),
                axis.title = element_text(face = "bold"),
                axis.text = element_text(colour = "black"),
                title = element_text(colour = "black",face = "bold"),
                legend.position = "none"))
  dev.off()
  
  
}
df_WAIC$Improves=as.numeric(df_WAIC$WAIC_prop<df_WAIC$WAIC_standard)
table(df_WAIC$Improves)
View(df_WAIC)

# WAIC comparison (Figure 3 of the paper)
table(df_WAIC$Improves,df_WAIC$kappa_choice)
df_WAIC_plot=c()
for (i in 1:nrow(df_WAIC)){
  df_WAIC_plot=rbind(df_WAIC_plot,data.frame(kappa_choice=df_WAIC$kappa_choice[i],
                                             Replicate=df_WAIC$Replicate[i],
                                             WAIC=-df_WAIC$WAIC_prop[i]+df_WAIC$WAIC_standard[i]))
}
df_WAIC_plot$kappa_choice[df_WAIC_plot$kappa_choice=="A"]="1"
df_WAIC_plot$kappa_choice[df_WAIC_plot$kappa_choice=="B"]="2"
df_WAIC_plot$kappa_choice[df_WAIC_plot$kappa_choice=="C"]="3"

postscript(paste0("Figures/WAIC_boxplots.eps"),family = "Helvetica")
ggplot(df_WAIC_plot, aes(x=kappa_choice,y=WAIC))+
  geom_boxplot(fill=brewer.pal(6,"BuPu")[4])+
  theme_bw()+
  ggtitle("Model comparison across scenarios")+
  xlab("Scenario")+ylab("WAIC variation")+
  theme(text = element_text(size = 23),
        axis.title = element_text(face = "bold"),
        # axis.text = element_blank(),
        # axis.ticks = element_blank(),
        title = element_text(colour = "black",face = "bold"),
        legend.position = "right")
dev.off()
