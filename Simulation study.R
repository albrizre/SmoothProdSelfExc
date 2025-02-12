library(nimble)
library(sp)
library(rgdal)
library(rgeos)
library(sf)
library(spdep)
library(maptools)
library(spatstat)
library(MASS)
library(npreg)

# Set the working directory
# setwd("...")

# Load codes for the standard and the proposed self-exciting model. There are two versions for each one
source("Model codes/standard_self.R")
source("Model codes/standard_self2.R")
source("Model codes/smoothprod_self.R")
source("Model codes/smoothprod_self2.R")

# Create window for simulation
study_area=square(r = 1500)
study_area=as(study_area,"SpatialPolygons")

###################################################################################
# Simulating point patterns according to Algorithm 1 in the paper and model fitting
###################################################################################

# Choose a kappa(x1,x2,t) function for the offspring parameter (see Algorithm 1 in the paper)
# Here we include the three functions employed for the simulation study described in the paper
kappa_choice="A"
if (kappa_choice=="A"){kappa_function <- function(x,y,t) (1/3000)*(x+y)*exp(-t/100)}
if (kappa_choice=="B"){kappa_function <- function(x,y,t) 0.75*as.numeric(y>=1000)*as.numeric(t>=150)}
if (kappa_choice=="C"){kappa_function <- function(x,y,t) 1*exp(-(sqrt((x-750)^2+(y-750)^2))/500)*exp(-abs(t-150)/100)}

# Choose values for lambda0, sigma0 and rho0 (see Algorithm 1 in the paper)
sigma0=50
rho0=1
lambda0_sim=2e-07
Tmax=365

# Choose model type ("standard" or "smoothprod", "standard" needs to be fit first for each specific pattern)
model="standard"

#######################################################################################
# Simulation of the spatio-temporal point pattern starts (see Algorithm 1 in the paper)
#######################################################################################

# Choose replicate seed (for the paper, seed was varied from 1 to 10 in each case)
seed=1
set.seed(seed)

pattern=c()
pattern_back=rpoispp(lambda = lambda0_sim*Tmax, win = as.owin(study_area)) # background events
times_sim=runif(pattern_back$n,0,Tmax)
pattern=rbind(pattern,cbind(pattern_back$x,pattern_back$y,times_sim))
offspring=rpois(pattern_back$n,kappa_function(pattern_back$x,pattern_back$y,times_sim))
condition_while=T
l = 0
pattern=cbind(pattern,0)
colnames(pattern)=c("x","y","t","Gen")
pattern_aux=pattern
while (condition_while){
  print(l)
  l=l+1
  if (!is.null(nrow(pattern_aux))){
    for (i in 1:nrow(pattern_aux)){
      if (offspring[i]>0){
        times_sim_offspring=pattern_aux[i,"t"]+rexp(offspring[i],rho0)
        aux_mvrnorm=mvrnorm(1,c(pattern_aux[i,"x"],pattern_aux[i,"y"]),diag(sigma0^2,2))
        x_sim_offspring=aux_mvrnorm[1]
        y_sim_offspring=aux_mvrnorm[2]
        pattern=rbind(pattern,cbind(x_sim_offspring,y_sim_offspring,times_sim_offspring,l))
      }
    }
  } else {
    if (offspring>0){
      times_sim_offspring=pattern_aux[3]+rexp(offspring,rho0)
      aux_mvrnorm=mvrnorm(1,c(pattern_aux[1],pattern_aux[2]),diag(sigma0^2,2))
      x_sim_offspring=aux_mvrnorm[1]
      y_sim_offspring=aux_mvrnorm[2]
      pattern=rbind(pattern,cbind(x_sim_offspring,y_sim_offspring,times_sim_offspring,l))
    }
  }
  if (sum(pattern[,"Gen"]==l)>1){
    pattern_aux=pattern[pattern[,"Gen"]==l,]
    offspring=rpois(nrow(pattern_aux),kappa_function(pattern_aux[,"x"],pattern_aux[,"y"],pattern_aux[,"t"]))
  } else if (sum(pattern[,"Gen"]==l)==1){
    pattern_aux=pattern[pattern[,"Gen"]==l,]
    offspring=rpois(1,kappa_function(pattern_aux[1],pattern_aux[2],pattern_aux[3]))
  } else {
    condition_while=F
  }
}
pattern=pattern[order(pattern[,3]),]
pattern=pattern[pattern[,3]<=Tmax,]
rownames(pattern)=1:nrow(pattern)

# Filter events and grid cells
events_sp=data.frame(x=pattern[,1],y=pattern[,2],Id=1:nrow(pattern))
coordinates(events_sp)=~x+y
proj4string(events_sp)=proj4string(study_area)
events_study_area=events_sp

# Time events
time_events=pattern[,3]

# Disaggregate grid
aux_grid=as(st_make_grid(as(study_area,"sf"), cellsize = 75),"Spatial")
grid_square=SpatialPolygonsDataFrame(aux_grid,data=data.frame(IdS=1:length(aux_grid)),match.ID = F)

# Partition of the time interval
time_points=seq(0,Tmax,1)

#################################################################################################
# Variable construction for model fitting (computation of distances across events, splines, etc.)
#################################################################################################

n_previous_aux=c()
for (i in 1:length(time_points)){
  if (length(which(time_events<time_points[i]))==0){
    n_previous_aux=c(n_previous_aux,0)
  } else {
    n_previous_aux=c(n_previous_aux,which(time_events<time_points[i])[length(which(time_events<time_points[i]))])
  }
}

# Identify time_points depending on n_previous
previous_0=which(n_previous_aux==0)
previous_1=which(n_previous_aux==1)
previous_2=which(n_previous_aux>=2)
N_0=length(events_sp)+length(previous_0)
if (length(previous_1)>0){
  N_1=N_0+length(previous_1)
} else {
  N_1=NULL
}

constants <- list(
  
  N = length(events_sp)+length(time_points),
  N_events = length(events_sp),
  N_integral = length(time_points),
  w=c(rep(0,length(events_sp)),rep(1,length(time_points))), # 1 is the length of the time interval (in days)
  I=c(rep(1,length(events_sp)),rep(0,length(time_points))),
  n_previous=c(0:(length(events_sp)-1),n_previous_aux),
  time_events = time_events,
  time_points = time_points,
  time_all = c(time_events,time_points),
  C = 100000000,
  pi = pi,
  lambda0_base = length(events_sp)/(sum(raster::area(grid_square))*Tmax),
  S = length(grid_square),
  N_0 = N_0,
  N_1 = N_1,
  Tmax = Tmax,
  area_cell = unique(gArea(grid_square,byid = T))
  
)

# Temporal distances (events-events)
matrix_distances=matrix(Inf,nrow=constants$N_events,constants$N_events)
for (i in 1:constants$N_events){
  if (constants$n_previous[i]>=1){
    matrix_distances[i,1:constants$n_previous[i]]=constants$time_all[i]-constants$time_events[1:constants$n_previous[i]]
  }
}
# Temporal distances (events-time points)
constants$matrix_distances=matrix_distances
matrix_distances_events_time_points=matrix(Inf,nrow=constants$N_events,constants$N_integral)
for (i in 1:constants$N_events){
  matrix_distances_events_time_points[i,]=abs(constants$time_all[i]-time_points)
}
constants$matrix_distances_events_time_points=matrix_distances_events_time_points

# Spatial distances (events-grid cells)
constants$matrix_spatial_distances_events_grid=gDistance(grid_square,events_sp,byid = T)
rownames(constants$matrix_spatial_distances_events_grid)=NULL
colnames(constants$matrix_spatial_distances_events_grid)=NULL
# Spatial distances (events-events)
constants$matrix_spatial_distances_events=gDistance(events_sp,events_sp,byid = T)
constants$matrix_spatial_distances_events[upper.tri(constants$matrix_spatial_distances_events)]=Inf
diag(constants$matrix_spatial_distances_events)=Inf

# Spline basis for b_S (spatial correction factor)
centroids_grid=gCentroid(grid_square,byid = T)
x=centroids_grid@coords[,1]
y=centroids_grid@coords[,2]
input=cbind(x,y)
k=6
knots=expand.grid(quantile(x,1:(k-1)/k),quantile(y,1:(k-1)/k))
X_splines=basis.tps(input, knots, m = 2, rk = F, intercept = FALSE, ridge = FALSE) 
X_splines=cbind(1,scale(X_splines[,3:ncol(X_splines)])) # scale
constants$nbasis=ncol(X_splines)
constants$X_splines=X_splines
events_grid=raster::intersect(events_sp,grid_square)
constants$IdS_events=apply(constants$matrix_spatial_distances_events_grid,1,which.min)

# Spline basis for b_T (temporal correction factor)
kT=20
knotsT=quantile(time_points,1:(kT-1)/kT)
X_splinesT=as.matrix(splines::bs(0:Tmax, knots = knotsT))
X_splinesT=cbind(1,scale(X_splinesT[,1:ncol(X_splinesT)])) # scale
constants$nbasisT=ncol(X_splinesT)
constants$X_splinesT=X_splinesT

# Fit model with NIMBLE
set.seed(12345)
data_aux <- list(zeros = rep(0,constants$N))
inits <- function() list(lambda0 = constants$lambda0_base,
                         phi_T = 25,
                         phi_S = 50,
                         theta = 0.01,
                         alpha = rep(0,constants$nbasis),
                         alphaT = rep(0,constants$nbasisT))

Sys.time()
print(paste0("N_1",N_1))
if (!is.null(N_1) & model=="standard"){code=standard_self2}
if (is.null(N_1) & model=="standard"){code=standard_self}
if (!is.null(N_1) & model=="smoothprod"){code=smoothprod_self2}
if (is.null(N_1) & model=="smoothprod"){code=smoothprod_self}
# Prior information for the proposed model is based on the estimates of the standard model (you need to fit the standard model first)
if (model!="standard"){
  if (file.exists(paste0("Models/self_sim_revised_kappa_",kappa_choice,"_","standard","_",seed,".rda"))){
    load(paste0("Models/self_sim_revised_kappa_",kappa_choice,"_","standard","_",seed,".rda"))
    constants$mean_theta=mean(mcmc.output$samples[,"theta"])
    constants$sd_theta=sd(mcmc.output$samples[,"theta"])
    constants$mean_phi_S=mean(mcmc.output$samples[,"phi_S"])
    constants$sd_phi_S=sd(mcmc.output$samples[,"phi_S"])
    constants$mean_phi_T=mean(mcmc.output$samples[,"phi_T"])
    constants$sd_phi_T=sd(mcmc.output$samples[,"phi_T"])
    params_monitor=c("lambda0","phi_T","phi_S","theta","b_S","b_T")
  }
  else {
    print("You need to fit the standard self-exciting model first. The estimates of this model are used for constructing the priors in the proposed model with smooth space-time varying productivity parameter")
  }
} else {
  params_monitor=c("lambda0","phi_T","phi_S","theta")
}
Sys.time()
mcmc.output <- nimbleMCMC(code, data = data_aux, inits = inits, constants = constants,
                          monitors = params_monitor, thin = 10,
                          niter = 10000, nburnin = 5000, nchains = 1, # small number of iterations for testing
                          summary = TRUE, WAIC = TRUE)
Sys.time()
save(mcmc.output,file=paste0("Models/self_sim_revised_kappa_",kappa_choice,"_",model,"_",seed,".rda"))
  



