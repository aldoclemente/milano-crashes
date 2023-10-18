## milano crashes analysis -----------------------------------------------------

if(!require(pacman)) install.packages("pacman")
pacman::p_load("rstudioapi")

# setting working directory
setwd(dirname(getActiveDocumentContext()$path))

source("packages.R")
source("LN-simulations/utils/utils.R")
source("LN-simulations/utils/plot.R")

# loading data -----------------------------------------------------------------
district <- "Lambrate"
foldername <- paste0("data/", district,"/")
if(!dir.exists(foldername)) warning("Cannot run data analysis. Try running preprocessing.R")
folder.img <- paste0(foldername, "imgs/")

load(paste0(foldername,"data_raw.RData"))

# loading road network(mesh) and preprocessed data
load(paste0(foldername, "mesh.RData"))
load(paste0(foldername, "data_raw.RData"))

# normalizing domain -----------------------------------------------------------
mesh_norm <- normalize_mesh_unit(mesh)

# locations
locs <- cbind(crashes$longitudine, crashes$latitudine)

x.min = mesh_norm$x.min
x.max = mesh_norm$x.max
y.min = mesh_norm$y.min
y.max = mesh_norm$y.max

x.norm = (locs[,1] - x.min)/(x.max-x.min)
y.norm = (locs[,2] - y.min)/(y.max-y.min)

locs = cbind(x.norm, y.norm)
mesh_coarse <- mesh_norm$mesh

locs <- projection.points.1.5D(mesh_coarse, locs)
colnames(locs) <- c("x", "y")
locs <- as.data.frame(locs)

# plot(mesh) returns a ggplot object! (see LN-simulations/utils/plot.R)
plot(mesh_coarse) + geom_point(data=locs, aes(x=x, y=y), col="red")

# refining mesh
delta = 0.035
mesh = refine.mesh.1.5D(mesh_coarse, delta=delta)
nrow(mesh$nodes)
FEMbasis = create.FEM.basis(mesh)

# Create regions ---------------------------------------------------------------

# set_region
nregion = nrow(mesh$nodes)
centroids = data.frame(x=mesh$nodes[1:nregion,1],
                        y=mesh$nodes[1:nregion,2])

plot(mesh) + geom_point(data=locs, aes(x=x, y=y), col="red") +
             geom_point(data=centroids, aes(x=x, y=y), col="green4")
  
nedges = nrow(mesh$edges)

source("utils.R")
data_to_edge = is_inside(locs, mesh)
range(data_to_edge)
edge_to_region <- aggregate_edges(centroids, mesh) 

mask_= rep(0, times=nregion)
mask_[ which( !(1:nregion %in% unique(edge_to_region)) == T) ] = 1
while(sum(mask_)){
  centroids = centroids[-which(mask_==1),]
  edge_to_region <- aggregate_edges(centroids, mesh)
  nregion = nrow(centroids)
  
  mask_= rep(0, times=nregion)
  mask_[ which( !(1:nregion %in% unique(edge_to_region)) == T) ] = 1
}

plot(mesh) + geom_point(data=locs, aes(x=x, y=y), col="red") +
             geom_point(data=centroids, aes(x=x, y=y), col="green4")

incidence_matrix = matrix(0, nrow=nregion, ncol=nedges)
for(i in 1:nedges){
  incidence_matrix[ edge_to_region[i],i] = 1
}

# CHECKs on  REGIONs -----------------------------------------------------------
mask_= rep(0, times=nregion)
mask_[ which( !(1:nregion %in% unique(edge_to_region)) == T) ] = 1
sum(mask_) # :)

mask_= rep(0, times=nregion)
mask_[which( apply(incidence_matrix, MARGIN = 1,sum) == 0) ] = 1
sum(mask_)

# Response ---------------------------------------------------------------------
response = rep(0, times= nregion)
for( i in 1:nrow(data_to_edge)){
  response[edge_to_region[data_to_edge[i]]]=response[edge_to_region[data_to_edge[i]]] + 1
}
range(response)

# Plot REGIONs -----------------------------------------------------------------
#
#
#

# GSR-PDE ----------------------------------------------------------------------
lambda <- 10^seq(from=-5, to=0, by=0.25)
GSR_PDE <- smooth.FEM(incidence_matrix=incidence_matrix, observations=response,
                      FEMbasis=FEMbasis, family = "poisson", lambda=lambda, 
                      lambda.selection.criterion = "grid",
                      lambda.selection.lossfunction = "GCV",
                      DOF.evaluation = "exact")

lambda_opt <- GSR_PDE$optimization$lambda_position
plot(log10(lambda), GSR_PDE$optimization$GCV_vector, xlab="log10(lambda)", ylab="GCV")

plot(FEM(GSR_PDE$fit.FEM$coeff[,lambda_opt], FEMbasis),
     linewidth=1)  + scale_color_viridis()

# GSR-PDE with covariates ------------------------------------------------------  
# playing with data_to_edge and edge_to_region
COVARIATES <- cbind(crashes$INTERSEZIONE_NONINTERSEZIONE, 
                    crashes$FONDO_STRADALE, crashes$TIPO_STRADA)
COVARIATES <- as.matrix(COVARIATES)

# Design Matrix
X <- matrix(0,nrow=nregion,ncol=ncol(COVARIATES))

# Filling design matrix: replacing raw categories with the average response value of the category
for(i in 1:nrow(crashes)){
  for(j in 1:ncol(COVARIATES)){
    X[edge_to_region[data_to_edge[i]],j] <- X[edge_to_region[data_to_edge[i]],j] + COVARIATES[i,j]/response[edge_to_region[data_to_edge[i]]]
  }
}

# GSR-PDE with covariates
lambda <- 10^seq(from=-5, to=0, by=0.25)
GSR_PDE <- smooth.FEM(incidence_matrix=incidence_matrix, observations=response, covariates = X,
                      FEMbasis=FEMbasis, family = "poisson", lambda=lambda, 
                      lambda.selection.criterion = "grid",
                      lambda.selection.lossfunction = "GCV",
                      DOF.evaluation = "exact")

lambda_opt <- GSR_PDE$optimization$lambda_position
plot(log10(lambda), GSR_PDE$optimization$GCV_vector, xlab="log10(lambda)", ylab="GCV")
plot(FEM(GSR_PDE$fit.FEM$coeff[,lambda_opt], FEMbasis),
     linewidth=1)  + scale_color_viridis()


