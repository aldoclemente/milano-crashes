if(!require(spatstat)) install.packages("spatstat")
pacman::p_load("spatstat")

aggregate_edges <- function(centroids, mesh){
  LN = as.linnet(mesh)
  lpp_centroids <- spatstat.linnet::lpp(X= ppp(centroids[,1], y=centroids[,2], window= LN$window), 
                                        L= LN)
  
  nedges = nrow(mesh$edges)
  mid_points = matrix(0, nrow=nedges, ncol=2)
  for(e in 1:nedges){  
    mid_points[e,1] =  ( mesh$nodes[ mesh$edges[e,1],1] + mesh$nodes[ mesh$edges[e,2],1] )/2
    mid_points[e,2] =  ( mesh$nodes[ mesh$edges[e,1],2] + mesh$nodes[ mesh$edges[e,2],2] )/2
  }
  
  lpp_midpoints <- spatstat.linnet::lpp(X = ppp(x=mid_points[,1], y= mid_points[,2], window= LN$window),
                                        L = LN)
  
  network_dist = spatstat.linnet::crossdist.lpp(lpp_midpoints, lpp_centroids)
  edges_to_region = matrix(0, nrow = nedges, ncol=1)
  
  for(e in 1:nedges){
    idxs = which( network_dist[e,] == min(network_dist[e,]))
    if(length(idxs) == 1) 
      edges_to_region[e] = idxs
    else
      edges_to_region[e] = idxs[sample(1:length(idxs), size=1)]
  }
  
  return(edges_to_region)
}

is_inside = function(points_, mesh){
  result = matrix(0, nrow=nrow(points_),ncol=1)
  for(i in 1:nrow(points_)){
    for(j in 1:nrow(mesh$edges)){
      node1 = mesh$nodes[mesh$edges[j,1],] 
      node2 = mesh$nodes[mesh$edges[j,2],]
      # x 
      if(abs(node2[2]-node1[2])< 10^5*.Machine$double.eps & 
         abs(points_[i,2] - node2[2])< 10^5*.Machine$double.eps){
        if((node1[1]<node2[1] & points_[i,1] <= node2[1] & points_[i,1] >= node1[1]) | 
           (node2[1]<node1[1] & points_[i,1] <= node1[1] & points_[i,1] >= node2[1]) ){
          result[i] = j
          break
        }
      }else if(abs(node2[1]-node1[1])< 10^5*.Machine$double.eps & 
               abs(points_[i,1] - node2[1])< 10^5*.Machine$double.eps){
        if((node1[2]<node2[2] & points_[i,2] <= node2[2] & points_[i,2] >= node1[2]) | 
           (node2[2]<node1[2] & points_[i,2] <= node1[2] & points_[i,2] >= node2[2]) ){
          result[i] = j
          break
        }
      }else if(abs((points_[i,2] - node1[2])/(node2[2] - node1[2]) - 
                   (points_[i,1] - node1[1])/(node2[1] - node1[1])) < 10^5 * .Machine$double.eps){
        if( (node1[1]<node2[1] & node1[1] <=points_[i,1] & points_[i,1] <= node2[1]) |
            (node2[1]<node1[1] & node2[1] <=points_[i,1] & points_[i,1] <= node1[1] )){
          result[i] = j
          break
        }
      }  
    }
  }
  return(result)
}

# set_region <- function(centroids, mesh, LN = chicago){
#   
#   lpp_centroids <- spatstat.linnet::lpp(X= ppp(centroids[,1], y=centroids[,2], window= LN$domain$window), 
#                                         L= LN$domain)
#   
#   nedges = nrow(mesh$edges)
#   mid_points = matrix(0, nrow=nedges, ncol=2)
#   for(e in 1:nedges){  
#     mid_points[e,1] =  ( mesh$nodes[ mesh$edges[e,1],1] + mesh$nodes[ mesh$edges[e,2],1] )/2
#     mid_points[e,2] =  ( mesh$nodes[ mesh$edges[e,1],2] + mesh$nodes[ mesh$edges[e,2],2] )/2
#   }
#   
#   lpp_midpoints <- spatstat.linnet::lpp(X = ppp(x=mid_points[,1], y= mid_points[,2], window= LN$domain$window),
#                                         L = LN$domain)
#   
#   network_dist = spatstat.linnet::crossdist.lpp(lpp_midpoints, lpp_centroids)
#   edges_to_region = matrix(0, nrow = nedges, ncol=1)
#   
#   for(e in 1:nedges){
#     idxs = which( network_dist[e,] == min(network_dist[e,]))
#     if(length(idxs) == 1) 
#       edges_to_region[e] = idxs
#     else
#       edges_to_region[e] = idxs[sample(1:length(idxs), size=1)]
#   }
#   
#   return(edges_to_region)
# }

# plot_region <-function(lines_to_region, 
#                        response,
#                        LN = chicago, 
#                        mesh = fdaPDE::create.mesh.1.5D(nodes=cbind(LN$domain$vertices$x, LN$domain$vertices$y),
#                                                        edges=cbind(LN$domain$from, LN$domain$to)),
#                        nregion = NULL){
#   
#   set.seed(NULL)
#   set.seed(314156)
#   
#   if(is.null(nregion)){
#     nregion = max(lines_to_region)
#     colors_ = jet.col(nregion)
#     colors_ = sample(colors_, nregion)
#     
#     nedges = nrow(mesh$edges)
#     plot(mesh, pch=".")
#     for(e in 1:nregion){
#       mask_ = which(lines_to_region == e)
#       segments(mesh$nodes[mesh$edges[mask_,1],1], mesh$nodes[mesh$edges[mask_,1],2],
#                mesh$nodes[mesh$edges[mask_,2],1],  mesh$nodes[mesh$edges[mask_,2],2],
#                col = colors_[e], lwd=4.5)
#       
#     }
#     legend("topright", legend=response, col=colors_, lwd=7)
#   }else{
#     colors_ = jet.col(nregion)
#     colors_ = sample(colors_, nregion)
#     
#     nedges = nrow(mesh$edges)
#     plot(mesh, pch=".")
#     mask_ = which(lines_to_region == nregion)
#     segments(mesh$nodes[mesh$edges[mask_,1],1], mesh$nodes[mesh$edges[mask_,1],2],
#              mesh$nodes[mesh$edges[mask_,2],1],  mesh$nodes[mesh$edges[mask_,2],2],
#              col = "red", lwd=4.5)
#     legend("topright", legend=response[nregion], col="red", lwd=7)
#     
#   }
# }
# 
# plot_region_gradient_color <- function(mesh,
#                                        response,
#                                        lines_to_region,
#                                        palette = jet.col,
#                                        line.size= 1,
#                                        title.size = 26){
#   
#   x=vector(mode="double")
#   y=vector(mode="double")
#   grp.nodes=vector(mode="integer")
#   
#   num_edges= dim(mesh$edges)[1]
#   for(e in 1:num_edges){
#     x = append(x, c(mesh$nodes[ mesh$edges[e,1], 1], mesh$nodes[ mesh$edges[e,2], 1]))
#     y = append(y, c(mesh$nodes[ mesh$edges[e,1], 2], mesh$nodes[ mesh$edges[e,2], 2]))
#     grp.nodes = append(grp.nodes, rep(e,times=2))
#   }
#   
#   MyTheme <- theme(
#     axis.text = element_text(size=(title.size-2)),
#     axis.title = element_text(size=title.size),
#     title = element_text(size=title.size),
#     legend.text = element_text(size=(title.size-6)),
#     legend.key.size = unit(1,"cm"),
#     legend.position = "right")
#   
#   num.col = diff(range(response))
#   p = palette(n=num.col)
#   color.min=min(response) 
#   color.max=max(response)
#   
#   mu = vector(mode="double", length= length(lines_to_region)*2)
#   for(i in 1:length(lines_to_region)){
#     mu[(2*i-1)] = response[lines_to_region[i]]
#     mu[(2*i)] = response[lines_to_region[i]]
#   }
#   
#   ratio = diff(range(mesh$nodes[,1]))/diff(range(mesh$nodes[,2]))
#   
#   data=data.frame(x,y,grp.nodes)
#   
#   ggplot(data=NULL) + 
#     geom_point(data=data, aes(x=x,y=y,group=grp.nodes), 
#                alpha=0.0) + 
#     geom_line(data=data, aes(x=x,y=y,group=grp.nodes, color= mu), 
#               size=line.size) +
#     labs(x="",y="",color="") + 
#     scale_color_gradientn(colours=p, limits = c(color.min, color.max))+ 
#     coord_fixed(ratio=ratio) + 
#     theme_void() +
#     MyTheme +   
#     theme(plot.title = element_text(hjust=0.5),
#           legend.title = element_blank(),
#           axis.title = element_blank(),
#           axis.text.x=element_blank(),
#           axis.text.y=element_blank(),
#           legend.key.size = unit(1,"cm"),
#           legend.key.width = unit(0.5,"cm"))
#   
# }
# 
# plot_nodes<-function(mesh){
#   
#   pdf("nodes.pdf")
#   for(i in 1:nrow(mesh$nodes)){
#     print(plot(mesh,pch=".", main=paste("node ",i,sep="")))
#     print(points(mesh$nodes[i,1],mesh$nodes[i,2], pch=16, col="red3"))
#   }
#   dev.off()
#   }
