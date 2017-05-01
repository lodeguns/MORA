library(GGally)
library(sda)
library(ggplot2)
library(network)
library(graph)
library(igraph)


build.random.graph <- function(names.g, seed.val, M, p.val)
{
  set.seed(seed.val)
  g1 <- randomGraph(names.g, 1:M, p=p.val)
  g1.igraph <- igraph.from.graphNEL(g1)
  return(g1.igraph)
}

lt <- function(x,...) {x[lower.tri(x,...)]}

calc.mean.shortest.path <- function(graph.ig, dir=TRUE)
{
  
  return(mean_distance(graph.ig, directed = TRUE, unconnected = TRUE))
  
}

calc.mean.cluster.coef <- function(g.nel)
{
  Ceach <- clusteringCoefficient(g.nel, selfLoops=TRUE)
  str <- paste(Ceach,collapse="-")
  #print(paste("CLustering coef::", str, sep=" "))
  Ceach <- Ceach *2 #because it is oriented
  Chat.g1 <- sum(Ceach)/(length(Ceach))
  return(Chat.g1)
}

conv.igraph.to.grapNEL <- function(g.ig)
{ #nota altri modi non sono supportati.
  g2.mat = as.matrix(as_adj(g.ig, attr = NULL))
  g2.mat <- (g2.mat +t(g2.mat )) > pi/4
  rownames(g2.mat ) <- colnames(g2.mat ) <- V(g.ig)$name
  g2.NEL <- as(g2.mat,"graphNEL")
  return(g2.NEL)
  
}

plot.alt.graph <- function(g){
  plot.igraph(g,vertex.label=names(count),
              layout=layout.circle, 
              vertex.label.color="black",
              edge.color="black",
              edge.width=count, 
              edge.arrow.size=0.5,
              edge.curved=FALSE,
              edge.label = E(g)$ntewcos,
              edge.label.color ="blue",
              vertex.size=30,
              vertex.label.font=1,
              vertex.color="aliceblue")
}

path.gene.names <- gsub("eco:", "" , V(list.kegg.path.igraph[[1]])$name)

##################################################################################################
####
####   MORA (multi-omic relational adjacencies)
##################################################################################################
get.pattern.adjacent.influences <- function(pattern, delta, psi, g){
  #delta     <- 2                           #minimum distance delta=2, 2 element of the pattern in neighborhood 1
  #the adjacent ones.
  couple    <- length(pattern) - (delta-2) #max number of couple at this distance: x-0-1-x-1-0-x-1-0-x
  dir.edges <- "all"                       #outgoing edges:out, ingoing edges:in
  #psi       <- 2:(ceiling(calc.mean.shortest.path(g)))    #max number of vertexes in a path: O-O-O---O 
  #from 0 to the average pathlength+1.
  count <- c(rep(0, length(pattern)))
  names(count)<- as.character((pattern))
  count <- c(count,count[1:max(psi)-1])
  
  pattern <- c(pattern, pattern[1:max(psi)-1])
  
  
  #psi       <- 2
  
  for(y in 1:length(psi)){
    for(i in 1:couple){
      p.from <- pattern[i]
      # therefore the i-th element is at distance i+(delta-1)
      p.to   <- pattern[i + (delta-1)] 
     # print(paste("", p.from, " ", p.to))
      patty=NULL
      if(distances(g, v = p.from , to = p.to , mode = dir.edges)[[1]] != Inf)
      {tryCatch( patty <-get.shortest.paths(g, p.from , p.to , mode = dir.edges,
                                            weights = NULL, 
                                            output=c("vpath", "epath", "both")) ,
                 warning = function(w) {print(paste("warning... "));},
                 error   = function(e) {print(paste("no path. "));  }
      )}
      # Vertexes Sequence
      # exstract the sequence of vertex from the path,
      # check if the path is of length delta-1
      # if it is of delta-1 save the occurences with +1 in the vector
      # psi = delta
      if(length(patty)!=0){
        for( k in 1:length(patty$vpath)){
          if(length(V(g)[patty$vpath[[k]]]$name)==psi[[y]]){
            pr = paste("Pattern elements delta = ", delta-1,
                       ".  Vertex psi neighbouring = ", psi[[y]]-1, "edges.", sep =" ")
            # pr = paste(pr, " Shortest path's end nodes: ", names(count)[i], " and ", names(count)[i + (delta-1)], sep=" ")
            #print(pr)
            namenodes<- V(g)[patty$vpath[[k]]]$name
            for( z in 1:length(namenodes)){
              zx<-grep(paste("\\b",namenodes[z], "\\b", sep=""), names(count))
              print(paste("", namenodes))
              count[c(zx)] <- (count[c(zx)]+(1/(psi[[y]]-1)))
              print(paste("", count[c(zx)]))
            }
          }
        }
      } 
    }
  } 
  
  count<- count[!duplicated(names(count))]   #stable
  
  return(count)
}

# 
# library(igraph)
# from = c("1", "2", "3","4","5","6","7","10","9","11","12","13","8")
# df.g <- data.frame(
#   from = c("1","1", "2", "2","3","4","5","5","7","7","7","7","7","1",
#            "13","8","12","11","6","9", "11","6","6","6", "8","12",
#            "11","6","9","7","1","1","3","3","3","3", "6","12","7","6", "7", "8", "9") , 
#   to =   c("13","2", "4", "11", "1","3","7","8","5","8","10","10","1",
#            "7","1","1","3","3","3","3", "6","12","7","9","4","10",
#            "13","11","7","13","8","12","11","6","9", "11","6","6","6","7", "8", "9", "11")
# )
# 
# g <- graph_from_data_frame(df.g, directed=TRUE, vertices=from)
# print(g, e=TRUE, v=TRUE)
# calc.mean.cluster.coef(conv.igraph.to.grapNEL(g))
# calc.mean.shortest.path(g)
# 
# ## Ordered pattern
# pattern = c("3", "1", "2", "4","5","6","7","8","10","11","13","9", "12")
# 
# count<-get.pattern.adjacent.influences(pattern, delta=2, psi=2:(ceiling(calc.mean.shortest.path(g))), g)
# 
# ## Form an array of zero values for the pattern
# #foresterno
# m=as_adjacency_matrix(g)
# m=as.matrix(m)
# g1 <- network(m, directed=TRUE)
# g1 %v% "structural influence" = ifelse(count >= median(count), "more adjacent", "less adjacet")
# g1 %v% "adjacency weight" = count
# #netg <- as.network(ig)
# #ig2 <- as.igraph(netg)
# ggnet2(g1, label = TRUE, label.alpha = 0.95,  arrow.size = 7, 
#        arrow.gap = 0.03, color = "structural influence", 
#        palette = "Set2",  label.color = "black", label.size = 3,
#        size = "adjacency weight",  edge.color = c("color", "grey50"))

###########################################################################################
###########################################################################################
###########################################################################################
# library(igraph)
# from = c("1", "2", "3","4","5","6","7","10","9","11","12","13","8")
# df.g <- data.frame(
#   from = c("1","1", "2", "2","3","4","5","5","7","7","7","7","7","1",
#            "13","8","12","11","6","9", "11","6","6","6", "8","12",
#            "11","6","9","7","1","1","3","3","3","3", "6","12","7","6", "7", "8", "9") , 
#   to =   c("13","2", "4", "11", "1","3","7","8","5","8","10","10","1",
#            "7","1","1","3","3","3","3", "6","12","7","9","4","10",
#            "13","11","7","13","8","12","11","6","9", "11","6","6","6","7", "8", "9", "11")
# )
# 
# global.graph <- graph_from_data_frame(df.g, directed=TRUE, vertices=from)
# print(global.graph, e=TRUE, v=TRUE)
# ##plot
# m=as_adjacency_matrix(global.graph)
# m=as.matrix(m)
# g1 <- network(m, directed=TRUE)
# g1 %v% "sub.network" = ifelse(names(global.pattern) %in% names(pathway_pattern), 0.8, 0.2)
# g1 %v% "multi-omic alternation" = global.pattern[order(as.numeric(names(global.pattern)))]
# #netg <- as.network(ig)
# #ig2 <- as.igraph(netg)
# ggnet2(g1, label = TRUE, label.alpha = 0.95,  arrow.size = 7, 
#        arrow.gap = 0.03, color = "multi-omic alternation", 
#        palette = "Set1",  label.color = "black", label.size = 3,  alpha = "sub.network", 
#        #  size = "sub.network",
#        edge.color = c("color", "grey50")
#        , legend.position = "bottom")
# 
# 
# ## Input data.
# ## Global ordered pattern
# global.pattern      <- c("3", "1", "2", "4","5","6","7","8","10","11","13","9","12")
# pattern_val         <- c( 1,   0,   1,   0,  1,  1,  0,  1,   1,   1,  0 ,  1,  1  )
# names(pattern_val)  <- global.pattern
# global.pattern      <- pattern_val
# ## Nodes of pathway (individuation sub-graph)
# pathway.gene.list   <- c("1","2","3","9","7","10","12")


get.multi.omic.alternation <- function(global.pattern, pathway.gene.list, global.graph, N=2){
  ## Multi-omic array and multi-omic sub-graph.
  pathway_pattern     <- global.pattern[pathway.gene.list ]
  pathway_pat_g       <- induced.subgraph(graph=global.graph, vids=pathway.gene.list )
  ##Average path length, clustering coefficient for the sub-graph
  cc.path    <- calc.mean.cluster.coef(conv.igraph.to.grapNEL(pathway_pat_g ))
  apl.path     <- calc.mean.shortest.path(pathway_pat_g )
  ## Compute alternations
  alternation <- compute.alternations(pathway_pattern, N)
  #alternation <- compute.alternations(global.pattern, N)
  ## Compute deviations with the best possible score.
  deviations.frame <- compute.deviations(alternation, global.graph, global.pattern,  N=2)
  ## Compute these extensions
  compute.extension.by.deviations(deviations.frame, global.pattern, pathway_pattern)
  ## Compute alternation on the network
  network.alternations.m(pathway_pat_g, pathway_pattern) 
  
}




compute.alternations <- function(pathway_pattern,N=2)
{
  N           <- 2
  l           <- length(pathway_pattern)
  divisor     <- c()
  alternation <- c()
  for(i in 1:(l-1))
  { divisor <- c(divisor, paste(names(pathway_pattern)[i], names(pathway_pattern)[i+1], sep="|"))
  alternation <- c(alternation , abs(pathway_pattern[i] - pathway_pattern[i+1]))
  }
  
  names(alternation) <- divisor
  alternation
}

#path extensions
compute.deviations <- function(alternation, global.graph, global.pattern.net,  N=2){

  deviations         <- list()
  deviations.frame   <- data.frame()
  names.deviation <- c()
  for(j in 1:length(alternation))
  {
    if(alternation[j]==0)
    { 
      fromto <- unlist(strsplit( names(alternation[j]), "|", fixed=TRUE))
      from   <- fromto[1]
      to     <- fromto[2]
      
      s.p <- all_shortest_paths(global.graph , from, to = to, mode = c("out", "all", "in"),
                                weights = NULL)
      if(length(s.p$res)!=0)
        for(i in 1:length(s.p$res))
        { 
          if(length(s.p$res[[i]]$name)!=2){
            # print(paste("--(",s.p$res[[i]]$name,")--" ))
            names.altern     <- names(global.pattern.net[s.p$res[[i]]$name])
            altern.sub       <- as.vector(global.pattern.net[s.p$res[[i]]$name])
            names(altern.sub ) <-  names.altern 
            deviations[[length(deviations)+1]] <- alternation.scores(j,i, altern.sub, N)
            #print(paste("-------", j, "--- ", i ))
            #print(deviations.frame)
            df   <- as.data.frame(deviations[[length(deviations)]])
            nn.c <-  as.data.frame(deviations[[length(deviations)]]$alter, 
                                   stringsAsFactors=TRUE)
            names.deviation <- c(names.deviation, rownames(nn.c))
            deviations.frame <- rbind(deviations.frame,  as.data.frame(df))
          }
        }
      
      
    }
    
  }
  
  deviations.frame$alter.name <- names.deviation
  return(deviations.frame)
}

compute.extension.by.deviations <- function(deviations.frame, global.pattern, local.pattern)
{ 
  
  
  if(nrow(deviations.frame)!=0){
  dev.list <- unique(deviations.frame$j)
  
  for(j in 1:length(dev.list))
  {
    dev.el   <- deviations.frame[deviations.frame$j == dev.list[j],]
    dev.el   <- dev.el[with(dev.el, order(sim.sc, decreasing = TRUE)), ] #high similarity
    dev.el.i <- dev.el[dev.el$i == dev.el$i[1],]
    dev.el.u <- unique(unlist(strsplit( dev.el.i$alter.name, "|", fixed=TRUE)))
    from     <- dev.el.u[[1]]
    to       <- dev.el.u[[length(dev.el.u)]]
    p.v      <- global.pattern[dev.el.u]
    from.g   <- grep(paste("^",from, "\\b", sep=""), names(local.pattern))
    to.g     <- grep(paste("^",to  , "\\b", sep=""), names(local.pattern))
    local.pattern_ls <- local.pattern[c(1: from.g)]
    local.pattern_rx <- local.pattern[c(to.g:length(local.pattern))]
    names(p.v) <- paste("dev[",names(p.v), "]", sep="")
    p.v         <- p.v[-1]
    p.v         <- p.v[1:length(p.v)-1]
    p.t                <- c(local.pattern_ls, p.v, local.pattern_rx)
    p.t                <- p.t[!is.na(p.t)]
    local.pattern      <- p.t
  }
  }
  return(local.pattern)
}

alternation.scores <- function(j=0,k=0, altern, N=2)
{
  l.sc   <- length(altern)-1  ##l-1
  divisor     <- c()
  alter       <- c()
  for(i in 1:(l.sc))
  { divisor <- c(divisor, paste(names(altern)[i], names(altern)[i+1], sep="|"))
  alter   <- c(alter  , abs(altern[i] - altern[i+1]))
  }
  names(alter) <- divisor
  w.s.id  <- 1/N
  a.s <- sum(alter * w.s.id) 
  w.s <- a.s/l.sc
  dis.sc <- (w.s.id - w.s) * N
  sim.sc <- 1 - dis.sc
  
  list(j=j, i=k, alter=alter, sim.sc = sim.sc , dis.sc = dis.sc, 
       abs.sc = w.s , rel.sc=a.s, id.sc= w.s.id, l.sc =l.sc)
}


#pathway_pattern        <- effective_pattern_path  
network.alternations.m <- function(pathway_pattern, global.network)
{
  names.pathway.net          <- names(pathway_pattern)
  pathway_network            <- induced.subgraph(graph=global.network, vids=names.pathway.net )
  
  E.M <- length(E(pathway_network))
  V.N <- length(V(pathway_network))
  delta.p <- E.M/(V.N*(V.N-1))
  
  n1 <- length(pathway_pattern[pathway_pattern==1])
  n0 <- length(pathway_pattern[pathway_pattern==0])
  
  E.alter     <- n1*n0*delta.p
  #E.no_alter  <- (((n1*(n1-1))/2)*delta.p + ((n0*(n0-1))/2)*delta.p )
  E.property1 <- ((n1*(n1-1))/2)*delta.p
  E.property0 <- ((n0*(n0-1))/2)*delta.p 
  
  edges <- E(pathway_network)
  count_alter     <- 0
  count_not_alter_p1 <- 0
  count_not_alter_p0 <- 0
  for(e in 1:length(edges)){
    get.edg <- paste(V(pathway_network)$name[get.edges(pathway_network, E(pathway_network)[e])])
    dyad    <- pathway_pattern[get.edg]
    
    if(alternation.scores(altern=dyad)$alter == 1)
    {
      count_alter <- count_alter +1    #m10,m01
    } else
    {
      if(sum(dyad) == 2)
        count_not_alter_p1 <- count_not_alter_p1 +1  #m00,m11
      
      if(sum(dyad)==  0)
        count_not_alter_p0 <- count_not_alter_p0 +1
    }
    
  }
  
  alter.index <- count_alter/E.alter
  not.alter.index <- ((count_not_alter_p1/E.property1)  + (count_not_alter_p0/E.property0))/2
  count_not_alter <- count_not_alter_p0 + count_not_alter_p1
  
  alter.test <- ifelse(alter.index >= 1, T, F)
  no.alter.test <- ifelse(not.alter.index  >= 1, T, F)
  E.no_alter <- E.property1 + E.property0
  list(alter.test= alter.test,
       no.alter.test = no.alter.test,
       alter.index = alter.index,
       no.alter.index= not.alter.index, 
       m.alter = count_alter,
       m.no.alter =count_not_alter, E.alter = E.alter, E.no_alter=E.no_alter, delta.p = delta.p, E.property0=E.property0,E.property1=E.property1)
}

###########################################################################################
###########################################################################################
### v   BEST CASE  : fully alternating on array,    fully alternating on network
###    Worst cases: fully alternating on array,    low alternating on network
###    Worst cases: fully alternating on network , low alternating on array.
###########################################################################################
# from = c("1", "2", "3","4","5")
# df.g <- data.frame(
#   from = c("1","2","3","4") , 
#   to =   c("2","3","4","5")
# )
# 
# ## Input data.
# ## Global ordered pattern
# global.pattern      <- c("1", "2", "3", "4","5")
# pattern_val         <- c( 1,   0,   1,   0,  1)
# names(pattern_val)  <- global.pattern
# global.pattern      <- pattern_val
# ## Nodes of pathway (individuation sub-graph)
# pathway.gene.list   <- c("1","2","3","4","5")
# pathway_pattern     <- global.pattern[pathway.gene.list ]
# pathway_pat_g       <- induced.subgraph(graph=global.graph, vids=pathway.gene.list )
# 
# global.graph <- graph_from_data_frame(df.g, directed=TRUE, vertices=from)
# print(global.graph, e=TRUE, v=TRUE)
# ##plot
# m=as_adjacency_matrix(global.graph)
# m=as.matrix(m)
# g1 <- network(m, directed=TRUE)
# g1 %v% "sub.network" = ifelse(names(global.pattern) %in% names(pathway_pattern), 0.8, 0.2)
# g1 %v% "multi-omic alternation" = global.pattern[order(as.numeric(names(global.pattern)))]
# p1 <- ggnet2(g1, label = TRUE, label.alpha = 0.95,  arrow.size = 7, 
#              arrow.gap = 0.03, color = "multi-omic alternation", 
#              palette = "Set1",  label.color = "black", label.size = 3,  alpha = "sub.network", 
#              #  size = "sub.network",
#              edge.color = c("color", "grey50")
#              , legend.position = "bottom")
# 
# cc.path    <- calc.mean.cluster.coef(conv.igraph.to.grapNEL(pathway_pat_g ))
# apl.path     <- calc.mean.shortest.path(pathway_pat_g )
# ## Compute alternations
# alternation <- compute.alternations(pathway_pattern, N)
# #alternation <- compute.alternations(global.pattern, N)
# ## Compute deviations with the best possible score.
# #deviations.frame <- compute.deviations(alternation, global.graph, global.pattern,  N=2)
# ## Compute these extensions
# #compute.extension.by.deviations(deviations.frame, global.pattern, pathway_pattern)
# ## Compute alternation on the network
# network.alternations.m(pathway_pat_g, pathway_pattern) 
# alternation.scores(0,0, pathway_pattern, 2)
# 
# nn.p <- names(global.pattern)
# count<-get.pattern.adjacent.influences( nn.p, delta=2, 
#                                         psi=2:(ceiling(calc.mean.shortest.path(global.graph))),
#                                         global.graph)
# 
# ## Form an array of zero values for the pattern
# #foresterno
# m=as_adjacency_matrix(global.graph)
# m=as.matrix(m)
# g1 <- network(m, directed=TRUE)
# g1 %v% "structural influence" = ifelse(count >= median(count), "more adjacent", "less adjacet")
# g1 %v% "adjacency weight" = count
# p2 <- ggnet2(g1, label = TRUE, label.alpha = 0.95,  arrow.size = 7, 
#              arrow.gap = 0.03, color = "structural influence", 
#              palette = "Set2",  label.color = "black", label.size = 3,
#              size = "adjacency weight",  edge.color = c("color", "grey50"))
# 
# 
# multiplot(p1,p2, cols = 2)
# 
# ###########################################################################################
# ###########################################################################################
# ###    BEST CASE  : fully alternating on array,    fully alternating on network
# ###   v Worst cases: fully alternating on network , low alternating on array.
# ###    Worst cases: fully alternating on array,    low alternating on network
# ###########################################################################################
# from = c("1", "2", "3","4","5")
# df.g <- data.frame(
#   from = c("1","2","3","4") , 
#   to =   c("2","3","4","5")
# )
# 
# ## Input data.
# ## Global ordered pattern
# global.pattern      <- c("1", "2", "3", "4","5")
# pattern_val         <- c( 1,   0,   1,   0,  1)
# names(pattern_val)  <- global.pattern
# global.pattern      <- pattern_val
# ## Nodes of pathway (individuation sub-graph)
# pathway.gene.list   <- c("1","3","5","2","4")
# pathway_pattern     <- global.pattern[pathway.gene.list ]
# pathway_pat_g       <- induced.subgraph(graph=global.graph, vids=pathway.gene.list )
# 
# global.graph <- graph_from_data_frame(df.g, directed=TRUE, vertices=from)
# print(global.graph, e=TRUE, v=TRUE)
# ##plot
# m=as_adjacency_matrix(global.graph)
# m=as.matrix(m)
# g1 <- network(m, directed=TRUE)
# g1 %v% "sub.network" = ifelse(names(global.pattern) %in% names(pathway_pattern), 0.8, 0.2)
# g1 %v% "multi-omic alternation" = global.pattern[order(as.numeric(names(global.pattern)))]
# p1<-ggnet2(g1, label = TRUE, label.alpha = 0.95,  arrow.size = 7, 
#            arrow.gap = 0.03, color = "multi-omic alternation", 
#            palette = "Set1",  label.color = "black", label.size = 3,  alpha = "sub.network", 
#            #  size = "sub.network",
#            edge.color = c("color", "grey50")
#            , legend.position = "bottom")
# 
# cc.path    <- calc.mean.cluster.coef(conv.igraph.to.grapNEL(pathway_pat_g ))
# apl.path     <- calc.mean.shortest.path(pathway_pat_g )
# ## Compute alternations
# alternation <- compute.alternations(pathway_pattern, N)
# #alternation <- compute.alternations(global.pattern, N)
# ## Compute deviations with the best possible score.
# #deviations.frame <- compute.deviations(alternation, global.graph, global.pattern,  N=2)
# ## Compute these extensions
# #compute.extension.by.deviations(deviations.frame, global.pattern, pathway_pattern)
# ## Compute alternation on the network
# network.alternations.m(pathway_pat_g, pathway_pattern) 
# alternation.scores(0,0, pathway_pattern, 2)
# 
# nn.p <- names(global.pattern)
# count<-get.pattern.adjacent.influences( nn.p, delta=2, 
#                                         psi=2:(ceiling(calc.mean.shortest.path(global.graph))),
#                                         global.graph)
# 
# ## Form an array of zero values for the pattern
# #foresterno
# m=as_adjacency_matrix(global.graph)
# m=as.matrix(m)
# g1 <- network(m, directed=TRUE)
# g1 %v% "structural influence" = ifelse(count >= median(count), "more adjacent", "less adjacet")
# g1 %v% "adjacency weight" = count
# #netg <- as.network(ig)
# #ig2 <- as.igraph(netg)
# p2<-ggnet2(g1, label = TRUE, label.alpha = 0.95,  arrow.size = 7, 
#            arrow.gap = 0.03, color = "structural influence", 
#            palette = "Set2",  label.color = "black", label.size = 3,
#            size = "adjacency weight",  edge.color = c("color", "grey50"))
# 
# 
# multiplot(p1,p6, cols = 2)
# 
# ###########################################################################################
# ###########################################################################################
# ###    BEST CASE  : fully alternating on array,    fully alternating on network
# ###    Worst cases: fully alternating on network , low alternating on array.
# ###   v Worst cases: fully alternating on array,    low alternating on network
# ###########################################################################################
# from = c("1", "2", "3","4","5")
# df.g <- data.frame(
#   from = c("1","2","3","4") , 
#   to =   c("2","3","4","5")
# )
# 
# ## Input data.
# ## Global ordered pattern
# global.pattern      <- c("1", "2", "3", "4","5")
# pattern_val         <- c( 1,   1,   1,   0,  0)
# names(pattern_val)  <- global.pattern
# global.pattern      <- pattern_val
# ## Nodes of pathway (individuation sub-graph)
# pathway.gene.list   <- c("1","4","2","5","3")
# pathway_pattern     <- global.pattern[pathway.gene.list ]
# pathway_pat_g       <- induced.subgraph(graph=global.graph, vids=pathway.gene.list )
# 
# global.graph <- graph_from_data_frame(df.g, directed=TRUE, vertices=from)
# print(global.graph, e=TRUE, v=TRUE)
# ##plot
# m=as_adjacency_matrix(global.graph)
# m=as.matrix(m)
# g1 <- network(m, directed=TRUE)
# g1 %v% "sub.network" = ifelse(names(global.pattern) %in% names(pathway_pattern), 0.8, 0.2)
# g1 %v% "multi-omic alternation" = global.pattern[order(as.numeric(names(global.pattern)))]
# p1<-ggnet2(g1, label = TRUE, label.alpha = 0.95,  arrow.size = 7, 
#            arrow.gap = 0.03, color = "multi-omic alternation", 
#            palette = "Set1",  label.color = "black", label.size = 3,  alpha = "sub.network", 
#            #  size = "sub.network",
#            edge.color = c("color", "grey50")
#            , legend.position = "bottom")
# 
# cc.path    <- calc.mean.cluster.coef(conv.igraph.to.grapNEL(pathway_pat_g ))
# apl.path     <- calc.mean.shortest.path(pathway_pat_g )
# ## Compute alternations
# alternation <- compute.alternations(pathway_pattern, N)
# #alternation <- compute.alternations(global.pattern, N)
# ## Compute deviations with the best possible score.
# #deviations.frame <- compute.deviations(alternation, global.graph, global.pattern,  N=2)
# ## Compute these extensions
# #compute.extension.by.deviations(deviations.frame, global.pattern, pathway_pattern)
# ## Compute alternation on the network
# network.alternations.m(pathway_pat_g, pathway_pattern) 
# alternation.scores(0,0, pathway_pattern, 2)
# 
# nn.p <- names(global.pattern)
# count<-get.pattern.adjacent.influences( nn.p, delta=2, 
#                                         psi=2:(ceiling(calc.mean.shortest.path(global.graph))),
#                                         global.graph)
# 
# ## Form an array of zero values for the pattern
# #foresterno
# m=as_adjacency_matrix(global.graph)
# m=as.matrix(m)
# g1 <- network(m, directed=TRUE)
# g1 %v% "structural influence" = ifelse(count >= median(count), "more adjacent", "less adjacet")
# g1 %v% "adjacency weight" = count**2
# #netg <- as.network(ig)
# #ig2 <- as.igraph(netg)
# p2<-ggnet2(g1, label = TRUE, label.alpha = 0.95,  arrow.size = 7, 
#            arrow.gap = 0.03, color = "structural influence", 
#            palette = "Set2",  label.color = "black", label.size = 3,
#            size = "adjacency weight",  edge.color = c("color", "grey50"))
# 
# multiplot(p1,p2, cols = 2)

###################################################################################




















# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}