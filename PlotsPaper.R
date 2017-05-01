###################################################################################
## Plot alternations with path extensions or without, influences etc
##  
###################################################################################
##
##  Plot influences
##
library(GGally)
library(sda)
library(ggplot2)
library(network)
library(graph)
library(igraph)

from = c("3", "1", "2", "4","5","6","7","8","10","11","13","9", "12")
df.g <- data.frame(
  from = c("1","1", "2", "2","3","4","5","5","7","7","7","7","7","1",
           "13","8","12","11","6","9", "11","6","6","6", "8","12",
           "11","6","9","7","1","1","3","3","3","3", "6","12","7","6", "7", "8", "9") ,
  to =   c("13","2", "4", "11", "1","3","7","8","5","8","10","10","1",
           "7","1","1","3","3","3","3", "6","12","7","9","4","10",
           "13","11","7","13","8","12","11","6","9", "11","6","6","6","7", "8", "9", "11")
)

g <- graph_from_data_frame(df.g, directed=TRUE, vertices=from)
print(g, e=TRUE, v=TRUE)
calc.mean.cluster.coef(conv.igraph.to.grapNEL(g))
calc.mean.shortest.path(g)

## Ordered pattern
pattern = c("3", "1", "2", "4","5","6","7","8","10","11","13","9", "12")

count<-get.pattern.adjacent.influences(pattern, delta=2, psi=2:(ceiling(calc.mean.shortest.path(g))), g)

## Form an array of zero values for the pattern
#foresterno
m=as_adjacency_matrix(g)
m=as.matrix(m)
g1 <- network(m, directed=TRUE)
g1 %v% "structural influence" = ifelse(as.numeric(count) >= median(as.numeric(count)), "more adjacent", "less adjacet")
g1 %v% "adjacency weight" = as.numeric(count)
#netg <- as.network(ig)
#ig2 <- as.igraph(netg)
g.plot.1 <- ggnet2(g1, label = TRUE, label.alpha = 0.95,  arrow.size = 7,
       arrow.gap = 0.03, color = "structural influence",
       palette = "Set2",  label.color = "black", label.size = 3,
       size = "adjacency weight",  edge.color = c("color", "grey50"))




####################################################################################
## Case base significance




from = c("4","5","2", "3", "1")
df.g <- data.frame(
  from = c("1","2","3","3","4", "5") ,
  to =   c("2","1","1","5","3", "3")
)

g <- graph_from_data_frame(df.g, directed=TRUE, vertices=from)
print(g, e=TRUE, v=TRUE)
calc.mean.cluster.coef(conv.igraph.to.grapNEL(g))
calc.mean.shortest.path(g)

## Ordered pattern
pattern = c("4","5","2", "3", "1")

count<-get.pattern.adjacent.influences(pattern, delta=2, psi=2:(ceiling(calc.mean.shortest.path(g))), g)

## Form an array of zero values for the pattern
#foresterno
m=as_adjacency_matrix(g)
m=as.matrix(m)
g1 <- network(m, directed=TRUE)
g1 %v% "structural influence" = ifelse(as.numeric(count) >= median(as.numeric(count)), "more adjacent", "less adjacet")
g1 %v% "adjacency weight" = as.numeric(count)
#netg <- as.network(ig)
#ig2 <- as.igraph(netg)
ggnet2(g1, label = TRUE, label.alpha = 0.95,  arrow.size = 7,
                   arrow.gap = 0.03, color = "structural influence",
                   palette = "Set2",  label.color = "black", label.size = 3,
                   size = "adjacency weight",  edge.color = c("color", "grey50"))








#####################################################################################
## Network alternations and path extensions
#####################################################################################
load("./gene_assoc.Rdata")
load("./list.kegg.path.igraph.Rdata")
load("./ecocyc.kegg.igraph.Rdata")
load("./computations.steady.state.Rdata")
load("./eco.df.paths.Rdata")
load("./faith1.Rdata")
global.pattern.names      <- sub("eco:", "", V(ecocyc.kegg.igraph)$name)
global.pattern            <- gene_assoc[global.pattern.names,]
global.pattern            <- global.pattern[!is.na(global.pattern$gene_loc),]
names.global.net          <- paste("eco:", global.pattern$gene_loc, sep="")
global.network            <- induced.subgraph(graph=ecocyc.kegg.igraph, vids=names.global.net )
global.pattern_names      <- as.vector(global.pattern$gene_loc)
global.pattern.net        <- as.vector(global.pattern$list.caipa)
names(global.pattern.net) <- global.pattern_names 
global.network            <- set.vertex.attribute(global.network , "name", value=global.pattern_names)
APL                       <- ceiling(calc.mean.shortest.path(global.network ))
ll1                       <- steady.state[[1]]


g  <- induced.subgraph(graph=ecocyc.kegg.igraph, vids= paste("eco:",names(ll1$multi.omic.pattern), sep=""))
V(g)$name <- b.to.sym.orderd(gene_assoc, gsub("eco:","",V(g)$name))
pattern = b.to.sym.orderd(gene_assoc, names(ll1$multi.omic.pattern))
count   = ll1$multi.omic.df[4,]
names(count) <- b.to.sym.orderd(gene_assoc, names(count))
m=as_adjacency_matrix(g)
m=as.matrix(m)
g1 <- network(m, directed=TRUE)
g1 %v% "structural influence" = ifelse(count >= median(count), "more adjacent", "less adjacet")
g1 %v% "adjacency weight" = count
p1 <- ggnet2(g1, label = TRUE, label.alpha = 0.95,  arrow.size = 7, 
       arrow.gap = 0.03, color = "structural influence", 
       palette = "Set2",  label.color = "black", label.size = 3,
       size = "adjacency weight",  edge.color = c("color", "grey50"))

#######################################################################################
## Multi-omics and path extensions
#######################################################################################
#"mhpA" "galU" "trpA" "gatZ" "tktB" "tktA" "garK" "qorB"
glb  <- induced.subgraph(graph = ecocyc.kegg.igraph, 
                         vids  = paste("eco:",names(ll1$multi.omic.pattern.c), sep=""))
V(glb)$name <- b.to.sym.orderd(gene_assoc, gsub("eco:","",V(glb)$name))
not.deviated <- b.to.sym.orderd(gene_assoc, gsub("eco:","",names(ll1$multi.omic.pattern)))
m=as_adjacency_matrix(glb)
m=as.matrix(m)
g1 <- network(m, directed=TRUE)
g1 %v% "deviations" = ifelse(V(glb)$name %in% not.deviated, 0.3, 0.8) #deviatons = extensions
g1 %v% "multi-omic alternation" = ll1$multi.omic.pattern.c
g1 %v% "structural influence" = ifelse(ll1$multi.omic.df.dev[4,] >= median(ll1$multi.omic.df.dev[4,]), "more adjacent", "less adjacet")

#netg <- as.network(ig)
#ig2 <- as.igraph(netg)
p3 <- ggnet2(g1, label = TRUE, label.alpha = 0.95,  arrow.size = 3, 
       arrow.gap = 0.03, color = "multi-omic alternation", 
       palette = "Set1",  label.color = "black", label.size = 3,  alpha = "deviations", 
       #size = "deviations",
       size = "structural influence",
      # edge.alpha = "multi-omic alternation",
       edge.color = c("color", "orange")
       , legend.position = "bottom")
p3 <- p3 + guides(color = FALSE, size = FALSE)
glb  <- induced.subgraph(graph = ecocyc.kegg.igraph, 
                         vids  = paste("eco:",names(ll1$multi.omic.pattern), sep=""))
V(glb)$name <- b.to.sym.orderd(gene_assoc, gsub("eco:","",V(glb)$name))
m=as_adjacency_matrix(glb)
m=as.matrix(m)
g1 <- network(m, directed=TRUE)
g1 %v% "deviations" = ifelse(names(ll1$multi.omic.pattern) 
                             %in% names(ll1$multi.omic.pattern), 0.3, 0.8)
g1 %v% "multi-omic alternation" = ll1$multi.omic.pattern
g1 %v% "structural influence" = ifelse(ll1$multi.omic.df[4,] >= median(ll1$multi.omic.df[4,]), "more adjacent", "less adjacet")
#netg <- as.network(ig)
#ig2 <- as.igraph(netg)
p4 <- ggnet2(g1, label = TRUE, label.alpha = 0.95,  arrow.size = 3, 
             arrow.gap = 0.03, color = "multi-omic alternation", 
             palette = "Set1",  label.color = "black", label.size = 3,  alpha = "deviations", 
             size = "structural influence",
             #size = "deviations",
             # edge.alpha = "multi-omic alternation",
             edge.color = c("color", "orange")
             , legend.position = "bottom")
p4 <- p4 + guides(color = FALSE, size = FALSE)

#allied or rival ? Ora ll1 è lo steady state.
length(names(ll1$multi.omic.pattern)) - length(names(ll1$multi.omic.pattern.c))
b.to.sym.orderd(gene_assoc, names(ll1$multi.omic.pattern.c[names(ll1$multi.omic.pattern.c) %in% names(ll1$multi.omic.pattern) == FALSE]))
ll1$array.alternations$sim.sc
ll1$network.alternations$alter.test
ll1$network.alternations$no.alter.test
ll1$network.alternations$alter.index
ll1$network.alternations$no.alter.index
sum(ll1$multi.omic.df["weight",])/length(ll1$multi.omic.df["weight",])
sum(ll1$omic.amount.df["CAI",])/length(ll1$omic.amount.df["CAI",])
sum(ll1$omic.amount.df["pv",])/length(ll1$omic.amount.df["pv",])

ll1$array.alternations.d$sim.sc
ll1$network.alternations.d$alter.test
ll1$network.alternations.d$no.alter.test
ll1$network.alternations.d$alter.index
ll1$network.alternations.d$no.alter.index
sum(ll1$multi.omic.df.dev["weight",])/length(ll1$multi.omic.df.dev["weight",])
sum(ll1$omic.amount.df.dev["CAI",])/length(ll1$omic.amount.df.dev["CAI",])
sum(ll1$omic.amount.df.dev["pv",])/length(ll1$omic.amount.df.dev["pv",])

multiplot(p4,p1,p3,p2, cols=2)
multiplot(p4,p3, cols=2)
####################################################################################
###    Plot array alternations 
####################################################################################

library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
get(load("./computations.steady.state.Rdata"))
get(load("./eco.df.paths.Rdata"))

for(i in 1:10)
{
  s <- paste("./suzuki", i, ".Rdata", sep="")
  get(load(s))
}

for(i in 1:59)
{
  s <- paste("./faith", i, ".Rdata", sep="")
  get(load(s))
  
}

get(load("./suzukietall.ma.Rdata"))
names.exp <- c()
for(j in 1:length(suzukietall.ma))
{
  names.exp <- c(names.exp, suzukietall.ma[[j]][[1]][[3]])
}

get(load("./norfloxacin.ma.Rdata"))
names.exp2 <- c()
for(j in 1:length(norfloxacin.ma))
{
  names.exp2 <- c(names.exp2, norfloxacin.ma[[j]][[1]][[3]])
}

list.steady.state <- list(steady.state)



strand.sensus <- function(gene_assoc, gene_list)
{
  path.ext         <- gene_assoc[gene_list,]
  path.ext         <- path.ext[!is.na(path.ext$gene_loc),]
  path.ext         <- path.ext[with(path.ext, order(orf_from)), ]
  return(as.vector(path.ext$direction))
}



make.line.ss <- function(omic, list.given, fixed_path, id, x, dir, pos){
    if(length(list.given[[1]][[fixed_path]])>1)
    {
      Y <- scale(x)
      X <- c(1:length(x))
    }
    return(data.frame(X=X,Y=Y,id=id, omic=omic, dir=dir, pos=pos))
  }
  


list.treats  <- list(faith1,faith2,faith3,faith4,faith5,faith6,faith7,faith8,faith9,faith10,
                     faith11,faith12,faith13,faith14,faith15,faith16,faith17,faith18,faith19,faith20,
                     faith21,faith22,faith23,faith24,faith25,faith26,faith27,faith28,faith29,faith30,
                     faith31,faith32,faith33,faith34,faith35,faith36,faith37,faith38,faith39,faith40,
                     faith41,faith42,faith43,faith44,faith45,faith46,faith47,faith48,faith49,faith50,
                     faith51,faith52,faith53,faith54,faith55,faith56,faith57,faith58,faith59)




build.alternations.plot <- function(fixed_path, list.steady.state, list.treats)
{
  nn  <- names(list.treats[[1]][[fixed_path]]$omic.amount.df[1,])
  cai.to <- list.steady.state[[1]][[fixed_path]]$omic.amount.df[1,nn]
  pa.to  <- list.steady.state[[1]][[fixed_path]]$omic.amount.df[2,nn]
  sym.to <- b.to.sym.orderd(gene_assoc, names(cai.to))
  dir.to <- strand.sensus(gene_assoc, names(cai.to))
  dir.to <- gsub("plus", "+", dir.to)
  dir.to <- gsub("minus", "-", dir.to)
  sc.mov <- (scale(as.numeric(pa.to)) + scale(as.numeric(cai.to)))/2
  df.std<- make.line.ss("mov_std", list.steady.state, fixed_path, sym.to, as.vector(cai.to), dir.to, -1)
  df <- data.frame()
  df<-rbind(df,df.std)

  for( x in 1:length(list.treats))
  {
    lb <-  paste("mov:",x, sep = "")
    pa.to  <- list.treats[[x]][[fixed_path]]$omic.amount.df[1,]
    cai.to <- list.treats[[x]][[fixed_path]]$omic.amount.df[2,]
    sc.mov <- (scale(as.numeric(pa.to)) + scale(as.numeric(cai.to)))/2
    sym.to <- b.to.sym.orderd(gene_assoc, names(pa.to))
    dir.to <- strand.sensus(gene_assoc, names(pa.to))
    dir.to <- gsub("plus", "+", dir.to)
    dir.to <- gsub("minus", "-", dir.to)
    df.pa  <- make.line.ss(lb, list.steady.state, fixed_path, sym.to, as.vector(sc.mov), dir.to, x)
    
    
    #df1<-makeLine(label, path.temp$genesym, path.temp[,label], path.temp$senso, x)
    
    df<-rbind(df,df.pa)
    
  }
  
  return(data.frame(group = df$omic, x = df$X, y=df$Y, n=df$pos, label = df$id))
  
}

dataalt <- build.alternations.plot(1, list.steady.state, list.treats)
operons <- list.treats[[1]][[1]]$multi.omic.df[2,] 
names(operons) <- b.to.sym.orderd(gene_assoc, names(operons))
operons.flag <- unique(operons[which(operons > 0)])
rain<-terrain.colors(length(operons.flag))
rain<-colorRampPalette(c("darkorange", "darkblue"))(length(operons.flag)) 

complex <- list.treats[[1]][[1]]$multi.omic.df[3,] 
names(complex) <- b.to.sym.orderd(gene_assoc, names(complex))
complex.flag <- unique(complex[which(complex > 0)])
rain2<-terrain.colors(length(complex.flag))
rain2<-colorRampPalette(c("black", "darkred"))(length(complex.flag)) 


#p1 <- qplot(x=x, y=y, data=dataalt, colour=group, group="mov_std", label=label, alpha=0.3) +
p1  <- qplot(x=x, y=y, data=dataalt, colour=group, group="mov_std", label=label, alpha=0.3) +
  #scale_y_continuous(paste("Multi-omics variations on alternation (59 experiments)")) +
  scale_x_discrete(name="",limits=dataalt[dataalt[,4] ==1,]$label) +
  scale_y_continuous(name="")+
  scale_shape_identity() +
  geom_point(data=dataalt, mapping=aes(x=x, y=y, shape=95,size=1, alpha=0.8 )) +
  geom_blank(mapping = NULL, data = NULL, stat = "identity", position = "identity") +
  theme(
 #   axis.title.x = element_text(face="bold", colour="#990000", size=8),
    axis.text.x  = element_text(angle=90, vjust=0.5, size=7)) +
  geom_line(aes(group="mov_std", alpha=0.7))+
  #theme(legend.position="bottom", legend.title=element_blank())+
  #ggtitle("E.coli K-12 MG165 Glycolysis pathway")+
  scale_fill_hue(l=max(dataalt[,2])) + 
  theme(legend.position="none")


rect <- data.frame(xmin=-Inf, xmax=+Inf, ymin=3.5, ymax=4)
p1 = p1 + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                    fill="#CC79A7",
                    alpha=0.1,
                    linetype = 0,
                    inherit.aes = FALSE)
rect <- data.frame(xmin=-Inf, xmax=+Inf, ymin=3, ymax=3.5)
p1 = p1 + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                    alpha=0.1,fill="#0072B2",
                    linetype = 0,
                    inherit.aes = FALSE)


for(j in 1:length(operons.flag)){
names.op <- c(names(operons[which(operons == operons.flag[j])]))
for(i in 1:length(names.op)){
  p1 = p1 + annotation_custom(grob = textGrob(label= "*" , 
                                    gp = gpar(fontsize = 27, col = rain[j])),  
                    xmin = which(names(operons) == names.op[[i]]), 
                    xmax = which(names(operons) == names.op[[i]]), 
                    ymin = 3.1, 
                    ymax = 3.1)
}
}

for(j in 1:length(complex.flag)){
  names.cm <- c(names(complex[which(complex == complex.flag[j])]))
  for(i in 1:length(names.cm)){
    p1 = p1 + annotation_custom(grob = textGrob(label= "*" , 
                                                gp = gpar(fontsize = 27, col = rain2[j])),  
                                xmin = which(names(complex) == names.cm[[i]]), 
                                xmax = which(names(complex) == names.cm[[i]]), 
                                ymin = 3.6, 
                                ymax = 3.6)
  }
}


p1 


build.alternations.plot.dev <- function(fixed_path, list.steady.state, list.treats)
{
  nn     <- names(list.treats[[1]][[fixed_path]]$omic.amount.df.dev[1,])
  cai.to <- list.steady.state[[1]][[fixed_path]]$omic.amount.df.dev[1,nn]
  pa.to  <- list.steady.state[[1]][[fixed_path]]$omic.amount.df.dev[2,nn]
  sym.to <- gene_assoc[names(cai.to),]$gene_sym
  dir.to <- strand.sensus(gene_assoc, names(cai.to))
  dir.to <- gsub("plus", "+", dir.to)
  dir.to <- gsub("minus", "-", dir.to)
  sc.mov <- (scale(as.numeric(pa.to)) + scale(as.numeric(cai.to)))/2
  df.std <- make.line.ss("mov_std", list.steady.state, fixed_path, sym.to, as.vector(cai.to), dir.to, -1)
  df.std$id <- make.unique(as.character(df.std$id))
  df <- data.frame()
  df<-rbind(df,df.std)
  
  for( x in 1:length(list.treats))
  {
    lb <-  paste("mov:",x, sep = "")
    pa.to  <- list.treats[[x]][[fixed_path]]$omic.amount.df.dev[1,]
    cai.to <- list.treats[[x]][[fixed_path]]$omic.amount.df.dev[2,]
    sc.mov <- (scale(as.numeric(pa.to)) + scale(as.numeric(cai.to)))/2
    sym.to <- gene_assoc[names(cai.to),]$gene_sym
    dir.to <- strand.sensus(gene_assoc, names(pa.to))
    dir.to <- gsub("plus", "+", dir.to)
    dir.to <- gsub("minus", "-", dir.to)
    df.pa  <- make.line.ss(lb, list.steady.state, fixed_path, sym.to, as.vector(sc.mov), dir.to, x)
    df.pa$id <- make.unique(as.character(df.pa$id))
    
    #df1<-makeLine(label, path.temp$genesym, path.temp[,label], path.temp$senso, x)
    
    df<-rbind(df,df.pa)
    
  }
  
  return(data.frame(group = df$omic, x = df$X, y=df$Y, n=df$pos, label = df$id))
  
}

dataalt <- build.alternations.plot.dev(1, list.steady.state, list.treats)
operons <- list.treats[[1]][[1]]$multi.omic.df.dev[2,] 
names(operons) <-gene_assoc[names(operons),]$gene_sym
names(operons) <- make.unique(as.character(names(operons)))
operons.flag <- unique(operons[which(operons > 0)])
rain<-terrain.colors(length(operons.flag))
rain<-colorRampPalette(c("darkorange", "darkblue"))(length(operons.flag)) 

complex <- list.treats[[1]][[1]]$multi.omic.df.dev[3,]
names(complex) <- gene_assoc[names(complex),]$gene_sym
names(complex) <- make.unique(as.character(names(complex)))
complex.flag <- unique(complex[which(complex > 0)])
rain2<-terrain.colors(length(complex.flag))
rain2<-colorRampPalette(c("black", "darkred"))(length(complex.flag)) 
dataalt$label <- gsub(".1","[r]",dataalt$label)
dataalt$label <- gsub(".2","[r ]",dataalt$label)
p2 <- qplot(x=x, y=y, data=dataalt, colour=group, group="mov_std", label=label, alpha=0.3) +
  #scale_y_continuous(paste("Multi-omics variations on alternation (59 experiments)")) +
  scale_x_discrete(name="",limits=dataalt[dataalt[,4] ==1,]$label) +
  scale_y_continuous(name="")+
  scale_shape_identity() +
  geom_point(data=dataalt, mapping=aes(x=x, y=y, shape=95,size=1, alpha=0.8 )) +
  geom_blank(mapping = NULL, data = NULL, stat = "identity", position = "identity") +
  theme(
    #   axis.title.x = element_text(face="bold", colour="#990000", size=8),
    axis.text.x  = element_text(angle=90, vjust=0.5, size=7)) +
  geom_line(aes(group="mov_std", alpha=0.7))+
  #theme(legend.position="bottom", legend.title=element_blank())+
  #ggtitle("E.coli K-12 MG165 Glycolysis pathway")+
  scale_fill_hue(l=max(dataalt[,2])) + 
  theme(legend.position="none")

rect <- data.frame(xmin=-Inf, xmax=+Inf, ymin=3.5, ymax=4)
p2 = p2 + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                    fill="#CC79A7",
                    alpha=0.1,
                    linetype = 0,
                    inherit.aes = FALSE)
rect <- data.frame(xmin=-Inf, xmax=+Inf, ymin=3, ymax=3.5)
p2 = p2 + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                    alpha=0.1,fill="#0072B2",
                    linetype = 0,
                    inherit.aes = FALSE)


for(j in 1:length(operons.flag)){
  names.op <- c(names(operons[which(operons == operons.flag[j])]))
  for(i in 1:length(names.op)){
    p2 = p2 + annotation_custom(grob = textGrob(label= "*" , 
                                                gp = gpar(fontsize = 27, col = rain[j])),  
                                xmin = which(names(operons) == names.op[[i]]), 
                                xmax = which(names(operons) == names.op[[i]]), 
                                ymin = 3.1, 
                                ymax = 3.1)
  }
}

for(j in 1:length(complex.flag)){
  names.cm <- c(names(complex[which(complex == complex.flag[j])]))
  for(i in 1:length(names.cm)){
    p2 = p2 + annotation_custom(grob = textGrob(label= "*" , 
                                                gp = gpar(fontsize = 27, col = rain2[j])),  
                                xmin = which(names(complex) == names.cm[[i]]), 
                                xmax = which(names(complex) == names.cm[[i]]), 
                                ymin = 3.6, 
                                ymax = 3.6)
  }
}


p2 
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