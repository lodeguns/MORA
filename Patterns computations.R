library(GGally)
library(sda)
library(ggplot2)
library(network)
library(graph)
library(igraph)
library(KEGGREST)
library(progress)


kegg.pathws <- c("00010", "00020", "00030", "00040", "00051", "00052", "00053", "00061",
                 "00071", "00130", "00190", "00220", "00230", "00240", "00250", "00260",
                 "00261", "00270", "00280", "00290", "00300", "00310", "00330", "00340", 
                 "00350", "00360", "00362", "00380", "00400", "00410", "00440", "00450",
                 "00480", "00500", "00520", "00521", "00540", "00550", "00561", "00564",
                 "00620", "00627", "00630", "00640", "00650", "00660", "00670", "00680",
                 "00730", "00740", "00750", "00760", "00770", "00780", "00790", "00860",
                 "00900", "00910", "00970", "02010", "02020", "02030", "02040", "02060",
                 "03010", "03018", "03030", "03060", "03070", "03410", "03420", "03430",
                 "03440", "04122")

# eco.df.paths <- data.frame()
# for(i in 1:length(kegg.pathws))
# {
#   df<-as.data.frame(keggList(paste("path:eco",kegg.pathws[[i]], sep="")))
#   colnames(df) <- "name.path"
#   eco.df.paths<-rbind(eco.df.paths , df)
# }
# 
# 
# fname <- "./eco.df.paths.Rdata"
# save(eco.df.paths, file=fname)

load("./gene_assoc.Rdata")
load("./list.kegg.path.igraph.Rdata")
load("./ecocyc.kegg.igraph.Rdata")
load("./norfloxacin.treats.model1.Rdata")
load("./suzukietall.treats.model1.Rdata")
#get(load("./norfloxacin.treats.model1.approx.Rdata"))
#get(load("./suzukietall.treats.model1.approx.Rdata"))
load("./operons.matrix.Rdata")
load("./complex.matrix.Rdata")
load("./computations.steady.state.Rdata")
load("./eco.df.paths.Rdata")

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


get.mov.js  <- function(list.cai, list.pa)
{
  list.cai_hat    <- scale(list.cai, center = TRUE, scale=TRUE)
  list.pa_hat     <- scale(list.pa, center = TRUE, scale=TRUE)
  list.summed     <-   scale( ((list.cai_hat + list.pa_hat) / 2), center = TRUE, scale=TRUE )
  list.summed[list.summed>=0]   <- 1
  list.summed[list.summed<0]    <- 0
  
 return(as.vector(list.summed))
}

list.cai <- gene_assoc$CAI
list.pa  <- gene_assoc$pa_mc
list.cai_hat    <- scale(list.cai, center = TRUE, scale=TRUE)
list.pa_hat     <- scale(list.pa, center = TRUE, scale=TRUE)
list.summed     <-   scale( ((list.cai_hat + list.pa_hat) / 2), center = TRUE, scale=TRUE )
list.summed[list.summed>=0]   <- 1
list.summed[list.summed<0]    <- 0
gene_assoc$list.caipa <- list.summed




sym.to.b.orderd   <- function(gene_assoc, gene_sym)
{
  path.ext         <- gene_assoc[gene_assoc$gene_sym %in% gene_sym,]
  path.ext         <- path.ext[!is.na(path.ext$gene_loc),]
  path.ext         <- path.ext[with(path.ext, order(orf_from)), ]
  return(as.vector(path.ext$gene_loc))
}

b.to.sym.orderd   <- function(gene_assoc, gene_list)
{
  path.ext         <- gene_assoc[gene_list,]
  path.ext         <- path.ext[!is.na(path.ext$gene_loc),]
  path.ext         <- path.ext[with(path.ext, order(orf_from)), ]
  return(as.vector(path.ext$gene_sym))
}

get.the.patterns <- function(i, gene_assoc, list.kegg.path.igraph, global.network)
{
  #extract the graph
  pathway_clean_name <- gsub("eco:","", V(list.kegg.path.igraph[[i]])$name)
  pathway_clean <- set.vertex.attribute(list.kegg.path.igraph[[i]] , "name", value=pathway_clean_name)
  not.found <- setdiff(V(pathway_clean)$name, V(global.network)$name)
  graph_ext_name <- delete_vertices(pathway_clean, not.found)
  graph_ext <- induced.subgraph(graph=global.network, vids=V(graph_ext_name)$name)
  gene_list <- gsub("eco:","", V(graph_ext )$name)
  path.ext         <- gene_assoc[gene_list,]
  path.ext         <- path.ext[!is.na(path.ext$gene_loc),]
  path.ext         <- path.ext[with(path.ext, order(orf_from)), ]
  list.cai         <- path.ext$CAI
  list.pa          <- path.ext$pa_mc
  local_pattern   <- get.mov.js(list.cai, list.pa)
  names(local_pattern) <- path.ext$gene_loc 
  global_pattern   <- as.vector(path.ext$list.caipa)
  names(global_pattern) <- path.ext$gene_loc
  
  return(list(local_pattern_path=local_pattern, 
              global_pattern_path=global_pattern, 
              lost.nodes=  not.found, size.pattern = length(local_pattern), 
              list.cai = list.cai, list.pa = list.pa))
}
#treat <- treats2[,1]
get.global.pattern.treat <- function(treat, gene_assoc)
{ 
  gene_assoc_mod        <-  gene_assoc[gene_assoc$gene_sym %in% names(treat),]
  list.cai.glb          <-  as.vector(gene_assoc_mod$CAI)
  names(list.cai.glb)   <-  rownames(gene_assoc_mod)
  list.pa.glb           <-  as.vector(treat[gene_assoc_mod$gene_sym])
  names(list.pa.glb)    <-  rownames(gene_assoc_mod)
  list.pa.glb           <-  list.pa.glb[!is.na(list.pa.glb)]
  list.cai.glb          <-  list.cai.glb[names(list.pa.glb)]
  global_pattern        <-  get.mov.js(list.cai.glb , list.pa.glb )
  names(global_pattern) <-  names(list.cai.glb)
  return(global_pattern)
}

get.pattern.treatment <- function(treat, gene_assoc, pattern_path, global_pattern)
{   pattern_path1   <-  pattern_path[unique(names(pattern_path))]
    genes_df        <-  gene_assoc[names(pattern_path1),]
    gene_sym        <-  b.to.sym.orderd(gene_assoc, as.vector(genes_df$gene_loc))
    tts             <-  treat[gene_sym] 
    gene_lost       <-  sym.to.b.orderd(gene_assoc, gene_sym[which(is.na(tts))]) 
    gene_sym_lost   <-  gene_sym[which(is.na(tts))]
    tts             <-  tts[!is.na(tts)]
    n.n             <-  sym.to.b.orderd(gene_assoc, names(tts))
    list.cai        <-  gene_assoc[n.n,]$CAI
    list.pa         <-  as.vector(tts)
    
    local_pattern   <-  get.mov.js(list.cai, list.pa)
    names(local_pattern) <- n.n
    names(list.cai)      <- n.n
    names(list.pa)       <- n.n
    l.p             <- local_pattern[names(pattern_path)]
    l.p             <- l.p[!is.na(l.p)]
    
    list.cai.ext <- c()
    list.pa.ext  <- c()
    n.cai.ext    <- c()
    n.pa.ext     <- c()
    n.l.p        <- names(l.p)
    for(i in 1:length(l.p))
    {
      index1 <- grep(n.l.p[i], names(list.cai))
      cai.g  <- list.cai[[index1[[1]]]]
      n.cai.ext <- c(n.cai.ext  , n.l.p[i])
      list.cai.ext <- c(list.cai.ext, cai.g )
      
      index2 <- grep(n.l.p[i], names(list.pa))
      cai.g  <- list.pa[[index1[[1]]]]
      n.pa.ext <- c(n.pa.ext , n.l.p[i])
      list.pa.ext <- c(list.pa.ext, cai.g )
    }
    
    list.cai <- list.cai.ext
    names(list.cai) <- n.cai.ext
    list.pa  <- list.pa.ext
    names(list.pa) <- n.pa.ext
    
    global_pattern        <- global_pattern[names(l.p)]
    
    return(list(local_pattern_path=l.p, 
                global_pattern_path=global_pattern, 
                lost.nodes=  gene_lost, size.pattern = length(l.p), 
                list.cai = list.cai, list.pa = list.pa))
}

##omics in OR
compute.effective.pattern <- function(patterns)
{
  alternation_local  <- compute.alternations(patterns$local_pattern, N=2)
  alternation_global <- compute.alternations(patterns$global_pattern, N=2)
  ## Compute path extensions with the best possible score.
  effective_alternation <- alternation_global | alternation_local
  effective_alternation[effective_alternation == TRUE]  <- 1
  effective_alternation[effective_alternation == FALSE] <- 0
  
  effective_pattern <- c(rep(0,length(patterns$local_pattern)))
  names(effective_pattern) <- names(patterns$local_pattern)
  for(i in 1:length(effective_alternation))
  {
    genes <- unlist(strsplit(names(effective_alternation[i]), "|", fixed = TRUE))
    
    loc.alt <- abs(patterns$local_pattern[genes][[1]] - patterns$local_pattern[genes][[2]])
    glb.alt <-  abs(patterns$global_pattern[genes][[1]] -   patterns$global_pattern[genes][[2]])
    
    if(loc.alt == 1)
    {
      effective_pattern[names(patterns$local_pattern[genes][1])] <- patterns$local_pattern[genes][[1]]
      effective_pattern[names(patterns$local_pattern[genes][2])] <- patterns$local_pattern[genes][[2]]
    } else
    {
      
      effective_pattern[names(patterns$global_pattern[genes][1])] <- patterns$global_pattern[genes][[1]]
      effective_pattern[names(patterns$global_pattern[genes][2])] <- patterns$global_pattern[genes][[2]]
      
    }}
    return(effective_pattern)
}


get.the.operon<-function(gene_assoc, operons.matrix, gene){
  operon <- operons.matrix[!is.na(operons.matrix[, grep(gene, colnames(operons.matrix))]),]
  operon <- names(operon[!is.na(operon)])
  return(list(b.operon = operon, sym.operon =  b.to.sym.orderd(gene_assoc, operon)))
}

get.the.complex<-function(gene_assoc, complex.matrix, gene){
  complex <- complex.matrix[!is.na(complex.matrix[, grep(gene, colnames(complex.matrix))]),]
  complex <- names(complex[!is.na(complex)])
  return(list(b.complex = complex, sym.complex =  b.to.sym.orderd(gene_assoc, complex)))
}

enrich.with.operons <- function(gene_assoc, operons.matrix, pattern_path)
{ op <- 1
  op.c <- -1
  searched<- c()
  genes <-  names(pattern_path)
  operon.label <- rep(0, length(genes)) 
  names(operon.label) <- genes
  operon.path.list <- list()
  for(i in 1:length(genes))
  { gene <- genes[i]
    gto <- get.the.operon(gene_assoc, operons.matrix, gene)
    gg  <- grep(gene, searched)
   
    if(!is.null(gto$b.operon) && length(gg)==0)
    { #print(gto)
      if(length(operon.label[names(operon.label) %in% gto$b.operon ]) > 1)
        { operon.label[names(operon.label) %in% gto$b.operon ] <- op 
          op <- op +1 }
      else
      {
        operon.label[names(operon.label) %in% gto$b.operon ] <- op.c 
        op.c <- op.c -1
        }
      operon.path.list[[length(operon.path.list)+1]] <- gto
     searched <- c(searched, gto$b.operon)
    }
    
  }
 return(list(operon.adjacent = operon.label, all.operons.involved = operon.path.list))
}

enrich.with.complex <- function(gene_assoc, complex.matrix, pattern_path)
{ cp <- 1
  cp.c <- -1
  searched<- c()
  genes <-  names(pattern_path)
  complex.label <- rep(0, length(genes)) 
  names(complex.label) <- genes
  complex.path.list <- list()
  for(i in 1:length(genes))
  { gene <- genes[i]
    gto <- get.the.complex(gene_assoc, complex.matrix, gene)
    gg  <- grep(gene, searched)

if(!is.null(gto$b.complex) && length(gg)==0)
{ #print(gto)
  if(length(complex.label[names(complex.label) %in% gto$b.complex ]) > 1)
  { complex.label[names(complex.label) %in% gto$b.complex ] <- cp 
    cp <- cp +1 }
  else
  {
    complex.label[names(complex.label) %in% gto$b.complex ] <- cp.c 
    cp.c <- cp.c -1
  }
  complex.path.list[[length(complex.path.list)+1]] <- gto
  searched <- c(searched, gto$b.complex)
}

}
return(list(complex.adjacent = complex.label, all.complex.involved = complex.path.list))
}



clear.deviations <- function(clear.dev.pattern)
{ 
  nn  <- names(clear.dev.pattern)
  nn1 <- gsub("[^dev]", "", nn)
  nn1  <- gsub("dev", 1, nn1)
  nn1  <-  as.numeric(nn1)
  nn1[is.na(nn1)] <- 0
  nn <- gsub("dev\\[", "", nn)
  nn <- gsub("\\]", "", nn)
  names(clear.dev.pattern) <- nn
  return(list(clear.dev= clear.dev.pattern, dev = nn1))
}


compute.operon.compression <- function(rx)
{  #rx <- r2
  gg <- grep("dev", rownames(rx))
  devN = FALSE
  if(length(gg)!=0)
  { devN = TRUE
    dev <- rx["dev",]
    operon_t <- rx["operon", ]
    col_n   <- colnames(rx)
    col_comp   <- c()
    col_hidden <- c()
    oper_comp  <- c()
    for(i in 1:length(dev))
    {
      if(dev[[i]]==1)
      {
        oper_comp  <- c(oper_comp, operon_t[[i]])
        operon_t[[i]] = 0
        col_hidden <- c(col_hidden,  paste("hidden",i, sep=""))
        col_comp   <- c(col_comp, col_n[[i]])
        col_n[[i]] = paste("hidden",i, sep="")
        
      }
    }
    
    rx["operon", ] <- operon_t
    colnames(rx)   <- col_n 
  }

   
  temp <- data.frame()
  operons.to.compress <- rx[grep("operon", rownames(rx)),]
  distinct.op <- unique(operons.to.compress[operons.to.compress > 0 ])
  dim.flag=TRUE
  if(length(distinct.op)!=0)
    {
    distinct.op1 <- distinct.op
    for(j in 1:length(distinct.op))
    { o.t.c <- names(operons.to.compress[operons.to.compress == distinct.op[[j]]])
     if(length(o.t.c)!=1){
     compress.block <- rx[, o.t.c]
     n.c <- ncol(compress.block)
     if(is.null(n.c))
     { g1 <- grep(distinct.op[[j]], distinct.op1)
       if(length(g1)!=0)
        { distinct.op1 = distinct.op1[-g1[[1]]] 
       } }
     else{
     compress.block <- as.data.frame(compress.block[, !duplicated(colnames(compress.block))])
    
     n.r <- nrow(compress.block)
     if(is.null(n.r))
     { 
       g1 <- grep(distinct.op[[j]], distinct.op1)
       if(length(g1)!=0)
       { distinct.op1 = distinct.op1[-g1[[1]]] 
        }
     }}
     } else{
       g1 <- grep(distinct.op[[j]], distinct.op1)
        if(length(g1)!=0)
         { distinct.op1 = distinct.op1[-g1[[1]]] 
         }
     }}
    
   if(length(distinct.op1)!=0){
   distinct.op <- distinct.op1
   for(j in 1:length(distinct.op))
   { o.t.c <- names(operons.to.compress[operons.to.compress == distinct.op[[j]]])
     compress.block <- as.data.frame(rx[, o.t.c])
     if(ncol(compress.block)>1 && ncol(compress.block) != ncol(rx)){
     compress.block <- as.data.frame(compress.block[, !duplicated(colnames(compress.block))])
    #compress.block
   
     compress.col <- c()
     for(i in 1:nrow(compress.block))
     { tt.comp <- as.data.frame(table(compress.block[i,]))
       tt.comp <-  tt.comp[with( tt.comp, order(-Freq)), ]
       
       if(nrow(tt.comp) == 1)
       {     compress.col <- as.numeric(c(compress.col, paste(tt.comp[1,1])))
       } else
       {     
           
              compress.col <- as.numeric(c(compress.col, paste(tt.comp[1,1])))
       }
     }
     compress.col <- as.data.frame(compress.col)
     colnames(compress.col) <- paste("(", paste(colnames(compress.block), collapse ="-"), ")", sep="")
     rownames(compress.col) <- rownames(rx)
     gg <- grep(colnames(compress.block)[[1]], colnames(rx))
     gg.x <- grep(colnames(compress.block)[[length(colnames(compress.block))]],colnames(rx))
     #print(paste("-->",gg))
     if(gg[[1]]==1)
     {   rx = cbind(compress.col, rx[,  !(colnames(rx) %in% colnames(compress.block))])}
     
     if(gg.x[[1]] == length(colnames(rx)))
     {
         rx = cbind(rx[,  !(colnames(rx) %in% colnames(compress.block))], compress.col)
     }
     else
     {
       left.rx <- rx[,  colnames(rx)[1:gg[[1]]-1]]
       gg2 <- grep(colnames(compress.block)[[length(colnames(compress.block))]], colnames(rx))
       if(length(gg2)!=0)
       { limit <- gg2[[1]]
         l.c.n <- length(colnames(rx))
         #l.c.n <- abs(l.c.n - limit)
         plus.one <- gg2[[1]]+1
         #print(rx)
         
         nn.name  <- colnames(rx)[c(plus.one:l.c.n)]
         #print(nn.name)
         right.rx <- rx[,paste(nn.name)]
         
         rx <- cbind(left.rx, compress.col)
         rx <- cbind(rx, right.rx)
         
       } else { print("Something goes wrong.") }
     }
    }
   }
  }
  }
  if(devN)
    { if(!is.null(col_hidden)){
      for(i in 1:length(col_hidden))
      {
      gg <- grep(col_hidden[[i]], colnames(rx))
      if(length(gg)!=0){
      rx["operon",gg[[1]]]  <- oper_comp[[i]]
      colnames(rx)[gg[[1]]] <- col_comp[[i]]
      }
    }}
  } 
  return(rx)
  
}

#########################################################################################################
#### Steady state conditions for all the pathways.
#### Multi-omic patterns with operons compressions
#### Multi omic patterns with path extensions
#########################################################################################################

#pathway.under.examinations     <- 1
 steady.state.an <- function(pathway.under.examinations, eco.df.paths,
                             gene_assoc, list.kegg.path.igraph, global.network, global.pattern.net, 
                             operons.matrix, complex.matrix, APL=3)
 { 
   
   name.path <- paste(eco.df.paths[pathway.under.examinations,])
   path.code <- paste(rownames(eco.df.paths)[[pathway.under.examinations]])
   patterns  <- get.the.patterns(pathway.under.examinations, gene_assoc, list.kegg.path.igraph, global.network )
   #global and local effect influences
   effective_pattern_path         <- compute.effective.pattern(patterns)
   #pathway multi-omic alternations on arrays
   effective_alternation_path     <- alternation.scores(0,0, effective_pattern_path, 2)
   effective_altern_path          <- effective_alternation_path$alter
   #compute path extensions need global pattern
   deviations.frame               <- compute.deviations(effective_altern_path , 
                                                        global.network, global.pattern.net,  N=2)
#   #compute alternative path with path extensions
   effective_pattern_path_dev     <- compute.extension.by.deviations(deviations.frame, 
                                                                     global.pattern.net, effective_pattern_path )
   #compute alternations on arrays for the extended pathway
   effective_alternation_path_dev <- alternation.scores(0,0,effective_pattern_path_dev,2)
   #effective_alternation_path 
   #effective_alternation_path_dev 
   
  
   clear.dds                     <- clear.deviations(effective_pattern_path_dev)
   effective_pattern_path_dev.c  <- clear.dds$clear.dev
                             dds <- clear.dds$dev
                             
   #compute network alternations 
   net.alt     <- network.alternations.m(effective_pattern_path , global.network) 
   net.dev.alt <- network.alternations.m(effective_pattern_path_dev.c , global.network) 
   
   count.path     <- get.pattern.adjacent.influences(names(effective_pattern_path), delta=2, psi=2:APL, global.network)
   count.path.dev <- get.pattern.adjacent.influences(names(effective_pattern_path_dev.c), delta=2, psi=2:APL, global.network)
  
   count.path.dev.order <- c()
   count.path.dev.name  <- c()
    for(i in 1:length(names(effective_pattern_path_dev.c)))
    {
      count.path.dev.order <-c(count.path.dev.order, count.path.dev[names(effective_pattern_path_dev.c)[[i]]][[1]])
      count.path.dev.name  <-c(count.path.dev.name , names(effective_pattern_path_dev.c)[[i]])
    }
   names(count.path.dev.order) <- count.path.dev.name
   count.path.dev <- count.path.dev.order
   #compute operon compressions
   operons.list.path <- enrich.with.operons(gene_assoc, operons.matrix, effective_pattern_path )
   complex.list.path <- enrich.with.complex(gene_assoc, complex.matrix, effective_pattern_path )
   r1<-rbind(effective_pattern_path,operons.list.path$operon.adjacent)
   r1<-rbind(r1,complex.list.path$complex.adjacent)
   r1<-rbind(r1, count.path)
   r1<-rbind(r1, patterns$list.cai)
   r1<-rbind(r1, patterns$list.pa)
   rownames(r1) <- c("altern", "operon", "complex", "weight", "CAI", "pv")
   
   operons.list.path.d <- enrich.with.operons(gene_assoc, operons.matrix, effective_pattern_path_dev.c )
   complex.list.path.d <- enrich.with.complex(gene_assoc, complex.matrix, effective_pattern_path_dev.c  )
   
   list.cai <- gene_assoc[names(effective_pattern_path_dev.c),]$CAI
   list.pa  <- gene_assoc[names(effective_pattern_path_dev.c),]$pa_mc
   r2 <- rbind(effective_pattern_path_dev.c ,   operons.list.path.d$operon.adjacent)
   r2 <- rbind(r2, complex.list.path.d$complex.adjacent)
   r2 <- rbind(r2, count.path.dev)
   r2<-rbind(r2, list.cai)
   r2<-rbind(r2, list.pa)
   r2 <- rbind(r2, dds)
   
   rownames(r2) <- c("altern", "operon", "complex","weight", "CAI", "pv","dev" )
   
   op.compr      <- compute.operon.compression(r1)
   effective_pattern_path_op <- round(op.compr [1,],0)
   op.compr.dev  <- compute.operon.compression(r2)
   effective_pattern_path_op_dev <- round(op.compr.dev[1,],0)
   if(is.data.frame(effective_pattern_path_op)){
   e_p_p_o <- c()
   val.val.o <- c()
   for(u in 1:length(effective_pattern_path_op)){
     e_p_p_o <- c(e_p_p_o, gsub("\\)", "",gsub("\\(", "",gsub("-", "", names(effective_pattern_path_op)[[u]]))))
     val.val.o <- c(val.val.o, effective_pattern_path_op[1,u])
   }
   val.val.o <- as.numeric(val.val.o)
   names(val.val.o)<-  e_p_p_o
   op.alt       <- alternation.scores(0,0,val.val.o,2)
   } else
   {
     op.alt       <- alternation.scores(0,0,effective_pattern_path_op,2)
     val.val.o    <- effective_pattern_path_op
   }
   
   if(is.data.frame(effective_pattern_path_op_dev)){
      e_p_p_o_d <- c()
      val.val.o.d <- c()
      for(u in 1:length(effective_pattern_path_op_dev )){
        e_p_p_o_d <- c(e_p_p_o_d, gsub("\\)", "",gsub("\\(", "",gsub("-", "", names(effective_pattern_path_op_dev)[[u]]))))
        val.val.o.d <- c(val.val.o.d, effective_pattern_path_op_dev[1,u])
        }
       val.val.o.d <- as.numeric(val.val.o.d)
       names(val.val.o.d) <-  e_p_p_o_d
       op.compr.alt <- alternation.scores(0,0,val.val.o.d,2)
   } else { op.compr.alt <- alternation.scores(0,0,effective_pattern_path_op_dev,2)
            val.val.o.d <- effective_pattern_path_op_dev}
   
   
  return(list(name.path              =  name.path ,
              path.code              =  path.code ,
              multi.omic.pattern     =  effective_pattern_path ,
              multi.omic.pattern.d   =  effective_pattern_path_dev ,
              multi.omic.pattern.c   =  effective_pattern_path_dev.c,
              multi.omic.pattern.o   =  val.val.o,
              multi.omic.pattern.o.d =  val.val.o.d,
              array.alternations     =  effective_alternation_path, 
              array.alternations.d   =  effective_alternation_path_dev,
              array.alternations.o   =  op.alt,
              array.alternations.o.d =  op.compr.alt,
              network.alternations   =  net.alt,
              network.alternations.d =  net.dev.alt,
              operon.enric           =  operons.list.path,
              operon.enric.d         =  operons.list.path.d,
              complex.enrich         =  complex.list.path,
              complex.enrich.d       =  complex.list.path.d,
              multi.omic.df          =  round(r1[-c(5,6),],2),
              multi.omic.df.dev      =  round(r2[-c(5,6),],2),
              omic.amount.df         =  r1[c(5,6),],
              omic.amount.df.dev     =  r2[c(5,6),],
              multi.omic.df.compr    =  round(op.compr[-c(5,6),],2),
              multi.omic.df.compr.d  =  round(op.compr.dev[-c(5,6),],2),
              omic.amount.df.compr   =  op.compr[c(5,6),]            ,
              omic.amount.df.compr.d =  op.compr.dev[c(5,6),]
              
              ))
 }


########################################################################################################
### Suzuki et all 10 treatments.
###
########################################################################################################
#pattern_under_considerations <- i <- 1
#treatment.i <- 2 #primo esperimento
#treats1 = suzukietall.treats.model1
#treats2 = norfloxacin.treats.model1
#treats <- treats1
# res <- treatment.an(i, treats2, steady.state, treatment.i, 
#                     eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
#                     global.pattern.net, operons.matrix, complex.matrix, APL)


treatment.an <- function(pattern_under_considerations, treats, steady.state, treatment.i,
                         eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                         global.pattern.net, operons.matrix, complex.matrix, APL=3){
  name.path                            <- paste(eco.df.paths[pattern_under_considerations,])
  path.code                            <- paste(rownames(eco.df.paths)[[pattern_under_considerations]])
  ith.treatment                        <- treats[,treatment.i]
  global_pattern_treat                 <- get.global.pattern.treat(ith.treatment, gene_assoc)
  global.pattern.treat                 <- c(global_pattern_treat, global.pattern.net[setdiff(names(global.pattern.net), names(global_pattern_treat))])
  put.it                               <- steady.state[[pattern_under_considerations]]$multi.omic.pattern
  patterns.treat                       <- get.pattern.treatment(ith.treatment, gene_assoc, put.it , global.pattern.treat)
  effective_pattern_treat              <- compute.effective.pattern(patterns.treat)
  effective_pattern_treat              <- effective_pattern_treat[names(effective_pattern_treat) %in% names(global_pattern_treat)]
  effective_alternation_path_treat     <- alternation.scores(0,0, effective_pattern_treat  , 2)
  effective_altern_path_treat          <- effective_alternation_path_treat$alter
  put.it.dev <- steady.state[[pattern_under_considerations]]$multi.omic.pattern.c
  patterns.treat.dev                   <- get.pattern.treatment(ith.treatment, gene_assoc, put.it.dev , global.pattern.treat)
  effective_pattern_treat.dev          <- compute.effective.pattern(patterns.treat.dev)
  effective_pattern_path_dev_treat     <- effective_pattern_treat.dev[names(effective_pattern_treat.dev) %in% names(global_pattern_treat)]
  effective_alternation_path_dev_treat <- alternation.scores(0,0,effective_pattern_path_dev_treat  ,2)
  effective_pattern_path_dev.c         <- effective_pattern_path_dev_treat
  data.group <- effective_pattern_path_dev.c[names(patterns.treat.dev$local_pattern_path)]
  effective_pattern_path_dev.c <- data.group[!is.na(data.group)]
  dds                                  <- steady.state[[pattern_under_considerations]]$multi.omic.df.dev
  dds                                  <- dds[,names(effective_pattern_path_dev.c) ]
  net.alt              <- network.alternations.m(effective_pattern_treat, global.network) 
  net.dev.alt          <- network.alternations.m(effective_pattern_path_dev.c, global.network) 
  count.path           <- steady.state[[pattern_under_considerations]]$multi.omic.df[4,names(effective_pattern_treat)]
  count.path.dev       <- steady.state[[pattern_under_considerations]]$multi.omic.df.dev[4,names(effective_pattern_path_dev.c)]
  

  operons.list.path <- enrich.with.operons(gene_assoc, operons.matrix, effective_pattern_treat )
  complex.list.path <- enrich.with.complex(gene_assoc, complex.matrix, effective_pattern_treat )
  r1 <- rbind(effective_pattern_treat ,operons.list.path$operon.adjacent)
  r1 <- rbind(r1,complex.list.path$complex.adjacent)
  r1 <- rbind(r1, count.path)
  a        <- patterns.treat$list.pa
  names(a) <- names(patterns.treat$local_pattern_path)
  list.cai <- a[names(a) %in% names(effective_pattern_treat) ]
  b <- patterns.treat$list.cai
  names(b) <- names(patterns.treat$local_pattern_path)
  list.pa  <- b[names(b) %in% names(effective_pattern_treat) ]
  r1 <- rbind(r1, list.cai)
  r1 <- rbind(r1, list.pa )
  rownames(r1) <- c("altern", "operon", "complex", "weight", "CAI", "pv")
  

  operons.list.path.d <- enrich.with.operons(gene_assoc, operons.matrix,  effective_pattern_path_dev.c  )
  complex.list.path.d <- enrich.with.complex(gene_assoc, complex.matrix,  effective_pattern_path_dev.c  )

  r2 <- rbind(effective_pattern_path_dev.c  ,operons.list.path.d$operon.adjacent)
  r2 <- rbind(r2, complex.list.path.d$complex.adjacent)
  r2 <- rbind(r2, count.path.dev)
  a       <- patterns.treat.dev$list.pa
  names(a) <- names(patterns.treat.dev$local_pattern_path)
  list.cai <- a[names(a) %in% names(effective_pattern_path_dev.c ) ]
  b <- patterns.treat.dev$list.cai
  names(b) <- names(patterns.treat.dev$local_pattern_path)
  list.pa  <- b[names(b) %in% names(effective_pattern_path_dev.c ) ]
  r2 <- rbind(r2, list.cai)
  r2 <- rbind(r2, list.pa)
  r2 <- rbind(r2, dds[5,])
  
  rownames(r2) <- c("altern", "operon", "complex","weight", "CAI", "pv","dev" )
  
  op.compr      <- compute.operon.compression(r1)
  effective_pattern_path_op <- round(op.compr [1,],0)
  r2.not.dup    <- r2[,!duplicated(colnames(r2))]
  op.compr.dev  <- compute.operon.compression(r2.not.dup )
  effective_pattern_path_op_dev <- round(op.compr.dev[1,],0)
  if(is.data.frame(effective_pattern_path_op)){
    e_p_p_o <- c()
    val.val.o <- c()
    for(u in 1:length(effective_pattern_path_op)){
      e_p_p_o <- c(e_p_p_o, gsub("\\)", "",gsub("\\(", "",gsub("-", "", names(effective_pattern_path_op)[[u]]))))
      val.val.o <- c(val.val.o, effective_pattern_path_op[1,u])
    }
    val.val.o <- as.numeric(val.val.o)
    names(val.val.o)<-  e_p_p_o
    op.alt       <- alternation.scores(0,0,val.val.o,2)
  } else
  {
    op.alt       <- alternation.scores(0,0,effective_pattern_path_op,2)
    val.val.o    <- effective_pattern_path_op
  }
  
 
  
  if(is.data.frame(effective_pattern_path_op_dev)){
    e_p_p_o_d <- c()
    val.val.o.d <- c()
    for(u in 1:length(effective_pattern_path_op_dev )){
      e_p_p_o_d <- c(e_p_p_o_d, gsub("\\)", "",gsub("\\(", "",gsub("-", "", names(effective_pattern_path_op_dev)[[u]]))))
      val.val.o.d <- c(val.val.o.d, effective_pattern_path_op_dev[1,u])
      
    }
    val.val.o.d <- as.numeric(val.val.o.d)
    names(val.val.o.d) <-  e_p_p_o_d
    op.compr.alt <- alternation.scores(0,0,val.val.o.d,2)
  } else { op.compr.alt <- alternation.scores(0,0,effective_pattern_path_op_dev,2) 
           val.val.o.d  <- effective_pattern_path_op_dev}
  
  return(list(name.path              =  name.path ,
              path.code              =  path.code ,
              multi.omic.pattern     =  effective_pattern_treat ,
              multi.omic.pattern.d   =  effective_pattern_path_dev_treat,
              multi.omic.pattern.c   =  effective_pattern_path_dev.c,
              multi.omic.pattern.o   =  val.val.o,
              multi.omic.pattern.o.d =  val.val.o.d,
              array.alternations     =  effective_alternation_path_treat , 
              array.alternations.d   =  effective_alternation_path_dev_treat,
              array.alternations.o   =  op.alt,
              array.alternations.o.d =  op.compr.alt,
              network.alternations   =  net.alt,
              network.alternations.d =  net.dev.alt,
              operon.enric           =  operons.list.path,
              operon.enric.d         =  operons.list.path.d,
              complex.enrich         =  complex.list.path,
              complex.enrich.d       =  complex.list.path.d,
              multi.omic.df          =  round(r1[-c(5,6),],2),
              multi.omic.df.dev      =  round(r2[-c(5,6),],2),
              omic.amount.df         =  r1[c(5,6),],
              omic.amount.df.dev     =  r2[c(5,6),],
              multi.omic.df.compr    =  round(op.compr[-c(5,6),],2),
              multi.omic.df.compr.d  =  round(op.compr.dev[-c(5,6),],2),
              omic.amount.df.compr   =  op.compr[c(5,6),]            ,
              omic.amount.df.compr.d =  op.compr.dev[c(5,6),]
              
  ))
}

compute.ith.treatment<- function(treatment, treatment.i, steady.state,
                                 eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                 global.pattern.net, operons.matrix, complex.matrix, APL)
{ ll <- list()
pb <- progress_bar$new(
  format = " computing [:bar] :percent in :elapsed - :current on :total",
  total = length(eco.df.paths[,1]), clear = FALSE, width= 60)
#pb <- progress_bar$new(total = length(eco.df.paths[,1]))
for(i in 1:length(eco.df.paths[,1]))
{ #print(i)
  pb$tick()
  diffs.v <-  setdiff(V(list.kegg.path.igraph[[i]])$name, V(ecocyc.kegg.igraph)$name) 
  v.list   <- V(list.kegg.path.igraph[[i]])$name
  v.list   <- v.list[!(v.list %in% diffs.v)]
  if(length(v.list)!=0){
    if(length(E(induced.subgraph(graph=ecocyc.kegg.igraph, vids=v.list)))==0)
    {ll[[length(ll)+1]] <- paste(eco.df.paths[i,], "does not presents edges...", sep="")}
    else{
      ll[[length(ll)+1]] <- treatment.an(i, treatment, steady.state, treatment.i, 
                                         eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                         global.pattern.net, operons.matrix, complex.matrix, APL)
    }
  } else
  {
    ll[[length(ll)+1]] <- paste(eco.df.paths[i,], "does not presents edges...", sep="")
  }
}
treatment.list <- ll
return(treatment.list)
}
#######################################################################################
##  Steady state extractions  3 hours
#######################################################################################
#ll1 <- steady.state.an(1,eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
#                        global.pattern.net, operons.matrix, complex.matrix, APL)
# 

compute.steady.state <- function(eco.df.paths, gene_assoc, 
                                 list.kegg.path.igraph,
                                 global.network,
                                 global.pattern.net,
                                 operons.matrix, 
                                 complex.matrix, APL) {
 ll.i <- list()
 pb <- progress_bar$new(
   format = " computing [:bar] :percent in :elapsed - :current on :total",
   total = length(eco.df.paths[,1]), clear = FALSE, width= 60)
 for(i in 1:length(eco.df.paths[,1]))
 { #print(i)
   pb$tick()
    diffs.v <-  setdiff(V(list.kegg.path.igraph[[i]])$name, V(ecocyc.kegg.igraph)$name) 
    v.list   <- V(list.kegg.path.igraph[[i]])$name
    v.list   <- v.list[!(v.list %in% diffs.v)]
   if(length(v.list)!=0){
   if(length(E(induced.subgraph(graph=ecocyc.kegg.igraph, vids=v.list)))==0)
   {ll.i[[length(ll.i)+1]] <- paste(eco.df.paths[i,], "does not presents edges...", sep="")}
   else{
   ll.i[[length(ll.i)+1]] <- steady.state.an(i,eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                          global.pattern.net, operons.matrix, complex.matrix, APL)
   }
   } else
   {
     ll.i[[length(ll.i)+1]] <- paste(eco.df.paths[i,], "does not presents edges...", sep="")
   }
 }
 
 steady.state <- ll.i
 
 return(steady.state)
}

steady.state <- compute.steady.state(eco.df.paths, gene_assoc, 
                                     list.kegg.path.igraph,
                                     global.network,
                                     global.pattern.net,
                                     operons.matrix, 
                                     complex.matrix, APL)

 fname <- "./computations.steady.state.Rdata"
 save(steady.state, file=fname)

####################################################################################
##  
## Suzuki et All. 10 experiments with path extentions and oepron compressions.   7 min for exp       
##
####################################################################################

#pattern_under_considerations <- i <- 1

treatment = suzukietall.treats.model1
#treats2 = norfloxacin.treats.model1

treatment   <- suzukietall.treats.model1
treatment.i <- 1 #primo esperimento
suzuki1     <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                 eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                 global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./suzuki1.Rdata"
save(suzuki1, file=fname)


treatment.i <- 2 #primo esperimento
suzuki2     <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                     eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                     global.pattern.net, operons.matrix, complex.matrix, APL)


fname <- "./suzuki2.Rdata"
save(suzuki2, file=fname)


treatment.i <- 3 #primo esperimento
suzuki3     <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                     eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                     global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./suzuki3.Rdata"
save(suzuki3, file=fname)


treatment.i <- 4 #primo esperimento
suzuki4    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                     eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                     global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./suzuki4.Rdata"
save(suzuki4, file=fname)


treatment.i <- 5 #primo esperimento
suzuki5     <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                     eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                     global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./suzuki5.Rdata"
save(suzuki5, file=fname)


treatment.i <- 6 #primo esperimento
suzuki6     <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                     eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                     global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./suzuki6.Rdata"
save(suzuki6, file=fname)


treatment.i <- 7 #primo esperimento
suzuki7     <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                     eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                     global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./suzuki7.Rdata"
save(suzuki7, file=fname)


treatment.i <- 8 #primo esperimento
suzuki8     <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                     eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                     global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./suzuki8.Rdata"
save(suzuki8, file=fname)


treatment.i <- 9 #primo esperimento
suzuki9     <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                     eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                     global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./suzuki9.Rdata"
save(suzuki9, file=fname)


treatment.i <- 10 #primo esperimento
suzuki10     <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                     eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                     global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./suzuki10.Rdata"
save(suzuki10, file=fname)

####################################################################################
##  
## Suzuki et All. 10 experiments with path extensions and operon compressions.   7 min for exp       
##
####################################################################################

#pattern_under_considerations <- i <- 1

treatment = norfloxacin.treats.model1

treatment.i <- 1 #primo esperimento
faith1    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                     eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                     global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith1.Rdata"
save(faith1, file=fname)

treatment.i <- 2 #primo esperimento
faith2    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith2.Rdata"
save(faith2, file=fname)

treatment.i <- 3 #primo esperimento
faith3    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith3.Rdata"
save(faith3, file=fname)

treatment.i <- 4 #primo esperimento
faith4    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith4.Rdata"
save(faith4, file=fname)

treatment.i <- 5 #primo esperimento
faith5    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith5.Rdata"
save(faith5, file=fname)

treatment.i <- 6 #primo esperimento
faith6    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith6.Rdata"
save(faith6, file=fname)

treatment.i <- 7 #primo esperimento
faith7    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith7.Rdata"
save(faith7, file=fname)

treatment.i <- 8 #primo esperimento
faith8    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith8.Rdata"
save(faith8, file=fname)

treatment.i <- 9 #primo esperimento
faith9    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith9.Rdata"
save(faith9, file=fname)

treatment.i <- 10 #primo esperimento
faith10    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith10.Rdata"
save(faith10, file=fname)

treatment.i <- 11 #primo esperimento
faith11    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith11.Rdata"
save(faith11, file=fname)

treatment.i <- 12 #primo esperimento
faith12    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith12.Rdata"
save(faith12, file=fname)

treatment.i <- 13 #primo esperimento
faith13    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith13.Rdata"
save(faith13, file=fname)

treatment.i <- 14 #primo esperimento
faith14    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith14.Rdata"
save(faith14, file=fname)

treatment.i <- 15 #primo esperimento
faith15    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith15.Rdata"
save(faith15, file=fname)

treatment.i <- 16 #primo esperimento
faith16    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith16.Rdata"
save(faith16, file=fname)

treatment.i <- 17 #primo esperimento
faith17    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith17.Rdata"
save(faith17, file=fname)

treatment.i <- 18 #primo esperimento
faith18    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith18.Rdata"
save(faith18, file=fname)

treatment.i <- 19 #primo esperimento
faith19    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith19.Rdata"
save(faith19, file=fname)

treatment.i <- 20 #primo esperimento
faith20    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith20.Rdata"
save(faith20, file=fname)

treatment.i <- 21 #primo esperimento
faith21    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith21.Rdata"
save(faith21, file=fname)

treatment.i <- 22 #primo esperimento
faith22    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith22.Rdata"
save(faith22, file=fname)

treatment.i <- 23 #primo esperimento
faith23    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith23.Rdata"
save(faith23, file=fname)

treatment.i <- 24 #primo esperimento
faith24    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith24.Rdata"
save(faith24, file=fname)

treatment.i <- 25 #primo esperimento
faith25    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith25.Rdata"
save(faith25, file=fname)

treatment.i <- 26 #primo esperimento
faith26    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith26.Rdata"
save(faith26, file=fname)

treatment.i <- 27 #primo esperimento
faith27    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith27.Rdata"
save(faith27, file=fname)

treatment.i <- 28 #primo esperimento
faith28   <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith28.Rdata"
save(faith28, file=fname)

treatment.i <- 29 #primo esperimento
faith29    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith29.Rdata"
save(faith29, file=fname)

treatment.i <- 30 #primo esperimento
faith30    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith30.Rdata"
save(faith30, file=fname)

treatment.i <- 31 #primo esperimento
faith31    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith31.Rdata"
save(faith31, file=fname)

treatment.i <- 32 #primo esperimento
faith32    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith32.Rdata"
save(faith32, file=fname)

treatment.i <- 33 #primo esperimento
faith33    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith33.Rdata"
save(faith33, file=fname)

treatment.i <- 34 #primo esperimento
faith34    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith34.Rdata"
save(faith34, file=fname)

treatment.i <- 35 #primo esperimento
faith35    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith35.Rdata"
save(faith35, file=fname)

treatment.i <- 36 #primo esperimento
faith36    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith36.Rdata"
save(faith36, file=fname)

treatment.i <- 37 #primo esperimento
faith37    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith37.Rdata"
save(faith37, file=fname)

treatment.i <- 38 #primo esperimento
faith38    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith38.Rdata"
save(faith38, file=fname)

treatment.i <- 39 #primo esperimento
faith39    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith39.Rdata"
save(faith39, file=fname)

treatment.i <- 40 #primo esperimento
faith40     <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith40.Rdata"
save(faith40, file=fname)

treatment.i <- 41 #primo esperimento
faith41    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith41.Rdata"
save(faith41, file=fname)

treatment.i <- 42 #primo esperimento
faith42    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith42.Rdata"
save(faith42, file=fname)

treatment.i <- 43 #primo esperimento
faith43    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith43.Rdata"
save(faith43, file=fname)

treatment.i <- 44 #primo esperimento
faith44    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith44.Rdata"
save(faith44, file=fname)

treatment.i <- 45 #primo esperimento
faith45    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith45.Rdata"
save(faith45, file=fname)

treatment.i <- 46 #primo esperimento
faith46    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith46.Rdata"
save(faith46, file=fname)

treatment.i <- 47 #primo esperimento
faith47    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith47.Rdata"
save(faith47, file=fname)

treatment.i <- 48 #primo esperimento
faith48    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith48.Rdata"
save(faith48, file=fname)

treatment.i <- 49 #primo esperimento
faith49    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith49.Rdata"
save(faith49, file=fname)

treatment.i <- 50 #primo esperimento
faith50    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith50.Rdata"
save(faith50, file=fname)

treatment.i <- 51 #primo esperimento
faith51    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith51.Rdata"
save(faith51, file=fname)

treatment.i <- 52 #primo esperimento
faith52    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith52.Rdata"
save(faith52, file=fname)

treatment.i <- 53 #primo esperimento
faith53    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith53.Rdata"
save(faith53, file=fname)

treatment.i <- 54 #primo esperimento
faith54    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith54.Rdata"
save(faith54, file=fname)

treatment.i <- 55 #primo esperimento
faith55    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith55.Rdata"
save(faith55, file=fname)

treatment.i <- 56 #primo esperimento
faith56     <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith56.Rdata"
save(faith56, file=fname)

treatment.i <- 57 #primo esperimento
faith57    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith57.Rdata"
save(faith57, file=fname)

treatment.i <- 58 #primo esperimento
faith58    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith58.Rdata"
save(faith58, file=fname)

treatment.i <- 59 #primo esperimento
faith59    <- compute.ith.treatment(treatment, treatment.i, steady.state,
                                   eco.df.paths, gene_assoc, list.kegg.path.igraph, global.network, 
                                   global.pattern.net, operons.matrix, complex.matrix, APL)
fname <- "./faith59.Rdata"
save(faith59, file=fname)


update.network.alternations <- function(ll.pr)
{ 
  ll.pr$network.alternations   <- network.alternations.m(ll.pr$multi.omic.pattern, global.network)
  ll.pr$network.alternations.d <- network.alternations.m(ll.pr$multi.omic.pattern.c, global.network)
  return(ll.pr)
}

library(progress)
update.network.alter <- function(ss1){
pb <- progress_bar$new(
      format = " computing [:bar] :percent in :elapsed - :current on :total",
      total  = length(ss1), clear = FALSE, width= 60)
for(i in 1:length(ss1))
{ pb$tick()
  if(length(ss1[[i]]) > 1)
     ss1[[i]] <- update.network.alternations(ss1[[i]])
  
}
return(ss1)

}

#prova <- update.network.alter(suzuki1)

suzuki1 <- update.network.alter(suzuki1)
suzuki2 <- update.network.alter(suzuki2)
suzuki3 <- update.network.alter(suzuki3)
suzuki4 <- update.network.alter(suzuki4)
suzuki5 <- update.network.alter(suzuki5)
suzuki6 <- update.network.alter(suzuki6)
suzuki7 <- update.network.alter(suzuki7)
suzuki8 <- update.network.alter(suzuki8)
suzuki9 <- update.network.alter(suzuki9)
suzuki10 <- update.network.alter(suzuki10)
#####
fname <- "./suzuki1.Rdata"
save(suzuki1, file=fname)

fname <- "./suzuki2.Rdata"
save(suzuki2, file=fname)

fname <- "./suzuki3.Rdata"
save(suzuki3, file=fname)

fname <- "./suzuki4.Rdata"
save(suzuki4, file=fname)

fname <- "./suzuki5.Rdata"
save(suzuki5, file=fname)

fname <- "./suzuki6.Rdata"
save(suzuki6, file=fname)

fname <- "./suzuki7.Rdata"
save(suzuki7, file=fname)

fname <- "./suzuki8.Rdata"
save(suzuki8, file=fname)

fname <- "./suzuki9.Rdata"
save(suzuki9, file=fname)

fname <- "./suzuki10.Rdata"
save(suzuki10, file=fname)



faith1 <- update.network.alter(faith1)
faith2 <- update.network.alter(faith2)
faith3 <- update.network.alter(faith3)
faith4 <- update.network.alter(faith4)
faith5 <- update.network.alter(faith5)
faith6 <- update.network.alter(faith6)
faith7 <- update.network.alter(faith7)
faith8 <- update.network.alter(faith8)
faith9 <- update.network.alter(faith9)
faith10 <- update.network.alter(faith10)
faith11 <- update.network.alter(faith11)
faith12 <- update.network.alter(faith12)
faith13 <- update.network.alter(faith13)
faith14 <- update.network.alter(faith14)
faith15 <- update.network.alter(faith15)
faith16 <- update.network.alter(faith16)
faith17 <- update.network.alter(faith17)
faith18 <- update.network.alter(faith18)
faith19 <- update.network.alter(faith19)
faith20 <- update.network.alter(faith20)
faith21 <- update.network.alter(faith21)
faith22 <- update.network.alter(faith22)
faith23 <- update.network.alter(faith23)
faith24 <- update.network.alter(faith24)
faith25 <- update.network.alter(faith25)
faith26 <- update.network.alter(faith26)
faith27 <- update.network.alter(faith27)
faith28 <- update.network.alter(faith28)
faith29 <- update.network.alter(faith29)
faith30 <- update.network.alter(faith30)
faith31 <- update.network.alter(faith31)
faith32 <- update.network.alter(faith32)
faith33 <- update.network.alter(faith33)
faith34 <- update.network.alter(faith34)
faith35 <- update.network.alter(faith35)
faith36 <- update.network.alter(faith36)
faith37 <- update.network.alter(faith37)
faith38 <- update.network.alter(faith38)
faith39 <- update.network.alter(faith39)
faith40 <- update.network.alter(faith40)
faith41 <- update.network.alter(faith41)
faith42 <- update.network.alter(faith42)
faith43 <- update.network.alter(faith43)
faith44 <- update.network.alter(faith44)
faith45 <- update.network.alter(faith45)
faith46 <- update.network.alter(faith46)
faith47 <- update.network.alter(faith47)
faith48 <- update.network.alter(faith48)
faith49 <- update.network.alter(faith49)
faith50 <- update.network.alter(faith50)
faith51 <- update.network.alter(faith51)
faith52 <- update.network.alter(faith52)
faith53 <- update.network.alter(faith53)
faith54 <- update.network.alter(faith54)
faith55 <- update.network.alter(faith55)
faith56 <- update.network.alter(faith56)
faith57 <- update.network.alter(faith57)
faith58 <- update.network.alter(faith58)
faith59 <- update.network.alter(faith59)

fname <- "./faith1.Rdata"
save(faith1, file=fname)


fname <- "./faith2.Rdata"
save(faith2, file=fname)


fname <- "./faith3.Rdata"
save(faith3, file=fname)


fname <- "./faith4.Rdata"
save(faith4, file=fname)


fname <- "./faith5.Rdata"
save(faith5, file=fname)


fname <- "./faith6.Rdata"
save(faith6, file=fname)


fname <- "./faith7.Rdata"
save(faith7, file=fname)


fname <- "./faith8.Rdata"
save(faith8, file=fname)


fname <- "./faith9.Rdata"
save(faith9, file=fname)


fname <- "./faith10.Rdata"
save(faith10, file=fname)


fname <- "./faith11.Rdata"
save(faith11, file=fname)


fname <- "./faith12.Rdata"
save(faith12, file=fname)


fname <- "./faith13.Rdata"
save(faith13, file=fname)


fname <- "./faith14.Rdata"
save(faith14, file=fname)


fname <- "./faith15.Rdata"
save(faith15, file=fname)


fname <- "./faith16.Rdata"
save(faith16, file=fname)


fname <- "./faith17.Rdata"
save(faith17, file=fname)

fname <- "./faith18.Rdata"
save(faith18, file=fname)


fname <- "./faith19.Rdata"
save(faith19, file=fname)


fname <- "./faith20.Rdata"
save(faith20, file=fname)


fname <- "./faith21.Rdata"
save(faith21, file=fname)


fname <- "./faith22.Rdata"
save(faith22, file=fname)


fname <- "./faith23.Rdata"
save(faith23, file=fname)


fname <- "./faith24.Rdata"
save(faith24, file=fname)


fname <- "./faith25.Rdata"
save(faith25, file=fname)


fname <- "./faith26.Rdata"
save(faith26, file=fname)


fname <- "./faith27.Rdata"
save(faith27, file=fname)

fname <- "./faith28.Rdata"
save(faith28, file=fname)


fname <- "./faith29.Rdata"
save(faith29, file=fname)


fname <- "./faith30.Rdata"
save(faith30, file=fname)


fname <- "./faith31.Rdata"
save(faith31, file=fname)


fname <- "./faith32.Rdata"
save(faith32, file=fname)


fname <- "./faith33.Rdata"
save(faith33, file=fname)

fname <- "./faith34.Rdata"
save(faith34, file=fname)


fname <- "./faith35.Rdata"
save(faith35, file=fname)


fname <- "./faith36.Rdata"
save(faith36, file=fname)


fname <- "./faith37.Rdata"
save(faith37, file=fname)


fname <- "./faith38.Rdata"
save(faith38, file=fname)


fname <- "./faith39.Rdata"
save(faith39, file=fname)


fname <- "./faith40.Rdata"
save(faith40, file=fname)


fname <- "./faith41.Rdata"
save(faith41, file=fname)

fname <- "./faith42.Rdata"
save(faith42, file=fname)


fname <- "./faith43.Rdata"
save(faith43, file=fname)


fname <- "./faith44.Rdata"
save(faith44, file=fname)

fname <- "./faith45.Rdata"
save(faith45, file=fname)


fname <- "./faith46.Rdata"
save(faith46, file=fname)


fname <- "./faith47.Rdata"
save(faith47, file=fname)


fname <- "./faith48.Rdata"
save(faith48, file=fname)

fname <- "./faith49.Rdata"
save(faith49, file=fname)


fname <- "./faith50.Rdata"
save(faith50, file=fname)

fname <- "./faith51.Rdata"
save(faith51, file=fname)


fname <- "./faith52.Rdata"
save(faith52, file=fname)


fname <- "./faith53.Rdata"
save(faith53, file=fname)


fname <- "./faith54.Rdata"
save(faith54, file=fname)


fname <- "./faith55.Rdata"
save(faith55, file=fname)


fname <- "./faith56.Rdata"
save(faith56, file=fname)


fname <- "./faith57.Rdata"
save(faith57, file=fname)


fname <- "./faith58.Rdata"
save(faith58, file=fname)

fname <- "./faith59.Rdata"
save(faith59, file=fname)


steady.state <- update.network.alter(steady.state)
fname <- "./computations.steady.state.Rdata"
save(steady.state, file=fname)

