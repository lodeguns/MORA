# The content of interesting motifs and their emerging properties in multi omic metabolic networks

This repository contains the manuscript mentioned in the title, and associated code and data sets. Should you need help running our code, please contact us.

# Citation

Bardozzo F, Liò P, Tafliaferri R, “The content of interesting motifs and their emerging properties in multi omic metabolic networks” (2017).

# Abstract
** Background: **  The two important challenges in the analysis of molecular biology information are the following: (a) data (multi omic information) integration; (b) detecting patterns across large scale molecular networks. The two challenges are actually coupled as the integration of omic information may provide better means to detect multi 
omic patterns that could reveal multi scale or emerging properties at the phenotype levels.
                                
** Results: **  Here we address the problem of integrating various type of molecular information (a large collection of gene expression data, sequence data, codon usage and protein abundances) to analyse the Escherichia Coli metabolic response to treatments at the entire network level. Our algorithm, MORA (Multi-omic relational adjacency) is able to detect patterns which may represent metabolic network motifs at pathway and supra pathway level which could hint at some functional role. We provide description and insights on the algorithm by testing on good size database of responses to antibiotics.
Interestingly, some patterns reveal alternating multi-omics or position variation. Our framework is implemented in a software 
written in R which provides effective and friendly means to design intervention scenarios in the perspective of comparing and testing metabolic models, designing new pathways or the redesign of an existing metabolic pathway or the experimental validation of an in silico metabolic  model using a nearby species, require the information of how multi-omics data build up multi scale phenotypes.
                                
 ** Conclusions: ** The integration of multi-omic data reveals that bacterial multi-omic metabolic networks contain position dependent and multi-omic alternating patterns which could provide clues of long range correlation in the bacterial genome.
 
 ** Source Code **
Source for the MORA is in the main repository, and source used to generate the protein variation is in protein-v, source used to integrate the whole metabolic network from [KEGG](http://www.genome.jp/kegg/) and [EcoCyc](https://ecocyc.org/) is in netw-int. 

Here, for the impatient, is an implementation of MORA in [R](https://cran.r-project.org/)

##   MORA (multi-omic relational adjacencies)

 function(pattern, delta=2, psi, g){
  couple    <- length(pattern) - (delta-2) #max number of couple at this distance: x-0-1-x-1-0-x-1-0-x
  dir.edges <- "all"                       #outgoing edges:out, ingoing edges:in
  count <- c(rep(0, length(pattern)))
  names(count)<- as.character((pattern))
  count <- c(count,count[1:max(psi)-1])
  pattern <- c(pattern, pattern[1:max(psi)-1])
  for(y in 1:length(psi)){
    for(i in 1:couple){
      p.from <- pattern[i]
      p.to   <- pattern[i + (delta-1)] 
      
      patty=NULL
      if(distances(g, v = p.from , to = p.to , mode = dir.edges)[[1]] != Inf)
      {tryCatch( patty <-get.shortest.paths(g, p.from , p.to , mode = dir.edges,
                                            weights = NULL, 
                                            output=c("vpath", "epath", "both")) ,
                 warning = function(w) {print(paste("warning... "));},
                 error   = function(e) {print(paste("no path. "));  }
      )}
      if(length(patty)!=0){
        for( k in 1:length(patty$vpath)){
          if(length(V(g)[patty$vpath[[k]]]$name)==psi[[y]]){
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

** Licence **
Computer Laboratory (Cambridge, UK), NeuRoNe Lab (Salerno - IT)
