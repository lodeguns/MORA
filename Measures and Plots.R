
library(ggplot2)
get(load("./computations.steady.state.Rdata"))
get(load("./eco.df.paths.Rdata"))
for(i in 1:10)
{  s <- paste("./suzuki", i, ".Rdata", sep="")
   get(load(s))
}
for(i in 1:59)
{ s <- paste("./faith", i, ".Rdata", sep="")
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

names.exp <- c()
for(j in 1:length(suzukietall.ma))
{
  names.exp <- c(names.exp, suzukietall.ma[[j]][[1]][[3]])
}

list.data <- list(suzuki1, suzuki2, suzuki3, suzuki4,  suzuki5,
                  suzuki6,  suzuki7,  suzuki8,  suzuki9,  suzuki10,
                  faith1,faith2,faith3,faith4,faith5,faith6,faith7,faith8,faith9,faith10,
                  faith11,faith12,faith13,faith14,faith15,faith16,faith17,faith18,faith19,faith20,
                  faith21,faith22,faith23,faith24,faith25,faith26,faith27,faith28,faith29,faith30,
                  faith31,faith32,faith33,faith34,faith35,faith36,faith37,faith38,faith39,faith40,
                  faith41,faith42,faith43,faith44,faith45,faith46,faith47,faith48,faith49,faith50,
                  faith51,faith52,faith53,faith54,faith55,faith56,faith57,faith58,faith59)

get.scores.on.array <- function(j, list.scores) {
nn.path <- c()
aa.alt  <- c()
nt.alt  <- c()
nt.no.alt <- c()
nt.ifl  <- c()
jj      <- c()
for(i in 1:length(list.scores))
{
  if(length(list.scores[[i]])>1){
  nn.path <- c(nn.path, list.scores[[i]]$path.code[[1]])
  aa.alt  <- c(aa.alt,  list.scores[[i]]$array.alternations$sim.sc)
  nt.alt  <- c(nt.alt,  list.scores[[i]]$network.alternations$alter.index)
  nt.no.alt <- c(nt.no.alt, list.scores[[i]]$network.alternations$no.alter.index)
  nt.ifl  <- c(nt.ifl,  median(list.scores[[i]]$multi.omic.df[4,]))
  jj      <- c(jj, j)
  }
}

return(list(name.path=nn.path, array.alt= aa.alt, network.alt= nt.alt, network.no.alt = nt.no.alt, influences =nt.ifl, j=jj))
}

get.scores.deviations <- function(j, list.scores) {
  nn.path <- c()
  aa.alt  <- c()
  nt.alt  <- c()
  nt.no.alt <- c()
  nt.ifl  <- c()
  jj      <- c()
  for(i in 1:length(list.scores))
  {
    if(length(list.scores[[i]])>1){
      nn.path <- c(nn.path, list.scores[[i]]$path.code[[1]])
      aa.alt  <- c(aa.alt,  list.scores[[i]]$array.alternations.d$sim.sc)
      nt.alt  <- c(nt.alt,  list.scores[[i]]$network.alternations.d$alter.index)
      nt.no.alt <- c(nt.no.alt, list.scores[[i]]$network.alternations.d$no.alter.index) 
      nt.ifl  <- c(nt.ifl,  median(list.scores[[i]]$multi.omic.df.dev[4,]))
      jj      <- c(jj, j)
      }
  }
  
  return(list(name.path=nn.path, array.alt= aa.alt, network.alt= nt.alt, network.no.alt=nt.no.alt, influences =nt.ifl, j=jj))
}


get.scores.operons.compression <- function(j, list.scores) {
  nn.path <- c()
  aa.alt  <- c()
  nt.alt  <- c()
  nt.no.alt <- c()
  nt.ifl  <- c()
  jj      <- c()
  for(i in 1:length(list.scores))
  {
    if(length(list.scores[[i]])>1){
      nn.path <- c(nn.path, list.scores[[i]]$path.code[[1]])
      aa.alt  <- c(aa.alt,  list.scores[[i]]$array.alternations.o$sim.sc)
      nt.alt  <- c(nt.alt,  list.scores[[i]]$network.alternations$alter.index)
      nt.no.alt <- c(nt.no.alt, list.scores[[i]]$network.alternations$no.alter.index)
      nt.ifl  <- c(nt.ifl,  median(as.numeric(list.scores[[i]]$multi.omic.df.compr[4,])))
      jj      <- c(jj, j)
    }
  }
  
  return(list(name.path=nn.path, array.alt= aa.alt, network.alt= nt.alt, network.no.alt=nt.no.alt, influences=nt.ifl, j=jj))
}


get.scores.operons.deviations <- function(j, list.scores) {
  nn.path <- c()
  aa.alt  <- c()
  nt.alt  <- c()
  nt.no.alt <- c()
  nt.ifl  <- c()
  jj      <- c()
  for(i in 1:length(list.scores))
  {
    if(length(list.scores[[i]])>1){
      nn.path <- c(nn.path, list.scores[[i]]$path.code[[1]])
      aa.alt  <- c(aa.alt,  list.scores[[i]]$array.alternations.o.d$sim.sc)
      nt.alt  <- c(nt.alt,  list.scores[[i]]$network.alternations.d$alter.index)
      nt.no.alt <- c(nt.no.alt, list.scores[[i]]$network.alternations.d$no.alter.index)
      nt.ifl  <- c(nt.ifl,  median(as.numeric(list.scores[[i]]$multi.omic.df.compr.d[4,])))
      jj      <- c(jj, j)
    }
  }
  
  return(list(name.path=nn.path, array.alt= aa.alt, network.alt= nt.alt, network.no.alt=nt.no.alt, influences=nt.ifl, j=jj))
}



###########################################################################
##
##   Table different modality 
##   Standard conditions and treatments
##
###########################################################################

### Treatments pathway
obtain.df.treatments_standard <- function(j, list.data)
{ ss1 <- get.scores.on.array(j, list.data[[1]])
med.val.array.alt <- ss1$array.alt
med.val.network.alt <- ss1$network.alt
med.val.network.no.alt <-  ss1$network.no.alt
med.val.influences  <- ss1$influences
for(i in 2:length(list.data))
{
  ss2 <- get.scores.on.array(j, list.data[[i]])
  med.val.array.alt <- (med.val.array.alt +ss2$array.alt)/2
  med.val.network.alt <- (med.val.network.alt +ss2$network.alt)/2
  med.val.network.no.alt <- ( med.val.network.no.alt +ss2$network.no.alt)/2
  med.val.influences  <- (med.val.influences + ss2$influences)/2
  
}

df1 <- cbind(med.val.array.alt, med.val.network.alt)
df1 <- cbind(df1, med.val.network.no.alt)
df1 <- cbind(df1, med.val.influences)
df1 <- cbind(df1, ss1$j)
rownames(df1) <- ss1$name.path 
colnames(df1) <- c("array", "network","network.no", "influences","std") 
return(df1)
}

## treatments pathway with path extensions
obtain.df.treatments_deviations <- function(j, list.data)
{ ss1 <- get.scores.deviations(j, list.data[[1]])
med.val.array.alt <- ss1$array.alt
med.val.network.alt <- ss1$network.alt
med.val.network.no.alt <-  ss1$network.no.alt
med.val.influences  <- ss1$influences
for(i in 2:length(list.data))
{
  ss2 <- get.scores.deviations(j, list.data[[i]])
  med.val.array.alt <- (med.val.array.alt + ss2$array.alt)/2
  med.val.network.alt <- (med.val.network.alt +ss2$network.alt)/2
  med.val.network.no.alt <- ( med.val.network.no.alt +ss2$network.no.alt)/2
  med.val.influences  <- (med.val.influences + ss2$influences)/2
  
}

df1 <- cbind(med.val.array.alt, med.val.network.alt)
df1 <- cbind(df1, med.val.network.no.alt)
df1 <- cbind(df1, med.val.influences)
df1 <- cbind(df1, ss1$j)
rownames(df1) <- ss1$name.path 
colnames(df1) <- c("array", "network","network.no", "influences","std") 
return(df1)
}

## treatments pathway with operon compression
obtain.df.treatments_operon <- function(j, list.data)
{ ss1 <- get.scores.operons.compression(j, list.data[[1]])
med.val.array.alt <- ss1$array.alt
med.val.network.alt <- ss1$network.alt
med.val.network.no.alt <-  ss1$network.no.alt
med.val.influences  <- ss1$influences
for(i in 2:length(list.data))
{
  ss2 <- get.scores.operons.compression(j, list.data[[i]])
  med.val.array.alt <- (med.val.array.alt +ss2$array.alt)/2
  med.val.network.alt <- (med.val.network.alt +ss2$network.alt)/2
  med.val.network.no.alt <- ( med.val.network.no.alt +ss2$network.no.alt)/2
  med.val.influences  <- (med.val.influences + ss2$influences)/2
  
}

df1 <- cbind(med.val.array.alt, med.val.network.alt)
df1 <- cbind(df1, med.val.network.no.alt)
df1 <- cbind(df1, med.val.influences)
df1 <- cbind(df1, ss1$j)
rownames(df1) <- ss1$name.path 
colnames(df1) <- c("array", "network","network.no", "influences","std") 
return(df1)
}
## treatments pathway with operon compression and path extensions
obtain.df.treatments_operon_deviations <- function(j, list.data)
{ ss1 <- get.scores.operons.deviations(j, list.data[[1]])
med.val.array.alt <- ss1$array.alt
med.val.network.alt <- ss1$network.alt
med.val.network.no.alt <-  ss1$network.no.alt
med.val.influences  <- ss1$influences
for(i in 2:length(list.data))
{
  ss2 <- get.scores.operons.deviations(j, list.data[[i]])
  med.val.array.alt <- (med.val.array.alt +ss2$array.alt)/2
  med.val.network.alt <- (med.val.network.alt +ss2$network.alt)/2
  med.val.network.no.alt <- ( med.val.network.no.alt +ss2$network.no.alt)/2
  med.val.influences  <- (med.val.influences + ss2$influences)/2
  
}

df1 <- cbind(med.val.array.alt, med.val.network.alt)
df1 <- cbind(df1, med.val.network.no.alt)
df1 <- cbind(df1, med.val.influences)
df1 <- cbind(df1, ss1$j)
rownames(df1) <- ss1$name.path 
colnames(df1) <- c("array", "network","network.no", "influences","std") 
return(df1)
}



obtain.df.sst <- function(ss1)
{df1 <- cbind(ss1$array.alt, ss1$network.alt)
 df1 <- cbind(df1, ss1$network.no.alt)
 df1 <- cbind(df1,ss1$influences)
 df1 <- cbind(df1,ss1$j)
 rownames(df1) <- ss1$name.path 
 colnames(df1) <- c("array", "network","network.no", "influences","std") 
 return(df1)}

table0 <- get.scores.on.array(0, steady.state)
table00 <- obtain.df.sst(table0)
table1 <- get.scores.deviations(1, steady.state)
table11 <- obtain.df.sst(table1)
table2 <- get.scores.operons.compression(2, steady.state)
table22 <- obtain.df.sst(table2)
table3 <- get.scores.operons.deviations(3, steady.state)
table33 <- obtain.df.sst(table3)

names.df       <- gsub(" - Escherichia coli K-12 MG1655", "", eco.df.paths[rownames(table00),1])

table44 <- obtain.df.treatments_standard(4, list.data)
table55 <- obtain.df.treatments_deviations(5, list.data)
table66 <- obtain.df.treatments_operon(6, list.data)
table77 <- obtain.df.treatments_operon_deviations(7, list.data)

#Supplementary Tables

library(xtable)
options(xtable.floating = TRUE)
options(xtable.timestamp = "")
tli.table <- xtable(data.frame(names.df, table00[, c(1:4)], table11[,c(1:4)], table22[,c(1:4)], table33[,c(1:4)]))
digits(tli.table) <- 2



library(xtable)
options(xtable.floating = TRUE)
options(xtable.timestamp = "")
tli.table <- xtable(data.frame(names.df, table44[, c(1:4)], table55[,c(1:4)], table66[,c(1:4)], table77[,c(1:4)]))
digits(tli.table) <- 2


################################################################################
##  Plot 1
##  Steady State - Pathways without deviations

get.the.plot.density <- function(df)
{
  p1 <- ggplot(as.data.frame(df), aes(x=network.no, y=network, z=array, color=std, label=rownames(df))) + 
    geom_point(aes(group=influences))+
    geom_point(aes(group=influences, alpha = influences, size = ifelse( influences > 1.5, 0.4, 0.3)))+
    # geom_smooth(method=lm, aes(fill=array), linetype="blank", alpha=0.1)+
    # geom_density_2d(alpha=.5) +
    geom_text(aes(label =  gsub("path:eco","", rownames(df))),
              parse = TRUE, check_overlap = TRUE, vjust = 0.3, nudge_x=0, nudge_y = 0.05) +
    theme(legend.position="none")
  
  rect <- data.frame(xmin=-Inf, xmax=+Inf, ymin=1, ymax=1.5)
  p1 = p1 + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                      fill="#CC79A7",
                      alpha=0.07,
                      linetype = 0,
                      inherit.aes = FALSE)
  
  rect <- data.frame(xmin=1, xmax=5, ymin=1, ymax=+Inf)
  p1 = p1 + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                      fill="#0072B2",
                      alpha=0.07,
                      linetype = 0,
                      inherit.aes = FALSE)
  
  p1 = p1 + theme_bw() +theme(panel.grid.major = element_line(colour = "gray80"),
                              panel.grid.minor = element_blank(),
                              legend.position="none")
  
  return(p1)
}

get.magnitude.sim.plot <- function(df)
{
  p2 <- ggplot(as.data.frame(df), aes(x=array, y=network, color=std, shape = as.factor(std), label=rownames(df))) + 
    geom_point(aes(group=influences))+
    geom_point(aes(group=influences, alpha = influences, size = ifelse( influences > 1.5, 0.4, 0.3)))+
    geom_smooth(method=lm, aes(fill=std), linetype="blank", alpha=0.1)+
    geom_density_2d(alpha=.5) +
    geom_text(aes(label = gsub("path:eco","",rownames(df))),
              parse = TRUE, check_overlap = TRUE, vjust = 0.3, nudge_y = 0.05) +
    theme(legend.position="none")
  
  p2 <- p2 + theme_bw() +theme(panel.grid.major = element_line(colour = "gray80"),
                               panel.grid.minor = element_blank(),
                               legend.position="none")
  return(p2)
}

## Steady State - Pathways and Pathways with path extensions

table0  <- get.scores.on.array(0, steady.state)
table00 <- obtain.df.sst(table0)
table1  <- get.scores.deviations(1, steady.state)
table11 <- obtain.df.sst(table1)
## Dyadic vs anti-dyadic
p1 <- get.the.plot.density(as.data.frame(rbind(table00,table11))) +
  scale_x_continuous(  limits = c(-0.21,7.5) )+
  scale_y_continuous(  limits = c(0.1,3.7),  name="" )
## Dyadic vs simalarity - allied competitor
p2 <- get.magnitude.sim.plot(as.data.frame(rbind(table00,table11)))+
  scale_x_continuous(  limits = c(0.23,1) )+
  scale_y_continuous(  limits = c(0.1,3.7),  name="")

multiplot(p1,p2, cols=2)
#blackdots
median(table00[,1])  #sim
sd(table00[,1])

median(table00[,2]) #anti-dyadic
sd(table00[,2])

median(table00[,3]) # dyadic
sd(table00[,3])

median(table00[,4])
sd(table00[,4])
#bluedots
median(table11[,1]) #sim
sd(table11[,1])

median(table11[,2]) #anti
sd(table11[,2])

median(table11[,3]) #dyad
sd(table11[,3])

median(table11[,4])
sd(table11[,4])

table2 <- get.scores.operons.compression(2, steady.state)
table22 <- obtain.df.sst(table2)
median(table22[,1]) #sim
sd(table22[,1])

median(table22[,2])
sd(table22[,2])

median(table22[,3])
sd(table22[,3])

median(table22[,4])
sd(table22[,4])
table3 <- get.scores.operons.deviations(3, steady.state)
table33 <- obtain.df.sst(table3)
median(table33[,1]) #sim
sd(table33[,1])

median(table33[,2])
sd(table33[,2])

median(table33[,3])
sd(table33[,3])

median(table33[,4])
sd(table33[,4])


## n competitor and allied sigma > 0.7

a <- table00[table00[,1] > 0.6,]

length(a[a[,2] > 1,1]) #18

a <- table11[table11[,1] > 0.6,]

length(a[a[,2] > 1,1]) #47

a <- table22[table22[,1] > 0.6,]

length(a[a[,2] > 1,1]) #14

a <- table33[table33[,1] > 0.6,]

length(a[a[,2] > 1,1]) #43

#################################################################
###
###  Treatment average vals versus control --> Standard 
###

table0  <- get.scores.on.array(0, steady.state)
table00 <- obtain.df.sst(table0)
table44 <- obtain.df.treatments_standard(4, list.data)
## Dyadic vs simalarity - allied competitor
p3 <- get.magnitude.sim.plot(as.data.frame(rbind(table44,table00)))
median(table44[,1]) #sim
sd(table44[,1])

median(table44[,2])
sd(table44[,2])

median(table44[,3])
sd(table44[,3])

median(table44[,4])
sd(table44[,4])

table1  <- get.scores.deviations(1, steady.state)
table11 <- obtain.df.sst(table1)
table55 <- obtain.df.treatments_deviations(5, list.data)
p3 <- get.magnitude.sim.plot(as.data.frame(rbind(table11,table55)))

median(table11[,1]) #sim
sd(table11[,1])
median(table55[,1]) #sim
sd(table55[,1])

median(table11[,2]) #dyd
sd(table11[,2])
median(table55[,2])
sd(table55[,2])

median(table11[,3]) #dyd
sd(table11[,3])
median(table55[,3])
sd(table55[,3])

a <- table55[table55[,1] > 0.6,]

length(a[a[,2] > 1,1]) #47-->25

table2 <- get.scores.operons.compression(2, steady.state)
table22 <- obtain.df.sst(table2)
table66 <- obtain.df.treatments_operon(6, list.data)

median(table22[,1]) #sim
sd(table22[,1])
median(table66[,1]) #sim
sd(table66[,1])

median(table11[,2]) #dyd
sd(table11[,2])
median(table55[,2])
sd(table55[,2])

median(table11[,3]) #dyd
sd(table11[,3])
median(table55[,3])
sd(table55[,3])


table3 <- get.scores.operons.deviations(3, steady.state)
table33 <- obtain.df.sst(table3)
table77 <- obtain.df.treatments_operon_deviations(7, list.data)

median(table33[,1]) #sim
sd(table33[,1])
median(table77[,1]) #sim
sd(table77[,1])

median(table33[,2]) #dyd
sd(table33[,2])
median(table77[,2])
sd(table77[,2])

median(table33[,3]) #dyd
sd(table33[,3])
median(table77[,3])
sd(table77[,3])

#ok




#######################################################################
###
### Stubs
ss1 <- get.scores.on.array(1, list.data[[1]])
med.val.array.alt <- ss1$array.alt
med.val.network.alt <- ss1$network.alt
med.val.influences  <- ss1$influences
for(i in 2:length(list.data))
{
  ss2 <- get.scores.on.array(1, list.data[[i]])
  med.val.array.alt <- (med.val.array.alt +ss2$array.alt)/2
  med.val.network.alt <- (med.val.network.alt +ss2$network.alt)/2
  med.val.influences  <- (med.val.influences + ss2$influences)/2
  
}

ss2 <- get.scores.on.array(0, steady.state)

df2 <- cbind(ss2$array.alt, ss2$network.alt)
df2 <- cbind(df2,ss2$influences)
df2 <- cbind(df2,ss2$j)
rownames(df2) <- ss2$name.path 
colnames(df2) <- c("array", "network", "influences","std")

df1 <- cbind(med.val.array.alt, med.val.network.alt)
df1 <- cbind(df1, med.val.influences)
df1 <- cbind(df1, ss1$j)
rownames(df1) <- ss1$name.path 
colnames(df1) <- c("array", "network", "influences","std")

df <- rbind(df1,df2)
p6<-ggplot(as.data.frame(df), aes(x=array, y=network, color=std, shape = as.factor(std), 
                                  label=gsub("path:","",rownames(df)))) + 
  geom_point(aes(group=influences))+
  geom_point(aes(group=influences, alpha = influences, size = ifelse( influences > mean(influences), 2, 1)))+
  geom_smooth(method=lm, aes(fill=std), linetype="blank", alpha=0.1)+
  geom_density_2d(alpha=.5)+
  geom_text(aes(label = gsub("path:eco","",rownames(df))),
            parse = TRUE, check_overlap = TRUE, vjust = 0.3, nudge_y = 0.05) +
  theme(legend.position="none")

p6 <- p6 +theme(
  panel.background = element_rect(fill = "gray97"),
  panel.grid.minor = element_line(linetype = "dotted"),
  legend.position = "none")+
  scale_y_continuous( name="")

#http://www.genome.jp/dbget-bin/www_bget?eco00627
#http://www.genome.jp/dbget-bin/www_bget?eco03430
#http://www.genome.jp/dbget-bin/www_bget?eco03410

#################################################################
###
###  Treatment average vals versus control --> Standard with devitions
###
### Steady state with path extensions
ss2 <- get.scores.deviations(0, steady.state)

df2 <- cbind(ss2$array.alt, ss2$network.alt)
df2 <- cbind(df2,ss2$influences)
df2 <- cbind(df2,ss2$j)
rownames(df2) <- ss2$name.path 
colnames(df2) <- c("array", "network", "influences","std")

### Treatments with path extensions
ss1 <- get.scores.deviations(1, list.data[[1]])
med.val.array.alt <- ss1$array.alt
med.val.network.alt <- ss1$network.alt
med.val.influences  <- ss1$influences
for(i in 2:length(list.data))
{
  ss2 <- get.scores.deviations(1, list.data[[i]])
  med.val.array.alt <- (med.val.array.alt +ss2$array.alt)
  med.val.network.alt <- (med.val.network.alt +ss2$network.alt)
  med.val.influences  <- (med.val.influences + ss2$influences)
  
}

med.val.array.alt <-  med.val.array.alt/length(list.data)
med.val.network.alt <- med.val.network.alt/length(list.data)
med.val.influences <- med.val.network.alt/length(list.data)

df1 <- cbind(med.val.array.alt, med.val.network.alt)
df1 <- cbind(df1, med.val.influences)
df1 <- cbind(df1, ss1$j)
rownames(df1) <- ss1$name.path 
colnames(df1) <- c("array", "network", "influences","std")

df <- rbind(df1, df2)
p7<-ggplot(as.data.frame(df), aes(x=array, y=network, color=std, shape = as.factor(std), label=rownames(df))) + 
  geom_point(aes(group=influences))+
  geom_point(aes(group=influences, alpha = influences, size = ifelse( influences > mean(influences), 1, 2)))+
  geom_smooth(method=lm, aes(fill=std), linetype="blank", alpha=0.1)+
  geom_density_2d(alpha=.5)+
  geom_text(aes(label = gsub("path:eco","",rownames(df))),
            parse = TRUE, check_overlap = TRUE, vjust = 0.3, nudge_y = 0.05) +
  theme(legend.position="none")

p7 <- p7 +theme(
  panel.background = element_rect(fill = "gray97"),
  panel.grid.minor = element_line(linetype = "dotted"),
  legend.position = "none") +
  scale_y_continuous( name="")

#################################################################
###
###  Treatment average vals versus control --> Standard with operon compression
###
ss1 <- get.scores.operons.compression(1, list.data[[1]])
med.val.array.alt <- ss1$array.alt
med.val.network.alt <- ss1$network.alt
med.val.influences  <- ss1$influences
for(i in 2:length(list.data))
{
  ss2 <- get.scores.operons.compression(i, list.data[[4]])
  med.val.array.alt <- (med.val.array.alt +ss2$array.alt)/2
  med.val.network.alt <- (med.val.network.alt +ss2$network.alt)/2
  med.val.influences  <- (med.val.influences + ss2$influences)/2
  
}

ss2 <- get.scores.operons.compression(2, steady.state)

df2 <- cbind(ss2$array.alt, ss2$network.alt)
df2 <- cbind(df2,ss2$influences)
df2 <- cbind(df2,ss2$j)
rownames(df2) <- ss2$name.path 
colnames(df2) <- c("array", "network", "influences","std")

df1 <- cbind(med.val.array.alt, med.val.network.alt)
df1 <- cbind(df1, med.val.influences)
df1 <- cbind(df1, ss1$j)
rownames(df1) <- ss1$name.path 
colnames(df1) <- c("array", "network", "influences","std")

df <- rbind(df1, df2)
p8 <- ggplot(as.data.frame(df), aes(x=array, y=network, color=std, shape = as.factor(std), label=rownames(df))) + 
  geom_point(aes(group=influences))+
  geom_point(aes(group=influences, alpha = influences, size = ifelse( influences > mean(influences), 1, 2)))+
  geom_smooth(method=lm, aes(fill=std), linetype="blank", alpha=0.1)+
  geom_density_2d(alpha=.5)+
  geom_text(aes(label = gsub("path:eco","",rownames(df))),
            parse = TRUE, check_overlap = TRUE, vjust = 0.3, nudge_y = 0.05) +
  theme(legend.position="none")

p8 <- p8 +theme(
  panel.background = element_rect(fill = "gray97"),
  panel.grid.minor = element_line(linetype = "dotted"),
  legend.position = "none") +
  scale_y_continuous( name="")
p8
#################################################################
###
###  Treatment average vals versus control --> Standard operon compression path extensions
###
ss1 <- get.scores.operons.deviations(1, list.data[[1]])
med.val.array.alt <- ss1$array.alt
med.val.network.alt <- ss1$network.alt
med.val.influences  <- ss1$influences
for(i in 2:length(list.data))
{
  ss2 <- get.scores.operons.deviations(i, list.data[[4]])
  med.val.array.alt <- (med.val.array.alt +ss2$array.alt)/2
  med.val.network.alt <- (med.val.network.alt +ss2$network.alt)/2
  med.val.influences  <- (med.val.influences + ss2$influences)/2
  
}

ss2 <- get.scores.operons.deviations(2, steady.state)

df2 <- cbind(ss2$array.alt, ss2$network.alt)
df2 <- cbind(df2,ss2$influences)
df2 <- cbind(df2,ss2$j)
rownames(df2) <- ss2$name.path 
colnames(df2) <- c("array", "network", "influences","std")

df1 <- cbind(med.val.array.alt, med.val.network.alt)
df1 <- cbind(df1, med.val.influences)
df1 <- cbind(df1, ss1$j)
rownames(df1) <- ss1$name.path 
colnames(df1) <- c("array", "network", "influences","std")

df <- rbind(df1, df2)
p9 <- ggplot(as.data.frame(df), aes(x=array, y=network, color=std, shape = as.factor(std), label=rownames(df))) + 
  geom_point(aes(group=influences))+
  geom_point(aes(group=influences, alpha = influences, size = ifelse( influences > mean(influences), 1, 2)))+
  geom_smooth(method=lm, aes(fill=std), linetype="blank", alpha=0.1)+
  geom_density_2d(alpha=.5)+
  geom_text(aes(label = gsub("path:eco","",rownames(df))),
            parse = TRUE, check_overlap = TRUE, vjust = 0.3, nudge_y = 0.05) +
  theme(legend.position="none")

p9 <- p9 +theme(
  panel.background = element_rect(fill = "gray97"),
  panel.grid.minor = element_line(linetype = "dotted"),
  legend.position = "none") +
  scale_y_continuous( name="")
p9

