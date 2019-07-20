#Network analysis 

# Load required packages for bipartite
library(permute)
library(lattice)
library(vegan)
## This is vegan 2.4-4
library(statnet.common)
## 
## Attaching package: 'statnet.common'
## The following object is masked from 'package:base':
## 
##     order
library(network)
## network: Classes for Relational Data
## Version 1.13.0 created on 2015-08-31.
## copyright (c) 2005, Carter T. Butts, University of California-Irvine
##                     Mark S. Handcock, University of California -- Los Angeles
##                     David R. Hunter, Penn State University
##                     Martina Morris, University of Washington
##                     Skye Bender-deMoll, University of Washington
##  For citation information, type citation("network").
##  Type help("network-package") to get started.
library(sna)
## sna: Tools for Social Network Analysis
## Version 2.4 created on 2016-07-23.
## copyright (c) 2005, Carter T. Butts, University of California-Irvine
##  For citation information, type citation("sna").
##  Type help(package="sna") to get started.
# Load main package
library(bipartite)
##  This is bipartite 2.08
##  For latest changes see versionlog in  ?"bipartite-package".
##  For citation see: citation("bipartite").
##  Have a nice time plotting and analysing two-mode networks.
## 
## Attaching package: 'bipartite'
## The following object is masked from 'package:vegan':
## 
##     nullmodel
# Read all the interaction network matrices into a list

# other packages

library(reshape2)
library(igraph)
library(networkD3)


#setwd("C:/Users/?caro Menezes/Desktop/sampling_completeness/Appendix S1")
setwd("~/Dropbox/MS_Icaro_Ilheus/MS_Icaro_UESC/data")


###Loading the data of the municipality of Una

lencois=read.csv("67_2.csv", head=T, row.names=1, sep=";", na.strings = "NA")
cariri<-read.csv("24_1.csv", head=T,row.names=1,sep=";", na.strings = "NA")
santamaria<-read.csv("33_2.csv", head=T,row.names=1,sep=";", na.strings = "NA")
santahelena<-read.csv("35_2.csv", head=T,row.names=1,sep=";", na.strings = "NA")
rebio<-read.csv("Rebio.csv", head=T,row.names=1,sep=";", na.strings = "NA")
alianca<-read.csv("50_2.csv", head=T,row.names=1,sep=";", na.strings = "NA")
juerana<-read.csv("51_1.csv", head=T,row.names=1,sep=";", na.strings = "NA")
ceplac<-read.csv("78_4.csv", head=T,row.names=1,sep=";", na.strings = "NA")
novaangelica<-read.csv("78_1.csv", head=T,row.names=1,sep=";", na.strings = "NA")
farao<-read.csv("89_3.csv", head=T,row.names=1,sep=";")

###Loading the data of the municipality of Belmonte
taquara<-read.csv("112_2.csv", head=T,row.names=1,sep=";", na.strings = "NA")
otavio<-read.csv("116_4.csv", head=T,row.names=1,sep=";", na.strings = "NA")
fazrenato<-read.csv("118_1.csv", head=T,row.names=1,sep=";", na.strings = "NA")
sempreviva<-read.csv("120_1.csv", head=T,row.names=1,sep=";", na.strings = "NA")
veracel1<-read.csv("131_1.csv", head=T,row.names=1,sep=";", na.strings = "NA")
veracel2<-read.csv("171_1.csv", head=T,row.names=1,sep=";", na.strings = "NA")
ouroverde1<-read.csv("187_1.csv", head=T,row.names=1,sep=";", na.strings = "NA")
ouroverde2<-read.csv("187_3.csv", head=T,row.names=1,sep=";", na.strings = "NA")
rampineli<-read.csv("46_1.csv", head=T,row.names=1,sep=";", na.strings = "NA")

# Read all the interaction network matrices into a list
webs <- list(lencois,cariri,santamaria,santahelena,rebio,
             alianca,juerana,ceplac,novaangelica,farao,
             taquara,otavio,fazrenato,sempreviva,veracel1,
             veracel2, ouroverde1,ouroverde2,rampineli)

#### I THINK IS EASIER TO HAVE THE NAMES OF THE SITES INSTEAD OF THE CODES BECAUSE 
# IT IS EASIER WHEN THINKING ABOUT THE RESULTS - BUT THAT IS UP TO YOU, OF COURSE!
webs.names <- c("67_1","24_1","33_2","35_2","Rebio",
                "50_1","51_1","78_4","78_1","89_3",
                "112_2","116_4","118_1","120_1","131_1",
                "171_1","187_1","187_3","46_1") 
names(webs) <- webs.names

# View data (interaction matrix, for generating a network, or a web)
lapply(webs, head, n = 2L)# Only display the first two rows in the dataset

# Visualize the observed networks from the datasets
visweb(lencois) 
plotweb(lencois, text.rot=90, col.low = "green", col.high = "blue")



#connectance
net.metrics.connec <- lapply(webs, networklevel, index = 'connectance') 
net.metrics.connec

net.metrics.connec1<-sapply (webs, FUN=networklevel, index="connectance")
net.metrics.connec1

#weighted connectance
net.metrics.wconnec <- lapply(webs, networklevel, index = 'weighted connectance') 
net.metrics.wconnec

### I PREFER TO USE 'SAPPLY' FOR VISUALIZING THE RESULTS IN A EASIER WAY
net.metrics.wconnec1<-sapply (webs, FUN=networklevel, index="weighted connectance")
net.metrics.wconnec1

######### BOTH CONNECTANCE ARE HIGHLY CORRELATED, SO I WOULD USE THE WEIGHTED ONE
plot(net.metrics.connec1,net.metrics.wconnec1)


#weighted NODF
net.metrics.nest <- sapply(webs, networklevel, index = 'weighted NODF') 
net.metrics.nest

#links per species
net.metrics.links <- sapply(webs, networklevel, index = 'links per species') 
net.metrics.links

#linkage density
net.metrics.linkage <- sapply(webs, networklevel, index = 'linkage density') 
net.metrics.linkage

#Shannon diversity
net.metrics.shd <- lapply(webs, networklevel, index = 'Shannon diversity') 
net.metrics.shd

#H2 specialization
net.metrics.spec <- sapply(webs, networklevel, index = 'H2') 
net.metrics.spec

#interaction evenness
net.metrics.evenn <- sapply(webs, networklevel, index = 'interaction evenness') 
net.metrics.evenn


#Modularity

m_lencois<-computeModules(lencois,method="Beckett")
m_lencois
m_cariri<-computeModules(cariri,method="Beckett")
m_cariri
m_santamaria<-computeModules(santamaria,method="Beckett")
m_santamaria
m_santahelena<-computeModules(santahelena,method="Beckett")
m_santahelena
m_rebio<-computeModules(rebio,method="Beckett")
m_rebio
m_alianca<-computeModules(alianca,method="Beckett")
m_alianca
m_juerana<-computeModules(juerana,method="Beckett")
m_juerana
m_ceplac<-computeModules(ceplac,method="Beckett")
m_ceplac
m_novaangelica<-computeModules(novaangelica,method="Beckett")
m_novaangelica
m_farao<-computeModules(farao,method="Beckett")
m_farao
m_taquara<-computeModules(taquara,method="Beckett")
m_taquara
m_otavio<-computeModules(otavio,method="Beckett")
m_otavio
md_fazrenato<-computeModules(fazrenato, method= "Beckett")
md_fazrenato #Erro
#??computeModules
m_sempreviva<-computeModules(sempreviva,method="Beckett")
m_sempreviva
m_veracel1<-computeModules(veracel1,method="Beckett")
m_veracel1
m_veracel2<-computeModules(veracel2,method="Beckett")
m_veracel2
m_ouroverde1<-computeModules(ouroverde1,method="Beckett")
m_ouroverde1
m_ouroverde2<-computeModules(ouroverde2,method="Beckett")
m_ouroverde2
m_rampineli<-computeModules(rampineli,method="Beckett")
m_rampineli

####### easier code for Modularity multiple communities


# Read all the interaction network matrices into a list - NOTE I REMOVE OUROVERDE1
webs1 <- list(lencois,cariri,santamaria,santahelena,rebio,alianca,juerana,ceplac,novaangelica,farao,
             taquara,otavio,fazrenato,sempreviva,veracel1, veracel2, #ouroverde1,
             ouroverde2,rampineli)
webs1.names <- c("67_1","24_1","33_2","35_2","Rebio","50_1","51_1","78_4","78_1","89_3",
                "112_2","116_4","118_1","120_1","131_1","171_1",#"187_1"
                "187_3","46_1") 
names(webs1) <- webs1.names

mod<-sapply(webs1,computeModules)
mod[[1]]
mod_comm<-c(mod[[1]]@likelihood,mod[[2]]@likelihood,
            mod[[3]]@likelihood,mod[[4]]@likelihood,mod[[5]]@likelihood,mod[[6]]@likelihood,mod[[7]]@likelihood,mod[[8]]@likelihood,mod[[9]]@likelihood,mod[[10]]@likelihood,mod[[11]]@likelihood,mod[[12]]@likelihood,mod[[13]]@likelihood,
            mod[[14]]@likelihood,mod[[15]]@likelihood,mod[[16]]@likelihood, mod[[17]]@likelihood,
            mod[[18]]@likelihood)

write.csv(mod_comm,"results_modularity.csv")



##################
## weighted connectance, weighted nestedness, modularity, links per species, 
#linkage density, interaction evenness

results_metrics<-cbind(net.metrics.wconnec1,net.metrics.nest,
                   net.metrics.links,net.metrics.linkage,
                   net.metrics.evenn)
write.csv(results_metrics,"results_metrics.csv")

save.image("Icaro_networks.RData")
