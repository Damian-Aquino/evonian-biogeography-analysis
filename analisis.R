

library(data.table)
library(plyr)
library(tidyr)
library(stringr)
library(betapart)
library(geosphere)
library(vegan)
library(cluster)
source("matriz_fun.R") #carga nuestra funcion para hacer las matrices de las listas de la PBDB
source("cJ1.R") #funcion para extrapolar la diversidad

#datos de Alemania
path1 = (paste0("https://paleobiodb.org/data1.2/occs/list.csv?", 
                "base_name=Brachiopoda&",
                "taxon_reso=lump_genus&",
                "interval=Lochkovian,Emsian&",  # bottom, top
                "cc=DE&",  # Alemania
                "show=full,strat,lith,acconly&",
                "occs_created_before=2025-07-31" # stable data set
))

#d_pbdb_DE = fread(path1, na.strings="", stringsAsFactors=TRUE)  # read from API
#write.csv(d_pbdb_DE, file="d_pbdb_DE.csv")
d_pbdb_DE =read.table(file="d_pbdb_DE.csv", header=TRUE, sep=",")

d_pbdb_DE$realm = "OW"
d_pbdb_DE$region = "CO"
d_pbdb_DE$basin = "DE"

mat_DE = matricear(d_pbdb_DE)



#datos de Sudafrica
path2 = (paste0("https://paleobiodb.org/data1.2/occs/list.csv?", 
                "base_name=Brachiopoda&",
                "taxon_reso=lump_genus&",
                "interval=Lochkovian,Emsian&",  # bottom, top
                "cc=ZA&",  # Sudafrica
                "show=full,strat,lith,acconly&",
                "occs_created_before=2025-07-31" # stable data set
))

#d_pbdb_ZA = fread(path2, na.strings="", stringsAsFactors=TRUE)  # read from API
#write.csv(d_pbdb_ZA, file="d_pbdb_ZA.csv")
d_pbdb_ZA = read.table(file="d_pbdb_ZA.csv", header=TRUE, sep=",")

d_pbdb_ZA$realm = "MK"
d_pbdb_ZA$region = "MK"
d_pbdb_ZA$basin = "ZA"

mat_ZA = matricear(d_pbdb_ZA)


#datos de China:Guangxi
path3 = (paste0("https://paleobiodb.org/data1.2/occs/list.csv?", 
                "base_name=Brachiopoda&",
                "taxon_reso=lump_genus&",
                "interval=Lochkovian,Emsian&",  # bottom, top
                "cc=CN&state=Guangxi&",  # China
                "show=full,strat,lith,acconly&",
                "occs_created_before=2025-07-31" # stable data set
))

#d_pbdb_CN = fread(path3, na.strings="", stringsAsFactors=TRUE)  # read from API
#write.csv(d_pbdb_CN, file="d_pbdb_CN.csv")
d_pbdb_CN = read.table(file="d_pbdb_CN.csv", header=TRUE, sep=",")

d_pbdb_CN$realm = "OW"
d_pbdb_CN$region = "SC"
d_pbdb_CN$basin = "CN"

mat_CN = matricear(d_pbdb_CN)


#datos de Bolivia
path4 = (paste0("https://paleobiodb.org/data1.2/occs/list.csv?", 
                "base_name=Brachiopoda&",
                "taxon_reso=lump_genus&",
                "interval=Lochkovian,Emsian&",  # bottom, top
                "lngmin=-68&lngmax=-61&latmin=-20&latmax=-18&",  # Bolivia ANDINA
                "show=full,strat,lith,acconly&",
                "occs_created_before=2025-07-31" # stable data set
))

#d_pbdb_BO = fread(path4, na.strings="", stringsAsFactors=TRUE)  # read from API
#write.csv(d_pbdb_BO, file="d_pbdb_BO.csv")
d_pbdb_BO = read.table(file="d_pbdb_BO.csv", header=TRUE, sep=",")

d_pbdb_BO$realm = "MK"
d_pbdb_BO$region = "MK"
d_pbdb_BO$basin = "BO"

mat_BO = matricear(d_pbdb_BO)


#datos de Francia-España
path5 = (paste0("https://paleobiodb.org/data1.2/occs/list.csv?", 
                "base_name=Brachiopoda&",
                "taxon_reso=lump_genus&",
                "interval=Lochkovian,Emsian&",  # bottom, top
                "lngmin=-7.614&lngmax=7.52&latmin=40.88&latmax=44.66&",  # ~Pirineos
                "show=full,strat,lith,acconly&",
                "occs_created_before=2025-07-31" # stable data set
))

#d_pbdb_PI = fread(path5, na.strings="", stringsAsFactors=TRUE)  # read from API
#write.csv(d_pbdb_PI, file="d_pbdb_PI.csv")
d_pbdb_PI = read.table(file="d_pbdb_PI.csv", header=TRUE, sep=",")

d_pbdb_PI$realm = "OW"
d_pbdb_PI$region = "RB"
d_pbdb_PI$basin = "PI"

mat_PI = matricear(d_pbdb_PI)


#datos de EEUU:Apalaches
path6 = (paste0("https://paleobiodb.org/data1.2/occs/list.csv?", 
                "base_name=Brachiopoda&",
                "taxon_reso=lump_genus&",
                "interval=Lochkovian,Emsian&",  # bottom, top
                "lngmin=-80&lngmax=-72.4&latmin=38.74&latmax=45.1&",  # ~Apalaches: Estados de Nueva York y Pennsylvania
                "show=full,strat,lith,acconly&",
                "occs_created_before=2025-07-31" # stable data set
))

#d_pbdb_AP = fread(path6, na.strings="", stringsAsFactors=TRUE)  # read from API
#write.csv(d_pbdb_AP, file="d_pbdb_AP.csv")
d_pbdb_AP = read.table(file="d_pbdb_AP.csv", header=TRUE, sep=",")

d_pbdb_AP$realm = "EA"
d_pbdb_AP$region = "AP"
d_pbdb_AP$basin = "AP"

mat_AP = matricear(d_pbdb_AP)



#datos de EEUU:Nevada
path7 = (paste0("https://paleobiodb.org/data1.2/occs/list.csv?", 
                "base_name=Brachiopoda&",
                "taxon_reso=lump_genus&",
                "interval=Lochkovian,Emsian&",  # bottom, top
                "cc=US&state=Nevada&",  # Nevada
                "show=full,strat,lith,acconly&",
                "occs_created_before=2025-07-31" # stable data set
))

#d_pbdb_NV = fread(path7, na.strings="", stringsAsFactors=TRUE)  # read from API
#write.csv(d_pbdb_NV, file="d_pbdb_NV.csv")
d_pbdb_NV = read.table(file="d_pbdb_NV.csv", header=TRUE, sep=",")

d_pbdb_NV$realm = "EA"
d_pbdb_NV$region = "NV" #CORDILLERAN
d_pbdb_NV$basin = "NV"
mat_NV = matricear(d_pbdb_NV)


#datos de Republica Checa
path8 = (paste0("https://paleobiodb.org/data1.2/occs/list.csv?", 
                "base_name=Brachiopoda&",
                "taxon_reso=lump_genus&",
                "interval=Lochkovian,Emsian&",  # bottom, top
                "cc=CZ&",  # Checa
                "show=full,strat,lith,acconly&",
                "occs_created_before=2025-07-31" # stable data set
))
#d_pbdb_CZ = fread(path8, na.strings="", stringsAsFactors=TRUE)  # read from API
#write.csv(d_pbdb_CZ, file="d_pbdb_CZ.csv")
d_pbdb_CZ = read.table(file="d_pbdb_CZ.csv", header=TRUE, sep=",")

d_pbdb_CZ$realm = "OW"
d_pbdb_CZ$region = "RB"
d_pbdb_CZ$basin = "CZ"

mat_CZ = matricear(d_pbdb_CZ)

#datos de Brasil-PARANA
path9 = (paste0("https://paleobiodb.org/data1.2/occs/list.csv?", 
                "base_name=Brachiopoda&",
                "taxon_reso=lump_genus&",
                "interval=Lochkovian,Emsian&",  # bottom, top
                "lngmin=-59&lngmax=-48&latmin=-27&latmax=-13&",  # Brasil
                "show=full,strat,lith,acconly&",
                "occs_created_before=2025-07-31" # stable data set
))

#d_pbdb_BR = fread(path9, na.strings="", stringsAsFactors=TRUE)  # read from API
#write.csv(d_pbdb_BR, file="d_pbdb_BR.csv")
d_pbdb_BR = read.table(file="d_pbdb_BR.csv", header=TRUE, sep=",")

d_pbdb_BR$realm = "MK"
d_pbdb_BR$region = "MK"
d_pbdb_BR$basin = "BR"

mat_BR = matricear(d_pbdb_BR)

#Datos de Colombia
path10 = (paste0("https://paleobiodb.org/data1.2/occs/list.csv?", 
                 "base_name=Brachiopoda&",
                 "taxon_reso=lump_genus&",
                 "interval=Lochkovian,Emsian&",  # bottom, top
                 "cc=CO&",  # cOLOMBIA
                 "show=full,strat,lith,acconly&",
                 "occs_created_before=2025-07-31" # stable data set
))

#d_pbdb_CO = fread(path10, na.strings="", stringsAsFactors=TRUE)  # read from API
#write.csv(d_pbdb_CO, file="d_pbdb_CO.csv")
d_pbdb_CO = read.table(file="d_pbdb_CO.csv", header=TRUE, sep=",")

d_pbdb_CO$realm = "EA"
d_pbdb_CO$region = "VC"  #region VENEZUELA-COLOMBIA
d_pbdb_CO$basin = "CO"# cOLOMBIA

mat_CO = matricear(d_pbdb_CO)

#Datos de Venezuela
path11 = (paste0("https://paleobiodb.org/data1.2/occs/list.csv?", 
                 "base_name=Brachiopoda&",
                 "taxon_reso=lump_genus&",
                 "interval=Lochkovian,Emsian&",  # bottom, top
                 "cc=VE&",  # Venezuela
                 "show=full,strat,lith,acconly&",
                 "occs_created_before=2025-07-31" # stable data set
))

#d_pbdb_VE = fread(path11, na.strings="", stringsAsFactors=TRUE)  # read from API
#write.csv(d_pbdb_VE, file="d_pbdb_VE.csv")
d_pbdb_VE = read.table(file="d_pbdb_VE.csv", header=TRUE, sep=",")

d_pbdb_VE$realm = "EA"
d_pbdb_VE$region = "VC" #region COLOMBIA-VENEZUELA
d_pbdb_VE$basin = "VE"


mat_VE = matricear(d_pbdb_VE)

#datos Colombia+Venezuela
d_pbdb_VC = rbind(d_pbdb_VE,d_pbdb_CO)
d_pbdb_VC$basin = "VC" #venezuela colombia como basin
mat_VC = matricear(d_pbdb_VC)


#Datos de Marruecos
path12 = (paste0("https://paleobiodb.org/data1.2/occs/list.csv?", 
                 "base_name=Brachiopoda&",
                 "taxon_reso=lump_genus&",
                 "interval=Lochkovian,Emsian&",  # bottom, top
                 "cc=MA&",  # Marruecos
                 "show=full,strat,lith,acconly&",
                 "occs_created_before=2025-07-31" # stable data set
))

#d_pbdb_MA = fread(path12, na.strings="", stringsAsFactors=TRUE)  # read from API
#write.csv(d_pbdb_MA, file="d_pbdb_MA.csv")
d_pbdb_MA = read.table(file="d_pbdb_MA.csv", header=TRUE, sep=",")

d_pbdb_MA$realm = "OW"
d_pbdb_MA$region = "RB"
d_pbdb_MA$basin = "MA"
mat_MA = matricear(d_pbdb_MA)

#Datos de Argelia
path13 = (paste0("https://paleobiodb.org/data1.2/occs/list.csv?", 
                 "base_name=Brachiopoda&",
                 "taxon_reso=lump_genus&",
                 "interval=Lochkovian,Emsian&",  # bottom, top
                 "cc=DZ&",  # Argelia
                 "show=full,strat,lith,acconly&",
                 "occs_created_before=2025-07-31" # stable data set
))

#d_pbdb_DZ = fread(path13, na.strings="", stringsAsFactors=TRUE)  # read from API
#write.csv(d_pbdb_DZ, file="d_pbdb_DZ.csv")
d_pbdb_DZ = read.table(file="d_pbdb_DZ.csv", header=TRUE, sep=",")

d_pbdb_DZ$realm = "OW"
d_pbdb_DZ$region = "RB"
d_pbdb_DZ$basin = "DZ"
mat_DZ = matricear(d_pbdb_DZ)

#datos Marruecos+Argelia
d_pbdb_MD = rbind(d_pbdb_MA,d_pbdb_DZ)
d_pbdb_MD$basin = "MD" #antiatlas como cuenca
mat_MD = matricear(d_pbdb_MA)


#datos de Australia
path14 = (paste0("https://paleobiodb.org/data1.2/occs/list.csv?", 
                 "base_name=Brachiopoda&",
                 "taxon_reso=lump_genus&",
                 "interval=Lochkovian,Emsian&",  # bottom, top
                 "lngmin=142&lngmax=153&latmin=-29&latmax=-38&",  # Australia
                 "show=full,strat,lith,acconly&",
                 "occs_created_before=2025-07-31" # stable data set
))

#d_pbdb_AU = fread(path14, na.strings="", stringsAsFactors=TRUE)  # read from API
#write.csv(d_pbdb_AU, file="d_pbdb_AU.csv")
d_pbdb_AU = read.table(file="d_pbdb_AU.csv", header=TRUE, sep=",")

d_pbdb_AU$realm = "OW"
d_pbdb_AU$region = "AU"
d_pbdb_AU$basin = "AU"
mat_AU = matricear(d_pbdb_AU)

#datos de Amazonia
path15 = (paste0("https://paleobiodb.org/data1.2/occs/list.csv?", 
                 "base_name=Brachiopoda&",
                 "taxon_reso=lump_genus&",
                 "interval=Lochkovian,Emsian&",  # bottom, top
                 "cc=BR&",  # Brasil
                 "lngmin=-55.1&lngmax=-41.4&latmin=-7.9&latmax=-0.4&", 
                 "show=full,strat,lith,acconly&",
                 "occs_created_before=2025-07-31" # stable data set
))

#d_pbdb_AM = fread(path15, na.strings="", stringsAsFactors=TRUE)  # read from API
#write.csv(d_pbdb_AM, file="d_pbdb_AM.csv")
d_pbdb_AM = read.table(file="d_pbdb_AM.csv", header=TRUE, sep=",")

d_pbdb_AM$realm = "EA"
d_pbdb_AM$region = "VC"
d_pbdb_AM$basin = "AM"
mat_AM = matricear(d_pbdb_AM)

####################################################################################################################

#Datos de Venezuela+Colombia y Marruecos+Argelia
datos_pbdb = rbind(d_pbdb_AP, d_pbdb_CN, d_pbdb_PI, d_pbdb_DE, d_pbdb_CZ, d_pbdb_NV, d_pbdb_VC, d_pbdb_MD, d_pbdb_AU, d_pbdb_BO, d_pbdb_ZA, d_pbdb_BR, d_pbdb_AM) 

localidades=c("AP", "CN", "PI", "DE", "CZ", "NV", "VC", "MD", "AU", "BO", "ZA", "BR" , "AM")
loc_realm = c("EA", "OW", "OW", "OW", "OW", "EA", "EA", "OW", "OW", "MK", "MK", "MK", "EA")

realm_labels=unique(datos_pbdb$realm)
region_labels=unique(datos_pbdb$region)
basin_labels=unique(datos_pbdb$basin)
collect_labels=unique(datos_pbdb$collection_no)

mat_EA = matricear(datos_pbdb[datos_pbdb$realm == "EA",])
mat_VM = matricear(datos_pbdb[datos_pbdb$realm == "OW",])
mat_MK = matricear(datos_pbdb[datos_pbdb$realm == "MK",])

mat_EA_basin = matricear.region(datos_pbdb[datos_pbdb$realm == "EA",], "basin")
mat_VM_basin = matricear.region(datos_pbdb[datos_pbdb$realm == "OW",], "basin")
mat_MK_basin = matricear.region(datos_pbdb[datos_pbdb$realm == "MK",], "basin")

mat_realm = matricear.region(datos_pbdb, "realm")
mat_region = matricear.region(datos_pbdb, "region")
mat_basin = matricear.region(datos_pbdb, "basin")
mat_col = matricear.region(datos_pbdb, "collection_no")
mat_col2 = matricear.region(datos_pbdb[!is.na(datos_pbdb$paleolat),], "collection_no") #sin las muestras sin coordenadas
col_localidad = unique(datos_pbdb[!is.na(datos_pbdb$paleolat),c(5,129:131)])

mat_col3=cbind(mat_col, unique(datos_pbdb[,c("collection_no", "realm", "region", "basin")])[,2:4])


lat_coll = unique(datos_pbdb[,c("collection_no", "paleolng", "paleolat")])

cor.test(lat_coll$paleolat, rowSums(mat_col), method="kendall")


#COMPARACIONES DE A PARES

data.mat_region <- beta.pair(mat_region)
data.mat_realm <- beta.pair(mat_realm)
data.mat_basin <- beta.pair(mat_basin)
data.mat_collection <-beta.pair(mat_col)
data.mat_col <- beta.pair(mat_col[rowSums(mat_col)>3,])
data.mat_col2 <- beta.pair(mat_col2) #sin muestras sin coordenadas

write.csv(datos_pbdb, "ocs.braquiopodos.scv")
###################################################################################################################


#OBJETOS BETAPART

mat_BO.core <- betapart.core (mat_BO)
mat_DE.core <- betapart.core (mat_DE)
mat_CN.core <- betapart.core (mat_CN)
mat_ZA.core <- betapart.core (mat_ZA)
mat_PI.core <- betapart.core (mat_PI)
mat_AP.core <- betapart.core (mat_AP)
mat_NV.core <- betapart.core (mat_NV)
mat_CZ.core <- betapart.core (mat_CZ)
mat_BR.core <- betapart.core (mat_BR)
mat_MD.core <- betapart.core (mat_MD)
mat_VC.core <- betapart.core (mat_VC)
mat_AU.core <- betapart.core (mat_AU)
mat_AM.core <- betapart.core (mat_AM)

mat_EA.core <- betapart.core (mat_EA)
mat_VM.core <- betapart.core (mat_VM)
mat_MK.core <- betapart.core (mat_MK)

mat_EA_basin.core <- betapart.core (mat_EA_basin)
mat_VM_basin.core <- betapart.core (mat_VM_basin)
mat_MK_basin.core <- betapart.core (mat_MK_basin)


#MEDIDA DE SITIOS MULTIPLES

data.mat_BO.multi <- beta.multi(mat_BO.core)
data.mat_CN.multi <- beta.multi(mat_CN.core)
data.mat_DE.multi <- beta.multi(mat_DE.core)
data.mat_ZA.multi <- beta.multi(mat_ZA.core)
data.mat_PI.multi <- beta.multi(mat_PI.core)
data.mat_AP.multi <- beta.multi(mat_AP.core)
data.mat_NV.multi <- beta.multi(mat_NV.core)
data.mat_CZ.multi <- beta.multi(mat_CZ.core)
data.mat_BR.multi <- beta.multi(mat_BR.core)
data.mat_MD.multi <- beta.multi(mat_MD.core)
data.mat_VC.multi <- beta.multi(mat_VC.core)
data.mat_AU.multi <- beta.multi(mat_AU.core)
data.mat_AM.multi <- beta.multi(mat_AM.core)

data.mat_EA.multi <- beta.multi(mat_EA.core)
data.mat_VM.multi <- beta.multi(mat_VM.core)
data.mat_MK.multi <- beta.multi(mat_MK.core)

data.mat_EA_basin.multi <- beta.multi(mat_EA_basin.core)
data.mat_VM_basin.multi <- beta.multi(mat_VM_basin.core)
data.mat_MK_basin.multi <- beta.multi(mat_MK_basin.core)

#NUESTREOS A TRAVES DE SITIOS IGUALES

data.mat_BO.samp <- beta.sample(mat_BO.core, sites = 8, samples = 1000)
data.mat_CN.samp <- beta.sample(mat_CN.core, sites = 8, samples = 1000)
data.mat_DE.samp <- beta.sample(mat_DE.core, sites = 8, samples = 1000)
data.mat_ZA.samp <- beta.sample(mat_ZA.core, sites = 8, samples = 1000)
data.mat_PI.samp <- beta.sample(mat_PI.core, sites = 8, samples = 1000)
data.mat_AP.samp <- beta.sample(mat_AP.core, sites = 8, samples = 1000)
data.mat_NV.samp <- beta.sample(mat_NV.core, sites = 8, samples = 1000)
data.mat_CZ.samp <- beta.sample(mat_CZ.core, sites = 8, samples = 1000)
data.mat_BR.samp <- beta.sample(mat_BR.core, sites = 8, samples = 1000)
data.mat_MD.samp <- beta.sample(mat_MD.core, sites = 8, samples = 1000)
data.mat_VC.samp <- beta.sample(mat_VC.core, sites = 8, samples = 1000)
data.mat_AU.samp <- beta.sample(mat_AU.core, sites = 8, samples = 1000)
data.mat_AM.samp <- beta.sample(mat_AM.core, sites = 8, samples = 1000)


data.mat_EA.samp <- beta.sample(mat_EA.core, sites = 50, samples = 10000)
data.mat_VM.samp <- beta.sample(mat_VM.core, sites = 50, samples = 10000)
data.mat_MK.samp <- beta.sample(mat_MK.core, sites = 50, samples = 10000)

data.mat_EA_basin.samp <- beta.sample(mat_EA_basin.core, sites = 3, samples = 100)
data.mat_VM_basin.samp <- beta.sample(mat_VM_basin.core, sites = 3, samples = 100)
data.mat_MK_basin.samp <- beta.sample(mat_MK_basin.core, sites = 3, samples = 100)

#PLOTEO DE LA DISTRIBUCION DE LOS COMPONENTES

dist.BO <- data.mat_BO.samp$sampled.values
dist.CN <- data.mat_CN.samp$sampled.values
dist.DE <- data.mat_DE.samp$sampled.values
dist.ZA <- data.mat_ZA.samp$sampled.values
dist.PI <- data.mat_PI.samp$sampled.values
dist.AP <- data.mat_AP.samp$sampled.values
dist.NV <- data.mat_NV.samp$sampled.values
dist.CZ <- data.mat_CZ.samp$sampled.values
dist.BR <- data.mat_BR.samp$sampled.values
dist.MD <- data.mat_MD.samp$sampled.values
dist.VC <- data.mat_VC.samp$sampled.values
dist.AU <- data.mat_AU.samp$sampled.values
dist.AM <- data.mat_AM.samp$sampled.values

dist.EA <- data.mat_EA.samp$sampled.values
dist.VM <- data.mat_VM.samp$sampled.values
dist.MK <- data.mat_MK.samp$sampled.values

dist.EA_basin <- data.mat_EA_basin.samp$sampled.values
dist.VM_basin <- data.mat_VM_basin.samp$sampled.values
dist.MK_basin <- data.mat_MK_basin.samp$sampled.values


#######################################################################################################################
#SOR, SIM y SNE vs paleolatitudes

#latitudes promedio de realm
median_lat_realm=c()
quantile_lat_realm = NULL
for (i in 1:length(realm_labels)) {
  median_lat_realm[i] = median(datos_pbdb$paleolat[datos_pbdb$realm==realm_labels[i]], na.rm=TRUE)
  quantile_lat_realm = rbind(quantile_lat_realm, quantile(datos_pbdb$paleolat[datos_pbdb$realm==realm_labels[i]], na.rm=TRUE, probs=c(0.25,0.75)))
}
names(median_lat_realm)=realm_labels
rownames(quantile_lat_realm) = realm_labels

#latitudes promedio de region
median_lat_region=c()
quantile_lat_region = NULL
for (i in 1:length(region_labels)) {
  median_lat_region[i] = median(datos_pbdb$paleolat[datos_pbdb$region==region_labels[i]], na.rm=TRUE)
  quantile_lat_region = rbind(quantile_lat_region, quantile(datos_pbdb$paleolat[datos_pbdb$region==region_labels[i]], na.rm=TRUE, probs=c(0.25,0.75)))
}
names(median_lat_region)=region_labels
rownames(quantile_lat_region) = region_labels

#latitudes promedio de basin
median_lat_basin=c()
quantile_lat_basin = NULL
for (i in 1:length(basin_labels)) {
  median_lat_basin[i] = median(datos_pbdb$paleolat[datos_pbdb$basin==basin_labels[i]], na.rm=TRUE)
  quantile_lat_basin = rbind(quantile_lat_basin, quantile(datos_pbdb$paleolat[datos_pbdb$basin==basin_labels[i]], na.rm=TRUE, probs=c(0.25,0.75)))
}
names(median_lat_basin)=basin_labels
rownames(quantile_lat_basin) = basin_labels

#por regiones
median_SOR = c(median(dist.AP$beta.SOR), 
               median(dist.CN$beta.SOR), 
               median(dist.PI$beta.SOR), 
               median(dist.DE$beta.SOR), 
               median(dist.CZ$beta.SOR), 
               median(dist.NV$beta.SOR), 
               median(dist.VC$beta.SOR),
               median(dist.AU$beta.SOR),
               median(dist.MD$beta.SOR),  
               median(dist.BO$beta.SOR), 
               median(dist.ZA$beta.SOR),
               median(dist.BR$beta.SOR),
               median(dist.AM$beta.SOR))
names(median_SOR) = localidades

median_SIM = c(median(dist.AP$beta.SIM), 
               median(dist.CN$beta.SIM), 
               median(dist.PI$beta.SIM), 
               median(dist.DE$beta.SIM), 
               median(dist.CZ$beta.SIM),
               median(dist.NV$beta.SIM), 
               median(dist.VC$beta.SIM), 
               median(dist.AU$beta.SIM), 
               median(dist.MD$beta.SIM),
               median(dist.BO$beta.SIM), 
               median(dist.ZA$beta.SIM),
               median(dist.BR$beta.SIM),
               median(dist.AM$beta.SIM))
names(median_SIM) = localidades

median_SNE = c(median(dist.AP$beta.SNE), 
               median(dist.CN$beta.SNE), 
               median(dist.PI$beta.SNE), 
               median(dist.DE$beta.SNE), 
               median(dist.CZ$beta.SNE),
               median(dist.NV$beta.SNE), 
               median(dist.VC$beta.SNE),
               median(dist.AU$beta.SNE),
               median(dist.MD$beta.SNE), 
               median(dist.BO$beta.SNE), 
               median(dist.ZA$beta.SNE),
               median(dist.BR$beta.SNE),
               median(dist.AM$beta.SNE))
names(median_SNE) = localidades

#distancias latitudinales y longitudinales
paleocoord = stats::aggregate(datos_pbdb[,c("paleolng", "paleolat")], by=list(datos_pbdb$basin),FUN= median)
ordenar = match (data.frame(localidades)$localidades, paleocoord$Group.1)
paleocoord = paleocoord[ordenar,]
rownames(paleocoord) = paleocoord[,1]
paleocoord = data.frame(paleocoord[,2:3])


#distancias lineales
dist_lineales_basin = distm(paleocoord)
colnames(dist_lineales_basin) = localidades
rownames(dist_lineales_basin) = localidades
dist_lineales_basin = as.dist(dist_lineales_basin)


median_lat = c(median(d_pbdb_AP$paleolat, na.rm=TRUE), 
               median(d_pbdb_CN$paleolat, na.rm=TRUE),
               median(d_pbdb_PI$paleolat, na.rm=TRUE), 
               median(d_pbdb_DE$paleolat, na.rm=TRUE), 
               median(d_pbdb_CZ$paleolat, na.rm=TRUE), 
               median(d_pbdb_NV$paleolat, na.rm=TRUE),
               median(d_pbdb_VC$paleolat, na.rm=TRUE),
               median(d_pbdb_MD$paleolat, na.rm=TRUE), 
               median(d_pbdb_AU$paleolat, na.rm=TRUE),
               median(d_pbdb_BO$paleolat, na.rm=TRUE), 
               median(d_pbdb_ZA$paleolat, na.rm=TRUE),
               median(d_pbdb_BR$paleolat, na.rm=TRUE),
               median(d_pbdb_AM$paleolat, na.rm=TRUE))
names(median_lat) = localidades

median_lng = c(median(d_pbdb_AP$paleolng, na.rm=TRUE), 
               median(d_pbdb_CN$paleolng, na.rm=TRUE), 
               median(d_pbdb_PI$paleolng, na.rm=TRUE), 
               median(d_pbdb_DE$paleolng, na.rm=TRUE), 
               median(d_pbdb_CZ$paleolng, na.rm=TRUE), 
               median(d_pbdb_NV$paleolng, na.rm=TRUE), 
               median(d_pbdb_VC$paleolng, na.rm=TRUE), 
               median(d_pbdb_MD$paleolng, na.rm=TRUE),
               median(d_pbdb_AU$paleolng, na.rm=TRUE),
               median(d_pbdb_BO$paleolng, na.rm=TRUE), 
               median(d_pbdb_ZA$paleolng, na.rm=TRUE),
               median(d_pbdb_BR$paleolng, na.rm=TRUE),
               median(d_pbdb_AM$paleolng, na.rm=TRUE))
names(median_lng) = localidades


#Diversidad beta según reino y latitud 
#armo tablas de SOR, SIM y SNE
Sor_realm = cbind(c(mean(dist.EA$beta.SOR), mean(dist.VM$beta.SOR), mean(dist.MK$beta.SOR)), c(min(dist.EA$beta.SOR), min(dist.VM$beta.SOR), min(dist.MK$beta.SOR)), c(max(dist.EA$beta.SOR), max(dist.VM$beta.SOR), max(dist.MK$beta.SOR)))
Sor_realm = data.frame(Sor_realm)
colnames(Sor_realm) = c("Sor", "lo_Sor", "hi_Sor")
Sor_realm$realm = realm_labels

Sim_realm = cbind(c(mean(dist.EA$beta.SIM), mean(dist.VM$beta.SIM), mean(dist.MK$beta.SIM)), c(min(dist.EA$beta.SIM), min(dist.VM$beta.SIM), min(dist.MK$beta.SIM)), c(max(dist.EA$beta.SIM), max(dist.VM$beta.SIM), max(dist.MK$beta.SIM)))
Sim_realm = data.frame(Sim_realm)
colnames(Sim_realm) = c("Sim", "lo_Sim", "hi_Sim")
Sim_realm$realm = realm_labels

Sne_realm = cbind(c(mean(dist.EA$beta.SNE), mean(dist.VM$beta.SNE), mean(dist.MK$beta.SNE)), c(min(dist.EA$beta.SNE), min(dist.VM$beta.SNE), min(dist.MK$beta.SNE)), c(max(dist.EA$beta.SNE), max(dist.VM$beta.SNE), max(dist.MK$beta.SNE)))
Sne_realm = data.frame(Sne_realm)
colnames(Sne_realm) = c("Sne", "lo_Sne", "hi_Sne")
Sne_realm$realm = realm_labels



#por colecciones
coord_col = unique(datos_pbdb[,c("collection_no","paleolng","paleolat")])
coord_col = coord_col[rowSums(mat_col)>3,]
dist_lineales_col = as.dist(distm(coord_col[!is.na(coord_col$paleolat),2:3]))
dist_lat_col = dist(coord_col$paleolat[!is.na(coord_col$paleolat)])

#por placas tectónicas
geoplates=unique(datos_pbdb$geoplate)
coord_geoplate=NULL
for (i in 1:length(geoplates)) {
  temp_lng = median(datos_pbdb$paleolng[datos_pbdb$geoplate == geoplates [i]], na.rm=TRUE)
  temp_lat = median(datos_pbdb$paleolat[datos_pbdb$geoplate == geoplates [i]], na.rm=TRUE)    
  coord_geoplate = rbind(coord_geoplate, c(geoplates[i], temp_lng, temp_lat))
}
coord_geoplate = data.frame(coord_geoplate)
colnames(coord_geoplate) = c("geoplate", "paleolng", "paleolat")
coord_geoplate$geoplate = geoplates 
coord_geoplate = coord_geoplate[!is.na(coord_geoplate$paleolat),]

dist_lineales_geo = as.dist(distm(coord_geoplate[,2:3]))

mat_geoplate = matricear.region(datos_pbdb[!is.na(datos_pbdb$paleolat),], "geoplate")
data.mat_geoplate <- beta.pair(mat_geoplate)


###############################################################################################
#Test de Mantel distancias vs similitudes

dist_lat <- dist(median_lat)
dist_lineal <- dist_lineales_basin / 1000


mantel_sor_lat <- mantel(dist_lat, data.mat_basin$beta.sor)
mantel_sim_lat <- mantel(dist_lat, data.mat_basin$beta.sim)
mantel_sne_lat <- mantel(dist_lat, data.mat_basin$beta.sne)

mantel_sor_lineal <- mantel(dist_lineal, data.mat_basin$beta.sor)
mantel_sim_lineal <- mantel(dist_lineal, data.mat_basin$beta.sim)
mantel_sne_lineal <- mantel(dist_lineal, data.mat_basin$beta.sne)

#tabla
library(dplyr)

# Crear tabla resumen
tabla_mantel <- data.frame(
  Componente = c("beta.sor", "beta.sim", "beta.sne",
                 "beta.sor", "beta.sim", "beta.sne"),
  Distancia = c("Latitud", "Latitud", "Latitud",
                "Lineal (km)", "Lineal (km)", "Lineal (km)"),
  r_Mantel = c(mantel_sor_lat$statistic,
               mantel_sim_lat$statistic,
               mantel_sne_lat$statistic,
               mantel_sor_lineal$statistic,
               mantel_sim_lineal$statistic,
               mantel_sne_lineal$statistic),
  p_value = c(mantel_sor_lat$signif,
              mantel_sim_lat$signif,
              mantel_sne_lat$signif,
              mantel_sor_lineal$signif,
              mantel_sim_lineal$signif,
              mantel_sne_lineal$signif)
)

tabla_mantel

###############################################################################################
#ANÁLISIS MULTIVARIADOS
#PAM------------------------
pam_sim_basin <- pam(data.mat_basin$beta.sim, k = 3, diss = TRUE)
pam_anosim <- anosim(data.mat_basin$beta.sim, pam_sim_basin$cluster, permutations=4999)

pairwise_anosim <- NULL
comp <- NULL
for (i in 1:3) {
  temp_id <- pam_sim_basin$cluster != i 
  temp_dist <- as.dist(as.matrix(data.mat_basin$beta.sim)[temp_id,temp_id])
  temp <- anosim(temp_dist, pam_sim_basin$cluster[temp_id], permutations=4999)
  pairwise_anosim <- rbind(pairwise_anosim, c(temp$statistic, temp$signif))
  comp_temp = unique(loc_realm[temp_id])
  comp = cbind(comp, paste0(paste0(comp_temp[1],"-",comp_temp[2])))
}

pairwise_anosim <- data.frame(pairwise_anosim)
rownames(pairwise_anosim) <- comp
colnames(pairwise_anosim) <- c("R", "p-val")


#PAM con bootstrap
pam.sim <- function(matriz.orig, k, boot=100, rare, anosim.test=FALSE, cluster.method) {
  out.cluster <- NULL
  out.anosim <- NULL
  pamSim <- NULL
  rare <- ceiling(rare)
  for (i in 1:boot) {
    matriz <- matriz.orig[, sample(c(1:ncol(matriz.orig)), rare, replace=FALSE)]
    
    while( (sum(colSums(matriz) == 0) > 0) | (sum(rowSums(matriz) == 0) > 0) ) {
      matriz <- matriz[rowSums(matriz) > 0, colSums(matriz) > 0] 
    }
    
    disSim <- beta.pair(matriz)$beta.sim
    if(cluster.method=="pam") pamSim <- pam(disSim, k, diss=TRUE)
    if(cluster.method=="cluster") pamSim$cluster <- cutree(hclust(disSim, "ward.D2"), k)
    if (anosim.test) {
      anosim.temp <- anosim(disSim, grouping = pamSim$cluster)
      out.anosim <- rbind(out.anosim, c(anosim.temp$statistic, anosim.temp$signif))
    }
    
    out.cluster <- rbind(out.cluster, pamSim$cluster)
  }
  
  out <- list(data.frame(out.cluster), data.frame(out.anosim))
  names(out) <- c("clusters", "anosim")
  return(out)
}

pam_sim_boot <- pam.sim(mat_basin, k=3, boot=1000, rare=(ncol(mat_basin)*(2/3)), anosim.test=FALSE, cluster.method="pam")
pam_prob <- apply(pam_sim_boot[[1]], MARGIN=2, FUN=function(x) (c(sum(x==1), sum(x==2), sum(x==3))/length(x)) )  


#cantidad de clusters
anosim.R = NULL
for (i in 2:10) {
  pam.temp <-  pam(data.mat_basin$beta.sim, k = i, diss = TRUE)
  clust.temp <- cutree(hclust(data.mat_basin$beta.sim, "ward.D2"), k=i)
  anosim.pam.temp <- anosim(data.mat_basin$beta.sim, grouping = pam.temp$cluster)
  anosim.clust.temp <- anosim(data.mat_basin$beta.sim, grouping = clust.temp)
  anosim.R <- rbind(anosim.R, c(i, anosim.pam.temp$statistic, anosim.clust.temp$statistic))
}



#NMDS-------------------------------
mat_basin.nmds <- metaMDS(mat_basin)
site_scores_nmds <- scores(mat_basin.nmds, display = "sites")
site_scores <- scores(mat_basin.ca, display = "sites", scaling = 2)


###############################################################################################

###Plots###


#Map------------------------------------------------------------------------------------------

coastlines <- rgplates::reconstruct("coastlines", age = 410, model = "MERDITH2021")
edge <- mapedge()

epsg <- "ESRI:54030"

coastsRob <- sf::st_transform(coastlines, crs = epsg)
edgeRob   <- sf::st_transform(edge, crs = epsg)

# -------------------------
# 2. todas las ocurrencias
# -------------------------
df_occ <- data.frame(
  basin = datos_pbdb$basin,
  lon = as.numeric(datos_pbdb[, 26]),
  lat = as.numeric(datos_pbdb[, 27])
)

df_occ <- df_occ[complete.cases(df_occ), ]

# reconstrucción paleogeográfica
occ_paleocoord <- rgplates::reconstruct(
  as.matrix(df_occ[, c("lon", "lat")]),
  age = 410,
  model = "MERDITH2021"
)

# conservar solo filas válidas
ok <- complete.cases(occ_paleocoord)
df_occ <- df_occ[ok, ]
occ_paleocoord <- occ_paleocoord[ok, ]

# colores por reino
df_occ$color <- "black"
df_occ$color[df_occ$basin %in% c("BO", "ZA", "BR")] <- "darkgreen"
df_occ$color[df_occ$basin %in% c("AP", "NV", "VC", "AM")] <- "red"

# convertir a sf y proyectar
occ_sf <- sf::st_as_sf(
  data.frame(
    basin = df_occ$basin,
    color = df_occ$color,
    lon = occ_paleocoord[, 1],
    lat = occ_paleocoord[, 2]
  ),
  coords = c("lon", "lat"),
  crs = 4326
)

occRob <- sf::st_transform(occ_sf, crs = epsg)

# -------------------------
# 3. medianas por cuenca
# -------------------------
cuencas <- unique(datos_pbdb$basin)

df_medianas <- data.frame(
  basin = cuencas,
  lon = median_lng,
  lat = median_lat
)

df_medianas <- df_medianas[complete.cases(df_medianas), ]

# colores por reino
df_medianas$color <- "black"
df_medianas$color[df_medianas$basin %in% c("BO", "ZA", "BR")] <- "darkgreen"
df_medianas$color[df_medianas$basin %in% c("AP", "NV", "VC", "AM")] <- "red"

# convertir a sf y proyectar
med_sf <- sf::st_as_sf(df_medianas, coords = c("lon", "lat"), crs = 4326)
medRob <- sf::st_transform(med_sf, crs = epsg)

coords_med <- sf::st_coordinates(medRob)

# -------------------------
# 4. gráfico
# -------------------------
plot(edgeRob, col = "#1A6BB0", border = "gray30")
plot(coastsRob, border = NA, col = "gray90", add = TRUE)

# todas las ocurrencias
plot(sf::st_geometry(occRob),
     col = occRob$color,
     pch = 16,
     cex = 0.45,
     add = TRUE)

# medianas por cuenca
plot(sf::st_geometry(medRob),
     col = df_medianas$color,
     pch = 19,
     cex = 1.2,
     add = TRUE)

# etiquetas de las cuencas
text(coords_med[, 1], coords_med[, 2],
     labels = df_medianas$basin,
     pos = 4,
     offset = 0.35,
     cex = 0.8,
     col = "black")

box_colors <- c("#548B54", "#EE3B3B", "gray42")   # verde MK, rojo EA, gris OW
group_levels <- c("MK", "EA", "OW")

# ----------- Beta Sørensen -------------
beta_values <- c(dist.MK$beta.SOR, dist.EA$beta.SOR, dist.VM$beta.SOR)
groups <- factor(rep(group_levels,
                     times = c(length(dist.MK$beta.SOR),
                               length(dist.EA$beta.SOR),
                               length(dist.VM$beta.SOR))),
                 levels = group_levels)

par(bg = "gray95", mar = c(5, 5, 4, 1), las = 1, bty = "o", tck = -0.015,
    cex.axis = 0.9, cex.lab = 1.2, cex.main = 1.4,
    col.axis = "gray20", col.lab = "gray10", col.main = "black")

plot(1, type = "n", xlim = c(0.5, 3.5), ylim = range(beta_values),
     xlab = "Realm",
     ylab = expression(paste("Sørensen dissimilarity ", beta[SOR], )),
     main = expression(paste("Internal structure of each realm: ", beta[SOR])),
     axes = FALSE)
rect(0.5, par("usr")[3], 3.5, par("usr")[4], col = "gray90", border = NA)
grid(nx = NULL, ny = NULL, col = "white", lty = 1, lwd = 1)
boxplot(beta_values ~ groups,
        col = box_colors,
        border = "gray14",
        medcol = "gray14",
        boxwex = 0.4, medlwd = 2, whisklwd = 1, staplewex = 0.4,
        outpch = 16, outcex = 0.6, outcol = "gray40",
        add = TRUE, axes = TRUE, frame = FALSE)

# ----------- Beta Simpson -------------
beta_values <- c(dist.MK$beta.SIM, dist.EA$beta.SIM, dist.VM$beta.SIM)
groups <- factor(rep(group_levels,
                     times = c(length(dist.MK$beta.SIM),
                               length(dist.EA$beta.SIM),
                               length(dist.VM$beta.SIM))),
                 levels = group_levels)

par(bg = "gray95", mar = c(5, 5, 4, 1), las = 1, bty = "o", tck = -0.015,
    cex.axis = 0.9, cex.lab = 1.2, cex.main = 1.4,
    col.axis = "gray20", col.lab = "gray10", col.main = "black")

plot(1, type = "n", xlim = c(0.5, 3.5), ylim = range(beta_values),
     xlab = "Realm",
     ylab = expression(paste("Simpson dissimilarity ", beta[SIM], )),
     main = expression(paste("Internal structure of each realm: ", beta[SIM])),
     axes = FALSE)
rect(0.5, par("usr")[3], 3.5, par("usr")[4], col = "gray90", border = NA)
grid(nx = NULL, ny = NULL, col = "white", lty = 1, lwd = 1)

boxplot(beta_values ~ groups,
        col = box_colors,
        border = "gray14",
        medcol = "gray14",
        boxwex = 0.4, medlwd = 2, whisklwd = 1, staplewex = 0.4,
        outpch = 16, outcex = 0.6, outcol = "gray40",
        add = TRUE, axes = TRUE, frame = FALSE)


# ----------- Beta Nestedness -------------
jpeg("FIGURE 12.jpg", width = 1200, height = 800, res = 150)
beta_values <- c(dist.MK$beta.SNE, dist.EA$beta.SNE, dist.VM$beta.SNE)
groups <- factor(rep(group_levels,
                     times = c(length(dist.MK$beta.SNE),
                               length(dist.EA$beta.SNE),
                               length(dist.VM$beta.SNE))),
                 levels = group_levels)

par(bg = "gray95", mar = c(5, 5, 4, 1), las = 1, bty = "o", tck = -0.015,
    cex.axis = 0.9, cex.lab = 1.2, cex.main = 1.4,
    col.axis = "gray20", col.lab = "gray10", col.main = "black")

plot(1, type = "n", xlim = c(0.5, 3.5), ylim = range(beta_values),
     xlab = "Realm",
     ylab = expression(paste("Nestedness dissimilarity ", beta[SNE], )),
     main = expression(paste("Internal structure of each realm: ", beta[SNE])),
     axes = FALSE)

rect(0.5, par("usr")[3], 3.5, par("usr")[4], col = "gray90", border = NA)
grid(nx = NULL, ny = NULL, col = "white", lty = 1, lwd = 1)

boxplot(beta_values ~ groups,
        col = box_colors,
        border = "gray14",
        medcol = "gray14",
        boxwex = 0.4, medlwd = 2, whisklwd = 1, staplewex = 0.4,
        outpch = 16, outcex = 0.6, outcol = "gray40",
        add = TRUE, axes = TRUE, frame = FALSE)


# Paleta de colores
colores_realm <- c("#EE3B3B", "gray42", "#548B54")
par(mfrow=c(1,3), pty="s", mar=c(4.5,4.5,3,2), bg = "gray90", bty = "o")

# --- Sørensen ---
plot(median_lat_realm, c(Sor_realm$Sor), ylim=c(0.94, max(Sor_realm$hi_Sor)), xlim=c(-90, 0),
     pch=19, cex=1.5, col=colores_realm, 
     xlab="paleolatitude (º)", ylab="Sørensen")
for (i in 1:nrow(quantile_lat_realm)) {
  segments(x0=quantile_lat_realm[i,1], x1=quantile_lat_realm[i,2], 
           y0=Sor_realm$Sor[i], y1=Sor_realm$Sor[i], col=colores_realm[i])
  segments(x0=median_lat_realm[i], x1=median_lat_realm[i], 
           y0=Sor_realm$lo_Sor[i], y1=Sor_realm$hi_Sor[i], col=colores_realm[i])
}
points(median_lat_realm, Sor_realm$Sor, pch=19, cex=1.5, col=colores_realm)
mtext("A", side=3, line=1, cex=2, at=-91)

# --- Simpson ---
plot(median_lat_realm, c(Sim_realm$Sim), ylim=c(0.88, max(Sim_realm$hi_Sim)), xlim=c(-90, 0),
     pch=17, cex=1.5, col=colores_realm,
     xlab="paleolatitude (º)", ylab="Simpson", main="Realms")
for (i in 1:nrow(quantile_lat_realm)) {
  segments(x0=quantile_lat_realm[i,1], x1=quantile_lat_realm[i,2], 
           y0=Sim_realm$Sim[i], y1=Sim_realm$Sim[i], col=colores_realm[i])
  segments(x0=median_lat_realm[i], x1=median_lat_realm[i], 
           y0=Sim_realm$lo_Sim[i], y1=Sim_realm$hi_Sim[i], col=colores_realm[i])
}
points(median_lat_realm, Sim_realm$Sim, pch=17, cex=1.5, col=colores_realm)
mtext("B", side=3, line=1, cex=2, at=-91)

# --- Anidamiento ---
plot(median_lat_realm, c(Sne_realm$Sne), ylim=c(0.01, max(Sne_realm$hi_Sne)), xlim=c(-90, 0),
     pch=17, cex=1.5, col=colores_realm, 
     xlab="paleolatitude (º)", ylab="Nestedness")
for (i in 1:nrow(quantile_lat_realm)) {
  segments(x0=quantile_lat_realm[i,1], x1=quantile_lat_realm[i,2], 
           y0=Sne_realm$Sne[i], y1=Sne_realm$Sne[i], col=colores_realm[i])
  segments(x0=median_lat_realm[i], x1=median_lat_realm[i], 
           y0=Sne_realm$lo_Sne[i], y1=Sne_realm$hi_Sne[i], col=colores_realm[i])
}
points(median_lat_realm, Sne_realm$Sne, pch=17, cex=1.5, col=colores_realm)
mtext("C", side=3, line=1, cex=2, at=-91)


#indices vs latitud y dist. lineales

par(mfrow = c(1, 2),
    bg = "gray95",
    mar = c(4.5, 5, 4, 3),
    las = 1,
    bty = "n",
    tck = -0.015,
    cex.axis = 0.95,
    cex.lab = 1.2,
    col.axis = "gray20",
    col.lab = "gray10",
    font.lab = 2,
    pty = "s")  


# --- PANEL A: Disimilitud vs ΔLatitud ---
x_vals <- dist(median_lat)

plot(x_vals, data.mat_basin$beta.sor,
     type = "p", pch = 21, col = "black", bg = "black",
     xlim = c(0, 90), ylim = c(0, 1),
     ylab = expression(beta[SOR]~"/"~beta[SIM]~"/"~beta[SNE]),
     xlab = expression(Delta~Latitude),
     cex = 0.8,
     xaxt = "n", yaxt = "n")

axis(1, at = pretty(x_vals), lwd = 1, col = NA, col.ticks = "gray40", col.axis = "gray20")
axis(3, at = pretty(x_vals), lwd = 1, col = NA, col.ticks = "gray40", labels = FALSE)
axis(2, at = seq(0, 1, 0.2), lwd = 1, col = NA, col.ticks = "gray40", col.axis = "gray20")
axis(4, at = seq(0, 1, 0.2), lwd = 1, col = NA, col.ticks = "gray40", col.axis = "gray20", labels = FALSE)
abline(h = seq(0, 1, 0.2), col = "gray80", lty = "dotted")
abline(v = pretty(x_vals), col = "gray80", lty = "dotted")

points(x_vals, data.mat_basin$beta.sim,
       pch = 22, col = "red", bg = "pink", cex = 0.8)
points(x_vals, data.mat_basin$beta.sne,
       pch = 24, col = "darkgreen", bg = "lightgreen", cex = 0.8)


legend("topright",
       inset = c(-0.02, -0.02),  # ← la aleja un poco del gráfico
       legend = c(expression(beta[SOR]), expression(beta[SIM]), expression(beta[SNE])),
       col = c("black", "red", "darkgreen"),
       pt.bg = c("black", "pink", "lightgreen"),
       pch = c(21, 22, 24),
       pt.cex = 0.9,
       bty = "o",                # con borde
       box.lwd = 1.2,            # grosor del borde
       box.col = "gray30",       # color del borde
       bg = adjustcolor("white", alpha.f = 0.8),  # fondo semitransparente
       text.font = 2,            # texto en negrita
       cex = 0.9,
       text.col = "gray20")

mtext("A", side = 3, line = 1.2, adj = 0.5, cex = 2, font = 2)


# --- PANEL B: Disimilitud vs Distancia lineal ---
x_vals <- dist_lineales_basin / 1000

plot(x_vals, data.mat_basin$beta.sor,
     type = "p", pch = 21, col = "black", bg = "black",
     ylim = c(0, 1),
     xlab = expression(Delta~Distance~"(km)"),
     ylab = expression(beta[SOR]~"/"~beta[SIM]~"/"~beta[SNE]),
     cex = 0.9,
     xaxt = "n", yaxt = "n")

axis(1, at = pretty(x_vals), lwd = 1, col = NA, col.ticks = "gray40", col.axis = "gray20")
axis(3, at = pretty(x_vals), lwd = 1, col = NA, col.ticks = "gray40", labels = FALSE)
axis(2, at = seq(0, 1, 0.2), lwd = 1, col = NA, col.ticks = "gray40", col.axis = "gray20")
axis(4, at = seq(0, 1, 0.2), lwd = 1, col = NA, col.ticks = "gray40", col.axis = "gray20", labels = FALSE)
abline(h = seq(0, 1, 0.2), col = "gray80", lty = "dotted")
abline(v = pretty(x_vals), col = "gray80", lty = "dotted")

points(x_vals, data.mat_basin$beta.sim,
       pch = 22, col = "red", bg = "pink", cex = 0.9)
points(x_vals, data.mat_basin$beta.sne,
       pch = 24, col = "darkgreen", bg = "lightgreen", cex = 0.9)

# Leyenda arriba a la derecha, igual estilo
legend("topright",
       inset = c(-0.02, -0.02),
       legend = c(expression(beta[SOR]), expression(beta[SIM]), expression(beta[SNE])),
       col = c("black", "red", "darkgreen"),
       pt.bg = c("black", "pink", "lightgreen"),
       pch = c(21, 22, 24),
       pt.cex = 0.9,
       bty = "o",
       box.lwd = 1.2,
       box.col = "gray30",
       bg = adjustcolor("white", alpha.f = 0.8),
       text.font = 2,
       cex = 0.9,
       text.col = "gray20")

mtext("B", side = 3, line = 1.2, adj = 0.5, cex = 2, font = 2)


#NMDS
point_col <- "black"
# Datos eje 1 vs latitud
x_lat <- median_lat
y_nmds1 <- mat_basin.nmds$points[,1]

layout(matrix(c(1, 2), nrow = 1))

par(bg = "gray95",
    mar = c(4.5, 4.5, 3, 1),
    las = 1,
    bty = "o",
    tck = -0.015,
    col.axis = "gray20",
    col.lab = "gray10",
    cex.axis = 1,
    cex.lab = 1.1)

# -------------------
# Panel A: NMDS
# -------------------
plot(site_scores_nmds,
     type = "n",
     xlab = "NMDS Axis 1",
     ylab = "NMDS Axis 2",
     main = paste("NMDS (stress =", round(mat_basin.nmds$stress, 3), ")"),
     asp = 1)

grid(col = "white", lwd = 1)

points(site_scores_nmds,
       pch = 21,
       bg = point_col,
       col = "gray45",
       cex = 1.5)

text(site_scores_nmds,
     labels = localidades,
     pos = 3,
     cex = 0.75,
     col = "gray20")

mtext("A", side = 3, line = 1, adj = 0, cex = 1.3)

# -------------------
# Panel B: NMDS Axis 1 vs Paleolatitude
# -------------------
plot(x_lat, y_nmds1,
     pch = 21,
     bg = point_col,
     col = "gray45",
     cex = 1.5,
     xlab = "Paleolatitude (°S)",
     ylab = "NMDS Axis 1",
     main = "NMDS1 vs Paleolatitude")

grid(col = "white", lwd = 1)


text(x_lat, y_nmds1,
     labels = localidades,
     pos = 3,
     cex = 0.75,
     col = "gray20")

mtext("B", side = 3, line = 1, adj = 0, cex = 1.3)

#PAM----------------------------------------------------------------------------------------------------------------------------------
jpeg("FIGURE_9.jpg", width = 1000, height = 900, res = 150)
reinos <- c( "Easterm Americas", 
             "Old Word",
             "Malvinoxhosian")

colores <- c("red3","black", "darkgreen")
colores_cluster <- colores[pam_sim_basin$clustering]

par(bg = "gray95", mar = c(4.5, 5, 3, 2), las = 1, bty = "o", 
    col.axis = "gray30", col.lab = "gray20", cex.axis = 0.9, cex.lab = 1.1)
plot(x = median_lng, y = median_lat,
     col = colores_cluster,
     ylim = c(-90, 90), xlim = c(-180, 180),
     pch = 21, bg = colores_cluster, cex = 1.8, lwd = 1.2,
     text(median_lng, median_lat, labels = localidades, pos = 3, cex = 0.8, col = "gray20"),
     xlab = "Paleolongitude (°)", ylab = "Paleolatitude (°)",
     main = "Clustering PAM(βsim)")
grid(nx = NA, ny = NULL, col = "white", lty = 1)
legend("bottomleft",
       legend = reinos,
       pt.bg = colores,
       pch = 21, pt.cex = 1.8, bty = "n", cex = 0.9)

#ANOSIM-----------------------------------------------------------------------------------------------------------------------------------
plot(anosim.R[,1], anosim.R[,2], type="l")
lines(anosim.R[,1], anosim.R[,3], col=3)

# Parámetros gráficos
par(bg = "gray95", mar = c(5,5,3,2), las = 1, bty = "n",
    col.axis = "gray30", col.lab = "gray20",
    cex.axis = 1, cex.lab = 1.2)

# Gráfico base
plot(anosim.R[,1], anosim.R[,2],
     type = "l",
     lwd = 3,
     col = "#1b9e77",
     xlab = "Permutations",
     ylab = "ANOSIM R statistic",
     main = "ANOSIM R across Permutations")

grid(nx = NA, ny = NULL, col = "white", lwd = 1)

# Segunda línea
lines(anosim.R[,1], anosim.R[,3],
      col = "#d95f02",
      lwd = 3,
      lty = 2)


# Leyenda
legend("topright",
       legend = c("Observed R", "Expected / Permuted R"),
       col = c("#1b9e77", "#d95f02"),
       lwd = 3,
       lty = c(1,2),
       bty = "n",
       cex = 1)

