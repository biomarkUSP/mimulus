#'----------------------------------------------------------------------------------------#
#' Authors: Rafael M. Yassue
#'          Ana Letycia B. Garcia
#'          Nathalia S. Silva
#'          Germano Costa Neto
#'----------------------------------------------------------------------------------------#

library(onemap)
dados <-read_mapmaker( file = "m_feb06.raw")

dados$geno
str(dados)
class(dados)

#'----------------------------------------------------------------------------------------#
# Segregation type visualization
#'----------------------------------------------------------------------------------------#

plot_by_segreg_type(dados)

#'----------------------------------------------------------------------------------------#
# Sugested LOD score estimation
#'----------------------------------------------------------------------------------------#

LOD_sug <- suggest_lod(dados) 

#'----------------------------------------------------------------------------------------#
# Map function
#'----------------------------------------------------------------------------------------#

set_map_fun(type="kosambi")

#'----------------------------------------------------------------------------------------#
# Segregation test
#'----------------------------------------------------------------------------------------#

f2_test <-test_segregation(x = dados) 
class(f2_test)
print(f2_test)

plot(f2_test)

#'----------------------------------------------------------------------------------------#
# Checking for reduntant markers
#'----------------------------------------------------------------------------------------#

bins <- find_bins(dados) 


# First attempt to assembly linkage groups
tp <- rf_2pts(dados, LOD=6, max.rf=.35)
Seque <- make_seq(tp, arg="all")
LGall <- group(Seque)

# Second attempt 
tp <- rf_2pts(dados, LOD=6, max.rf=.5)
Seque <- make_seq(tp, arg="all")
LGall <- group(Seque)

# Third attempt
tp <- rf_2pts(dados, LOD=7, max.rf=.35)
Seque <- make_seq(tp, arg="all")
LGall <- group(Seque)

# Fourth attempt
tp <- rf_2pts(dados, LOD=6, max.rf=.25)
Seque <- make_seq(tp, arg="all")
LGall <- group(Seque)

## If we would like remove distorted markers according to segregation test
#remo<-select_segreg(f2_test, distorted = TRUE)#é preciso manter estes marcadores
#a <- match(remo, colnames(dados$geno))
#f2<-drop_marker(Seque, a)
#'----------------------------------------------------------------------------------------#

# A diagnostics to remove markers 
# A loop to create and order linkage groups
grupos <- list()

tp <- rf_2pts(dados2, LOD=6, max.rf=.25)
Seque <- make_seq(tp, arg="all")
LGall <- group(f2)
rm(grupos)
grupos <- list()
for (i in 1:14){ 
  grupos[[i]] <- make_seq(LGall,i) 
  grupos[[i]]  <- order_seq(input.seq= grupos[[i]], n.init = 4)
  grupos[[i]]<- make_seq(grupos[[i]] , "force")
}

#'----------------------------------------------------------------------------------------#
# Heatmap to check distorted markers and take decisions ad hoc 
#'----------------------------------------------------------------------------------------#

rf_graph_table(grupos[[1]] , main ="LG1", inter=F) 
rf_graph_table(grupos[[2]] , main ="LG2", inter=F ) 
rf_graph_table(grupos[[3]] , main ="LG3", inter=F ) 
rf_graph_table(grupos[[4]] , main ="LG4", inter=F ) 
rf_graph_table(grupos[[5]] , main ="LG5", inter=F ) 
rf_graph_table(grupos[[6]] , main ="LG6", inter=F ) 
rf_graph_table(grupos[[7]] , main ="LG7", inter=F ) 
rf_graph_table(grupos[[8]] , main ="LG8", inter=F ) 
rf_graph_table(grupos[[9]] , main ="LG9", inter=F ) 
rf_graph_table(grupos[[10]] , main ="LG10", inter=F ) 
rf_graph_table(grupos[[11]] , main ="LG11", inter=F ) 
rf_graph_table(grupos[[12]] , main ="LG12", inter=F ) 
rf_graph_table(grupos[[13]] , main ="LG13", inter=F ) 
rf_graph_table(grupos[[14]] , main ="LG14", inter=F ) 

# Remove distorted marks after visual analysis of heatmap
tp <- rf_2pts(dados2, LOD=7, max.rf=.27)
Seque <- make_seq(tp, arg="all")

f2<-drop_marker(Seque, c("212","165","30","93","229","64","152","258","143","240","150",
                         "213","312","204","6","242","182","146","126", "276",
                         "53","145","81", "140", "44", "17", "405","376", "12",
                         "91" ))

LGall <- group(f2)

for (i in 1:14){ 
  grupos[[i]] <- make_seq(LGall,i) 
  grupos[[i]]  <- order_seq(input.seq= grupos[[i]], n.init = 4)
  grupos[[i]]<- make_seq(grupos[[i]] , "force")
}

# Final ordering
for (i in 1:14){ 
  cat("grupo de ligacao", i, "\n")
ripple_seq(grupos[[i]], ws = 5, LOD = 3)
}

#'----------------------------------------------------------------------------------------#
# Final heatmaps
#'----------------------------------------------------------------------------------------#

g1<-rf_graph_table(grupos[[1]] , main ="LG1", inter=F) 
g2<-rf_graph_table(grupos[[2]] , main ="LG2", inter=F ) 
g3<-rf_graph_table(grupos[[3]] , main ="LG3", inter=F ) 
g4<-rf_graph_table(grupos[[4]] , main ="LG4", inter=F ) 
g5<-rf_graph_table(grupos[[5]] , main ="LG5", inter=F ) 
g6<-rf_graph_table(grupos[[6]] , main ="LG6", inter=F ) 
g7<-rf_graph_table(grupos[[7]] , main ="LG7", inter=F ) 
g8<-rf_graph_table(grupos[[8]] , main ="LG8", inter=F ) 
g9<-rf_graph_table(grupos[[9]] , main ="LG9", inter=F ) 
g10<-rf_graph_table(grupos[[10]] , main ="LG10", inter=F ) 
g11<-rf_graph_table(grupos[[11]] , main ="LG11", inter=F ) 
g12<-rf_graph_table(grupos[[12]] , main ="LG12", inter=F ) 
g13<-rf_graph_table(grupos[[13]] , main ="LG13", inter=F ) 
g14<-rf_graph_table(grupos[[14]] , main ="LG14", inter=F ) 

# Final linkage map
maps<- list(grupos[[1]], grupos[[2]], grupos[[3]], grupos[[4]], grupos[[5]], grupos[[6]]
            , grupos[[7]], grupos[[8]], grupos[[9]], grupos[[10]], grupos[[11]], grupos[[12]], grupos[[13]], grupos[[14]])

map<-draw_map(maps, names=T, grid=T, cex.mrk=0.7)

library(ggpubr)

figurafinal<- ggarrange(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,
  nrow = 5,
  ncol = 3)

write_map(maps, file="test_sum.txt")

png(filename = "figurafinal2.png", width = 20, height = 30, units="in", res=450)
figurafinal
dev.off()
maps<- list(grupos[[1]], grupos[[2]], grupos[[3]], grupos[[4]], grupos[[5]], grupos[[6]]
            , grupos[[7]], grupos[[8]], grupos[[9]], grupos[[10]], grupos[[11]], grupos[[12]], grupos[[13]], grupos[[14]])

map <-draw_map(maps, names=T, grid=T, cex.mrk=0.7)

#'----------------------------------------------------------------------------------------#
# EDF