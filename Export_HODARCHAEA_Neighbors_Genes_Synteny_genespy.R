# ------------------------------------------------------------------------------
# Script Name: Export_Neighbors_fasta.R
# Author: Antoine Hocher
# Date: 2023-10-10
# Description: Function to export protein sequences of neighbors from a list of fasta sequences, a gff folder and a proteome fasta
# Version: 1.0.0
# 
# Copyright (C) 2023 Antoine Hocher
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# ------------------------------------------------------------------------------







ProteomeFile="~/ASGARD/GENOMES/PROTEINS/AnnotatedGenomes/230930_Concat_ASGARD.faa"

#Protein names have been reset for compatibility using Fomart_Prodigal_and_Genbank_GFF_to_Genespy

GFFFolder="~/ASGARD/GENOMES/GFF/HODARCHAEA/ALL_MOD/"
ProteinOfInterestFile="~/ASGARD/HMMSEARCH/PF00808_vs-Asgard_hits_ordered_o__Hodarchaeales.fa"
ProteinList=getName(read.fasta(ProteinOfInterestFile))

setwd(GFFFolder)
GFFFiles=list.files(pattern = ".gff")
ExportProteinList=c()

for(i in ProteinList){
  print(i)
  Assembly=strsplit(i,"\\$")[[1]][1]
  Protein=strsplit(i,"\\$")[[1]][2]
  
  #The names had to be modified to be compatible with genespy
  ProteinNoDot=gsub("\\.1","",Protein)
  GFF=readGFF(GFFFiles[grep(Assembly,GFFFiles)])
  GFFCDS=GFF[which(GFF$type=="CDS"),]
  
  PoI_Coord=grep(ProteinNoDot,GFFCDS$ID)
  
  PoIGFF=GFFCDS[PoI_Coord,]
  Sequences_to_Export=GFFCDS[c((PoI_Coord-5):(PoI_Coord+5)),]
  Sequences_to_Export=Sequences_to_Export[which(as.character(Sequences_to_Export$seqid)==as.character(PoIGFF$seqid)),]$ID
  Sequences_to_Export=paste(Assembly,Sequences_to_Export,sep="$")
  ExportProteinList=c(ExportProteinList,Sequences_to_Export)
}

OriginalNames=gsub("\\-",".1-",ExportProteinList)

#Retrieve protein sequences
Proteome=read.fasta(ProteomeFile)
ProteinsToExport=Proteome[which(getName(Proteome)%in%OriginalNames)]

ExporFileName=gsub(".fa","_neighboringGenes.fa",ProteinOfInterestFile)
write.fasta(ProteinsToExport,getName(ProteinsToExport),ExporFileName)



#Cluster the sequences
NeighborhoodFolder="~/ASGARD/SEQUENCES/Gene_Neighboorhoods"
Commandline=paste0("mmseqs easy-cluster ",ExporFileName," ",NeighborhoodFolder," ",NeighborhoodFolder," --min-seq-id 0.1 -c 0.5 --cov-mode 0 -s 7")
system(Commandline)



#Retrieve the clusters
Clusterstable=read.table("~/ASGARD/SEQUENCES/Gene_Neighboorhoods/Gene_Neighboorhoods_cluster.tsv",header=F,sep="\t")
names(Clusterstable)=c("Representative","names")
ClusterCount=table(Clusterstable$Representative)
ClusterCount=ClusterCount[order(-ClusterCount)[1:50]]
ClusterFinal=Clusterstable[which(Clusterstable$Representative%in%names(ClusterCount)),]
ClusterFinal$ClusterNum=NA
for(i in 1:50){
  ClusterFinal[which(ClusterFinal$Representative%in%names(ClusterCount)[i]),]$ClusterNum=i
  
}

ClusterFinal$ProteinNames=unlist(lapply(ClusterFinal$names,function(x) strsplit(x,"\\$")[[1]][2]))

ClusterFinal$ProteinNames=gsub(".1-","-",ClusterFinal$ProteinNames)

#Modify all GFF to add the colors
for(i in GFFFiles){
  GFF=readGFF(i)
  GFF$gene=NA
  Selecta=ClusterFinal$ProteinNames[which(ClusterFinal$ProteinNames%in%GFF$ID)]
  for(j in Selecta){
  GFF[which(GFF$ID==j),]$gene=ClusterFinal[which(ClusterFinal$ProteinNames==j),]$ClusterNum
  GFF[which(GFF$ID==j),]$Name=ClusterFinal[which(ClusterFinal$ProteinNames==j),]$ClusterNum
  
  
   }
  export.gff3(GFF,paste("~/ASGARD/GENOMES/GFF/HODARCHAEA/ALL_MOD_Clusters_HistonesSeqs/",i,sep=""))
}



#Externally we have to create a database in genespy with the modified GFFs





#Creating colors for later : 
library(RColorBrewer)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

#Generate a colorfile
List=c(1:100)
List[seq(1,100,by=2)]=paste0(">",col_vector[1:50])

List[seq(2,101,by=2)]=1:50

write.table(List,"~/ASGARD/ITOL/Genespy_Clusters_color.txt",row.names = F,col.names = F,sep="\t",quote=F)





#At last we modify the names of tree leaves as genespy has modified them

HodarchaealHF=read.fasta("~/ASGARD/HMMSEARCH/PF00808_vs-Asgard_hits_ordered_o__Hodarchaeales.fa")
HodarchaealHFNames=getName(HodarchaealHF)


iTOLgenespy=read.delim("~/ASGARD/ITOL/HF_Hodarchaea_itol_output.txt",header=F)
iTOLgenespyMod=iTOLgenespy
for(i in 5:dim(iTOLgenespy)[1]){
  print(i)
  String=iTOLgenespy[i,1]
  FullName=strsplit(String,"\\,")[[1]][1]
  ProteinName=strsplit(FullName,"\\_")[[1]][3]
  
  #For protein names from prodigal -those contain a "-" we need to further modifiy the names due to name changes iimposed by genespy
  if(length(grep("\\-",ProteinName))>0){
    ProteinName=gsub("\\-",".1-",ProteinName)
  }
  
  CorrectedName=HodarchaealHFNames[grep(ProteinName,HodarchaealHFNames)]
  
  #iTol also modifies names and refuses $ so we replace them by "_ :) 
  
  CorrectedName2=gsub("\\$","_",CorrectedName)
  
  
  String2=gsub(FullName,CorrectedName2,String)
  iTOLgenespyMod[i,1]=String2
}
write.table(iTOLgenespyMod,"~/ASGARD/ITOL/HF_Hodarchaea_itol_Mod.txt",row.names = F,col.names = F,sep="\t",quote=F)

