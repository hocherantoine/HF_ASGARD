# ------------------------------------------------------------------------------
# Script Name: ASGARDS_HF.R
# Author: Antoine Hocher
# Date: 2023-10-16
# Description: Retrieve genomes and Histones from asgards
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

######################################
#PART I : Retrieve ASGARD genomeslist from gtdb and then retrieve genomes and proteomes from NCBI
# And reformat them
#######################


GTDB_Genomes=read.table("~/ASGARD/230929_Asgards_gtdb-search.tsv",header=T,sep="\t")

Sys.setenv(PATH="/Users/ahocher/anaconda3/bin:/Users/ahocher/anaconda3/condabin:/usr/local/bin:/System/Cryptexes/App/usr/bin:/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/ncbi/blast/bin")






#Re compute all proteomes using prodigal
setwd("~/ASGARD/GENOMES/PROTEINS/PRODIGAL/")


FastaGenomeDir="~/ASGARD/GENOMES/DNA/"
ProteomeDir="~/ASGARD/GENOMES/PROTEINS/"


GenomeFiles=dir(FastaGenomeDir,pattern = "*.fna")
ProteomesFiles=dir(ProteomeDir,pattern = "*.faa")


AssemblyNb=unlist(lapply(GenomeFiles,function(x) strsplit(x,"_")[[1]][2]))
ProtAssemblyNb=unlist(lapply(ProteomesFiles,function(x) strsplit(x,"_")[[1]][2]))

ToComputeProdigal=AssemblyNb[-which(AssemblyNb%in%ProtAssemblyNb)]
setwd(FastaGenomeDir)
for(i in ToComputeProdigal){
  GenomeOfInterest=GenomeFiles[grep(pattern = i,GenomeFiles)]
  ProtName=sub(".fna","_prodigal.faa",GenomeOfInterest)
  GFFName=sub(".fna","_prodigal.gff",GenomeOfInterest)
  CommandProdigal=paste("/usr/local/bin/prodigal -i ",GenomeOfInterest," -f gff ","-o ",GFFName," -a ",ProtName,sep="")
  setwd(FastaGenomeDir)
  system(CommandProdigal)
}



#Rename proteins for proteomes computed with prodigal : 
library(seqinr)

setwd(FastaGenomeDir)
for(i in ToComputeProdigal){

    GenomeOfInterest=GenomeFiles[grep(pattern = i,GenomeFiles)]
    ProtName=sub(".fna","_prodigal.faa",GenomeOfInterest)
    
    A=read.fasta(ProtName,strip.desc = F)
    
    #Export directly within proteomes directory
    write.fasta(A,names = unlist(lapply(getName(A),function(x) paste(strsplit(x,"_")[[1]][1:2],collapse="-"))),file.out = paste(ProteomeDir,ProtName,sep="") )} 






#Re-annotate all proteins to include their assembly id within the protein id (useful for downstream)
setwd(ProteomeDir)

subDirName <- "AnnotatedGenomes"
mainDir=ProteomeDir
if (file.exists(subDirName)){
  setwd(file.path(mainDir))
} else {
  dir.create(file.path(mainDir, subDirName))
  setwd(file.path(mainDir))
}

subDir=paste(mainDir,subDirName,"/",sep="")



library(seqinr)


#List of all genome 

FastaList=dir(pattern = ".faa")
for (i in 1:length(FastaList)) {
  print(i)
  FastaName=FastaList[i]
  FastaNoext=sub('\\.faa$', '', FastaName) 
  
    print(FastaName)
    AssemblyName=paste(strsplit(FastaName,"_")[[1]][1],strsplit(FastaName,"_")[[1]][2],sep="_")
    
    A=seqinr::read.fasta(FastaName)
      
    write.fasta(sequences = A,names =     paste(AssemblyName,"$",getName(A),sep="")
,file.out = paste(subDir,FastaNoext,".faa",sep=""))}
      

#Renaming the Vienna SP b35 genome (L. ossiferum)
#simply because it wasn't i my initial list (of gtdb)

  FastaName="~/ASGARD/GENOMES/SP_B35/GCA_025839675.1_ASM2583967v1_protein.faa"
  
  FastaNoext=sub('\\.faa$', '', FastaName) 
  
  print(FastaName)
  AssemblyName=paste(strsplit(basename(FastaName),"_")[[1]][1],strsplit(basename(FastaName),"_")[[1]][2],sep="_")
  
  A=seqinr::read.fasta(FastaName)
  
  write.fasta(sequences = A,names =     paste(AssemblyName,"$",getName(A),sep=""),file.out = paste(subDir,basename(FastaNoext),".faa",sep=""))

#Once this is done, I concatenate all the proteomes into one big fasta
# 230930_Concat_ASGARD.faa

  
  
  
  
  
  #PART II GETTING THE HISTONES HITS
  
  
  
#Look for histones using hmm search option cut-ga when using pfam models

  #All those lines have to be run in the terminal ( or using system() in R)
#CBFD
system("/Users/ahocher/Documents/Softwares/hmmer-3.4/src/hmmsearch --cut_ga -A ~/ASGARD/HMMSEARCH/PF00808_vs-Asgard_ali.fa --tblout ~/ASGARD/HMMSEARCH/PF00808_vs-Asgard.txt ~/ASGARD/HMM/PF00808.hmm ~/ASGARD/GENOMES/PROTEINS/AnnotatedGenomes/230930_Concat_ASGARD.faa  > ~/ASGARD/HMMSEARCH/PF00808_vs-Asgard_full.txt")

#Histone core
/Users/ahocher/Documents/Softwares/hmmer-3.4/src/hmmsearch --cut_ga -A ~/ASGARD/HMMSEARCH/PF00125_vs-Asgard_ali.fa --tblout ~/ASGARD/HMMSEARCH/PF00125_vs-Asgard.txt ~/ASGARD/HMM/PF00125.hmm ~/ASGARD/GENOMES/PROTEINS/AnnotatedGenomes/230930_Concat_ASGARD.faa  > ~/ASGARD/HMMSEARCH/PF00125_vs-Asgard_full.txt 

#DUF1931
/Users/ahocher/Documents/Softwares/hmmer-3.4/src/hmmsearch --cut_ga -A ~/ASGARD/HMMSEARCH/PF09123_vs-Asgard_ali.fa --tblout ~/ASGARD/HMMSEARCH/PF09123_vs-Asgard.txt ~/ASGARD/HMM/PF09123.hmm ~/ASGARD/GENOMES/PROTEINS/AnnotatedGenomes/230930_Concat_ASGARD.faa  > ~/ASGARD/HMMSEARCH/PF09123_vs-Asgard_full.txt 
  

#Alba
/Users/ahocher/Documents/Softwares/hmmer-3.4/src/hmmsearch --cut_ga -A ~/ASGARD/HMMSEARCH/PF01918_vs-Asgard_ali.fa --tblout ~/ASGARD/HMMSEARCH/PF01918_vs-Asgard.txt ~/ASGARD/HMM/PF01918.hmm ~/ASGARD/GENOMES/PROTEINS/AnnotatedGenomes/230930_Concat_ASGARD.faa  > ~/ASGARD/HMMSEARCH/PF01918_vs-Asgard_full.txt 

  
  
  #Using H2A, using Nick Irvin HMM models (https://doi.org/10.1101/2023.09.20.558576)
  /Users/ahocher/Documents/Softwares/hmmer-3.4/src/hmmsearch --tblout ~/ASGARD/HMMSEARCH/H2A_vs_ASGARD.txt ~/ASGARD/EXTERNAL_DATA/Figshare_data_N_Irvin_Bioarkiv/HMMs/H2A.hmm ~/ASGARD/GENOMES/PROTEINS/AnnotatedGenomes/230930_Concat_ASGARD.faa  > ~/ASGARD/HMMSEARCH/H2A_vs_ASGARD_full.txt 

  /Users/ahocher/Documents/Softwares/hmmer-3.4/src/hmmsearch --tblout ~/ASGARD/HMMSEARCH/H2B_vs_ASGARD.txt ~/ASGARD/EXTERNAL_DATA/Figshare_data_N_Irvin_Bioarkiv/HMMs/H2B.hmm ~/ASGARD/GENOMES/PROTEINS/AnnotatedGenomes/230930_Concat_ASGARD.faa  > ~/ASGARD/HMMSEARCH/H2B_vs_ASGARD_full.txt 

  /Users/ahocher/Documents/Softwares/hmmer-3.4/src/hmmsearch --tblout ~/ASGARD/HMMSEARCH/H3_vs_ASGARD.txt ~/ASGARD/EXTERNAL_DATA/Figshare_data_N_Irvin_Bioarkiv/HMMs/H3.hmm ~/ASGARD/GENOMES/PROTEINS/AnnotatedGenomes/230930_Concat_ASGARD.faa  > ~/ASGARD/HMMSEARCH/H3_vs_ASGARD_full.txt 

  /Users/ahocher/Documents/Softwares/hmmer-3.4/src/hmmsearch -A ~/ASGARD/HMMSEARCH/H4_vs_ASGARD_HMMali.fa --tblout ~/ASGARD/HMMSEARCH/H4_vs_ASGARD.txt ~/ASGARD/EXTERNAL_DATA/Figshare_data_N_Irvin_Bioarkiv/HMMs/H4.hmm ~/ASGARD/GENOMES/PROTEINS/AnnotatedGenomes/230930_Concat_ASGARD.faa  > ~/ASGARD/HMMSEARCH/H4_vs_ASGARD_full.txt
  
  
  
  
  
  
  
  
  
  

#Reformat alignement
/Users/ahocher/Documents/Softwares/hmmer-3.4/easel/miniapps/esl-reformat a2m ~/ASGARD/HMMSEARCH/PF00808_vs-Asgard_ali.fa > ~/ASGARD/HMMSEARCH/PF00808_vs-Asgard_aligned_HMM_Domains.fa

/Users/ahocher/Documents/Softwares/hmmer-3.4/easel/miniapps/esl-reformat a2m ~/ASGARD/HMMSEARCH/PF09123_vs-Asgard_ali.fa > ~/ASGARD/HMMSEARCH/PF09123_vs-Asgard_aligned_HMM_Domains.fa
  
  /Users/ahocher/Documents/Softwares/hmmer-3.4/easel/miniapps/esl-reformat a2m ~/ASGARD/HMMSEARCH/H4_vs_ASGARD_HMMali.fa > ~/ASGARD/HMMSEARCH/H4_vs_ASGARD_HMM_aligned.fa
#Rtrieve the ssequences
#index
/Users/ahocher/Documents/Softwares/hmmer-3.4/easel/miniapps/esl-sfetch --index ~/ASGARD/GENOMES/PROTEINS/AnnotatedGenomes/230930_Concat_ASGARD.faa

#OUtput the proteins from the Hmm res file

grep -v "^#" ~/ASGARD/HMMSEARCH/PF00808_vs-Asgard.txt | awk '{print $1}' | /Users/ahocher/Documents/Softwares/hmmer-3.4/easel/miniapps/esl-sfetch -f ~/ASGARD/GENOMES/PROTEINS/AnnotatedGenomes/230930_Concat_ASGARD.faa - > ~/ASGARD/HMMSEARCH/PF00808_vs-Asgard_hits.fa

grep -v "^#" ~/ASGARD/HMMSEARCH/PF00125_vs-Asgard.txt | awk '{print $1}' | /Users/ahocher/Documents/Softwares/hmmer-3.4/easel/miniapps/esl-sfetch -f ~/ASGARD/GENOMES/PROTEINS/AnnotatedGenomes/230930_Concat_ASGARD.faa - > ~/ASGARD/HMMSEARCH/PF00125_vs-Asgard_hits.fa


grep -v "^#" ~/ASGARD/HMMSEARCH/PF09123_vs-Asgard.txt | awk '{print $1}' | /Users/ahocher/Documents/Softwares/hmmer-3.4/easel/miniapps/esl-sfetch -f ~/ASGARD/GENOMES/PROTEINS/AnnotatedGenomes/230930_Concat_ASGARD.faa - > ~/ASGARD/HMMSEARCH/PF09123_vs-Asgard_hits.fa


grep -v "^#" ~/ASGARD/HMMSEARCH/PF01918_vs-Asgard.txt | awk '{print $1}' | /Users/ahocher/Documents/Softwares/hmmer-3.4/easel/miniapps/esl-sfetch -f ~/ASGARD/GENOMES/PROTEINS/AnnotatedGenomes/230930_Concat_ASGARD.faa - > ~/ASGARD/HMMSEARCH/PF01918_vs-Asgard_hits.fa


#Nick's HMM models
grep -v "^#" ~/ASGARD/HMMSEARCH/H2A_vs_ASGARD.txt | awk '{print $1}' | /Users/ahocher/Documents/Softwares/hmmer-3.4/easel/miniapps/esl-sfetch -f ~/ASGARD/GENOMES/PROTEINS/AnnotatedGenomes/230930_Concat_ASGARD.faa - > ~/ASGARD/HMMSEARCH/H2A_vs-Asgard_hits.fa

grep -v "^#" ~/ASGARD/HMMSEARCH/H2B_vs_ASGARD.txt | awk '{print $1}' | /Users/ahocher/Documents/Softwares/hmmer-3.4/easel/miniapps/esl-sfetch -f ~/ASGARD/GENOMES/PROTEINS/AnnotatedGenomes/230930_Concat_ASGARD.faa - > ~/ASGARD/HMMSEARCH/H2B_vs-Asgard_hits.fa

grep -v "^#" ~/ASGARD/HMMSEARCH/H3_vs_ASGARD.txt | awk '{print $1}' | /Users/ahocher/Documents/Softwares/hmmer-3.4/easel/miniapps/esl-sfetch -f ~/ASGARD/GENOMES/PROTEINS/AnnotatedGenomes/230930_Concat_ASGARD.faa - > ~/ASGARD/HMMSEARCH/H3_vs-Asgard_hits.fa

grep -v "^#" ~/ASGARD/HMMSEARCH/H4_vs_ASGARD.txt | awk '{print $1}' | /Users/ahocher/Documents/Softwares/hmmer-3.4/easel/miniapps/esl-sfetch -f ~/ASGARD/GENOMES/PROTEINS/AnnotatedGenomes/230930_Concat_ASGARD.faa - > ~/ASGARD/HMMSEARCH/H4_vs-Asgard_hits.fa


#Re order by size: 
library(ape)

HistoneHits=read.fasta("~/ASGARD/HMMSEARCH/PF00808_vs-Asgard_hits.fa")

HistoneHits.r=HistoneHits[order(getLength(HistoneHits))]

write.fasta(HistoneHits.r,getName(HistoneHits.r),"~/ASGARD/HMMSEARCH/PF00808_vs-Asgard_hits_orderedLength.fa")

#re order by name
HistoneHits.r2=HistoneHits[order(getName(HistoneHits))]

write.fasta(HistoneHits.r2,getName(HistoneHits.r2),"~/ASGARD/HMMSEARCH/PF00808_vs-Asgard_hits_orderedName.fa")

#Check the overlap between the two hmm models from pfamA
HistoneHitsCore=read.fasta("~/ASGARD/HMMSEARCH/PF00125_vs-Asgard_hits.fa")

HistoneHitsCore.r=HistoneHitsCore[order(getLength(HistoneHitsCore))]

write.fasta(HistoneHitsCore.r,getName(HistoneHitsCore.r),"~/ASGARD/HMMSEARCH/PF00125_vs-Asgard_hits_orderedLength.fa")


#Overlap ?
table(getName(HistoneHitsCore)%in%getName(HistoneHits))
getName(HistoneHitsCore)[which(getName(HistoneHitsCore)%in%getName(HistoneHits)==F)]


#Check the overlap between the different HMM model made using eukaryotic histones

H2A=read.fasta("~/ASGARD/HMMSEARCH/H2A_vs-Asgard_hits.fa")
H2B=read.fasta("~/ASGARD/HMMSEARCH/H2B_vs-Asgard_hits.fa")
H3=read.fasta("~/ASGARD/HMMSEARCH/H3_vs-Asgard_hits.fa")
H4=read.fasta("~/ASGARD/HMMSEARCH/H4_vs-Asgard_hits.fa")

AllHF_eukHmm=c(H2A,H2B,H3,H4)
table(getName(HistoneHits)%in%c(getName(AllHF_eukHmm)))
table(c(getName(AllHF_eukHmm))%in%getName(HistoneHits))

NotinPfamA=unique(AllHF_eukHmm[which(!(getName(AllHF_eukHmm))%in%getName(HistoneHits))])

NotinPfamA=NotinPfamA[order(getLength(NotinPfamA))]
write.fasta(getSequence(NotinPfamA),getName(NotinPfamA),"~/ASGARD/HMMSEARCH/ASGARD_detected_w_EukHMM_not_PF808.fa")

H4_NotPFAMA=unique(H4[which(!(getName(H4))%in%getName(HistoneHits))])
H4_NotPFAMA=H4_NotPFAMA[order(getLength(H4_NotPFAMA))]
write.fasta(getSequence(H4_NotPFAMA),getName(H4_NotPFAMA),"~/ASGARD/HMMSEARCH/ASGARD_detected_with_H4_EukHMM_not_PF808.fa")


H4_NotPFAMA=unique(H4[which(!(getName(H4))%in%getName(HistoneHits))])
write.fasta(getSequence(H4_NotPFAMA),getName(H4_NotPFAMA),"~/ASGARD/HMMSEARCH/ASGARD_detected_with_H4_EukHMM_not_PF808.fa")

#After visual inspection and looking at the structures i keep only seauences <120aa 
ToExport=NotinPfamA[which(getLength(NotinPfamA)<120)]
write.fasta(getSequence(ToExport),getName(ToExport),"~/ASGARD/HMMSEARCH/ASGARD_detected_with_H4_EukHMM_not_PF808_inf120.fa")




#Concatenate all different histone fold hits:

cat ~/ASGARD/HMMSEARCH/ASGARD_detected_with_H4_EukHMM_not_PF808_inf120.fa ~/ASGARD/HMMSEARCH/PF00808_vs-Asgard_hits.fa ~/ASGARD/HMMSEARCH/PF09123_vs-Asgard_hits.fa > ~/ASGARD/HMMSEARCH/All_HF_Asgard.fa


#Remove special characters
#Remove duplicated hits from different HMMs

system(paste("tr -d '*' < ","~/ASGARD/HMMSEARCH/All_HF_Asgard.fa"," > ","~/ASGARD/HMMSEARCH/All_HF_Asgard_nostar.fa",sep=""))

ASGARDS_HF=read.fasta("~/ASGARD/HMMSEARCH/All_HF_Asgard_nostar.fa")
ASGARDS_HF.d=ASGARDS_HF[which(duplicated(ASGARDS_HF)==F)]
write.fasta(getSequence(ASGARDS_HF.d),getName(ASGARDS_HF.d),"~/ASGARD/HMMSEARCH/All_HF_Asgard_nodup.fa")

#Manually added HmfA, B h2A,H2B, h3, h4 from cerevisiae for outgroups
#Align
system("mafft-linsi --thread 16 --reorder ~/ASGARD/HMMSEARCH/All_HF_Asgard_nodup_HmfAB_HFcere.fa > ~/ASGARD/HMMSEARCH/All_HF_Asgard_nodup_HmfAB_HFcere_aligned.fa")


#Now this line is for a larger alignement that include euk histones and viral histones
#this alingment is made on the domains that match the hmm models so not the whole seq
#Collect those proteins to put them together for an HMM alignment output (-A from hmmsearch) only for larger tree building
HFASGARDs=getName(read.fasta("~/ASGARD/HMMSEARCH/All_HF_Asgard_nodup.fa"))

PF00808=read.fasta("~/ASGARD/HMMSEARCH/PF00808_vs-Asgard_aligned_HMM_Domains.fa")

PF09123=read.fasta("~/ASGARD/HMMSEARCH/PF09123_vs-Asgard_aligned_HMM_Domains.fa")
H4=read.fasta("~/ASGARD/HMMSEARCH/H4_vs_ASGARD_HMM_aligned.fa")

AllHMM=c(PF00808,PF09123,H4)
OriginalName=unlist(lapply(getName(AllHMM),function(x) strsplit(x,"\\/")[[1]][1]))
AllHMM.r=AllHMM[which(OriginalName%in%HFASGARDs)]

write.fasta(getSequence(AllHMM.r),getName(AllHMM.r),"~/ASGARD/HMMSEARCH/All_HF_Asgard_nodup_HMM_AliOnly.fa")








#PART III : TREES AND ANNOTATIONS




#Export itol table for tree:
source("/Users/ahocher/Documents/Softwares/table2itol-master/table2itol.R")


#Create a pruned tree from the larger GTDB tree
ProteomeDir="~/ASGARD/GENOMES/PROTEINS/AnnotatedGenomes/"
ProteomesFiles=list.files(ProteomeDir,pattern  =".faa")

library(ape)
#Loading the tree :

  GTDBTree=read.tree(file = "~/ASGARD/TREES/GTDB/ar53.tree")
  
  GTDB_info=read.delim(file = "~/ASGARD/TREES/GTDB/ar53_metadata.tsv",header=T,sep="\t",stringsAsFactors = F)
  
  #Assemblies in our dataset:
  ProtAssemblyNb=unlist(lapply(ProteomesFiles,function(x) strsplit(x,"_")[[1]][2]))
  Assemblies=paste("GCA_",ProtAssemblyNb,sep="")
  
  GTDB_info$ncbi_genbank_assembly_accession
  
  for( i in 1:length(GTDBTree$tip.label)){
    print(i)
    GTDBTree$tip.label[i]=GTDB_info[which(GTDB_info$accession==GTDBTree$tip.label[i]),]$ncbi_genbank_assembly_accession
  }
  
  #Prune tree to only have species included in our dataset :
  AssemblyToDROP=GTDBTree$tip.label[which(!(GTDBTree$tip.label%in%Assemblies))]
  GTDBTreeR=drop.tip(GTDBTree,AssemblyToDROP)
  
  
  #Export the tree with it's new name : 
  write.tree(GTDBTreeR,file="~/ASGARD/TREES/GTDB/ar53_asgards_pruned.tree")
  


GTDB_info$Phylum=unlist(lapply(GTDB_info$gtdb_taxonomy,function(x) strsplit(x,";")[[1]][2]))
#Metadata of our assembly set: 
GTDB_info.asgard=GTDB_info[which(GTDB_info$Phylum=="p__Asgardarchaeota"),]


  
#Now exporting a hit table
HitTable=data.frame(TID=GTDB_info.asgard$ncbi_genbank_assembly_accession,Has_HF=0,GTDB_info.asgard$checkm_completeness,GTDB_info.asgard$genome_size,GTDB_info.asgard$coding_density)


HFAssemblies=unlist(lapply(HFASGARDs,function(x) strsplit(x,"\\$")[[1]][1]))

DoubletHF=unlist(lapply(getName(read.fasta("~/ASGARD/HMMSEARCH/Manual_Doublet_and_doubletLike.fa")),function(x) strsplit(x,"\\$")[[1]][1]))

Clust1=unlist(lapply(getName(read.fasta("~/ASGARD/HMMSEARCH/Manual_LongHF_Cluster1.fa")),function(x) strsplit(x,"\\$")[[1]][1]))

Clust2=unlist(lapply(getName(read.fasta("~/ASGARD/HMMSEARCH/Manual_LongHF_Cluster2.fa")),function(x) strsplit(x,"\\$")[[1]][1]))

Alba=unlist(lapply(getName(read.fasta("~/ASGARD/HMMSEARCH/PF01918_vs-Asgard_hits.fa")),function(x) strsplit(x,"\\$")[[1]][1]))


SumUpHF=table(HFAssemblies)
HitTable$Has_HF=ifelse(HitTable$TID%in%HFAssemblies,1,0)
HitTable$Has_Doublet=ifelse(HitTable$TID%in%DoubletHF,1,0)
HitTable$Has_Clust1=ifelse(HitTable$TID%in%Clust1,1,0)
HitTable$Has_Clust2=ifelse(HitTable$TID%in%Clust2,1,0)
HitTable$Has_Alba=ifelse(HitTable$TID%in%Alba,1,0)

#Counting the HF nb
HitTable$HFNb=0
for(i in 1:length(SumUpHF)){
  #Because not all gtdb genomes are in their tree we filter
  if(length(which(HitTable$TID==names(SumUpHF)[i]))>0){
    #Attributing the length
  HitTable[which(HitTable$TID==names(SumUpHF)[i]),]$HFNb=as.numeric(SumUpHF[i])}
}

write.table(HitTable,"~/ASGARD/ITOL/GTDB_Asgard_info_table_allHF.txt",row.names=F,quote=F,sep="\t")

setwd("~/ASGARD/ITOL/")
create_itol_files("~/ASGARD/ITOL/GTDB_Asgard_info_table_allHF.txt",identifier = "TID",separator = "\t",label="TID")



#Exporting information of clusters for hodarchaeal histones: 

GTDB_info.hod=read.table("~/ASGARD/GTDB_Hodarchaea.txt",header=T,sep="\t",quote="",stringsAsFactors = F)
ClusterInfo=read.table("~/ASGARD/SEQUENCES/HF_Hodarchaea_Clusters_from_tree.txt",header=T,sep="\t",quote="",stringsAsFactors = F)
ClusterInfo$Assemblies=paste("GCA_",unlist(lapply(ClusterInfo$ID,function(x) strsplit(x,"\\ ")[[1]][2])),sep="")

GTDB_info.hod[,paste("Hod_has_",names(table(ClusterInfo$Cluster)),sep='')]=0
for(i in names(table(ClusterInfo$Cluster))){
  Ass=ClusterInfo[which(ClusterInfo$Cluster==i),]$Assemblies
  GTDB_info.hod[which(GTDB_info.hod$ncbi_genbank_assembly_accession%in%Ass),paste0("Hod_has_",i)]=1
  
}


ExportTable=GTDB_info.hod[,which(names(GTDB_info.hod)%in%c(paste("Hod_has_",names(table(ClusterInfo$Cluster)),sep='')))]
ExportTable$TID=ExportTable$ncbi_genbank_assembly_accession

write.table(ExportTable,"~/ASGARD/ITOL/GTDB_Hodarchaea_info_table_allHF_clusters.txt",row.names=F,quote=F,sep="\t")

setwd("~/ASGARD/ITOL/")
create_itol_files("~/ASGARD/ITOL/GTDB_Hodarchaea_info_table_allHF_clusters.txt",identifier = "TID",separator = "\t",label="TID")



#Same fo lokiarchaea 
GTDB_info.loki=GTDB_info[grep("c__Lokiarchaeia",GTDB_info$gtdb_taxonomy),]


ClusterInfo=read.table("~/ASGARD/SEQUENCES/HF_Lokiarchaea_Clusters_from_tree.txt",header=T,sep="\t",quote="",stringsAsFactors = F)
ClusterInfo$Assemblies=paste("GCA_",unlist(lapply(ClusterInfo$ID,function(x) strsplit(x,"\\ ")[[1]][2])),sep="")

GTDB_info.loki[,paste("Loki_has_",names(table(ClusterInfo$Cluster)),sep='')]=0
for(i in names(table(ClusterInfo$Cluster))){
  Ass=ClusterInfo[which(ClusterInfo$Cluster==i),]$Assemblies
  GTDB_info.loki[which(GTDB_info.loki$ncbi_genbank_assembly_accession%in%Ass),paste0("Loki_has_",i)]=1
  
}


ExportTable=GTDB_info.loki[,which(names(GTDB_info.loki)%in%c(paste("Loki_has_",names(table(ClusterInfo$Cluster)),sep='')))]
ExportTable$TID=ExportTable$ncbi_genbank_assembly_accession

write.table(ExportTable,"~/ASGARD/ITOL/GTDB_Lokiarchaea_info_table_allHF_clusters.txt",row.names=F,quote=F,sep="\t")

setwd("~/ASGARD/ITOL/")
create_itol_files("~/ASGARD/ITOL/GTDB_Lokiarchaea_info_table_allHF_clusters.txt",identifier = "TID",separator = "\t",label="TID")



#Creating colors for later : 
library(RColorBrewer)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

BInaryFiles=list.files(pattern = "iTOL_binary")

#Here we change all shapes to circles and colors
for(j in 1:length(BInaryFiles)){
  i=BInaryFiles[j]
  A=read.delim2(i,header=F,sep="\t",stringsAsFactors = F)
  A[c(4,6,9),2]=col_vector[j]
  A[c(8,11),2]=2
  A[A[,2]==0,2] = -1
  write.table(A,file=i,quote = F,row.names = F,col.names = F,sep="\t")}


#######################################################
#After manually curating the paralogs, I Annotate the paralog back on the tree for verification
#######################################################

ClusterInfo=read.table("~/ASGARD/SEQUENCES/HF_Hodarchaea_Clusters_from_tree.txt",header=T,sep="\t",quote="",stringsAsFactors = F)
ClusterInfo[which(ClusterInfo$ID=="GCA 016839785.1 JAEOTQ010000142.1-7"),]

ClusterInfo$TID=ClusterInfo$ID

ClusterInfo$TID=gsub(" ","_", ClusterInfo$TID)

write.table(ClusterInfo[,2:3],"~/ASGARD/ITOL/GTDB_Hodarchaea_HF_Clusters.txt",row.names=F,quote=F,sep="\t")


setwd("~/ASGARD/ITOL/")
create_itol_files("~/ASGARD/ITOL/GTDB_Hodarchaea_HF_Clusters.txt",identifier = "TID",separator = "\t",label="TID")


#For Lokiarchaea:

ClusterInfo=read.table("~/ASGARD/SEQUENCES/HF_Lokiarchaea_Clusters_from_tree.txt",header=T,sep="\t",quote="",stringsAsFactors = F)

ClusterInfo$TID=ClusterInfo$ID

ClusterInfo$TID=gsub(" ","_", ClusterInfo$TID)
names(ClusterInfo)[2]="LokiCluster"
ClusterInfo$LokiCluster=paste("c_",ClusterInfo$LokiCluster,sep="")
write.table(ClusterInfo[,2:3],"~/ASGARD/ITOL/GTDB_Loki_HF_Clusters.txt",row.names=F,quote=F,sep="\t")


setwd("~/ASGARD/ITOL/")
create_itol_files("~/ASGARD/ITOL/GTDB_Loki_HF_Clusters.txt",identifier = "TID",separator = "\t",label="TID")










  




#Retrieve Hodarchaeal histones only
GTDB_info.asgard.Hod=GTDB_info[grep("o__Hodarchaeales",GTDB_info$gtdb_taxonomy),]


write.table(GTDB_info.asgard.Hod,file="~/ASGARD/GTDB_Hodarchaea.txt",row.names=F,sep="\t",quote=F)

#For iTol : move the right GFF

GFF_toMove=list.files("~/ASGARD/GENOMES/GFF/PRODIGAL_ASGARDS/MODS",pattern = ".gff")

GFF_toMove=gsub(".gff","",GFF_toMove)
GFF_toMoveList=GFF_toMove[which(GFF_toMove%in%GTDB_info.asgard.Hod$ncbi_genbank_assembly_accession)]
GFF_toMoveList=paste0(GFF_toMoveList,".gff")
GFF_toMoveList=paste0("~/ASGARD/GENOMES/GFF/PRODIGAL_ASGARDS/MODS/",GFF_toMoveList)

targetdir="~/ASGARD/GENOMES/GFF/HODARCHAEA/WITH_PRODIGAL"
file.copy(from=GFF_toMoveList, to=targetdir)


#Reload histones seq
Histones=read.fasta("~/ASGARD/HMMSEARCH/All_HF_Asgard_nodup.fa")
Assemblies=unlist(lapply(getName(Histones),function(x) strsplit(x,"\\$")[[1]][1]))
HistonesHods=Histones[which(Assemblies%in%GTDB_info.hod$ncbi_genbank_assembly_accession)]

HistonesHods=HistonesHods[order(getLength(HistonesHods))]
write.fasta(HistonesHods,getName(HistonesHods),"~/ASGARD/HMMSEARCH/All_HF_o__Hodarchaeales.fa")

HistonesHodsAssemblies=unlist(lapply(getName(HistonesHods),function(x) strsplit(x,"\\$")[[1]][1]))
mean(table(HistonesHodsAssemblies))
mean(getLength(HistonesHods)[which(getLength(HistonesHods)<130)])


#
GTDB_info.loki=GTDB_info.asgard[grep("c__Lokiarchaeia",GTDB_info.asgard$gtdb_taxonomy),]

#Reload histones seq
Histones=read.fasta("~/ASGARD/HMMSEARCH/All_HF_Asgard_nodup.fa")
Assemblies=unlist(lapply(getName(Histones),function(x) strsplit(x,"\\$")[[1]][1]))
Histones.Loki=Histones[which(Assemblies%in%c(GTDB_info.loki$ncbi_genbank_assembly_accession,"GCA_025839675.1"))]

Histones.Loki=Histones.Loki[order(getLength(Histones.Loki))]
write.fasta(Histones.Loki,getName(Histones.Loki),"~/ASGARD/HMMSEARCH/All_HF_Asgard_hits_ordered_c__Lokiarchaeia_and_Ossifereum.fa")





#
N_Tails=read.fasta("~/ASGARD/HMMSEARCH/All_ASGARD_N_Tails_AlignmentDefined.fa",seqtype = "AA")
N_Tails.o=N_Tails[order(getLength(N_Tails))]
write.fasta(getSequence(N_Tails.o),getName(N_Tails.o),"~/ASGARD/HMMSEARCH/All_ASGARD_N_Tails_AlignmentDefined_lengtho.fa")

#align hoarchaeal histones
mafft-linsi --reorder ~/ASGARD/HMMSEARCH/PF00808_vs-Asgard_hits_ordered_o__Hodarchaeales.fa > ~/ASGARD/HMMSEARCH/PF00808_vs-Asgard_hits_ordered_o__Hodarchaeales_aligned.fa

#Aligne Lokiarchaeal Histones
mafft-linsi --reorder ~/ASGARD/HMMSEARCH/PF00808_vs-Asgard_hits_ordered_c__Lokiarchaeia.fa > ~/ASGARD/HMMSEARCH/PF00808_vs-Asgard_hits_ordered_c__Lokiarchaeia_aligned_linsi.fa




#BUild a tree, auto model 
cd ~/ASGARD/TREES

#All the parameters have been happily copied from Nick :) 
#https://www.biorxiv.org/content/10.1101/2023.09.20.558576v1.full.pdf

#(run in the terminal)
/Users/ahocher/Documents/Softwares/iqtree-2.2.2.6-MacOSX/bin/iqtree2 -pers 0.2 -nstop 500 -alrt 1000 -nt AUTO -mem 8G -m TEST -s ~/ASGARD/HMMSEARCH/PF00808_vs-Asgard_hits_ordered_o__Hodarchaeales_aligned.fa

#Build a lokiarchaeal histones tree (run in the terminal)
cd ~/ASGARD/TREES/IQTREE/LOKI_HF
/Users/ahocher/Documents/Softwares/iqtree-2.2.2.6-MacOSX/bin/iqtree2 -pers 0.2 -nstop 500 -alrt 1000 -nt AUTO -mem 8G -m TEST -s ~/ASGARD/HMMSEARCH/PF00808_vs-Asgard_hits_ordered_c__Lokiarchaeia_aligned_linsi.fa











######PART III Retrive archaeal histone for comparison with the rest of asgards
##############
#Retrieving Archaeal HIstones
/Users/ahocher/Documents/Softwares/hmmer-3.4/src/hmmsearch --cut_ga -A ~/ASGARD/HMMSEARCH/PF00808_vs-ArchaeaEBMC2_ali.fa --tblout ~/ASGARD/HMMSEARCH/PF00808_vs-ArchaeaEBMC2.txt ~/ASGARD/HMM/PF00808.hmm ~/Final_analysis_archaeal_chromatin/SEQ_DB/DB_Archaea95_update012020/genome_assemblies_prot_fasta/ncbi-genomes-2020-05-14/AnnotatedGenomes/concatArchaeaEBMC2.faa  > ~/ASGARD/HMMSEARCH/PF00808_vs-ArchaeaEBMC2.txt 


#Change hmm alignement output format
/Users/ahocher/Documents/Softwares/hmmer-3.4/easel/miniapps/esl-reformat a2m ~/ASGARD/HMMSEARCH/PF00808_vs-ArchaeaEBMC2_ali.fa > ~/ASGARD/HMMSEARCH/PF00808_vs-ArchaeaEBMC2_ali_HMM_domain.fa


#Cluster
mmseqs easy-cluster ~/ASGARD/HMMSEARCH/PF00808_vs-ArchaeaEBMC2_ali_HMM_domain.fa ~/ASGARD/SEQUENCES/MMSEQS_ARCHAEA/MMSEQS_ARCHAEA ~/ASGARD/SEQUENCES/MMSEQS_ARCHAEA/MMSEQS_ARCHAEA  --min-seq-id 0.75 -c 0.6 --cov-mode 0 -s 1
  
  
#Rtrieve the ssequences
#index
/Users/ahocher/Documents/Softwares/hmmer-3.4/easel/miniapps/esl-sfetch --index ~/Final_analysis_archaeal_chromatin/SEQ_DB/DB_Archaea95_update012020/genome_assemblies_prot_fasta/ncbi-genomes-2020-05-14/AnnotatedGenomes/concatArchaeaEBMC2.faa

#OUtput from the Hmm res file

grep -v "^#" ~/ASGARD/HMMSEARCH/PF00808_vs-ArchaeaEBMC2.txt | awk '{print $1}' | /Users/ahocher/Documents/Softwares/hmmer-3.4/easel/miniapps/esl-sfetch -f ~/Final_analysis_archaeal_chromatin/SEQ_DB/DB_Archaea95_update012020/genome_assemblies_prot_fasta/ncbi-genomes-2020-05-14/AnnotatedGenomes/concatArchaeaEBMC2.faa - > ~/ASGARD/HMMSEARCH/PF00808_vs-ArchaeaEBMC2_hits.fa



#Aim here is to remove asgard histones and reorder by length
ArchaeaHF=read.fasta("~/ASGARD/HMMSEARCH/PF00808_vs-ArchaeaEBMC2_hits.fa")

ArchaeaHFAssemblies=unlist(lapply(getName(ArchaeaHF),function(x) strsplit(x,"\\$")[[1]][2]))

GTDB_info$Phylum=unlist(lapply(GTDB_info$gtdb_taxonomy,function(x) strsplit(x,"\\;")[[1]][2]))

#Only include seauences in gtdb and Remove asgards

ArchaeaHFinGTDB=ArchaeaHF[which(ArchaeaHFAssemblies%in%GTDB_info$ncbi_genbank_assembly_accession & !(ArchaeaHFAssemblies%in%GTDB_info[which(GTDB_info$Phylum%in%c("p__Asgardarchaeota","p__Halobacteriota")),]$ncbi_genbank_assembly_accession))]


ArchaeaHFinGTDB.o=ArchaeaHFinGTDB[order(getLength(ArchaeaHFinGTDB))]
#Export
write.fasta(getSequence(ArchaeaHFinGTDB.o),getName(ArchaeaHFinGTDB.o),"~/ASGARD/HMMSEARCH/PF00808_vs-ArchaeaEBMC2_hits_InGTDB_noHalo_noASGARDS.fa")

#Align 
mafft-linsi --thread 8 ~/ASGARD/HMMSEARCH/PF00808_vs-ArchaeaEBMC2_hits_InGTDB_noHalo_noASGARDS.fa > ~/ASGARD/HMMSEARCH/PF00808_vs-ArchaeaEBMC2_hits_InGTDB_noHalo_noASGARDS_ali.fa




#Exporting GTDB taxonomy for the tree

FastaPath="~/ASGARD/HMMSEARCH/All_HF_Asgard_Aligned.fa"
ExportTable=data.frame(TID=getName(read.fasta(FastaPath)))
ExportTable$Assemblies=unlist(lapply(ExportTable$TID,function(x) strsplit(x,"\\$")[[1]][1]))
GTDB_info$order=unlist(lapply(GTDB_info$gtdb_taxonomy,function(x) strsplit(x,"\\;")[[1]][4]))

GTDB_info$clade=unlist(lapply(GTDB_info$gtdb_taxonomy,function(x) strsplit(x,"\\;")[[1]][3]))
GTDB_info$family=unlist(lapply(GTDB_info$gtdb_taxonomy,function(x) strsplit(x,"\\;")[[1]][5]))

ExportTable.m=merge(ExportTable,GTDB_info,by.x="Assemblies",by.y="ncbi_genbank_assembly_accession")  

ExportTable.m$TID=gsub("\\$","_",ExportTable.m$TID)
write.table(ExportTable.m[,which(names(ExportTable.m)%in%c("TID","order","clade","family"))],"~/ASGARD/ITOL/HF_ASGARD_order.txt",row.names = F,col.names = T,sep="\t",quote=F)
source("/Users/ahocher/Dropbox/Scripts/table2itol-master/table2itol.R")

setwd("~/ASGARD/ITOL/")
create_itol_files("HF_ASGARD_order.txt",identifier = "TID")



#Families for loki clade

ExportTable.m=merge(ExportTable,GTDB_info,by.x="Assemblies",by.y="ncbi_genbank_assembly_accession")  

ExportTable.m$TID=gsub("\\$","_",ExportTable.m$TID)
write.table(ExportTable.m[which(ExportTable.m$clade=="c__Lokiarchaeia"),which(names(ExportTable.m)%in%c("TID","family"))],"~/ASGARD/ITOL/HF_ASGARD_Family_LokIOnly.txt",row.names = F,col.names = T,sep="\t",quote=F)
source("/Users/ahocher/Dropbox/Scripts/table2itol-master/table2itol.R")

setwd("~/ASGARD/ITOL/")
create_itol_files("~/ASGARD/ITOL/HF_ASGARD_Family_LokIOnly.txt",identifier = "TID")
