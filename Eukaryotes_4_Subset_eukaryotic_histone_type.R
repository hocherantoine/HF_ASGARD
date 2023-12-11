# ------------------------------------------------------------------------------
# Script Name: Subset_eukaryotic_histone_type.R
# Author: Antoine Hocher
# Date: 2023-10-16
# Description: split histones according to their hmmer scores
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



#Classify Using Nick Irvin HMM models from Nature micro
#Run in the terminal
hmmer-3.4/src/hmmsearch --tblout ~/ASGARD/HMMSEARCH/PF00125_vs-EUK_hits_H2A_HMM_Tblout.txt ~/ASGARD/EXTERNAL_DATA/Figshare_data_N_Irvin_Bioarkiv/HMMs/H2A.hmm ~/ASGARD/HMMSEARCH/PF00125_vs-EUK_hits_MMSEQS_no_duplicates_less386aa.fa  > ~/ASGARD/HMMSEARCH/PF00125_vs-EUK_hits_H2A_HMM__Full.txt

hmmer-3.4/src/hmmsearch --tblout ~/ASGARD/HMMSEARCH/PF00125_vs-EUK_hits_H2B_HMM_Tblout.txt ~/ASGARD/EXTERNAL_DATA/Figshare_data_N_Irvin_Bioarkiv/HMMs/H2B.hmm ~/ASGARD/HMMSEARCH/PF00125_vs-EUK_hits_MMSEQS_no_duplicates_less386aa.fa  > ~/ASGARD/HMMSEARCH/PF00125_vs-EUK_hits_H2B_HMM__Full.txt

hmmer-3.4/src/hmmsearch --tblout ~/ASGARD/HMMSEARCH/PF00125_vs-EUK_hits_H3_HMM_Tblout.txt ~/ASGARD/EXTERNAL_DATA/Figshare_data_N_Irvin_Bioarkiv/HMMs/H3.hmm ~/ASGARD/HMMSEARCH/PF00125_vs-EUK_hits_MMSEQS_no_duplicates_less386aa.fa  > ~/ASGARD/HMMSEARCH/PF00125_vs-EUK_hits_H3_HMM__Full.txt

hmmer-3.4/src/hmmsearch --tblout ~/ASGARD/HMMSEARCH/PF00125_vs-EUK_hits_H4_HMM_Tblout.txt ~/ASGARD/EXTERNAL_DATA/Figshare_data_N_Irvin_Bioarkiv/HMMs/H4.hmm ~/ASGARD/HMMSEARCH/PF00125_vs-EUK_hits_MMSEQS_no_duplicates_less386aa.fa  > ~/ASGARD/HMMSEARCH/PF00125_vs-EUK_hits_H4_HMM__Full.txt





AllHistones=read.fasta("~/ASGARD/HMMSEARCH/PF00125_vs-EUK_hits_MMSEQS_no_duplicates_less386aa.fa")
library(rhmmer)
H2A_Hits=read_tblout("~/ASGARD/HMMSEARCH/PF00125_vs-EUK_hits_H2A_HMM_Tblout.txt")
H2B_Hits=read_tblout("~/ASGARD/HMMSEARCH/PF00125_vs-EUK_hits_H2B_HMM_Tblout.txt")
H3_Hits=read_tblout("~/ASGARD/HMMSEARCH/PF00125_vs-EUK_hits_H3_HMM_Tblout.txt")
H4_Hits=read_tblout("~/ASGARD/HMMSEARCH/PF00125_vs-EUK_hits_H4_HMM_Tblout.txt")

AllHits=rbind(H2A_Hits,H2B_Hits,H3_Hits,H4_Hits)

#Filter on evalue

AllHits=AllHits[which(AllHits$sequence_evalue< 1e-10),]

AllHits$Attribution=NA
for(i in unique(AllHits$domain_name)){
Sub=AllHits[which(AllHits$domain_name==i),]
  TempAttribution=Sub[which.max(Sub$sequence_score),]$query_name
  TempAttribution=strsplit(TempAttribution,"\\.")[[1]][1]
AllHits[which(AllHits$domain_name==i),]$Attribution=TempAttribution
}

table(AllHits$Attribution)

#Export individual histones types
H2A=AllHistones[getName(AllHistones)%in%c(unique(AllHits[which(AllHits$Attribution=="H2A"),]$domain_name))]
H2B=AllHistones[getName(AllHistones)%in%c(unique(AllHits[which(AllHits$Attribution=="H2B"),]$domain_name))]
H3=AllHistones[getName(AllHistones)%in%c(unique(AllHits[which(AllHits$Attribution=="H3"),]$domain_name))]
H4=AllHistones[getName(AllHistones)%in%c(unique(AllHits[which(AllHits$Attribution=="H4"),]$domain_name))]

write.fasta(getSequence(H2A),getName(H2A),"~/ASGARD/HMMSEARCH/PF00125_vs-EUK_hits_H2A_HMM_Subset.fa")
write.fasta(getSequence(H2B),getName(H2B),"~/ASGARD/HMMSEARCH/PF00125_vs-EUK_hits_H2B_HMM_Subset.fa")
write.fasta(getSequence(H3),getName(H3),"~/ASGARD/HMMSEARCH/PF00125_vs-EUK_hits_H3_HMM_Subset.fa")
write.fasta(getSequence(H4),getName(H4),"~/ASGARD/HMMSEARCH/PF00125_vs-EUK_hits_H4_HMM_Subset.fa")

Input="~/ASGARD/HMMSEARCH/PF00125_vs-EUK_hits_H4_HMM_Subset.fa"


###############
#Align them all 
###############
Align=function(Input){
  Output=gsub(".fa","_mafftaligned.fa",Input)
  Commandline=paste0("mafft-linsi --thread 16 ",Input," > ",Output)
  system(Commandline)
}

Align("~/ASGARD/HMMSEARCH/PF00125_vs-EUK_hits_H2A_HMM_Subset.fa")
Align("~/ASGARD/HMMSEARCH/PF00125_vs-EUK_hits_H2B_HMM_Subset.fa")
Align("~/ASGARD/HMMSEARCH/PF00125_vs-EUK_hits_H3_HMM_Subset.fa")
Align("~/ASGARD/HMMSEARCH/PF00125_vs-EUK_hits_H4_HMM_Subset.fa")


########PART II : Same on the histones obtained using individual HMM models


AllHistones=read.fasta("~/ASGARD/HMMSEARCH/Individual_Euk_HF_Pooled_no_duplicates_NoPfAM_A_LengthOrdered_nostar.fa")
library(rhmmer)
H2A_Hits=read_tblout("~/ASGARD/HMMSEARCH/H2A_vs-eukprot_CLS.txt")
H2B_Hits=read_tblout("~/ASGARD/HMMSEARCH/H2B_vs-eukprot_CLS.txt")
H3_Hits=read_tblout("~/ASGARD/HMMSEARCH/H3_vs-eukprot_CLS.txt")
H4_Hits=read_tblout("~/ASGARD/HMMSEARCH/H4_vs-eukprot_CLS.txt")
NFYA_Hits=read_tblout("~/ASGARD/HMMSEARCH/PF00125_vs-EUK_hits_MMSEQS_no_duplicates_vsNFYC_tblout.txt")


AllHits=rbind(H2A_Hits,H2B_Hits,H3_Hits,H4_Hits,NFYA_Hits)

#Filter on evalue

AllHits=AllHits[which(AllHits$sequence_evalue< 1e-10),]

AllHits$Attribution=NA
for(i in unique(AllHits$domain_name)){
  Sub=AllHits[which(AllHits$domain_name==i),]
  TempAttribution=Sub[which.max(Sub$sequence_score),]$query_name
  TempAttribution=strsplit(TempAttribution,"\\.")[[1]][1]
  AllHits[which(AllHits$domain_name==i),]$Attribution=TempAttribution
}



##################

#Check for histones which have more than one non overlapping domain
#We will exclude them from the singlet histones
#Filter on evalue

Large=AllHits[which(AllHits$domain_name%in%getName(AllHistones[which(getLength(AllHistones)>120)])),]
Large=Large[which(duplicated(Large$domain_name)==F),]

table(Large$domain_number_env)

#Load Hmm alignemnts
setwd("~/ASGARD/HMMSEARCH/")
H2A_ali=read.fasta("H2A_vs-eukprot_CLS_HMM_ali.fa")
H2B_ali=read.fasta("H2B_vs-eukprot_CLS_HMM_ali.fa")
H3_ali=read.fasta("H3_vs-eukprot_CLS_HMM_ali.fa")
H4_ali=read.fasta("H4_vs-eukprot_CLS_HMM_ali.fa")
AllAliHmm=c(H2A_ali,H2B_ali,H3_ali,H4_ali)


AliStats=data.frame(ID=getName(AllAliHmm))
AliStats$BaseID=unlist(lapply(AliStats$ID,function(x) strsplit(x,"\\/")[[1]][1]))
AliStats$Coord=unlist(lapply(AliStats$ID,function(x) strsplit(x,"\\/")[[1]][2]))
AliStats$Start=as.numeric(unlist(lapply(AliStats$Coord,function(x) strsplit(x,"\\-")[[1]][1])))
AliStats$Stop=as.numeric(unlist(lapply(AliStats$Coord,function(x) strsplit(x,"\\-")[[1]][2])))
AliStats$Multiple=1

for(i in 1:dim(Large)[1]){
  NPoI=Large[i,]$domain_name
  if(dim(AliStats[which(AliStats$BaseID==NPoI),])[1]>1){
    MyProt=AliStats[which(AliStats$BaseID==NPoI),]
    
    MyProt=MyProt[order(MyProt$Start),]
    k=1
    for(i in 1:(dim(MyProt)[1]-1)){
      #Count non overlaping domains
      if(MyProt[i,]$Stop < MyProt[(i+1),]$Start){
        k=k+1}
    }
    AliStats[which(AliStats$BaseID==NPoI),]$Multiple=k
  }
}

DoubleHits=AliStats[which(AliStats$Multiple==2),]
DoubleHits$Hit=NA
#Annotate each domain: 
DoubleHits[which(DoubleHits$ID%in%getName(H2A_ali)),]$Hit="H2A"
DoubleHits[which(DoubleHits$ID%in%getName(H2B_ali)),]$Hit="H2B"
DoubleHits[which(DoubleHits$ID%in%getName(H3_ali)),]$Hit="H3"
DoubleHits[which(DoubleHits$ID%in%getName(H4_ali)),]$Hit="H4"
DoubleHits=DoubleHits[with(DoubleHits,order(BaseID,Start)),]
DoubleHits$Combination=NA
for(i in unique(DoubleHits$BaseID)){
  DoubleHits[which(DoubleHits$BaseID==i),]$Combination=paste(DoubleHits[which(DoubleHits$BaseID==i),]$Hit[1],DoubleHits[which(DoubleHits$BaseID==i),]$Hit[length(which(DoubleHits$BaseID==i))],sep="-")
  
  
}

for(i in names(table(DoubleHits$Combination))){
  DHits=AllHistones[which(getName(AllHistones)%in%DoubleHits[which(DoubleHits$Combination==i),]$BaseID)]
  
  write.fasta(getSequence(DHits),getName(DHits),paste("~/ASGARD/HMMSEARCH/Individual_Euk_HF_hits_PostProcessed_Doublets_",i,".fa",sep=""))
  
}

#Export All doublets
DHits=AllHistones[which(getName(AllHistones)%in%DoubleHits$BaseID)]

write.fasta(getSequence(DHits),getName(DHits),paste("~/ASGARD/HMMSEARCH/Individual_Euk_HF_hits_PostProcessed_Doublets.fa",sep=""))

#Export Triplets
TripleHits=AliStats[which(AliStats$Multiple==3),]
THits=AllHistones[which(getName(AllHistones)%in%TripleHits$BaseID)]
write.fasta(getSequence(THits),getName(THits),"~/ASGARD/HMMSEARCH/Individual_Euk_HF_hits_PostProcessed_Triplets.fa")

#Export Quadriplets !!!
QuadrupleHits=AliStats[which(AliStats$Multiple==4),]
QuadrupleHits=QuadrupleHits[order(QuadrupleHits$ID),]
QHits=AllHistones[which(getName(AllHistones)%in%QuadrupleHits$BaseID)]
write.fasta(getSequence(QHits),getName(THits),"~/ASGARD/HMMSEARCH/Individual_Euk_HF_hits_PostProcessed_Quadruplets.fa")



############
#Remove proteins with multiples Hits from the list of singlets>

AllHits=AllHits[-which(AllHits$domain_name%in%AliStats[which(AliStats$Multiple>1),]$BaseID),]

#Export individual histones types
H2A=AllHistones[getName(AllHistones)%in%c(unique(AllHits[which(AllHits$Attribution=="H2A"),]$domain_name))]
H2B=AllHistones[getName(AllHistones)%in%c(unique(AllHits[which(AllHits$Attribution=="H2B"),]$domain_name))]
H3=AllHistones[getName(AllHistones)%in%c(unique(AllHits[which(AllHits$Attribution=="H3"),]$domain_name))]
H4=AllHistones[getName(AllHistones)%in%c(unique(AllHits[which(AllHits$Attribution=="H4"),]$domain_name))]
NFYC=AllHistones[getName(AllHistones)%in%c(unique(AllHits[which(AllHits$Attribution=="Q13952_NFYC_HUMAN"),]$domain_name))]

length(getLength(H2A));length(getLength(H2B));length(getLength(H3));length(getLength(H4));length(getLength(NFYC))


#remove 5 % longuest seqeunces
H2AQt=quantile(getLength(H2A),probs=seq(0,1,by=0.05))[20]
H2BQt=quantile(getLength(H2B),probs=seq(0,1,by=0.05))[20]
H3Qt=quantile(getLength(H3),probs=seq(0,1,by=0.05))[20]
H4Qt=quantile(getLength(H4),probs=seq(0,1,by=0.05))[20]

H2A=H2A[which(getLength(H2A)<H2AQt)]
H2B=H2B[which(getLength(H2B)<H2BQt)]
H3=H3[which(getLength(H3)<H3Qt)]
H4=H4[which(getLength(H4)<H4Qt)]



write.fasta(getSequence(H2A),getName(H2A),"~/ASGARD/HMMSEARCH/Individual_Euk_HF_hits_PostProcessed_H2A_HMM_Subset.fa")
write.fasta(getSequence(H2B),getName(H2B),"~/ASGARD/HMMSEARCH/Individual_Euk_HF_hits_PostProcessed_H2B_HMM_Subset.fa")
write.fasta(getSequence(H3),getName(H3),"~/ASGARD/HMMSEARCH/Individual_Euk_HF_hits_PostProcessed_H3_HMM_Subset.fa")
write.fasta(getSequence(H4),getName(H4),"~/ASGARD/HMMSEARCH/Individual_Euk_HF_hits_PostProcessed_H4_HMM_Subset.fa")






###############
#Align them all 
###############
Align=function(Input){
  Output=gsub(".fa","_mafftaligned.fa",Input)
  Commandline=paste0("mafft-linsi --thread 16 ",Input," > ",Output)
  system(Commandline)
}





Align("~/ASGARD/HMMSEARCH/Individual_Euk_HF_hits_PostProcessed_H2A_HMM_Subset.fa")
Align("~/ASGARD/HMMSEARCH/Individual_Euk_HF_hits_PostProcessed_H2B_HMM_Subset.fa")
Align("~/ASGARD/HMMSEARCH/Individual_Euk_HF_hits_PostProcessed_H3_HMM_Subset.fa")
Align("~/ASGARD/HMMSEARCH/Individual_Euk_HF_hits_PostProcessed_H4_HMM_Subset.fa")


