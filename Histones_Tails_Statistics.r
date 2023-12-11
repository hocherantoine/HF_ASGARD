# ------------------------------------------------------------------------------
# Script Name: Tail_KmerCount.R
# Author: Antoine Hocher
# Date: 2023-10-16
# Description: Plots for Review on histones in asgard archaea
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

rm(list = ls())
library(seqinr)

setwd("~/ASGARD/HMMSEARCH/")
#Archaeal tails are defined based on HMM model hit boundaries

##########################
#ASGARDS, N-terminal tails
##########################
ASGARD_N_Ter_Tails=seqinr::read.fasta("~/ASGARD/HMMSEARCH/All_ASGARD_N_Tails_AlignmentDefined_lengtho_no_structured.fa",seqtype = "AA")


#Loki Tails
LokiHFNAmes=getName(read.fasta("~/ASGARD/HMMSEARCH/PF00808_vs-Asgard_hits_ordered_c__Lokiarchaeia.fa"))

#Hodarchaeal Tail
HodaHFNAmes=getName(read.fasta("~/ASGARD/HMMSEARCH/PF00808_vs-Asgard_hits_ordered_o__Hodarchaeales_aligned.fa"))


#######################
#Eukaryotic tails are defined based on alignment

#based on cerevisiae structures, only the most N-C terminal unstructured bits were split based on the corresponding alignement
#For N and ,Coords are respectively 
#H2A: 1-18, 99:..
#H2B:1-41,130:..
#H3:1-41;133:
#H4:1:31;94:
##########################
#Euk N-terminal tails
##########################
Euk_H2A_N_Ter_Tails=seqinr::read.fasta("Individual_Euk_HF_hits_PostProcessed_H2A_HMM_Subset_N_Tails.fa",seqtype = "AA")
Euk_H2A_N_Ter_Tails=Euk_H2A_N_Ter_Tails[order(getLength(Euk_H2A_N_Ter_Tails))]

Euk_H2B_N_Ter_Tails=seqinr::read.fasta("Individual_Euk_HF_hits_PostProcessed_H2B_HMM_Subset_N_Tails.fa",seqtype = "AA")
Euk_H2B_N_Ter_Tails=Euk_H2B_N_Ter_Tails[order(getLength(Euk_H2B_N_Ter_Tails))]

Euk_H3_N_Ter_Tails=seqinr::read.fasta("Individual_Euk_HF_hits_PostProcessed_H3_HMM_Subset_N_Tails.fa",seqtype = "AA")
Euk_H3_N_Ter_Tails=Euk_H3_N_Ter_Tails[order(getLength(Euk_H3_N_Ter_Tails))]

Euk_H4_N_Ter_Tails=seqinr::read.fasta("Individual_Euk_HF_hits_PostProcessed_H4_HMM_Subset_N_Tails.fa",seqtype = "AA")
Euk_H4_N_Ter_Tails=Euk_H4_N_Ter_Tails[order(getLength(Euk_H4_N_Ter_Tails))]

#just exporting the length ordered version to get a visual feeling of things
write.fasta(getSequence(Euk_H2A_N_Ter_Tails),getName(Euk_H2A_N_Ter_Tails),"Individual_Euk_HF_hits_PostProcessed_H2A_HMM_Subset_N_Tails_length_ordered.fa")
write.fasta(getSequence(Euk_H2B_N_Ter_Tails),getName(Euk_H2B_N_Ter_Tails),"Individual_Euk_HF_hits_PostProcessed_H2B_HMM_Subset_N_Tails_length_ordered.fa")
write.fasta(getSequence(Euk_H3_N_Ter_Tails),getName(Euk_H3_N_Ter_Tails),"Individual_Euk_HF_hits_PostProcessed_H3_HMM_Subset_N_Tails_length_ordered.fa")
write.fasta(getSequence(Euk_H4_N_Ter_Tails),getName(Euk_H4_N_Ter_Tails),"Individual_Euk_HF_hits_PostProcessed_H4_HMM_Subset_N_Tails_length_ordered.fa")


#pooling them all together
All_Euk_N_Tails=c(Euk_H2A_N_Ter_Tails,Euk_H2B_N_Ter_Tails,Euk_H3_N_Ter_Tails,Euk_H4_N_Ter_Tails)





#Safety check To check that we have most species in the CLS dataset of eukpro: 
Specie=unlist(lapply(getName(Euk_H4_N_Ter_Tails),function(x) strsplit(x,"\\$")[[1]][1]))
length(table(Specie))



#Start a table that will host various properties of histones
Tail_properties=data.frame(ID=getName(ASGARD_N_Ter_Tails),Length=getLength(ASGARD_N_Ter_Tails),Family="ASGARD",Segment="N")


Tail_properties=rbind(Tail_properties,data.frame(ID=getName(Euk_H2A_N_Ter_Tails),Length=getLength(Euk_H2A_N_Ter_Tails),Family="H2A",Segment="N"))
Tail_properties=rbind(Tail_properties,data.frame(ID=getName(Euk_H2B_N_Ter_Tails),Length=getLength(Euk_H2B_N_Ter_Tails),Family="H2B",Segment="N"))
Tail_properties=rbind(Tail_properties,data.frame(ID=getName(Euk_H3_N_Ter_Tails),Length=getLength(Euk_H3_N_Ter_Tails),Family="H3",Segment="N"))
Tail_properties=rbind(Tail_properties,data.frame(ID=getName(Euk_H4_N_Ter_Tails),Length=getLength(Euk_H4_N_Ter_Tails),Family="H4",Segment="N"))



#########################
#Now compute properties
#########################






AllNTails=c(ASGARD_N_Ter_Tails,Euk_H2A_N_Ter_Tails,Euk_H2B_N_Ter_Tails,Euk_H3_N_Ter_Tails,Euk_H4_N_Ter_Tails)


#Clip out the first methionine to exclude that from the stats:
#Here we show that as expected it's for the vast majority a methionine
FirstBase=unlist(lapply(getSequence(AllNTails), function(x) x[1]))
table(FirstBase)

which(FirstBase==">")
#Test
test=list(AllNTails[[100]],AllNTails[[101]],AllNTails[[102]])
lapply(getSequence(test), function(x) x[-1])
#Test with sequence of length 0 
test=list(AllNTails[[1]],AllNTails[[2]],AllNTails[[3]])
lapply(getSequence(test), function(x) x[-1])


ClippedSeq=lapply(getSequence(AllNTails), function(x) x[-1])

write.fasta(sequences = ClippedSeq,names = getName(AllNTails),"AllNTails_ASGARD_Euk_Clipped.fa",as.string = T)



#Important:
#To avoid a bug associated to sequences of length 0, one has to add a empty tab for all sequences of length 0 (either in the terminal or by using a software such as aliview)


#Reload Clipped sequences (using the file with added tabs for sequence of length 0)
AllNTails=seqinr::read.fasta("AllNTails_ASGARD_Euk_Clipped.fa",seqtype = "AA")

#Compute all kind of statistics using AAstats from seqinr
namesToAdd=names(unlist(AAstat(ASGARD_N_Ter_Tails[[100]],plot = F)))
Tail_properties[,namesToAdd]=NA

for(i in 1:dim(Tail_properties)[1]){
  print(i)
  if(Tail_properties[i,]$Length>5){
  Tail_properties[i,namesToAdd]=unlist(AAstat(getSequence(AllNTails[which(getName(AllNTails)==Tail_properties[i,]$ID)])[[1]],plot = F))
  }
}



#Normalise compo into proportion
Tail_properties[,5:25]=Tail_properties[,5:25]/Tail_properties$Length

names(Tail_properties)=gsub("Compo.","Prop.",names(Tail_properties))


#Compute some other metrics associated to IDR according to Zarin et al. 2019
Tail_properties$NetChargePerRes=(Tail_properties$Prop.R+Tail_properties$Prop.K)-(Tail_properties$Prop.D+Tail_properties$Prop.E)

Tail_properties$NetCharges=Tail_properties$Length*((Tail_properties$Prop.R+Tail_properties$Prop.K)-(Tail_properties$Prop.D+Tail_properties$Prop.E))

Tail_properties$RKratio=Tail_properties$Prop.R/Tail_properties$Prop.K
Tail_properties$EDratio=Tail_properties$Prop.E/Tail_properties$Prop.D


#Export that table
write.table(Tail_properties,"~/ASGARD/ANALYSIS/N_Terminal_Tails_Properties.txt",quote=F,row.names=F,sep="\t")

#Reload
Tail_properties=read.table("~/ASGARD/ANALYSIS/N_Terminal_Tails_Properties.txt",header=T,sep="\t",quote="",stringsAsFactors = F)
      



  
  #Load the protein alphabet
aa <- unique(unlist(getSequence(ASGARD_N_Ter_Tails)))
aa <- aa[aa != "*"]

#Make sure it's equal to 20:
length(aa)



# ========================================================================
# ========================================================================
#Compute K-mer counts
# ========================================================================
# ========================================================================

AllNTails.sub=AllNTails[which(getLength(AllNTails)>5)]

init=seqinr::count(AllNTails.sub[1],wordsize = 2,alphabet = aa)
BigList=lapply(AllNTails.sub,function(x) seqinr::count(x,wordsize = 2,alphabet = aa))
Dimers=as.data.frame(lapply(BigList,function(x) as.numeric(x)))
names(Dimers)=names(AllNTails.sub)
Dimers$ID=names(init)


head(Dimers[,1:3])

TailDimerSummary=data.frame(ID=Dimers$ID)

TailDimerSummary$LOKI=rowSums(Dimers[,which(names(Dimers)%in%LokiHFNAmes)])
TailDimerSummary$HOD=rowSums(Dimers[,which(names(Dimers)%in%HodaHFNAmes)])
TailDimerSummary$ASGARD=rowSums(Dimers[,which(names(Dimers)%in%getName(ASGARD_N_Ter_Tails))])
TailDimerSummary$EUK=rowSums(Dimers[,which(names(Dimers)%in%getName(All_Euk_N_Tails))])
TailDimerSummary$H2A=rowSums(Dimers[,which(names(Dimers)%in%getName(Euk_H2A_N_Ter_Tails))])
TailDimerSummary$H2B=rowSums(Dimers[,which(names(Dimers)%in%getName(Euk_H2B_N_Ter_Tails))])
TailDimerSummary$H3=rowSums(Dimers[,which(names(Dimers)%in%getName(Euk_H3_N_Ter_Tails))])
TailDimerSummary$H4=rowSums(Dimers[,which(names(Dimers)%in%getName(Euk_H4_N_Ter_Tails))])


TailDimerSummary$LOKI=100*TailDimerSummary$LOKI/sum(TailDimerSummary$LOKI)
TailDimerSummary$HOD=100*TailDimerSummary$HOD/sum(TailDimerSummary$HOD)
TailDimerSummary$ASGARD=100*TailDimerSummary$ASGARD/sum(TailDimerSummary$ASGARD)
TailDimerSummary$EUK=100*TailDimerSummary$EUK/sum(TailDimerSummary$EUK)
TailDimerSummary$H2A=100*TailDimerSummary$H2A/sum(TailDimerSummary$H2A)
TailDimerSummary$H2B=100*TailDimerSummary$H2B/sum(TailDimerSummary$H2B)
TailDimerSummary$H3=100*TailDimerSummary$H3/sum(TailDimerSummary$H3)
TailDimerSummary$H4=100*TailDimerSummary$H4/sum(TailDimerSummary$H4)


head(TailDimerSummary[order(-TailDimerSummary$ASGARD),],n=30)


write.table(TailDimerSummary,"~/ASGARD/ANALYSIS/N_Terminal_Tails_2merSummary.txt",quote=F,row.names=F,sep="\t")



# ========================================================================
#TriMers
# ========================================================================
AllNTails.sub=AllNTails[which(getLength(AllNTails)>5)]

init=seqinr::count(AllNTails.sub[1],wordsize = 3,alphabet = aa)
BigList=lapply(AllNTails.sub,function(x) seqinr::count(x,wordsize = 3,alphabet = aa))
Trimers=as.data.frame(lapply(BigList,function(x) as.numeric(x)))
names(Trimers)=names(AllNTails.sub)
Trimers$ID=names(init)



TailTrimerSummary=data.frame(ID=Trimers$ID)

TailTrimerSummary$LOKI=rowSums(Trimers[,which(names(Trimers)%in%LokiHFNAmes)])
TailTrimerSummary$HOD=rowSums(Trimers[,which(names(Trimers)%in%HodaHFNAmes)])
TailTrimerSummary$ASGARD=rowSums(Trimers[,which(names(Trimers)%in%getName(ASGARD_N_Ter_Tails))])
TailTrimerSummary$EUK=rowSums(Trimers[,which(names(Trimers)%in%getName(All_Euk_N_Tails))])
TailTrimerSummary$H2A=rowSums(Trimers[,which(names(Trimers)%in%getName(Euk_H2A_N_Ter_Tails))])
TailTrimerSummary$H2B=rowSums(Trimers[,which(names(Trimers)%in%getName(Euk_H2B_N_Ter_Tails))])
TailTrimerSummary$H3=rowSums(Trimers[,which(names(Trimers)%in%getName(Euk_H3_N_Ter_Tails))])
TailTrimerSummary$H4=rowSums(Trimers[,which(names(Trimers)%in%getName(Euk_H4_N_Ter_Tails))])


TailTrimerSummary$LOKI=100*TailTrimerSummary$LOKI/sum(TailTrimerSummary$LOKI)
TailTrimerSummary$HOD=100*TailTrimerSummary$HOD/sum(TailTrimerSummary$HOD)
TailTrimerSummary$ASGARD=100*TailTrimerSummary$ASGARD/sum(TailTrimerSummary$ASGARD)
TailTrimerSummary$EUK=100*TailTrimerSummary$EUK/sum(TailTrimerSummary$EUK)
TailTrimerSummary$H2A=100*TailTrimerSummary$H2A/sum(TailTrimerSummary$H2A)
TailTrimerSummary$H2B=100*TailTrimerSummary$H2B/sum(TailTrimerSummary$H2B)
TailTrimerSummary$H3=100*TailTrimerSummary$H3/sum(TailTrimerSummary$H3)
TailTrimerSummary$H4=100*TailTrimerSummary$H4/sum(TailTrimerSummary$H4)




TailTrimerSummary$Total=rowSums(TailTrimerSummary[,2:9])
TailTrimerSummary=TailTrimerSummary[-which(TailTrimerSummary$Total==0),]


write.table(TailTrimerSummary,"~/ASGARD/ANALYSIS/N_Terminal_Tails_3merSummary.txt",quote=F,row.names=F,sep="\t")




# ========================================================================
#TetraMers
# ========================================================================
AllNTails.sub=AllNTails[which(getLength(AllNTails)>5)]

init=seqinr::count(AllNTails.sub[1],wordsize = 4,alphabet = aa)
BigList=lapply(AllNTails.sub,function(x) seqinr::count(x,wordsize = 4,alphabet = aa))
Tetramers=as.data.frame(lapply(BigList,function(x) as.numeric(x)))
names(Tetramers)=names(AllNTails.sub)
Tetramers$ID=names(init)



TailTetramersummary=data.frame(ID=Tetramers$ID)

TailTetramersummary$LOKI=rowSums(Tetramers[,which(names(Tetramers)%in%LokiHFNAmes)])
TailTetramersummary$HOD=rowSums(Tetramers[,which(names(Tetramers)%in%HodaHFNAmes)])
TailTetramersummary$ASGARD=rowSums(Tetramers[,which(names(Tetramers)%in%getName(ASGARD_N_Ter_Tails))])
TailTetramersummary$EUK=rowSums(Tetramers[,which(names(Tetramers)%in%getName(All_Euk_N_Tails))])
TailTetramersummary$H2A=rowSums(Tetramers[,which(names(Tetramers)%in%getName(Euk_H2A_N_Ter_Tails))])
TailTetramersummary$H2B=rowSums(Tetramers[,which(names(Tetramers)%in%getName(Euk_H2B_N_Ter_Tails))])
TailTetramersummary$H3=rowSums(Tetramers[,which(names(Tetramers)%in%getName(Euk_H3_N_Ter_Tails))])
TailTetramersummary$H4=rowSums(Tetramers[,which(names(Tetramers)%in%getName(Euk_H4_N_Ter_Tails))])

TailTetramersummary$LOKI=100*TailTetramersummary$LOKI/sum(TailTetramersummary$LOKI)
TailTetramersummary$HOD=100*TailTetramersummary$HOD/sum(TailTetramersummary$HOD)

TailTetramersummary$ASGARD=100*TailTetramersummary$ASGARD/sum(TailTetramersummary$ASGARD)
TailTetramersummary$EUK=100*TailTetramersummary$EUK/sum(TailTetramersummary$EUK)
TailTetramersummary$H2A=100*TailTetramersummary$H2A/sum(TailTetramersummary$H2A)
TailTetramersummary$H2B=100*TailTetramersummary$H2B/sum(TailTetramersummary$H2B)
TailTetramersummary$H3=100*TailTetramersummary$H3/sum(TailTetramersummary$H3)
TailTetramersummary$H4=100*TailTetramersummary$H4/sum(TailTetramersummary$H4)



TailTetramersummary$Total=rowSums(TailTetramersummary[,2:9])
TailTetramersummary=TailTetramersummary[-which(TailTetramersummary$Total==0),]


write.table(TailTetramersummary,"~/ASGARD/ANALYSIS/N_Terminal_Tails_4merSummary.txt",quote=F,row.names=F,sep="\t")






#For Each Sequence count the number of unique Kmers, this is a proxy for complexity (non zero Dimer/tetramer)
head(Dimers)
DimerRepetiv=unlist(lapply(Dimers[,which(names(Dimers)%in%getName(AllNTails.sub))],function(x) length(which(x!=0))/sum(x)))

TetraRepetiv=unlist(lapply(Tetramers[,which(names(Tetramers)%in%getName(AllNTails.sub))],function(x) length(which(x!=0))/sum(x)))


Repetitiveness=data.frame(ID=names(DimerRepetiv),DimerRep=DimerRepetiv,TetraRep=TetraRepetiv)


Tail_properties.me=merge(Tail_properties,Repetitiveness,by="ID",all.x=T)

write.table(Tail_properties.me,"~/ASGARD/ANALYSIS/N_Terminal_Tails_Properties_w_Repetitiveness.txt",quote=F,row.names=F,sep="\t")

#Re-load 

Tail_properties.me=read.table("~/ASGARD/ANALYSIS/N_Terminal_Tails_Properties_w_Repetitiveness.txt",header=T,sep="\t",quote="",stringsAsFactors = F)


#Plot part 

DimerRepPlot=ggplot(data=Tail_properties.me)+geom_boxplot(aes(x=Family,group=Family,y=DimerRep))
ggsave(plot = DimerRepPlot,"~/ASGARD/PLOTS/KMERS/Dimermer_uniqueness_AllLength.pdf")

DimerRepPlot2=ggplot(data=Tail_properties.me[which(Tail_properties.me$Length>15),])+geom_boxplot(aes(x=Family,group=Family,y=DimerRep))
ggsave(plot = DimerRepPlot2,"~/ASGARD/PLOTS/KMERS/Dimermer_uniqueness_More15.pdf")


TetraRepPlot=ggplot(data=Tail_properties.me)+geom_boxplot(aes(x=Family,group=Family,y=TetraRep))
ggsave(plot = TetraRepPlot,"~/ASGARD/PLOTS/KMERS/Tetramer_uniqueness_AllLength.pdf")

TetraRepPlot2=ggplot(data=Tail_properties.me[which(Tail_properties.me$Length>15),])+geom_boxplot(aes(x=Family,group=Family,y=TetraRep))
ggsave(plot = TetraRepPlot2,"~/ASGARD/PLOTS/KMERS/Tetramer_uniqueness_More15.pdf")







#Just at the amino acid content what can we seee
library(reshape2)

Tail_properties.m=reshape2::melt(Tail_properties,id.vars =c("ID", "Family","Segment" ))

Tail_properties.m$value=as.numeric(Tail_properties.m$value)
LengthPI=ggplot(data=Tail_properties)+geom_point(aes(x=Length,y=Pi))+facet_wrap(~Family,scales = "free")
ggsave(plot  = LengthPI,"~/ASGARD/PLOTS/KMERS/All_Tails_length_vs_PI.pdf")


AllPropPlot=ggplot(data=Tail_properties.m)+geom_boxplot(aes(x=Family,y=as.numeric(value)),outlier.size = 0.5)+facet_wrap(~variable,scales = "free")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(plot  = AllPropPlot,"~/ASGARD/PLOTS/KMERS/All_Properties_Pannels.pdf",height=10,width = 10)


Tail_properties.m=reshape2::melt(Tail_properties[which(Tail_properties$Length>15),],id.vars =c("ID", "Family","Segment" ))
AllPropPlot=ggplot(data=Tail_properties.m)+geom_boxplot(aes(x=Family,y=as.numeric(value)),outlier.size = 0.5)+facet_wrap(~variable,scales = "free")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(plot  = AllPropPlot,"~/ASGARD/PLOTS/KMERS/All_Properties_Pannels_lengthsup15aa.pdf",height=10,width = 10)




#Individual properties ( only length include tail of length 0): 
library(ggthemes);library(ggpubr)
Length=ggplot(data=Tail_properties.me)+geom_boxplot(aes(x=Family,y=Length),outlier.size = 0.5,fill="grey90")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme_pubclean()+xlab("")+ylab("Tail length (a.a)")
ggsave(plot = Length,filename = "~/ASGARD/PLOTS/KMERS/BOXPLOTS/Tail_lengths.pdf",width=3,height=3)



Glycine=ggplot(data=Tail_properties.me)+geom_boxplot(aes(x=Family,y=100*Prop.G),outlier.size = 0.5,fill="grey90")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme_pubclean()+xlab("")+ylab("Glycine (%)")+ylim(c(0,100))
ggsave(plot = Glycine,filename = "~/ASGARD/PLOTS/KMERS/BOXPLOTS/Tail_Glycine.pdf",width=3,height=3)



Size=ggplot(data=Tail_properties.me)+geom_boxplot(aes(x=Family,y=100*(Prop.Small)),outlier.size = 0.5,fill="grey90")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme_pubclean()+xlab("")+ylab("Tiny / Small amino acids (%)")+ylim(c(0,100))
ggsave(plot = Size,filename = "~/ASGARD/PLOTS/KMERS/BOXPLOTS/Tail_Size.pdf",width=3,height=3)


PropKR=ggplot(data=Tail_properties.me)+geom_boxplot(aes(x=Family,y=100*(Prop.K+Prop.R)),outlier.size = 0.5,fill="grey90")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme_pubclean()+xlab("")+ylab("Lysine + Arginine (%) (a.a)")+ylim(c(0,100))
ggsave(plot = PropKR,filename = "~/ASGARD/PLOTS/KMERS/BOXPLOTS/Tail_PropKR.pdf",width=3,height=3)



PropH=ggplot(data=Tail_properties.me)+geom_boxplot(aes(x=Family,y=100*Prop.H),outlier.size = 0.5,fill="grey90")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme_pubclean()+xlab("")+ylab("Histidine (%) (a.a)")+ylim(c(0,100))
ggsave(plot = PropH,filename = "~/ASGARD/PLOTS/KMERS/BOXPLOTS/Tail_PropH.pdf",width=3,height=3)


Neg=ggplot(data=Tail_properties.me)+geom_boxplot(aes(x=Family,y=100*(Prop.D+Prop.E)),outlier.size = 0.5,fill="grey90")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme_pubclean()+xlab("")+ylab("Negatively Charged (%)")+ylim(c(0,100))
ggsave(plot = Neg,filename = "~/ASGARD/PLOTS/KMERS/BOXPLOTS/Tail_Neg.pdf",width=3,height=3)





#Dimer abundancy
MatPlot=as.matrix(TailDimerSummary[order(-TailDimerSummary$ASGARD)[1:20],2:9])
row.names(MatPlot)=TailDimerSummary[order(-TailDimerSummary$ASGARD)[1:20],]$ID
MatplotLong=reshape2::melt(MatPlot)


DimerHeatmapPLot=ggplot(data=MatplotLong,aes(x=Var2,y=Var1,fill=value))+geom_tile()+coord_equal()+scale_fill_steps(low = "white")+theme_pubr()+ylab("% of 2-mer")+theme_tufte()+xlab("")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))
ggsave(plot  = DimerHeatmapPLot,"~/ASGARD/PLOTS/KMERS/All_DimerHeatmapOrdered_ASGARD.pdf")

#Ordered according to H3
MatPlot=as.matrix(TailDimerSummary[order(-TailDimerSummary$H3)[1:20],2:9])
row.names(MatPlot)=TailDimerSummary[order(-TailDimerSummary$H3)[1:20],]$ID
MatplotLong=reshape2::melt(MatPlot)


DimerHeatmapPLotH3=ggplot(data=MatplotLong,aes(x=Var2,y=Var1,fill=value))+geom_tile()+coord_equal()+scale_fill_steps(low = "white")+theme_pubr()+ylab("% of 2-mer")+theme_tufte()+xlab("")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))
ggsave(plot  = DimerHeatmapPLotH3,"~/ASGARD/PLOTS/KMERS/All_DimerHeatmapOrdered_H3.pdf")



# ========================================================================
#KMER PLOTS PART
# ========================================================================
library(ggthemes)
library(reshape2)
library(ggpubr)

N=10
MatPlot=as.matrix(TailTetramerSummary[order(-TailTetramerSummary$ASGARD)[c(1:N)],2:9])
row.names(MatPlot)=TailTetramerSummary[order(-TailTetramerSummary$ASGARD)[c(1:N)],]$ID
MatplotLong=reshape2::melt(MatPlot)
plota=ggplot(MatplotLong[which(MatplotLong$Var2=="ASGARD"),])+geom_bar(aes(x=Var1,y=value),stat = "identity", width=0.7,fill="grey55")+ylim(c(0,3))+scale_y_reverse()+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of 5-mer")+coord_flip()+theme_tufte()+xlab("")+ggtitle("ASGARD")

plotb=ggplot(MatplotLong[which(MatplotLong$Var2=="EUK"),])+geom_bar(aes(x=Var1,y=value),stat = "identity", width=0.7,fill="grey55")+ylim(c(0,3))+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of 5-mer")+coord_flip()+theme_tufte()+xlab("")+ggtitle("EUK")
ggsave(plot = ggarrange(plota,plotb),filename = "~/ASGARD/PLOTS/KMERS/Top20_N_Tail_ASGARD_4mersvs_Euk.pdf")


#Hodarchaea vs Euk 
N=10
MatPlot=as.matrix(TailTetramerSummary[order(-TailTetramerSummary$HOD)[c(1:N)],2:9])
row.names(MatPlot)=TailTetramerSummary[order(-TailTetramerSummary$HOD)[c(1:N)],]$ID
MatplotLong=reshape2::melt(MatPlot)
plota=ggplot(MatplotLong[which(MatplotLong$Var2=="HOD"),])+geom_bar(aes(x=Var1,y=value),stat = "identity", width=0.7,fill="grey55")+ylim(c(0,3))+scale_y_reverse()+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of 5-mer")+coord_flip()+theme_tufte()+xlab("")+ggtitle("HOD")

plotb=ggplot(MatplotLong[which(MatplotLong$Var2=="EUK"),])+geom_bar(aes(x=Var1,y=value),stat = "identity", width=0.7,fill="grey55")+ylim(c(0,3))+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of 5-mer")+coord_flip()+theme_tufte()+xlab("")+ggtitle("EUK")
ggsave(plot = ggarrange(plota,plotb),filename = "~/ASGARD/PLOTS/KMERS/Top10_N_Tail_HOD_4mersvs_Euk.pdf")


#Hod vs h2a,h2b,h3,h4
plota=ggplot(MatplotLong[which(MatplotLong$Var2=="HOD"),])+geom_bar(aes(x=Var1,y=value),stat = "identity", width=0.7,fill="grey55")+scale_y_reverse()+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of 5-mer")+coord_flip()+theme_tufte()+xlab("")+ggtitle("HOD")+ylim(c(0,2))
plotb=ggplot(MatplotLong[which(MatplotLong$Var2=="H2A"),])+geom_bar(aes(x=Var1,y=value),stat = "identity", width=0.7,fill="grey55")+scale_y_reverse()+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of 5-mer")+coord_flip()+theme_tufte()+xlab("")+ggtitle("H2A")+ylim(c(0,2))
plotc=ggplot(MatplotLong[which(MatplotLong$Var2=="H2B"),])+geom_bar(aes(x=Var1,y=value),stat = "identity", width=0.7,fill="grey55")+scale_y_reverse()+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of 5-mer")+coord_flip()+theme_tufte()+xlab("")+ggtitle("H2B")+ylim(c(0,2))
plotd=ggplot(MatplotLong[which(MatplotLong$Var2=="H3"),])+geom_bar(aes(x=Var1,y=value),stat = "identity", width=0.7,fill="grey55")+scale_y_reverse()+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of 5-mer")+coord_flip()+theme_tufte()+xlab("")+ggtitle("H3")+ylim(c(0,2))
plote=ggplot(MatplotLong[which(MatplotLong$Var2=="H4"),])+geom_bar(aes(x=Var1,y=value),stat = "identity", width=0.7,fill="grey55")+scale_y_reverse()+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of 5-mer")+coord_flip()+theme_tufte()+xlab("")+ggtitle("H4")+ylim(c(0,2))

ggsave(plot = ggarrange(plota,plotb,plotc,plotd,plote,ncol=5),filename = "~/ASGARD/PLOTS/KMERS/Top10_N_Tail_HOD_4mersvs_HistoneTypes.pdf",height=3,width=8)







#Lokiarchaea vs Euk 
N=20
MatPlot=as.matrix(TailTetramersummary[order(-TailTetramersummary$LOKI)[c(1:N)],2:9])
row.names(MatPlot)=TailTetramersummary[order(-TailTetramersummary$LOKI)[c(1:N)],]$ID
MatplotLong=reshape2::melt(MatPlot)
plota=ggplot(MatplotLong[which(MatplotLong$Var2=="LOKI"),])+geom_bar(aes(x=Var1,y=value),stat = "identity", width=0.7,fill="grey55")+ylim(c(0,7.5))+scale_y_reverse()+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of 5-mer")+coord_flip()+theme_tufte()+xlab("")+ggtitle("LOKI")

plotb=ggplot(MatplotLong[which(MatplotLong$Var2=="EUK"),])+geom_bar(aes(x=Var1,y=value),stat = "identity", width=0.7,fill="grey55")+ylim(c(0,7.5))+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of 5-mer")+coord_flip()+theme_tufte()+xlab("")+ggtitle("EUK")
ggsave(plot = ggarrange(plota,plotb),filename = "~/ASGARD/PLOTS/KMERS/Top20_N_Tail_LOKI_4mersvs_Euk.pdf")


#LOKI vs h2a,h2b,h3,h4
plota=ggplot(MatplotLong[which(MatplotLong$Var2=="LOKI"),])+geom_bar(aes(x=Var1,y=value),stat = "identity", width=0.7,fill="grey55")+scale_y_reverse()+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of 5-mer")+coord_flip()+theme_tufte()+xlab("")+ggtitle("LOKI")+ylim(c(0,7.5))
plotb=ggplot(MatplotLong[which(MatplotLong$Var2=="H2A"),])+geom_bar(aes(x=Var1,y=value),stat = "identity", width=0.7,fill="grey55")+scale_y_reverse()+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of 5-mer")+coord_flip()+theme_tufte()+xlab("")+ggtitle("H2A")+ylim(c(0,7.5))
plotc=ggplot(MatplotLong[which(MatplotLong$Var2=="H2B"),])+geom_bar(aes(x=Var1,y=value),stat = "identity", width=0.7,fill="grey55")+scale_y_reverse()+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of 5-mer")+coord_flip()+theme_tufte()+xlab("")+ggtitle("H2B")+ylim(c(0,7.5))
plotd=ggplot(MatplotLong[which(MatplotLong$Var2=="H3"),])+geom_bar(aes(x=Var1,y=value),stat = "identity", width=0.7,fill="grey55")+scale_y_reverse()+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of 5-mer")+coord_flip()+theme_tufte()+xlab("")+ggtitle("H3")+ylim(c(0,7.5))
plote=ggplot(MatplotLong[which(MatplotLong$Var2=="H4"),])+geom_bar(aes(x=Var1,y=value),stat = "identity", width=0.7,fill="grey55")+scale_y_reverse()+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of 5-mer")+coord_flip()+theme_tufte()+xlab("")+ggtitle("H4")+ylim(c(0,7.5))

ggsave(plot = ggarrange(plota,plotb,plotc,plotd,plote,ncol=5),filename = "~/ASGARD/PLOTS/KMERS/Top20_N_Tail_LOKI_4mersvs_HistoneTypes.pdf")


#Clipped axis


#LOKI vs h2a,h2b,h3,h4
plota=ggplot(MatplotLong[which(MatplotLong$Var2=="LOKI"),])+geom_bar(aes(x=Var1,y=value),stat = "identity", width=0.7,fill="grey55")+scale_y_reverse()+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of 5-mer")+coord_flip()+theme_tufte()+xlab("")+ggtitle("LOKI")+ylim(c(0,3))
plotb=ggplot(MatplotLong[which(MatplotLong$Var2=="H2A"),])+geom_bar(aes(x=Var1,y=value),stat = "identity", width=0.7,fill="grey55")+scale_y_reverse()+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of 5-mer")+coord_flip()+theme_tufte()+xlab("")+ggtitle("H2A")+ylim(c(0,3))
plotc=ggplot(MatplotLong[which(MatplotLong$Var2=="H2B"),])+geom_bar(aes(x=Var1,y=value),stat = "identity", width=0.7,fill="grey55")+scale_y_reverse()+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of 5-mer")+coord_flip()+theme_tufte()+xlab("")+ggtitle("H2B")+ylim(c(0,3))
plotd=ggplot(MatplotLong[which(MatplotLong$Var2=="H3"),])+geom_bar(aes(x=Var1,y=value),stat = "identity", width=0.7,fill="grey55")+scale_y_reverse()+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of 5-mer")+coord_flip()+theme_tufte()+xlab("")+ggtitle("H3")+ylim(c(0,3))
plote=ggplot(MatplotLong[which(MatplotLong$Var2=="H4"),])+geom_bar(aes(x=Var1,y=value),stat = "identity", width=0.7,fill="grey55")+scale_y_reverse()+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of 5-mer")+coord_flip()+theme_tufte()+xlab("")+ggtitle("H4")+ylim(c(0,3))

ggsave(plot = ggarrange(plota,plotb,plotc,plotd,plote,ncol=5),filename = "~/ASGARD/PLOTS/KMERS/Top20_N_Tail_LOKI_4mersvs_HistoneTypesClippedaxis3.pdf")










#Now oredered on eukaryotic tails:
N=10
MatPlot=as.matrix(TailTetramerSummary[order(-TailTetramerSummary$EUK)[c(1:N)],2:9])
row.names(MatPlot)=TailTetramerSummary[order(-TailTetramerSummary$EUK)[c(1:N)],]$ID
MatplotLong=reshape2::melt(MatPlot)
plota=ggplot(MatplotLong[which(MatplotLong$Var2=="ASGARD"),])+geom_bar(aes(x=Var1,y=value),stat = "identity", width=0.7,fill="grey55")+ylim(c(0,3))+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of 5-mer")+coord_flip()+theme_tufte()+xlab("")+ggtitle("ASGARD")+scale_y_reverse()

plotb=ggplot(MatplotLong[which(MatplotLong$Var2=="EUK"),])+geom_bar(aes(x=Var1,y=value),stat = "identity", width=0.7,fill="grey55")+ylim(c(0,3))+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of 5-mer")+coord_flip()+theme_tufte()+xlab("")
ggsave(plot = ggarrange(plota,plotb),filename = "~/ASGARD/PLOTS/KMERS/Top20_N_Tail_Euk_4mersvs_ASGARD.pdf")




#Now Top5 ASGARD and Top5Euk
N=5
MatPlot=as.matrix(TailTetramerSummary[c(order(-TailTetramerSummary$EUK)[c(1:N)],rev(order(-TailTetramerSummary$ASGARD)[c(1:N)])),2:9])
row.names(MatPlot)=TailTetramerSummary[c(order(-TailTetramerSummary$EUK)[c(1:N)],rev(order(-TailTetramerSummary$ASGARD)[c(1:N)])),]$ID
MatplotLong=reshape2::melt(MatPlot)
plota=ggplot(MatplotLong[which(MatplotLong$Var2=="ASGARD"),])+geom_bar(aes(x=Var1,y=value),stat = "identity", width=0.7,fill="grey55")+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of 5-mer")+theme_tufte()+xlab("")+ggtitle("ASGARD")+coord_flip()+ylim(c(0,4))

plotb=ggplot(MatplotLong[which(MatplotLong$Var2=="EUK"),])+geom_bar(aes(x=Var1,y=value),stat = "identity", width=0.7,fill="grey55")+ylim(c(0,4))+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of 5-mer")+coord_flip()+theme_tufte()+xlab("")+ggtitle("EUK")

ggsave(plot = ggarrange(plota,plotb),filename = "~/ASGARD/PLOTS/KMERS/Top5_ASGARD_And_EUK_N_Tails.pdf",height=3,width=3)








#Ordred on H2A
N=20
MatPlot=as.matrix(TailTetramersummary[order(-TailTetramersummary$H2A)[c(1:N)],2:9])
row.names(MatPlot)=TailTetramersummary[order(-TailTetramersummary$H2A)[c(1:N)],]$ID
MatplotLong=reshape2::melt(MatPlot)
plota=ggplot(MatplotLong[which(MatplotLong$Var2=="ASGARD"),])+geom_bar(aes(x=Var1,y=value),stat = "identity", width=0.7,fill="grey55")+ylim(c(0,3))+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of 5-mer")+coord_flip()+theme_tufte()+xlab("")+ggtitle("ASGARD")

plotb=ggplot(MatplotLong[which(MatplotLong$Var2=="H2A"),])+geom_bar(aes(x=Var1,y=value),stat = "identity", width=0.7,fill="grey55")+ylim(c(0,3))+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of 5-mer")+coord_flip()+theme_tufte()+xlab("")+ggtitle("H2A")
ggsave(plot = ggarrange(plota,plotb),filename = "~/ASGARD/PLOTS/KMERS/Top20_N_Tail_H2A_4mersvs_ASGARD.pdf")


#Ordred on H2B
N=20
MatPlot=as.matrix(TailTetramersummary[order(-TailTetramersummary$H2B)[c(1:N)],2:9])
row.names(MatPlot)=TailTetramersummary[order(-TailTetramersummary$H2B)[c(1:N)],]$ID
MatplotLong=reshape2::melt(MatPlot)
plota=ggplot(MatplotLong[which(MatplotLong$Var2=="ASGARD"),])+geom_bar(aes(x=Var1,y=value),stat = "identity", width=0.7,fill="grey55")+ylim(c(0,3))+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of 5-mer")+coord_flip()+theme_tufte()+xlab("")+ggtitle("ASGARD")

plotb=ggplot(MatplotLong[which(MatplotLong$Var2=="H2B"),])+geom_bar(aes(x=Var1,y=value),stat = "identity", width=0.7,fill="grey55")+ylim(c(0,3))+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of 5-mer")+coord_flip()+theme_tufte()+xlab("")+ggtitle("H2B")
ggsave(plot = ggarrange(plota,plotb),filename = "~/ASGARD/PLOTS/KMERS/Top20_N_Tail_H2B_4mersvs_ASGARD.pdf")

#Ordred on H3
N=20
MatPlot=as.matrix(TailTetramersummary[order(-TailTetramersummary$H3)[c(1:N)],2:9])
row.names(MatPlot)=TailTetramersummary[order(-TailTetramersummary$H3)[c(1:N)],]$ID
MatplotLong=reshape2::melt(MatPlot)
plota=ggplot(MatplotLong[which(MatplotLong$Var2=="ASGARD"),])+geom_bar(aes(x=Var1,y=value),stat = "identity", width=0.7,fill="grey55")+ylim(c(0,3))+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of 5-mer")+coord_flip()+theme_tufte()+xlab("")+ggtitle("ASGARD")

plotb=ggplot(MatplotLong[which(MatplotLong$Var2=="H3"),])+geom_bar(aes(x=Var1,y=value),stat = "identity", width=0.7,fill="grey55")+ylim(c(0,3))+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of 5-mer")+coord_flip()+theme_tufte()+xlab("")+ggtitle("H3")
ggsave(plot = ggarrange(plota,plotb),filename = "~/ASGARD/PLOTS/KMERS/Top20_N_Tail_H3_4mersvs_ASGARD.pdf")

#Ordred on H4
N=20
MatPlot=as.matrix(TailTetramersummary[order(-TailTetramersummary$H4)[c(1:N)],2:9])
row.names(MatPlot)=TailTetramersummary[order(-TailTetramersummary$H4)[c(1:N)],]$ID
MatplotLong=reshape2::melt(MatPlot)

plota=ggplot(MatplotLong[which(MatplotLong$Var2=="ASGARD"),])+geom_bar(aes(x=Var1,y=value),stat = "identity", width=0.7,fill="grey55")+ylim(c(0,5))+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of 5-mer")+coord_flip()+theme_tufte()+xlab("")+ggtitle("ASGARD")

plotb=ggplot(MatplotLong[which(MatplotLong$Var2=="H4"),])+geom_bar(aes(x=Var1,y=value),stat = "identity", width=0.7,fill="grey55")+ylim(c(0,5))+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of 5-mer")+coord_flip()+theme_tufte()+xlab("")+ggtitle("H4")
ggsave(plot = ggarrange(plota,plotb),filename = "~/ASGARD/PLOTS/KMERS/Top20_N_Tail_H4_4mersvs_ASGARD.pdf")




####################################
#PCA PLOTS on K-mer CONTENT and biochemical ppties:


#Before computing those, we will remove 1 to the legnth of tail in tails properties to account for the fact taht we removed a methionine

Tail_properties[which(Tail_properties$Length>0),]$Length=Tail_properties[which(Tail_properties$Length>0),]$Length-1

#Creating colors for later : 
library(RColorBrewer)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

#Annotate Sub-families
Tail_properties$SubFamily=Tail_properties$Family
Tail_properties[which(Tail_properties$ID%in%HodaHFNAmes),]$SubFamily="Hodarchaea"
Tail_properties[which(Tail_properties$ID%in%LokiHFNAmes),]$SubFamily="Lokiarchaea"




####################################
##PCA on biochemical properties


#PCA according to http://www.sthda.com/english/wiki/fviz-pca-quick-principal-component-analysis-data-visualization-r-software-and-data-mining#fviz_pca_var-graph-of-variables
library("factoextra")


ToPCA=Tail_properties[which(Tail_properties$Length>5),]




ToPCAMat=as.matrix(ToPCA[,which(names(ToPCA)%in%c("Prop.Tiny","Prop.Small","Prop.Aliphatic","Prop.Non.polar","Prop.Polar","Prop.Charged","Prop.Basic","Prop.Acidic","NetChargePerRes"))])
row.names(ToPCAMat)=ToPCA$ID



PCA_PPTIES=prcomp(ToPCAMat, center = T, scale. = F)



p <- fviz_pca_ind(PCA_PPTIES, label="none", habillage=ToPCA$Family,
                  addEllipses=TRUE, ellipse.level=0.95)
p=p+scale_color_brewer(palette="Dark2") +
  theme_minimal()+theme(aspect.ratio = 1)


#Remove the Prop. for visualisation
row.names(PCA_PPTIES$rotation)=gsub("Prop.","",row.names(PCA_PPTIES$rotation))
#Contribution of each variable
g=fviz_pca_var(PCA_PPTIES, col.var="contrib",labelsize = 3)+scale_color_gradient2(low="white", mid="blue",high="red", midpoint=5) +theme_minimal()+theme(aspect.ratio = 1)



ggsave(plot = ggpubr::ggarrange(p,g),filename = "~/ASGARD/PLOTS/KMERS/PCA_individual_HIstones_AA_Biochemistry.pdf",width=8,height=8)


#Separating Loki and Hodarchaea
p <- fviz_pca_ind(PCA_PPTIES, label="none", habillage=ToPCA$SubFamily,
                  addEllipses=TRUE, ellipse.level=0.95)
p
p=p+scale_color_brewer(palette="Dark2") +
  theme_minimal()+theme(aspect.ratio = 1)


fviz_ellipses(PCA_PPTIES, habillage = as.factor(ToPCA$Family),
              graph.type = "ggplot", ggtheme = theme_pubr(),
             pointsize=0.7, geom="point" ,
              ellipse.level = 0.9, ellipse.type = c("norm"),legendposition="right")


#Contribution of each variable
g=fviz_pca_var(PCA_PPTIES, col.var="contrib",labelsize = 3)+scale_color_gradient2(low="white", mid="blue",high="red", midpoint=5) +theme_minimal()+theme(aspect.ratio = 1)



ggsave(plot = ggpubr::ggarrange(p,g),filename = "~/ASGARD/PLOTS/KMERS/PCA_individual_HIstones_AA_Biochemistry_w_Hod_Loki.pdf",width=8,height=8)

#On individual AA content

ToPCAMat=as.matrix(ToPCA[,c(6:25)])
row.names(ToPCAMat)=ToPCA$ID


PCA_PPTIES=prcomp(ToPCAMat)
summary(PCA_PPTIES)
PCA_Scores=as.data.frame(PCA_PPTIES$x)

PCA_Scores$Family=ToPCA$Family
PCA_Scores$Segment=ToPCA$Segment

p=fviz_pca_ind(PCA_PPTIES, label="none", habillage=ToPCA$Family,
             addEllipses=TRUE, ellipse.level=0.90,col.var = "black", #setas
             geom = "point",pointsize=0.5,mean.point = FALSE)+
  scale_shape_manual(values=c(rep(19,7)))+scale_color_manual(values=col_vector[c(1,2,11,3,24)])+scale_fill_manual(values=col_vector[c(1,2,11,3,24)])+theme_pubr()+theme(aspect.ratio = 1)




row.names(PCA_PPTIES$rotation)=gsub("Prop.","",row.names(PCA_PPTIES$rotation))

#Contribution of each variable
g=fviz_pca_var(PCA_PPTIES, col.var="contrib",labelsize = 3)+scale_color_gradient2(low="white", mid=col_vector[5],high=col_vector[22], midpoint=5) +theme_minimal()+theme(aspect.ratio = 1)

ggsave(plot = ggarrange(p,g),filename = "~/ASGARD/PLOTS/KMERS/PCA_individual_HIstones_AA_Compo.pdf",width=8,height=8)









#On Dimers, individual sequences
SequenceToKeep=Tail_properties[which(Tail_properties$Length>5),]$ID

ToPCA=Dimers[which(names(Dimers)%in%c(SequenceToKeep))]

#Normalise count into frequency
ToPCA=as.data.frame(t(apply(ToPCA,2, function(x) x/sum(x))))
names(ToPCA)=Dimers$ID
PCA_PPTIES=prcomp(ToPCA)

ToPCA$Family=unlist(lapply(row.names(ToPCA),function(x) Tail_properties[which(Tail_properties$ID==x),]$Family
))


summary(PCA_PPTIES)

p <- fviz_pca_ind(PCA_PPTIES, label="none", habillage=ToPCA$Family,
                  addEllipses=TRUE, ellipse.level=0.95)
p=p+scale_color_brewer(palette="Dark2") +
  theme_minimal()+theme(aspect.ratio = 1)

#Contribution of each variable
g=fviz_pca_var(PCA_PPTIES, col.var="contrib",labelsize = 3)+scale_color_gradient2(low="white", mid="blue",high="red", midpoint=5) +theme_minimal()+theme(aspect.ratio = 1)


ggsave(plot = ggarrange(p,g),filename = "~/ASGARD/PLOTS/KMERS/PCA_individual_HIstones_AA_Dimers.pdf",width=8,height=8)




#Using sparce PCA

SequenceToKeep=Tail_properties[which(Tail_properties$Length>5),]$ID

ToPCA=Dimers[which(names(Dimers)%in%c(SequenceToKeep))]

#Normalise count into frequency
ToPCA=as.data.frame(t(apply(ToPCA,2, function(x) x/sum(x))))
names(ToPCA)=Dimers$ID
PCA_PPTIES=sparsepca::spca(ToPCA,k=10,max_iter=300)
SCORES=as.data.frame(PCA_PPTIES$scores)
SCORES$Family=Tail_properties[which(Tail_properties$Length>5),]$Family
SPCADimers=ggplot(SCORES,aes(x=V1,y=V2,color=Family))+geom_point()

ggsave(plot = SPCADimers,filename = "~/ASGARD/PLOTS/KMERS/SPCA_individual_HIstones_AA_Dimers.pdf",width=8,height=8)

#Subfamillies of asgards
SCORES$Family=Tail_properties[which(Tail_properties$Length>5),]$SubFamily
SPCADimers=ggplot(SCORES,aes(x=V1,y=V2,color=Family))+geom_point()

ggsave(plot = SPCADimers,filename = "~/ASGARD/PLOTS/KMERS/SPCA_individual_HIstones_AA_Dimers_Subfamilies.pdf",width=8,height=8)



#on DIMERS summary
TailDimerSummary=read.table("~/ASGARD/ANALYSIS/N_Terminal_Tails_2merSummary.txt",header=T,sep="\t",quote="")

ToPCA=as.data.frame(t.data.frame(TailDimerSummary[,c(4,6:9)]))
names(ToPCA)=TailDimerSummary$ID
#There is a downstream bug where the DIMER NA is counted a non available 
#So we modifiy it ... 
names(ToPCA)[which(is.na(names(ToPCA)))]="NA_d"
PCADimer=prcomp(ToPCA)
ToPCA$Family=row.names(ToPCA)


p <- fviz_pca_ind(PCADimer,
                  addEllipses=TRUE, ellipse.level=0.95)
p=p+scale_color_brewer(palette="Dark2") +
  theme_minimal()+theme(aspect.ratio = 1)

#Contribution of each variable
g=fviz_pca_var(PCADimer, col.var="contrib",labelsize = 3)+scale_color_gradient2(low="white", mid="blue",high="red", midpoint=5) +theme_minimal()+theme(aspect.ratio = 1)



ggsave(plot = ggarrange(p,g),filename = "~/ASGARD/PLOTS/KMERS/PCA_individual_HIstones_AA_DimersSUMMARY.pdf",width=8,height=8)




#Trimers> Full and then Summary

if(length(which(rowSums(Trimers[,-which(names(Trimers)=="ID")])==0))>0){
Trimers=Trimers[-which(rowSums(Trimers[,-which(names(Trimers)=="ID")])==0),]}

SequenceToKeep=Tail_properties[which(Tail_properties$Length>5),]$ID

ToPCA=Trimers[,which(names(Trimers)%in%c(SequenceToKeep))]


#Normalise count into frequency
ToPCA=as.data.frame(t(apply(ToPCA,2, function(x) x/sum(x))))
names(ToPCA)=Trimers$ID
PCA_PPTIES=prcomp(ToPCA)

ToPCA$Family=unlist(lapply(row.names(ToPCA),function(x) Tail_properties[which(Tail_properties$ID==x),]$Family))




p <- fviz_pca_ind(PCA_PPTIES, label="none", habillage=ToPCA$Family,
                  addEllipses=TRUE, ellipse.level=0.95)
p=p+scale_color_brewer(palette="Dark2") +
  theme_minimal()+theme(aspect.ratio = 1)

#Contribution of each variable
g=fviz_pca_var(PCA_PPTIES, col.var="contrib",labelsize = 3)+scale_color_gradient2(low="white", mid="blue",high="red", midpoint=5) +theme_minimal()+theme(aspect.ratio = 1)




ggsave(plot = ggarrange(p,g),filename = "~/ASGARD/PLOTS/KMERS/PCA_individual_HIstones_AA_Trimers.pdf",width=8,height=8)


#using sparce PCA
#Normalise count into frequency

PCA_PPTIES=sparsepca::spca(ToPCA,k=10,max_iter=100)
SCORES=as.data.frame(PCA_PPTIES$scores)
SCORES$Family=Tail_properties[which(Tail_properties$Length>5),]$Family



SPCA_TRIMERS=ggplot(SCORES,aes(x=V1,y=V2,color=Family))+geom_point(size=0.7)+scale_color_manual(values = col_vector)


ggsave(plot = SPCA_TRIMERS,filename = "~/ASGARD/PLOTS/KMERS/SPCA_individual_HIstones_AA_Trimers.pdf",width=8,height=8)








#############
#On summary
#############

TailTrimerSummary=read.table("~/ASGARD/ANALYSIS/N_Terminal_Tails_3merSummary.txt",header=T,sep="\t",quote="")

ToPCA=as.data.frame(t.data.frame(TailTrimerSummary[,c(4,6:9)]))
names(ToPCA)=TailTrimerSummary$ID
#There is a downstream bug where the DIMER NA is counted a non available 
#So we modifiy it ... 
PCATrimer=prcomp(ToPCA)
ToPCA$Family=row.names(ToPCA)


p <- fviz_pca_ind(PCATrimer,
                  addEllipses=TRUE, ellipse.level=0.95)
p=p+scale_color_brewer(palette="Dark2") +
  theme_minimal()+theme(aspect.ratio = 1)

#Contribution of each variable
g=fviz_pca_var(PCATrimer, col.var="contrib",labelsize = 3)+scale_color_gradient2(low="white", mid="blue",high="red", midpoint=5) +theme_minimal()+theme(aspect.ratio = 1)


ggsave(plot = ggarrange(p,g),filename = "~/ASGARD/PLOTS/KMERS/PCA_individual_HIstones_AA_Trimers_SUMMARY.pdf",width=8,height=8)




#¢¢ TETRA mer 
TailTetramerSummary=read.table("~/ASGARD/ANALYSIS/N_Terminal_Tails_4merSummary.txt",header=T,sep="\t",quote="")

ToPCA=as.data.frame(t.data.frame(TailTetramerSummary[,c(4,6:9)]))
names(ToPCA)=TailTetramerSummary$ID
#There is a downstream bug where the DIMER NA is counted a non available 
#So we modifiy it ... 
PCATetramer=prcomp(ToPCA)
ToPCA$Family=row.names(ToPCA)


p <- fviz_pca_ind(PCATetramer,
                  addEllipses=TRUE, ellipse.level=0.95)
p=p+scale_color_brewer(palette="Dark2") +
  theme_minimal()+theme(aspect.ratio = 1)

#Contribution of each variable
g=fviz_pca_var(PCATetramer, col.var="contrib",labelsize = 3)+scale_color_gradient2(low="white", mid="blue",high="red", midpoint=5) +theme_minimal()+theme(aspect.ratio = 1)



ggsave(plot = ggarrange(p,g),filename = "~/ASGARD/PLOTS/KMERS/PCA_individual_HIstones_AA_TetraSUMMARY.pdf",width=8,height=8)




#Using sparce PCA on tetramers

SequenceToKeep=Tail_properties[which(Tail_properties$Length>3),]$ID

if(length(which(rowSums(Tetramers[,-which(names(Tetramers)=="ID")])==0))>0){
  Tetramers=Tetramers[-which(rowSums(Tetramers[,-which(names(Tetramers)=="ID")])==0),]}

ToPCA=Tetramers[which(names(Tetramers)%in%c(SequenceToKeep))]

#Normalise count into frequency
ToPCA=as.data.frame(t(apply(ToPCA,2, function(x) x/sum(x))))
names(ToPCA)=Tetramers$ID
PCA_PPTIES=sparsepca::spca(ToPCA,k=10,max_iter=100)
SCORES=as.data.frame(PCA_PPTIES$scores)
SCORES$Family=Tail_properties[which(Tail_properties$Length>5),]$Family
SPCATetramers=ggplot(SCORES,aes(x=V1,y=V2,color=Family))+geom_point()

ggsave(plot = SPCATetramers,filename = "~/ASGARD/PLOTS/KMERS/SPCA_individual_HIstones_AA_Tetramers.pdf",width=8,height=8)


######Classic PCA
ToPCA=Tetramers[,which(names(Tetramers)%in%c(SequenceToKeep))]


#Normalise count into frequency
ToPCA=as.data.frame(t(apply(ToPCA,2, function(x) x/sum(x))))
names(ToPCA)=Tetramers$ID
PCA_PPTIES=prcomp(ToPCA)

ToPCA$Family=unlist(lapply(row.names(ToPCA),function(x) Tail_properties[which(Tail_properties$ID==x),]$Family))




p <- fviz_pca_ind(PCA_PPTIES, label="none", habillage=ToPCA$Family,
                  addEllipses=TRUE, ellipse.level=0.95)
p=p+scale_color_brewer(palette="Dark2") +
  theme_minimal()+theme(aspect.ratio = 1)

#Contribution of each variable
g=fviz_pca_var(PCA_PPTIES, col.var="contrib",labelsize = 3)+scale_color_gradient2(low="white", mid="blue",high="red", midpoint=5) +theme_minimal()+theme(aspect.ratio = 1)




ggsave(plot = ggarrange(p,g),filename = "~/ASGARD/PLOTS/KMERS/PCA_individual_HIstones_AA_Tetramers.pdf",width=8,height=8)








###########################################################
##Reverse the sequences of tails for visualisation purposes

if(length(which(getLength(AllNTails)>0))>0){
AllNTails=AllNTails[which(getLength(AllNTails)>0)]}
AllNTails.Rev=lapply(getSequence(AllNTails),function(x) rev(x))
write.fasta(AllNTails.Rev,getName(AllNTails),"~/ASGARD/HMMSEARCH/All_N_Tails_RERVERSED.fa")


#export reversed tails per family
write.fasta(AllNTails.Rev[which(getName(AllNTails)%in%getName(ASGARD_N_Ter_Tails))], getName(AllNTails)[which(getName(AllNTails)%in%getName(ASGARD_N_Ter_Tails))],file.out = "~/ASGARD/HMMSEARCH/All_N_Tails_RERVERSED_ASGARD.fa")

#export reversed tails per family
write.fasta(AllNTails.Rev[which(getName(AllNTails)%in%getName(Euk_H2A_N_Ter_Tails))], getName(AllNTails)[which(getName(AllNTails)%in%getName(Euk_H2A_N_Ter_Tails))],file.out = "~/ASGARD/HMMSEARCH/All_N_Tails_RERVERSED_H2A.fa")
#export reversed tails per family
write.fasta(AllNTails.Rev[which(getName(AllNTails)%in%getName(Euk_H2B_N_Ter_Tails))], getName(AllNTails)[which(getName(AllNTails)%in%getName(Euk_H2B_N_Ter_Tails))],file.out = "~/ASGARD/HMMSEARCH/All_N_Tails_RERVERSED_H2B.fa")  
#export reversed tails per family
write.fasta(AllNTails.Rev[which(getName(AllNTails)%in%getName(Euk_H3_N_Ter_Tails))], getName(AllNTails)[which(getName(AllNTails)%in%getName(Euk_H3_N_Ter_Tails))],file.out = "~/ASGARD/HMMSEARCH/All_N_Tails_RERVERSED_H3.fa")  
write.fasta(AllNTails.Rev[which(getName(AllNTails)%in%getName(Euk_H4_N_Ter_Tails))], getName(AllNTails)[which(getName(AllNTails)%in%getName(Euk_H4_N_Ter_Tails))],file.out = "~/ASGARD/HMMSEARCH/All_N_Tails_RERVERSED_H4.fa")  

#Archaeal Non-asgard
ArchaeaTails=read.fasta("~/ASGARD/HMMSEARCH/PF00808_vs-ArchaeaEBMC2_hits_InGTDB_noHalo_noASGARDS_ali_curatedNTail_alidefined.fa",seqtype = "AA")
ArchaeaTails=ArchaeaTails[which(getLength(ArchaeaTails)>0)]
ArchaeaTails.Rev=lapply(getSequence(ArchaeaTails),function(x) rev(x))
write.fasta(ArchaeaTails.Rev,getName(ArchaeaTails),"~/ASGARD/HMMSEARCH/All_N_Tails_RERVERSED.fa")



###########################################################
#recode amino acid sequences into their boichemical properties
#https://www.cell.com/iscience/fulltext/S2589-0042(22)01866-1?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2589004222018661%3Fshowall%3Dtrue
# Mapping of amino acids to a single letter based on biochemical properties
aa_map <- list(
  
  'V' = 'V', # Hydrophobic
  'I' = 'V',
  'L' = 'V',
  'M' = 'V',

  'A' = 'S', #Small
  'P' = 'S',
  'S' = 'S',
  'G' = 'S',
  'T' = 'S',

  'K' = 'K',  #Positively charged
  'R' = 'K',
  'H' = 'K',
  
  'D' = 'E',  # Negatively charged
  'E' = 'E',
  
  'Q' = 'Q',  # Polar uncharged
  'N' = 'Q',
  
  
  
  'F' = 'F', #Aromatic
  'Y' = 'F',
  'W' = 'F',
  
  'C' = 'C',
  
  'X' = 'X' #NA residues
)



# Decode amino acid sequence based on biochemical properties
recode_sequence <- function(sequence) {
  recoded <- sapply(sequence, function(aa) aa_map[[aa]])
  return(recoded)
}



AllNTails.Rev.recoded=AllNTails

#Input must be a list of sequnces like you ge from seqinr getsequence()
RecodedSeq=lapply(AllNTails.Rev, function(x) recode_sequence(x))

  write.fasta(RecodedSeq,getName(AllNTails.Rev.recoded),file.out = "~/ASGARD/HMMSEARCH/All_N_Tails_RERVERSED_recoded.fa")
  

  
  #export recoded tails per family
  write.fasta(RecodedSeq[which(getName(AllNTails.Rev.recoded)%in%getName(ASGARD_N_Ter_Tails))], getName(AllNTails.Rev.recoded)[which(getName(AllNTails.Rev.recoded)%in%getName(ASGARD_N_Ter_Tails))],file.out = "~/ASGARD/HMMSEARCH/All_N_Tails_RERVERSED_recoded_ASGARD.fa")
  #export recoded tails per family
  write.fasta(RecodedSeq[which(getName(AllNTails.Rev.recoded)%in%getName(Euk_H2A_N_Ter_Tails))], getName(AllNTails.Rev.recoded)[which(getName(AllNTails.Rev.recoded)%in%getName(Euk_H2A_N_Ter_Tails))],file.out = "~/ASGARD/HMMSEARCH/All_N_Tails_RERVERSED_recoded_H2A.fa")
  #export recoded tails per family
  write.fasta(RecodedSeq[which(getName(AllNTails.Rev.recoded)%in%getName(Euk_H2B_N_Ter_Tails))], getName(AllNTails.Rev.recoded)[which(getName(AllNTails.Rev.recoded)%in%getName(Euk_H2B_N_Ter_Tails))],file.out = "~/ASGARD/HMMSEARCH/All_N_Tails_RERVERSED_recoded_H2B.fa")  
  #export recoded tails per family
  write.fasta(RecodedSeq[which(getName(AllNTails.Rev.recoded)%in%getName(Euk_H3_N_Ter_Tails))], getName(AllNTails.Rev.recoded)[which(getName(AllNTails.Rev.recoded)%in%getName(Euk_H3_N_Ter_Tails))],file.out = "~/ASGARD/HMMSEARCH/All_N_Tails_RERVERSED_recoded_H3.fa")  
  write.fasta(RecodedSeq[which(getName(AllNTails.Rev.recoded)%in%getName(Euk_H4_N_Ter_Tails))], getName(AllNTails.Rev.recoded)[which(getName(AllNTails.Rev.recoded)%in%getName(Euk_H4_N_Ter_Tails))],file.out = "~/ASGARD/HMMSEARCH/All_N_Tails_RERVERSED_recoded_H4.fa")  
  
#Recode and export archaea non asgard histones
  
  ArchaeaTails.Rev.Recode=ArchaeaTails
  ArchaeaTails.RevUp=toupper(ArchaeaTails.Rev)
  RecodedSeqArchaea=lapply(ArchaeaTails.Rev, function(x) recode_sequence(x))
  write.fasta(RecodedSeqArchaea,getName(ArchaeaTails.Rev.Recode),file.out = "~/ASGARD/HMMSEARCH/All_N_Tails_ARCHAEA_RERVERSED_recoded.fa")
  
  
  
  #Addendum: some stats on average length for the main text.
  boxplot(getLength(ArchaeaTails))
  mean(getLength(ArchaeaTails))
  mean(getLength(ASGARD_N_Ter_Tails))
  median(getLength(ArchaeaTails))
  median(getLength(ASGARD_N_Ter_Tails))
  
  
  

  
  