# ------------------------------------------------------------------------------
# Script Name: RemovePFAMHitsfromtails.r
# Author: Antoine Hocher
# Date: 2023-10-16
# Description: remove hits with a structure domain in the list of sequences that are on the sides of the core histone fold
#basically stuff like macroH2A.
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

library(rhmmer)
library(seqinr)



RemoveDomainsContainingSequences=function(FastaIn){
  
  Tableout=gsub(".fa","_tblout.txt",FastaIn)
  TableFout=gsub(".fa","_tbFullout.txt",FastaIn)
  SeqNoDomainOut=gsub(".fa","_NoPfAM_A.fa",FastaIn)
  SeqNoDomainOutNoStar=gsub(".fa","_NoPfAM_A_nostar.fa",FastaIn)
  
  SeqNoDomainOutOrdered=gsub(".fa","_NoPfAM_A_LengthOrdered.fa",FastaIn)
  SeqNoDomainOutOrderedNoStar=gsub(".fa","_NoPfAM_A_LengthOrdered_nostar.fa",FastaIn)
  
Commandline=paste("/Users/ahocher/Documents/Softwares/hmmer-3.4/src/hmmsearch --cut_ga --tblout ",Tableout," ","~/THERMOPHILY/HMM_MODELS/Pfam-A.hmm"," ",FastaIn,"  > ",TableFout,sep="")
system(Commandline)


Domains=read_tblout(Tableout)
FastaSeq=read.fasta(FastaIn)
if(length(Domains$domain_name)>0){
  
  #Remove domains which we know are affiliated to histone, we don't want to remove those
  if(length(which(Domains$query_name%in%c("Histone_H2A_C","CBFD_NFYB_HMF","Histone")))>0){
Domains=Domains[-which(Domains$query_name%in%c("Histone_H2A_C","CBFD_NFYB_HMF","Histone","CENP-T_C","TAF")),]}
  
  if(length(which(getName(FastaSeq)%in%Domains$domain_name))>0){
  FastaSeq.r=FastaSeq[-which(getName(FastaSeq)%in%Domains$domain_name)]
  FastaSeq.ro=FastaSeq.r[order(getLength(FastaSeq.r))]
  write.fasta(getSequence(FastaSeq.r),getName(FastaSeq.r),SeqNoDomainOut)
  write.fasta(getSequence(FastaSeq.ro),getName(FastaSeq.ro),SeqNoDomainOutOrdered)}
  system(paste("tr -d '*' < ",SeqNoDomainOut," > ",SeqNoDomainOutNoStar,sep=""))
  system(paste("tr -d '*' < ",SeqNoDomainOutOrdered," > ",SeqNoDomainOutOrderedNoStar,sep=""))
  
}
}



RemoveDomainsContainingSequences(FastaIn="~/ASGARD/HMMSEARCH/PF00125_vs-EUK_hits_MMSEQS_no_duplicates.fa")


RemoveDomainsContainingSequences(FastaIn="~/ASGARD/HMMSEARCH/Individual_Euk_HF_Pooled_no_duplicates.fa")



#Annotating NFYC, a clear outlier for which PFAM (A, B might be ok) is not perfomiong well
#We will use it to annotate histones later on
/Users/ahocher/Documents/Softwares/hmmer-3.4/src/jackhmmer -E 1e-10--tblout ~/ASGARD/HMMSEARCH/Histones_vs-EUK_hits_MMSEQS_no_duplicates_vsNFYC_tblout.txt ~/ASGARD/HMMSEARCH/Individual_Euk_HF_Pooled_no_duplicates_NoPfAM_A_nostar.fa > ~/ASGARD/HMMSEARCH/Histones_vs-EUK_hits_MMSEQS_no_duplicates_vsNFYC.txt

