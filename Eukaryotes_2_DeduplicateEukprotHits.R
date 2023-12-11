# ------------------------------------------------------------------------------
# Script Name: Deduplicate_EukprotHits.r
# Author: Antoine Hocher
# Date: 2023-10-16
# Description: remove hits that are duplicated in euk prot ( those are duplicated because they are computed from transcripts, we want the longest isoform)
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

library(seqinr)


deduplicate=function(InputFasta){
Histones=read.fasta(InputFasta)
Histones.nodup=Histones
NamesHistones=getName(Histones)
Specie=unlist(lapply(NamesHistones,function(x) strsplit(x,"\\$")[[1]][1]))

DedupList=c()


#Now using MMSEQS
for(i in unique(Specie)){
    #Get the name of duplicated sequences
    SpecieHistone=Histones[which(Specie==i)]
  
    setwd("~/ASGARD/HMMSEARCH/temp")
    write.fasta(getSequence(SpecieHistone),getName(SpecieHistone),"temp.fa")
    
    #Compute clusters using mmseqs2 (note the settings here only cluster almost identical sequences)
    Command=paste0("mmseqs easy-cluster ","temp.fa"," ","temp"," ","temp"," --min-seq-id 0.9 -c 0.9 --cov-mode 0 -s 1")
    system(Command)
    
    #Onyl retrieve representatives
    DedupList=c(DedupList,read.fasta("temp_rep_seq.fasta"))
    
}
OutputFasta=gsub(".fa","_no_duplicates.fa",InputFasta)
write.fasta(getSequence(DedupList),getName(DedupList),OutputFasta)}



deduplicate("~/ASGARD/HMMSEARCH/PF00125_vs-EUK_hits.fa")

deduplicate("~/ASGARD/HMMSEARCH/Individual_Euk_HF_Pooled.fa")

