# ------------------------------------------------------------------------------
# Script Name: Retrieve_eukaryotic_Histones.R
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


##############
#Retrieving Eukaryotic HIstones
#Using PFAM-A domains
/Users/ahocher/Documents/Softwares/hmmer-3.4/src/hmmsearch --cut_ga -A ~/ASGARD/HMMSEARCH/PF00125_vs-eukprot_CLS_ali.fa --tblout ~/ASGARD/HMMSEARCH/PF00125_vs-eukprot_CLS.txt ~/ASGARD/HMM/PF00125.hmm ~/CLK1/EukProt/SeqDB/Concat_eukprot_CLS.fa  > ~/ASGARD/HMMSEARCH/PF00125_vs-eukprot_CLS_full.txt 

/Users/ahocher/Documents/Softwares/hmmer-3.4/src/hmmsearch --cut_ga -A ~/ASGARD/HMMSEARCH/PF00808_vs-eukprot_CLS_ali.fa --tblout ~/ASGARD/HMMSEARCH/PF00808_vs-eukprot_CLS.txt ~/ASGARD/HMM/PF00808.hmm ~/CLK1/EukProt/SeqDB/Concat_eukprot_CLS.fa  > ~/ASGARD/HMMSEARCH/PF00808_vs-eukprot_CLS_full.txt 


#Using H2A, using Nick's HMM models
/Users/ahocher/Documents/Softwares/hmmer-3.4/src/hmmsearch -E 1e-10 -A ~/ASGARD/HMMSEARCH/H2A_vs-eukprot_CLS_HMM.fa --tblout ~/ASGARD/HMMSEARCH/H2A_vs-eukprot_CLS.txt ~/ASGARD/EXTERNAL_DATA/Figshare_data_N_Irvin_Bioarkiv/HMMs/H2A.hmm ~/CLK1/EukProt/SeqDB/Concat_eukprot_CLS.fa  > ~/ASGARD/HMMSEARCH/H2A_vs-eukprot_CLS_full.txt 
#Using H2B, using Nick's HMM models
/Users/ahocher/Documents/Softwares/hmmer-3.4/src/hmmsearch -E 1e-10 -A ~/ASGARD/HMMSEARCH/H2B_vs-eukprot_CLS_HMM.fa --tblout ~/ASGARD/HMMSEARCH/H2B_vs-eukprot_CLS.txt ~/ASGARD/EXTERNAL_DATA/Figshare_data_N_Irvin_Bioarkiv/HMMs/H2B.hmm ~/CLK1/EukProt/SeqDB/Concat_eukprot_CLS.fa  > ~/ASGARD/HMMSEARCH/H2B_vs-eukprot_CLS_full.txt 
#Using H3, using Nick's HMM models
/Users/ahocher/Documents/Softwares/hmmer-3.4/src/hmmsearch -E 1e-10 -A ~/ASGARD/HMMSEARCH/H3_vs-eukprot_CLS_HMM.fa --tblout ~/ASGARD/HMMSEARCH/H3_vs-eukprot_CLS.txt ~/ASGARD/EXTERNAL_DATA/Figshare_data_N_Irvin_Bioarkiv/HMMs/H3.hmm ~/CLK1/EukProt/SeqDB/Concat_eukprot_CLS.fa  > ~/ASGARD/HMMSEARCH/H3_vs-eukprot_CLS_full.txt 
#Using H4, using Nick's HMM models
/Users/ahocher/Documents/Softwares/hmmer-3.4/src/hmmsearch -E 1e-10 -A ~/ASGARD/HMMSEARCH/H4_vs-eukprot_CLS_HMM.fa --tblout ~/ASGARD/HMMSEARCH/H4_vs-eukprot_CLS.txt ~/ASGARD/EXTERNAL_DATA/Figshare_data_N_Irvin_Bioarkiv/HMMs/H4.hmm ~/CLK1/EukProt/SeqDB/Concat_eukprot_CLS.fa  > ~/ASGARD/HMMSEARCH/H4_vs-eukprot_CLS_full.txt


#Reformat alignement  of the HMM hit domain

#Change hmm alignement output format
/Users/ahocher/Documents/Softwares/hmmer-3.4/easel/miniapps/esl-reformat a2m ~/ASGARD/HMMSEARCH/H2A_vs-eukprot_CLS_HMM.fa > ~/ASGARD/HMMSEARCH/H2A_vs-eukprot_CLS_HMM_ali.fa
/Users/ahocher/Documents/Softwares/hmmer-3.4/easel/miniapps/esl-reformat a2m ~/ASGARD/HMMSEARCH/H2B_vs-eukprot_CLS_HMM.fa > ~/ASGARD/HMMSEARCH/H2B_vs-eukprot_CLS_HMM_ali.fa
/Users/ahocher/Documents/Softwares/hmmer-3.4/easel/miniapps/esl-reformat a2m ~/ASGARD/HMMSEARCH/H3_vs-eukprot_CLS_HMM.fa > ~/ASGARD/HMMSEARCH/H3_vs-eukprot_CLS_HMM_ali.fa
/Users/ahocher/Documents/Softwares/hmmer-3.4/easel/miniapps/esl-reformat a2m ~/ASGARD/HMMSEARCH/H4_vs-eukprot_CLS_HMM.fa > ~/ASGARD/HMMSEARCH/H4_vs-eukprot_CLS_HMM_ali.fa


  
#Retrieve full sequences
#index
/Users/ahocher/Documents/Softwares/hmmer-3.4/easel/miniapps/esl-sfetch --index ~/CLK1/EukProt/SeqDB/Concat_eukprot_CLS.fa

#OUtput sequence names from the Hmm res file
grep -v "^#" ~/ASGARD/HMMSEARCH/H2A_vs-eukprot_CLS.txt | awk '{print $1}' > ~/ASGARD/HMMSEARCH/H2A_vs-EUK_hits_Names.txt
grep -v "^#" ~/ASGARD/HMMSEARCH/H2B_vs-eukprot_CLS.txt | awk '{print $1}' > ~/ASGARD/HMMSEARCH/H2B_vs-EUK_hits_Names.txt
grep -v "^#" ~/ASGARD/HMMSEARCH/H3_vs-eukprot_CLS.txt | awk '{print $1}' > ~/ASGARD/HMMSEARCH/H3_vs-EUK_hits_Names.txt
grep -v "^#" ~/ASGARD/HMMSEARCH/H4_vs-eukprot_CLS.txt | awk '{print $1}' > ~/ASGARD/HMMSEARCH/H4_vs-EUK_hits_Names.txt

grep -v "^#" ~/ASGARD/HMMSEARCH/PF00125_vs-eukprot_CLS.txt | awk '{print $1}' > ~/ASGARD/HMMSEARCH/PF00125_vs-EUK_hits_Names.txt





#Command has to be run within another architecture
use_x86 

PathToGenomes="~/CLK1/EukProt/SeqDB/Concat_eukprot_CLS.fa"

SeqNames="~/ASGARD/HMMSEARCH/PF00125_vs-EUK_hits_Names.txt"
FastaNames=gsub("_Names.txt",".fa",SeqNames)

Command=paste("seqtk subseq ",dirname(PathToGenomes),"/",basename(PathToGenomes)," ",SeqNames," > ",FastaNames,sep="")
system(Command)

SeqNames="~/ASGARD/HMMSEARCH/H2A_vs-EUK_hits_Names.txt"
FastaNames=gsub("_Names.txt",".fa",SeqNames)

Command=paste("seqtk subseq ",dirname(PathToGenomes),"/",basename(PathToGenomes)," ",SeqNames," > ",FastaNames,sep="")
system(Command)

SeqNames="~/ASGARD/HMMSEARCH/H2B_vs-EUK_hits_Names.txt"
FastaNames=gsub("_Names.txt",".fa",SeqNames)

Command=paste("seqtk subseq ",dirname(PathToGenomes),"/",basename(PathToGenomes)," ",SeqNames," > ",FastaNames,sep="")
system(Command)

SeqNames="~/ASGARD/HMMSEARCH/H3_vs-EUK_hits_Names.txt"
FastaNames=gsub("_Names.txt",".fa",SeqNames)

Command=paste("seqtk subseq ",dirname(PathToGenomes),"/",basename(PathToGenomes)," ",SeqNames," > ",FastaNames,sep="")
system(Command)

SeqNames="~/ASGARD/HMMSEARCH/H4_vs-EUK_hits_Names.txt"
FastaNames=gsub("_Names.txt",".fa",SeqNames)

Command=paste("seqtk subseq ",dirname(PathToGenomes),"/",basename(PathToGenomes)," ",SeqNames," > ",FastaNames,sep="")
system(Command)


library(seqinr)
setwd("~/ASGARD/HMMSEARCH/")
#Agregate H2A,H2B,H3,H4 and remove duplicated names (Will re-attribute which type is which later on)
H2A=read.fasta("H2A_vs-EUK_hits.fa")
H2B=read.fasta("H2B_vs-EUK_hits.fa")
H3=read.fasta("H3_vs-EUK_hits.fa")
H4=read.fasta("H4_vs-EUK_hits.fa")

AllHF=c(H2A,H2B,H3,H4)
AllHF.r=AllHF[which(duplicated(getName(AllHF))==F)]
write.fasta(getSequence(AllHF.r),getName(AllHF.r),"Individual_Euk_HF_Pooled.fa")

