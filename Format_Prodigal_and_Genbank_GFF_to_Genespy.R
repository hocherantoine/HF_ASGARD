# ------------------------------------------------------------------------------
# Script Name: Format prodigal and genbank gff
# Author: Antoine Hocher
# Date: 2023-10-10
# Description: Function to Format Genkank and prodigal gff to be compatible with gene spy
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

library(rtracklayer)

ProdigalGFFs_files=list.files("~/ASGARD/GENOMES/GFF/PRODIGAL_ASGARDS")
setwd("~/ASGARD/GENOMES/GFF/PRODIGAL_ASGARDS")


for (i in ProdigalGFFs_files){
  MyGFF=readGFF(i)
  Assembly=paste("GCA_",strsplit(i,"_")[[1]][2],sep="")
  Temp=MyGFF$ID
  Temp2=paste(gsub("\\.1","",as.character(MyGFF$seqid)),unlist(lapply(MyGFF$ID,function(x) strsplit(x,"_")[[1]][2])),sep="-")
  
  MyGFFMod=MyGFF
  MyGFFMod$ID=Temp2
  MyGFFMod$locus_tag=Temp2
  MyGFFMod$protein_id=Temp2
  MyGFFMod$Name=Temp2
  MyGFFMod$ID=Temp2
  
  
  MyGFFMod2=MyGFFMod
  MyGFFMod2$type="gene"
  MyGFFMod3=rbind(MyGFFMod,MyGFFMod2)
  export.gff3(MyGFFMod3,paste("~/ASGARD/GENOMES/GFF/PRODIGAL_ASGARDS/MODS/",Assembly,".gff",sep=""))
  
}





#Also modifying genbank to remove the cds- in the ID


GenbankGFFs_files=list.files("~/ASGARD/GENOMES/GFF/HODARCHAEA/",pattern = ".gff")
setwd("~/ASGARD/GENOMES/GFF/HODARCHAEA/")




for (i in GenbankGFFs_files){
  MyGFF=rtracklayer::readGFF(i)
  Assembly=paste("GCA_",strsplit(i,"_")[[1]][2],sep="")
  Temp=MyGFF$ID
  Temp2=gsub("cds-","",as.character(Temp))
  
  MyGFFMod=MyGFF
  MyGFFMod$ID=Temp2
  MyGFFMod$locus_tag=Temp2
  MyGFFMod$protein_id=Temp2
  MyGFFMod$Name=Temp2
  MyGFFMod$ID=Temp2
  
  export.gff3(MyGFFMod,paste("~/ASGARD/GENOMES/GFF/HODARCHAEA/GENBANKMOD/",Assembly,".gff",sep=""))
  
}





#For Loki


#Also modifying genbank to remove the cds- in the ID
GenbankGFFs_files=list.files("~/ASGARD/GENOMES/GFF/LOKI/GENBANK_LOKI/",pattern = ".gff")
setwd("~/ASGARD/GENOMES/GFF/LOKI/GENBANK_LOKI/")




for (i in GenbankGFFs_files){
  MyGFF=rtracklayer::readGFF(i)
  Assembly=paste("GCA_",strsplit(i,"_")[[1]][2],sep="")
  Temp=MyGFF$ID
  Temp2=gsub("cds-","",as.character(Temp))
  
  MyGFFMod=MyGFF
  MyGFFMod$ID=Temp2
  MyGFFMod$locus_tag=Temp2
  MyGFFMod$protein_id=Temp2
  MyGFFMod$Name=Temp2
  MyGFFMod$ID=Temp2
  
  export.gff3(MyGFFMod,paste("~/ASGARD/GENOMES/GFF/LOKI/GENBANK_LOKI_MOD/",Assembly,".gff",sep=""))
  
}





#For All asgards



#Also modifying genbank to remove the cds- in the ID
GenbankGFFs_files=list.files("~/ASGARD/GENOMES/GFF/ALL_ASGARDS/",pattern = ".gff")
setwd("~/ASGARD/GENOMES/GFF/ALL_ASGARDS/")




for (i in GenbankGFFs_files){
  MyGFF=rtracklayer::readGFF(i)
  Assembly=paste("GCA_",strsplit(i,"_")[[1]][2],sep="")
  Temp=MyGFF$ID
  Temp2=gsub("cds-","",as.character(Temp))
  
  MyGFFMod=MyGFF
  MyGFFMod$ID=Temp2
  MyGFFMod$locus_tag=Temp2
  MyGFFMod$protein_id=Temp2
  MyGFFMod$Name=Temp2
  MyGFFMod$ID=Temp2
  
  export.gff3(MyGFFMod,paste("~/ASGARD/GENOMES/GFF/ALL_ASGARDS/MOD/",Assembly,".gff",sep=""))
  
}
