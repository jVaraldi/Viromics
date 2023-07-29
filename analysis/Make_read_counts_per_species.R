

grouping_info=read.table("/Users/jvaraldi/Documents/Manips/Viromics_Metagenomics/Snakemake/TABLES/summary_table.txt", h=T, sep="\t")
# ambidensovirus 2 is in fact a subset of ambidensovirus 1 (contig_18555 is included into contig_15192)
#grouping_info=grouping_info[-which(grouping_info$virus_name=="Ambidensovirus2 n=2"),]

# add genomic structure info
grouping_info$genome="RNA"
grouping_info$genome[which(grouping_info$virus_name %in% c("Parvoviridae_Pachy n=1", 
                                                           "Linvill_road_virus_D.sim n=2", 
                                                           "Vesantovirus_D.mel n=11", 
                                                           "Vesantovirus_D.sub n=10", 
                                                           "LbFV_L.b n=10",
                                                           "Parvoviridae2 n=1", 
                                                           "LhFV_L.h n=12"))] = "DNA"

n_viruses=length(unique(grouping_info$virus_name))
viruses=unique(grouping_info$virus_name)
  
# import data (number of reads mapped on each contig from each sample)
wga_reads_mapped=read.table(file = "../TABLES/wga_nreads_viruses_plus_unassigned.txt", h=T)
wta_reads_mapped=read.table(file = "../TABLES/wta_nreads_viruses_plus_unassigned.txt", h=T)
names(wga_reads_mapped)=sub(pattern = "_WGA", replacement = "", x = names(wga_reads_mapped))
names(wta_reads_mapped)=sub(pattern = "_WTA", replacement = "", x = names(wta_reads_mapped))
# merge tables wga and wta
all_reads_mapped=rbind(wga_reads_mapped,wta_reads_mapped)
dim(all_reads_mapped)

results=list()
for (i in 1:n_viruses){
  contigs=grouping_info$contig_name[grouping_info$virus_name==viruses[i]]
  subset=apply(X = all_reads_mapped[which(rownames(all_reads_mapped) %in% contigs), ], 2, sum)
  species=unlist(lapply(strsplit(names(subset), "_"), FUN = function(x){return(x[[1]])}))
  subset_table=data.frame(subset, species)
  # total number of reads per species
  res=tapply(subset_table$subset, subset_table$species, sum)
  # number of samples per species
  res_n_samples=tapply(subset_table$subset, subset_table$species, length)
  # normalize read number per number of samples
  res=res/res_n_samples
  
  res_n_samples=t(as.data.frame(res_n_samples))
  res=t(as.data.frame(res)) # transpose
  total=sum(res)
  res=cbind(viruses[i], res, total)
  colnames(res)=c("virus", paste0(colnames(res)[-c(1, 18)], " n=",res_n_samples), "total")
  results[[i]]=res
}

table_reads=Reduce(function(...) merge(..., all = TRUE),
       results)

table_reads_ordered=table_reads[order(table_reads$virus),]

write.table(table_reads_ordered, file = "../TABLES/table_reads_per_species_per_sample.txt", row.names = TRUE, col.names = TRUE, quote=FALSE, sep="\t")
#####################
# heatmaps
######################
library(ComplexHeatmap)
library(circlize)
col_fun = circlize::colorRamp2(c(0, 1), c("white", "black"))
col_fun(seq(-3, 3))

pdf("../figures/heatmap_n-reads-mapped.pdf")
t=table_reads_ordered[,-1]
t2=apply(t, 2, as.numeric)
rownames(t2)=table_reads_ordered$virus
genome_struct=factor(grouping_info$genome[match(rownames(t2), grouping_info$virus_name)])
levels(genome_struct)=c( "#009999","#000000")
t=t2[,-18]

# normalize by row
total_reads=apply(t, 1, max)
t3=t/total_reads

# add percent identity
# collect information on sequences (percentage identity...)
infos_names_wga=read.table("../TABLES/grouped_data_wga2.txt", h=T, sep="\t")
infos_names_wta=read.table("../TABLES/grouped_data_wta2.txt", h=T, sep="\t")
infos_names=rbind(infos_names_wta[,c("pident", "seqs_length" ,"ncontigs")], infos_names_wga[,c("pident", "seqs_length" ,"ncontigs")])
infos_names$pident=round(infos_names$pident*100,digits = 1) # round percentages
# create a new name including percent id and length
infos_names$newname=paste0(rownames(infos_names)," ",infos_names$pident , "% [",infos_names$seqs_length, "bp]")

rownames(t3)=infos_names$newname[match(rownames(t3), rownames(infos_names))]
t3=as.matrix(t3)
######


Heatmap(t3[,-c(3, 12,17)], row_names_gp  = gpar(fontsize = 7, col = as.character(genome_struct)), 
        col = col_fun, name="proportion of reads mapped\n      relative to the maximum",column_names_gp = gpar(fontsize=7),
        border = FALSE, rect_gp = gpar(col = "lightgrey"))

dev.off()


#####################
### another way of normalizing : by the number of isofemale lines
#####################

insects=c("C.amo", "D.hyd", "D.im", "D.kuntzei", "D.mel", "D.sim", 
          "D.sub", "D.sub.obs", "D.suz", "A.ruf.tab", "L.b", "L.h","Pachy", "Tricho")
n_isofemales=c(1,39,116, 8, 230, 51, 62, 19, 14, 9, 89, 108, 30, 23)
tab_n_iso=t(data.frame(n_isofemales, row.names = insects))
tab_n_iso=tab_n_iso[, sort(colnames(tab_n_iso))]

results=list()
for (i in 1:n_viruses){
  contigs=grouping_info$contig_name[grouping_info$virus_name==viruses[i]]
  subset=apply(X = all_reads_mapped[which(rownames(all_reads_mapped) %in% contigs), ], 2, sum)
  species=unlist(lapply(strsplit(names(subset), "_"), FUN = function(x){return(x[[1]])}))
  subset_table=data.frame(subset, species)
  # total number of reads per species
  res=tapply(subset_table$subset, subset_table$species, sum)
  # number of samples per species
  res_n_samples=tapply(subset_table$subset, subset_table$species, length)
  # remove controls
  res=res[-c(3,12)]
  res_n_samples=res_n_samples[-c(3,12)]
  # normalize read number per number of isofemale lines
  res=res/tab_n_iso
  
  res_n_samples=t(as.data.frame(res_n_samples))
  res=t(as.data.frame(res)) # transpose
  res=cbind(viruses[i], res)
  colnames(res)=c("virus", paste0(colnames(res)[-1], " n=",res_n_samples,"/",tab_n_iso))
  results[[i]]=res
}

table_reads=Reduce(function(...) merge(..., all = TRUE),
                   results)

table_reads_ordered=table_reads[order(table_reads$virus),]
table_reads_ordered[,-1]=apply(X = table_reads_ordered[,-1], 2, as.numeric)

write.table(table_reads_ordered, file = "../TABLES/table_reads_per_species_per_isofemale_line.txt", row.names = TRUE, col.names = TRUE, quote=FALSE, sep="\t")
#####################
# heatmaps
######################
library(ComplexHeatmap)
library(circlize)

pdf("../figures/heatmap_n-reads-mapped-per_isofemale_line.pdf")
t2=table_reads_ordered[,-1]
rownames(t2)=table_reads_ordered$virus
genome_struct=factor(grouping_info$genome[match(rownames(t2), grouping_info$virus_name)])
levels(genome_struct)=c( "#009999","#000000")

# normalize by row
#total_reads=apply(t2, 1, max)
#t3=t/total_reads

# add percent identity
# collect information on sequences (percentage identity...)
infos_names_wga=read.table("../TABLES/grouped_data_wga2.txt", h=T, sep="\t")
infos_names_wta=read.table("../TABLES/grouped_data_wta2.txt", h=T, sep="\t")
infos_names=rbind(infos_names_wta[,c("pident", "seqs_length" ,"ncontigs")], infos_names_wga[,c("pident", "seqs_length" ,"ncontigs")])
infos_names$pident=round(infos_names$pident*100,digits = 1) # round percentages
# create a new name including percent id and length
infos_names$newname=paste0(rownames(infos_names)," ",infos_names$pident , "% [",infos_names$seqs_length, "bp]")

rownames(t2)=infos_names$newname[match(rownames(t2), rownames(infos_names))]
t2=as.matrix(t2)
######

# scale by rows
t3=t2/apply(t2, 1, sum)
#scaled_mat = t(scale(t(t2)))

# define color scheme accordingly
col_fun = circlize::colorRamp2(c(min(t3), max(t3)), c("white", "black"))
col_fun(seq(-3, 3))

Heatmap(t3, row_names_gp  = gpar(fontsize = 7, col = as.character(genome_struct)), 
        col = col_fun, name="relative\n number\n of reads",column_names_gp = gpar(fontsize=7),
        border = FALSE, rect_gp = gpar(col = "lightgrey"))




dev.off()




#############################################
#### write sequences of each virus separately in a file
############################################

library(Biostrings)
seqs=readBStringSet("../sequences/assembly.fasta")

# modify one virus name that mess with the writing:
viruses[3]=sub(pattern = "/", replacement = "_", x = viruses[3])

for (i in 1:n_viruses){
  print(i)
  contigs_subset=grouping_info$contig_name[grouping_info$virus_name==viruses[i]]
  seqs_subset=seqs[contigs_subset]
  writeXStringSet(seqs_subset, filepath = paste0("../sequences/virus_genomes/",viruses[i], ".fa" ))
}






