suppressMessages(library(ggplot2))
suppressMessages(library(data.table))
#suppressMessages(library(reshape2))
#suppressMessages(library(plyr))
#suppressMessages(library(dplyr))

filter.bygene=function(position, start, end) {
  if (start <= position && position <= end) {
    
    return(position)
  }
  else {
    
    return(NA)
  }
  
}

#THIS FUNCTION FILTERS THE SUPPLIED DATA FRAME MAKING SURE THAT IT IS ON THE CORRECT CHROMOSOME AND FALLS WITHIN
#THE GENE BEING ANALYZED
filter.master=function(query, reference){
  
  chrom=as.character(unique(gene$CHROM))
  gene.start=reference[reference$FEATURE == "gene", 3]
  gene.end=reference[reference$FEATURE == "gene", 4]
  gene.transcripts=reference[reference$FEATURE == "transcript",]
  gene.transcripts=gene.transcripts[order(gene.transcripts$START),]
  
  query=query[query$CHROM == chrom,]
  position=query$START
  starts.filtered=na.omit(sapply(position, filter.bygene, start=gene.start, end=gene.end))
  query=query[query$START %in% starts.filtered,]
  return(query)
}
#READ IN ARGUMENTS
args=commandArgs(trailingOnly=TRUE)
inputtable = args[1]  
inputtid = args[2]
gene.name = args[3]


#THESE ARE THE NAMES OF THE PLOTS THAT GET GENERATED
file1 = paste("output/", gene.name, ".1.pdf", sep = "")
file2 = paste("output/", gene.name, ".2.pdf", sep = "")
file3 = paste("output/", gene.name, ".3.pdf", sep = "")
file4 = paste("output/", gene.name, ".output.table.txt", sep = "")


################### READ IN OUR DATA############################################################
recombine.v3=read.table("resources/ActiveCarrierMapMutations.sorted.bed")
colnames(recombine.v3)=c("CHROM","START","END","RECOMBINE_ID")
recombine.v3=recombine.v3[!duplicated(recombine.v3$START),]

#ASSIGN THE TRANSCRIPT ID'S 
tid=read.table(inputtid)
tid$V1=as.character(tid$V1)
ids=na.omit(unique(tid))

#gene = read.table("resources/TEMP.gene.table.txt")
gene=read.table(inputtable)
colnames(gene)=c("CHROM","FEATURE","START","END")
gene$TID = tid$V1
chrom=as.character(unique(gene$CHROM))

trusight = read.table("resources/TruSightOne.nochr.txt")
colnames(trusight) = c("CHROM","START","END","AMP_ID")
trusight = trusight[!duplicated(trusight$START),]

nextera = read.table("resources/nextera.expanded.nochr.txt")
colnames(nextera) = c("CHROM","START","END","NAME")
nextera = nextera[!duplicated(nextera$START),]

empty_df=data.frame(CHROM=0, START=0, END=0, NAME=0)
#SUBSET THE HOMO SAPIENS DATA BY GENETIC FEATURE

gene.justexon=gene[gene$FEATURE=="exon",]
gene.justexon$FEATURE="Exon"
gene.justCDS=gene[gene$FEATURE=="CDS",]
gene.justCDS[,"Y"]=1
gene.justCDS$FEATURE="Coding Sequence"
gene.justUTR=gene[gene$FEATURE=="UTR",]
gene.transcripts=gene[gene$FEATURE=="transcript",]

gene.start=gene[gene$FEATURE=="gene", 3]
gene.end=gene[gene$FEATURE=="gene", 4]

##### THE FOLLOWING IS USED USED TO GET THE NUMBER OF TID'S SO THEY 
##### CAN BE BROKEN UP AND DISPLAYED MORE EASILY IN THE 3RD PLOT OUTPUTTED

gene.justCDS.tid=as.character(gene.justCDS[!duplicated(gene.justCDS$TID),5])
gene.justCDS.num=as.numeric(length(gene.justCDS.tid))
CDSwidth=20/gene.justCDS.num
divisions=1/(gene.justCDS.num + 1)

for (i in 1:gene.justCDS.num) {
  TID=gene.justCDS.tid[i]
  Ynew = 0.5 + (i*divisions)
  CDS.subset = gene.justCDS[gene.justCDS$TID == TID,]
  rows=rownames(CDS.subset)
  gene.justCDS[rows,6] = Ynew
}

################## THIS IS USED FOR SCALING THE OUTPUT PLOTS ####################
numtranscripts1=as.numeric(nrow(gene.justexon[!duplicated(gene.justexon$TID),]))
numtranscripts2=as.numeric(nrow(gene.justCDS[!duplicated(gene.justCDS$TID),]))
numcol1=ceiling(numtranscripts1/4)
numcol2=ceiling(numtranscripts2/4)

################# FILTER OUR DATA SETS BY START POSITION #################
##### THE FOLLOWING CHECK THAT THERE ARE OCCURENCES OF THE GIVEN TECHNOLOGY IN THE GENE QUERIED ##

recombine.v3.df=filter.master(recombine.v3, gene)
if (nrow(recombine.v3.df) >= 1) { 
  recombine.v3.df[,"FEATURE"]="Recombine V3" 
  recombine.v3.df[,"END"]=recombine.v3.df$START + 4
} else {
recombine.v3.df=empty_df
recombine.v3.df$FEATURE="Recombine V3"
}

trusight.df=filter.master(trusight, gene)
if (nrow(trusight.df) >= 1) {
  trusight.df[,"FEATURE"]="Trusight"
} else {
  trusight.df=empty_df
  trusight.df$FEATURE="Trusight"
}

nextera.df=filter.master(nextera, gene)
if (nrow(nextera.df) >= 1) {
  nextera.df[,"FEATURE"]="Nextera"
} else {
 nextera.df=empty_df
 nextera.df$FEATURE="Nextera"
}



#TITLE OF OUR PLOTS
plottitle=paste(gene.name, " - Chrom ", chrom, sep = "")

############# DEFINE OUR COLOR SETS #########
colors1=c("dodgerblue", "limegreen", "firebrick1", "darkorange", "gold" )
colors2=c("dodgerblue", "darkslateblue", "limegreen", "firebrick1", "gold", "darkorange")

#THIS IS THE THEME USED FOR ALL GRAPHS MODIFY IT HERE TO CHANGE ALL OF THEM.
plot_theme = theme(panel.background=element_rect(fill="white"), legend.position ="right", legend.title=element_blank(),
                   strip.background=element_rect(fill="forestgreen"), strip.text.x=element_text(size=20, color="white", face="bold"),
                   axis.text.y=element_blank(), 
                   #axis.text.x=element_blank(), 
                   plot.background=element_rect(fill="gray90"), 
                   axis.title.x=element_text(size=30, color="forestgreen"), panel.grid.major.x=element_blank(), 
                   panel.grid.minor.x=element_blank(),
                   panel.grid.minor.y=element_line(size=1, color="gray90"), 
                   plot.title=element_text(size=30, face="bold", color="forestgreen"), axis.ticks.y = element_blank())

############ PLOT THE GRAPHS ########################
#THE THEME SECTION IS THE SAME FOR ALL OF THE PLOTS ENSURING A CONSISTENT LOOK

####### THIS PLOT IS BROKEN UP BY TRANSCRIPT ID ###########
plot2 = ggplot() +
  geom_segment(data=gene.justCDS, aes(x=START, y=1, xend=END, yend=1, color=FEATURE),size=20) +
  geom_segment(data=gene.justUTR, aes(x=START, y=1, xend=END, yend=1, color=FEATURE), size=20) +
  geom_segment(data=recombine.v3.df, aes(x=START, y=2, xend=END, yend=2, color=FEATURE), size=20) +
  geom_segment(data=trusight.df, aes(x=START, y=3, xend=END, yend=3, color=FEATURE), size=20) +
  geom_segment(data=nextera.df, aes(x=START, y=4, xend=END, yend=4, color=FEATURE), size=20) + 
  facet_wrap(~TID) +
  scale_color_manual(values=colors1) +
  plot_theme +
  scale_y_continuous(breaks=1:5) +
  labs(title = plottitle) + 
  xlab("Gene Position") + ylab("")

######## THIS PLOT DISPLAYS ALL OF THE CODING SEQUENCES ON TOP OF EACHOTHER  ######
plot3 = ggplot() +
  geom_segment(data=gene.justCDS, aes(x=START, y=1, xend=END, yend=1, color=FEATURE), size=20, alpha = 0.5) +
  geom_segment(data=gene.justCDS, aes(x=START, y=1, xend=START +5, yend=1), color="dodgerblue4", size=20) +
  geom_segment(data=gene.justUTR, aes(x=START, y=1, xend=END, yend=1, color=FEATURE), size=20) +
  geom_segment(data=recombine.v3.df, aes(x=START, y=2, xend=END, yend=2, color=FEATURE), size=20) +
  geom_segment(data=trusight.df, aes(x=START, y=3, xend=END, yend=3, color=FEATURE), size=20) +
  geom_segment(data=nextera.df, aes(x=START, y=4, xend=END, yend=4, color=FEATURE), size=20) + 
  scale_color_manual(values=colors1) +
  plot_theme + 
  scale_y_continuous(breaks=1:5) +
  labs(title = plottitle) + 
  xlab("Gene Position") + ylab("")

######## THIS PLOT DISPLAYS THE CODING SEQUENCES ON TOP OF EACHTOHER TO SHOW DIFFERENT REGIONS ###
plot4 = ggplot() +
  geom_segment(data=gene.justCDS, aes(x=START, y=Y, xend=END, yend=Y, color=FEATURE), size=CDSwidth) +
  geom_segment(data=gene.justUTR, aes(x=START, y=1, xend=END, yend=1, color=FEATURE), size=20) +
  geom_segment(data=recombine.v3.df, aes(x=START, y=2, xend=END, yend=2, color=FEATURE), size=20) +
  geom_segment(data=trusight.df, aes(x=START, y=3, xend=END, yend=3, color=FEATURE), size=20) +
  geom_segment(data=nextera.df, aes(x=START, y=4, xend=END, yend=4, color=FEATURE), size=20) + 
  scale_color_manual(values=colors1) +
  plot_theme + 
  scale_y_continuous(breaks=1:5) +
  labs(title = plottitle) + 
  xlab("Gene Position") + ylab("")


############EDIT OUR DATA FRAMES FOR LEGIBLE OUTPUT##########
gene.small=gene[c("CHROM", "START", "END", "FEATURE")]
recombine.v3.df.small=recombine.v3.df[c("CHROM", "START", "END", "FEATURE")]
nextera.df.small=nextera.df[c("CHROM", "START", "END", "FEATURE")]
trusight.df.small=trusight.df[c("CHROM", "START", "END", "FEATURE")]
output.df=rbind(gene.small, recombine.v3.df.small, nextera.df.small, trusight.df.small)

#############SAVE OUR FILES##########
write.table(output.df, file=file4, row.names=FALSE)
ggsave(file = file1, plot=plot2, width=10*numcol2, height=3*numtranscripts2, limitsize=F)
ggsave(file = file2, plot=plot3, width=10, height=7)
ggsave(file = file3, plot=plot4, width=10, height=7)



