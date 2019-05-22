suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,stringr))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,tibble))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,magrittr))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,assertthat))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,data.table))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,tidyverse))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,here))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,rtracklayer))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,GenomicFeatures))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,GenomicRanges))

args <- list(
  gtf='gencode.v30.annotation.cds.gtf',
  STARTCODTRIM=15,
  STOPCODTRIM=5,
  outtrimfile='gencode.v30.annotation.cds.trim_15_5.gtf'
)
args<-commandArgs(trail=T)


i='STARTCODTRIM'
for(i in names(args)) assign(i,args[[i]],envir=.GlobalEnv)

outdir<-dirname(outtrimfile)

#Ingola's reccomendation https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5582988/
#is to just exclude the first 15 and last 5 codons
#Pretty brutal to short orfs tho...

allanno<-here('pipeline/gencode.v30.annotation.gtf')%>%import
allcds<-allanno%>%subset(type=='CDS')%>%split(.,.$protein_id)
stopifnot(length(allcds)>0)
# testcds<-c(
#   allcds%>%subset(strand=='+')%>%split(.$protein_id)%>%.[666],
#   allcds%>%subset(strand=='-')%>%split(.$protein_id)%>%.[666]
# )%>%GRangesList

cdslens<-allcds%>%width%>%sum
isshort <- cdslens <= (3*(STARTCODTRIM+STOPCODTRIM))
isshort%>%table

newext<-paste0('.tooshort_',STARTCODTRIM,'_',STOPCODTRIM,'.gtf')
export(allcds[isshort]%>%unlist,here(outdir,gtf%>%basename%>%str_replace('.gtf$',newext)))

allcds <- allcds[!isshort]
cdslens <- cdslens[!isshort]

#define the trim zones
untrimcds<-GRanges(names(cdslens),IRanges(STARTCODTRIM*3+1,cdslens-(3*STOPCODTRIM)))%>%pmapFromTranscripts(x=.,allcds)%>%unlist%>%subset(hit)
untrimcds$trimmed_fp <- STARTCODTRIM
untrimcds$trimmed_tp <- STOPCODTRIM

trimstarts<-GRanges(names(cdslens),IRanges(1,STARTCODTRIM*3))%>%pmapFromTranscripts(x=.,allcds)%>%unlist%>%subset(hit)
trimstarts$fp_trim <- TRUE

trimends<-GRanges(names(cdslens),IRanges((STARTCODTRIM*3)+1,cdslens))%>%pmapFromTranscripts(x=.,allcds)%>%unlist%>%subset(hit)
trimends$tp_trim <- TRUE

newext<-paste0('.trim_',STARTCODTRIM,'_',STOPCODTRIM,'.gtf')
export(untrimcds,here(outdir,gtf%>%basename%>%str_replace('.gtf$',newext)))

newext<-paste0('.trimmed_fp_',STARTCODTRIM,'.gtf')
export(trimstarts,here(outdir,gtf%>%basename%>%str_replace('.gtf$',newext)))

newext<-paste0('.trimmed_tp_',STOPCODTRIM,'.gtf')
export(trimends,here(outdir,gtf%>%basename%>%str_replace('.gtf$',newext)))


