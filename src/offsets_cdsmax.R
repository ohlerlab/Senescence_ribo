# suppressMessages({library(svglite)})
suppressMessages({library(readr)})
suppressMessages({library(Biostrings)})
suppressMessages({library(Rsamtools)})
#suppressMessages({library(psd)})
suppressMessages({library(txtplot)})
suppressMessages({library(rtracklayer)})
suppressMessages({library(stringr)})
suppressMessages({library(data.table)})
suppressMessages({library(assertthat)})
suppressMessages({library(parallel)})
suppressMessages({library(dplyr)})
# suppressMessages({library(riboWaltz)})
suppressMessages({library(purrr)})
suppressMessages({library(here)})
suppressMessages({library(magrittr)})
suppressMessages({library(stringr)})
suppressMessages({library(tidyverse)})
suppressMessages({library(GenomicAlignments)})
suppressMessages({library(GenomicFeatures)})
suppressMessages({library(GenomicFiles)})
reduce <- GenomicRanges::reduce

MAPQTHRESH <- 50

mycache=memoise::cache_filesystem(here::here('R_cache'))

mymemoise <- function(f) if(!is.memoised(f)) memoise(f,cache=mycache) else f

fp <-function(gr)ifelse(strand(gr)=='-',end(gr),start(gr))
tp <-function(gr)ifelse(strand(gr)=='-',start(gr),end(gr))
strandshift<-function(gr,shift)shift(gr , ifelse(strand(gr)=='-',- shift,shift))

for(fname in lsf.str('package:dplyr')) assign(fname,get(fname,'package:dplyr'))

source(here('/src/R/Rprofile.R'))

argv <- c(
	bam = here('pipeline/star/data/riboseq_E1/riboseq_E1.bam'),
	gtf = here('pipeline/gencode.v30.annotation.gtf'),
	REF = here('pipeline/hg38.fa'),
	outfolder = 'ribo_cdsmax/riboseq_E1/'
)

argv[] <- commandArgs(trailing=TRUE)


for (nm in names(argv)) assign(nm,argv[[nm]])


bam %T>%{stopifnot(file.exists(.))}

#get exons
if(!exists('gtf_gr')) gtf_gr<-gtf%>%import
if(!is('exons','GRanges')) exons <- gtf_gr%>%subset(type=='exon')
if(!is('cds','GRanges')) cds <- gtf_gr%>%subset(type=='CDS')
if(!exists('startcods')) startcods <- gtf_gr%>%subset(type=='start_codon')
if(!exists('genes_gr')) genes_gr <- gtf_gr%>%subset(type=='gene')

# save.image('offsets_cdsmax.Rdata')
# stop()

#get counts over all cds
cdsstren<-mymemoise(bamCount)(bam,cds)
cds$count<-cdsstren
trcounts<-cds$count%>%split(cds$transcript_id)%>%map_dbl(sum)

compartments <- rep('nucl',length(seqlevels(cds)))%>%setNames(seqlevels(cds))
circs_in_data <- intersect(DEFAULT_CIRC_SEQS,names(compartments))
compartments[circs_in_data] <- circs_in_data

#gene tr relationshiop df
gtrdf<-exons%>%mcols%>%.[,c('gene_id','transcript_id')]%>%as.data.frame%>%distinct


startcodcount <- startcods%>%.$transcript_id%>%table

simpletrs<-startcodcount[startcodcount==1]%>%enframe%>%select(transcript_id=name)
ccds <- cds%>%mcols%>%as.data.frame%>%select(transcript_id,tag)%>%filter(tag=='CCDS')%>%distinct
stopifnot(nrow(ccds)>1e3)
#get the tr with the highest signal per gene
toptrs <- gtrdf%>%
	left_join(enframe(trcounts,'transcript_id','count'))%>%
	semi_join(simpletrs)%>%
	semi_join(ccds)%>%
	group_by(gene_id)%>%
	slice(which.max(count))%>%
	arrange(desc(count))%>%
	.$transcript_id%>%
	head(1e3)


#now take the top trs with only 1 start codon
topcds <- cds%>%
	subset(transcript_id %in% toptrs)%>%
	identity
	# head(1e3)

topstartcods<-startcods[match(topcds%>%.$transcript_id%>%unique,startcods$transcript_id)]

#get the exons for these
topcdsexons <- exons%>%subset(transcript_id %in% toptrs)

#mapped cds
topcdsmap<-topcds%>%pmapToTranscripts(topcdsexons%>%split(.,.$transcript_id)%>%.[topcds$transcript_id])%>%reduce
#also start codons
starcods_trmap<-topstartcods%>%mapToTranscripts(topcdsexons%>%split(.$transcript_id))
starcods_trmap%<>%setNames(as.character(seqnames(.)))


library(GenomicAlignments)
#get reads over them
topcdsreadsbak <- bam%>%BamFile(yieldSize=NA)%>%readGAlignments(param=ScanBamParam(which=topcds))
topcdsreads <- topcdsreadsbak

#reads as a gr
topcdsreadsgr<-topcdsreads%>%as("GRanges")
#get the shifts and store in hte gr
topcdsreadsgr$length<-qwidth(topcdsreads)
#now map these to transcript space

cdsread_trmap<-topcdsreadsgr%>%split(.,ceiling(seq_along(.)/ 5e3))%>%lapply(.%>%mapToTranscripts(topcdsexons%>%split(.$transcript_id)))%>%Reduce(f=c)

# testwidde<-cdsread_trmap%>%subset(width>40)%>%.[1])
cdsread_trmap$length<-topcdsreadsgr$length[cdsread_trmap$xHits]
#select only those which mapped cleanly onto transcripts
cdsread_trmap%<>%trim
cdsread_trmap<-cdsread_trmap%>%subset(width==length)

cdsread_trmap$phase <- ((start(cdsread_trmap) - start(starcods_trmap[seqnames(cdsread_trmap)])) %%3)

non3ttr<-coverage(cdsread_trmap)[topcdsmap]%>%lengths%>%map_dbl(`%%`,3)%>%.[.!=0]%>%names
stopifnot(length(non3ttr)==0)#all the cds should be multiples of 3

#assess entropy of the phase - should 
get_frame_entropy<-function(gr,topstartcodphases){
	stopifnot(all(seqnames(gr) %in% names(topstartcodphases)))
	adjstart <- start(gr) - topstartcodphases[as.character(gr@seqnames)]
	ptable<-adjstart %>% `%%`(3)%>%table
	list(ptable%>%{./sum(.)}%>%{-sum(.*log(.))},ptable/sum(ptable))
}


readsizes <- cdsread_trmap$length%>%table%>%{./sum(.)}%>%keep(. > 0.05)%>%names%>%as.numeric%>%setNames(.,.)
overalloffsets<-seq(6,max(readsizes),by=3)%>%setNames(.,.)

#add compartment to the mapped reads
trcomps<-data_frame(chr=as.character(topcdsexons@seqnames),transcript_id=topcdsexons$transcript_id)%>%mutate(compartment=compartments[chr])%>%
	distinct(transcript_id,compartment)%>%{setNames(CharacterList(as.list(.$compartment)),.$transcript_id)}
cdsread_trmap$compartment<-trcomps[cdsread_trmap@seqnames]%>%unlist

#cdsread_trmap%>%split(paste(cdsread_trmap$compartment,cdsread_trmap$length,cdsread_trmap$phase,sep=';'))

#Get the best offset for our score and 
offset_cds_scores<-lapply(overalloffsets,function(offset){
	lapply(unique(compartments),function(compartment_i){
		lapply(readsizes,function(length){
			lapply(0:2%>%setNames(.,.),function(phase_i){
				if((length - 6 - offset) < 0) return(data_frame(score=NA))
				cat('.')
				cdsread_trmap%>%
					subset(compartment==compartment_i)%>%
					subset(phase==phase_i)%>%
					subset(length==length)%>%
					resize(1,'start',ignore.strand=T)%>%
					shift(offset)%>%
					countOverlaps(topcdsmap)%>%
					`>`(0)%>%sum%>%data_frame(score=.)
			})%>%bind_rows(.id='phase')
		})%>%bind_rows(.id='length')
	})%>%bind_rows(.id='compartment')
})%>%bind_rows(.id='offset')

offset_cds_scores%<>%mutate_at(vars(everything()),as.numeric)
offset_cds_scores$compartment<-unique(compartments)[offset_cds_scores$compartment]

#We now have optimal offsets per readlength/phase
bestscores<-offset_cds_scores%>%group_by(length,phase)%>%slice(which.max(score))

#
outfolder%>%dir.create(showWarnings=F,rec=T)
bestscores%>%write_tsv(file.path(outfolder,'cdsmax_offsets.tsv'))


#Now also get counts for ALL the TRIMMED CDS, using these offsets
# trimmed_cds <- import('gencode.v30.annotation.trim_15_5.gtf')
# trimmed_cds <- import('gencode.v30.annotation.trim_15_5.gtf')

# trimmed_cds%<>%split(.,.$protein_id)




