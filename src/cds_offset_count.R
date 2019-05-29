

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

fp <-function(gr)ifelse(strand(gr)=='-',end(gr),start(gr))
tp <-function(gr)ifelse(strand(gr)=='-',start(gr),end(gr))
strandshift<-function(gr,shift)shift(gr , ifelse(strand(gr)=='-',- shift,shift))

for(fname in lsf.str('package:dplyr')) assign(fname,get(fname,'package:dplyr'))

source(here('/src/R/Rprofile.R'))

argv <- c(
	bam = here('pipeline/star/data/riboseq_E1/riboseq_E1.bam'),
	gtf = here('pipeline/gencode.v30.annotation.gtf'),
	REF = here('pipeline/hg38.fa'),
	offsetsfile = here('cdsmax/riboseq_E1/cdsmax_offsets.tsv'),
	outfolder = 'ribo_cdsmax/riboseq_E1/',
	STARTCODTRIM=15,
	STOPCODTRIM=5
)

argv[] <- commandArgs(trailing=TRUE)


for (nm in names(argv)) assign(nm,argv[[nm]])

STARTCODTRIM%<>%as.numeric
STOPCODTRIM%<>%as.numeric

bam %T>%{stopifnot(file.exists(.))}

#get exons
if(!exists('gtf_gr')) gtf_gr<-gtf%>%import
if(!is('exons','GRanges')) exons <- gtf_gr%>%subset(type=='exon')
if(!is('cds','GRanges')) cds <- gtf_gr%>%subset(type=='CDS')
if(!exists('startcods')) startcods <- gtf_gr%>%subset(type=='start_codon')
if(!exists('genes_gr')) genes_gr <- gtf_gr%>%subset(type=='gene')



#readGAlignments(bamob)%>%subset(njunc(.)>0)%>%head(1)%>%qnarrow(start=8,width=1)%>%as('GRanges')
cds%<>%split(.,.$transcript_id)




#with map reduce - faster?
riboparam<-ScanBamParam(scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE),mapqFilter=50)

bamob<-Rsamtools::BamFile(bam,yieldSize=5e4)

allcdstouse <- cds
allcdstouse%<>%setNames(.,NULL)%>%unlist%>%split(.,.$protein_id)

library(GenomicFiles)


i<-0

#So indeed CDS look unique per-gene
# cds%>%unlist%>%as.data.frame%>%select(gene_id,transcript_id,start)%>%group_by(gene_id,transcript_id,start)%>%tally%>%.$n%>%table





cdssplit<-cds
#iterate along 1000 cds at a time
cdssplit<-cdssplit%>%unlist%>%split(.,floor(as.numeric(factor(.$protein_id))/1000))

riboparam<-ScanBamParam(scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE),mapqFilter=50)
bamob<-Rsamtools::BamFile(bam,yieldSize=250e3)

file=bam
cdstouse<-cds%>%unlist%>%split(.,.$protein_id)

# protidcounts<-reduceByRange(iterate=TRUE,ranges=cdssplit,


save.image('offsets_cdsmax.Rdata')

protidcounts<-reduceByYield(X=bamob,parallel=FALSE,
	# files=bam%>%setNames(.,.),
	YIELD=function(file){reads=readGAlignmentsList(file,param=riboparam);reads},
	# MAP=function(cdstouse,file){
	MAP=function(reads){

		# riboparam<-ScanBamParam(scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE),mapqFilter=50,which=unlist(cdstouse))
		# reads=readGAlignmentsList(file,param=riboparam)
		# i<<-i+1
		# cdstouse%<>%split(.,.$protein_id)
		cat('.. ')


		jreads<-reads[unlist(njunc(reads)>0)]
		
		reads<-reads[njunc(reads)==0]%>%unlist
		
		
		cdstouse<-subsetByOverlaps(allcdstouse,reads)
		reads<-subsetByOverlaps(reads,allcdstouse)

		#first carry out mappign of the spliced reads
		#This is a bit convoluted - the more obvious approach with pmaptotrancripts
		#was somehow creating a vector too big for memory, so we chunk the junction reads, unlist
		#them intot heir start sites, map to transcripts, then resize them once mapped
		jmapped <- jreads%>%split(.,ceiling(seq_along(.)/ 5e3))%>%as.list%>%lapply(function(jreads){

			ov <- findOverlaps(jreads,cdstouse)
			ovenc<-encodeOverlaps(as(jreads[queryHits(ov)],"GRangesList"),cdstouse[subjectHits(ov)])
			ov<-ov[isCompatibleWithSplicing(ovenc)]
			jreads <- jreads[queryHits(ov)]
			qwidths<-qwidth(jreads)
			

			jmapped <- jreads %>%narrow(start=1,width=1) %>% unlist %>% as("GRanges")%>%pmapToTranscripts(cdstouse[subjectHits(ov)])
			end(jmapped) <- end(jmapped)-1+unlist(qwidths)
			jmapped
		})%>%Reduce(f=c)
		
		#now carry out mapping of unsliced reads
		mapped <- suppressWarnings({mapToTranscripts(as(reads,'GRanges'),cdstouse)})

		#and combine
		mapped <- append(mapped,jmapped)

		cdstousechrs<-seqnames(cdstouse)%>%runValue
		chrs<-unlist(cdstousechrs[match(seqnames(mapped),names(cdstouse))])%>%as.vector


		suppressWarnings({
			mapped$offset<-tibble(
			compartment=compartments[chrs],
			length=width(mapped),
			phase=start(mapped)%>%`%%`(3)
			)%>%left_join(bestscores,by=c('compartment','length','phase'))%>%.$offset
		})

		mapped%<>%subset(!is.na(offset))

		#turn the mapped reads into psites
		suppressWarnings({mapped%<>%resize(1,'start')%>%shift(mapped$offset)})

		#mapped[end(mapped)>seqlengths(mapped)[as.vector(seqnames(mapped))]]

		#use count only psites in designated area
		counts<-mapped%>%
			subset(start>(STARTCODTRIM*3))%>%
			subset(end < seqlengths(.)[as.vector(seqnames(.))] - (STARTCODTRIM*3))%>%
			seqnames%>%
			table
		#make it a vector not a table
		counts%<>%as.vector%>%setNames(names(counts))
		#mark things which are so short the trimming killid their counts
		tooshort <- seqlengths(mapped) < (((STARTCODTRIM*3))+(STOPCODTRIM*3))
		counts[tooshort] <- NA
		#output
		return(counts)
		
	})


protidcountsbak <- protidcounts

protidcounts %<>% lapply(enframe,'protein_id','count')%>%bind_rows%>%group_by(protein_id)%>%summarise(count=sum(count))

prot2gene_df<-as.data.frame(mcols(unlist(cds))[c('protein_id','gene_id','transcript_id')])%>%distinct

genecounts <- protidcounts%>%left_join(prot2gene_df,by='protein_id')%>%group_by(gene_id)%>%slice(which.max(count))

protidcounts%>%write_tsv(file.path(outfolder,'protein_id_counts.tsv'))
genecounts%>%write_tsv(file.path(outfolder,'gene_id_counts.tsv'))

# save.image('offsets_cdsmax.Rdata')
# load('offsets_cdsmax.Rdata')