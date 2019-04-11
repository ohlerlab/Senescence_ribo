# the qc step at the moment assumes that the data is single ended
#
shell.executable("bash")
shell.prefix("set -e  pipefail;")

def is_nonempty(file):
  assert os.stat(file).st_size
def is_over_size(file,n):
  assert os.stat(file).st_size > n

# user set parameter
TMPDIR = '$TMPDIR'
SCRIPTDIR = '../git/rna_seq/scripts'

# #reference genome
REF_orig = '../../genomes/hg19.fa'


#things this needs - all sorts of shit for the R scripts....
#samtools, bed tools, a bunhc of ucsc tools, picard, 

# import os
# def filebase(file): return()


# # used by star
# STARINDEX = "/fast/projects/cubit/0.12.0/static_data/precomputed/STAR/2.4.1d/GENCODE/M12/GRCm38/p5.chr_scaff/50/"
# # used by 
# GTF_orig = '/fast/projects/cubit/0.12.0/static_data/annotation/GENCODE/M12/GRCm38/gencode.vM12.annotation.gtf'
# GTF_cdsfilt = 'static_local/gencode.vM12.annotation_cdsfilt.gtf'
# CDSGTF = 'static_local/gencode.vM12.annotation.cds.gtf'

# # used by infer_experiment
# BED = 'static_local/gencode.vM12.annotation.bed'
# GFF = '/fast/projects/cubit/0.12.0/static_data/annotation/GENCODE/M12/GRCm38/gencode.vM12.annotation.gff3'

#for f in $(echo star/data/*/!(*star_transc*).bam ) ; do samtools view $f  | head -n 300 | cut -d$'\t' -f 6 | awk -v f="$f" '{print $0,f}'  ; done  > cigar_sample.txt

REF = 'my_'+os.path.splitext(os.path.split(REF_orig)[1])[0]+'.fa'

# used by infer_experiment
GTF_orig = '../annotation/gencode.v24lift37.annotation.gtf'
GFF_orig = '../annotation/gencode.v24lift37.annotation.gff3'

ANNOBASE = 'my_'+os.path.splitext(os.path.split(GTF_orig)[1])[0]

GFF = ANNOBASE+'.gff3'

BED = ANNOBASE+'.bed'

# used by 
GTF = ANNOBASE+'.gtf'
CDSGTF = ANNOBASE+'.cdsfilt.gtf'
#used to make indices
RNAFASTA = ANNOBASE+'.transcript.fa'
CDSFASTA = ANNOBASE+'.cds.fa'

# need to think on this.... we start with a genome and a gff3 file. We add the extra transcript to our gff3 file, (optionally eliminating all others) 
#then we create the transcript fasta file from our gff3 file, and the gtf file. Then, we use these to make salmon, and rsem/star indices

# used by qc
SeQC_GTF = ANNOBASE+'_SEQC_GTF'+'.gtf'
SeQC_REF = REF

SAMPLELINES = [line.strip().split(',') for line in open("sample_parameter.csv").readlines()]

#switch for testmode
if(config.get('test',0)): 
  print('\n\n--------testmode on -------------\n\n')
  SAMPLELINES = SAMPLELINES[0:2]  
  origSAMPLE = [ entry[SAMPLELINES[0].index('sample_id')] for entry in SAMPLELINES[1:]]
  SAMPLES = ['test']

else:
  SAMPLES = [ entry[SAMPLELINES[0].index('sample_id')] for entry in SAMPLELINES[1:]]



#get our info as dictionaries
def tab_to_dict(SAMPLELINES,valcol):
  valind = SAMPLELINES[0].index(valcol)
  vals   = [ entry[valind] for entry in SAMPLELINES[1:]]
  return dict(zip(SAMPLES,vals))

#sample - info dictionaries
LIBRARY_DICT          = tab_to_dict(SAMPLELINES,'library_layout')
READ_PATTERN_DICT     = tab_to_dict(SAMPLELINES,'read_pattern')
PROTOCOL_DICT         = tab_to_dict(SAMPLELINES, 'protocol')
FRAG_LENGTH_MEAN_DICT = tab_to_dict(SAMPLELINES, 'fragment_length_mean')
FRAG_LENGTH_SD_DICT   = tab_to_dict(SAMPLELINES, 'fragment_length_sd')

ASSAY_DICT            = tab_to_dict(SAMPLELINES, 'assay')
GTF_DICT              = {k: CDSGTF if 'ribo' in v else GTF for k, v in ASSAY_DICT.items()}

#for f in $(echo input/*); do for q in $( echo ${f}/* ); do echo $f $q; done; done | sed 's/input\///' > pipeline/sample_file.txt
SAMPLEFASTQLINES = [line.strip().split('\t') for line in open("sample_file.txt").readlines()]
FASTQS = [l[1] for l in SAMPLEFASTQLINES]
FASTQSAMPLES = [l[0] for l in SAMPLEFASTQLINES]
FASTQSAMPLEDICT = dict(zip(FASTQS,FASTQSAMPLES))
SAMPLEFASTQDICT = {v:[i for i in FASTQSAMPLEDICT.keys() if FASTQSAMPLEDICT[i] == v ] for k,v in FASTQSAMPLEDICT.items()}


assert set(FASTQSAMPLES) == set(SAMPLES)

#the group dict is structured differently, returns a list of samples
GROUP_DICT = tab_to_dict(SAMPLELINES,'group')
GROUP_SAMPLES = {}

for k, v in GROUP_DICT.items():
    GROUP_SAMPLES[v] = GROUP_SAMPLES.get(v, [])
    GROUP_SAMPLES[v].append(k)

GROUPS = list(GROUP_SAMPLES.keys())

#information on the strand of stuff
strands = ['pos','neg']
STRANDSYMS={strands[0]:'+',strands[1]:'-'}

#extensions for transcript and chromosome bigwigs
istransvals = ['.transcript','.chr']
#extensions used by STAR to denot the transcript/genomic bam
BEXTS={istransvals[0]:'.star_transcript',istransvals[1]:''}


# print('Samples are: ',SAMPLES)
# print('Groups are: ',GROUPS)


RIBO_TOTAL_DICT = dict(zip(
  list(filter(lambda x: 'ribo' in x,SAMPLES)),
  list(filter(lambda x: 'total' in x,SAMPLES))
))

GENEREGIONS = ['gene','cds','fputrs','tputrs']
# generegions = ['gene','cds','fputrs','tputrs','cds_tiles','fputr_tiles','tputr_tiles']

# TRNAs = ['gencode.vM12.tRNAs.gtf.gz']
TRNAs = ['tRNAs']

# READRANGES = ['25_30','1_26','27_28','29_100','1_300']
READRANGES = ['25_31','1_300']
# READRANGENUM = [[25,30],[1,26],[27,28],[29,100],[1,300]]
READRANGENUM = [[25,31],[1,300]]
READRANGEDICT = dict(zip(READRANGES,READRANGENUM))


assert set(FASTQSAMPLES) == set(SAMPLES)


ribosamples=list(filter(lambda s: ASSAY_DICT[s] == 'ribo',SAMPLES))
totalsamples=list(filter(lambda s: ASSAY_DICT[s] == 'total',SAMPLES))
riboms_total_file='/fast/work/groups/cubi/projects/2017-10-10_cortexomics/gdrive/cortexomics_ms_total/325_new_2_all_log2_LFQ_n7464.txt'
ms_spec_file='/fast/work/groups/cubi/projects/2017-10-10_cortexomics/gdrive/cortexomics_ms_cyt+80+Poly/proteinGroups.txt'

ribofastqs=list(filter(lambda fq: FASTQSAMPLEDICT[fq] in ribosamples,FASTQS))
totalfastqs=list(filter(lambda fq: FASTQSAMPLEDICT[fq] == totalsamples,FASTQS))


groupnames = ['OD5P','ONVC','OMM']
samplegroups = {
  'OD5P': list(filter(lambda s: 'OD5P' in s,ribosamples)),
  'ONVC': list(filter(lambda s: 'ONVC' in s,ribosamples)),
  'OMM': list(filter(lambda s: 'OMM' in s,ribosamples)),
  'T1014A': list(filter(lambda s: 'T1014A' in s,ribosamples)),
  'T1015A': list(filter(lambda s: 'T1015A' in s,ribosamples)),
  'T1185B': list(filter(lambda s: 'T1185B' in s,ribosamples)),
}

samplegroups.update(dict(zip(SAMPLES,[[s] for s in SAMPLES] )))
# print(samplegroups)


# print(SAMPLES)

rule all:
  input:
    FASTQS,
    expand("processed_reads/{sample}/.done", sample = SAMPLES),
    # expand("cutsequences/{sample}/cutseqs.txt.gz", sample = SAMPLES),
    # GTF,
    # CDSFASTA,
    #expand("rsem/data/{sample}/.done", sample = SAMPLES),
    # expand("fastqc/data/{sample}/.done", sample = SAMPLES),
    #    "fastqc/summary/fastqc_summary.tsv",
    expand("star/data/{sample}/.done", sample = SAMPLES),
    #expand("infer_experiment/data/{sample}/.done", sample = SAMPLES),
    expand("qc/data/{sample}/.done", sample = SAMPLES),
    ("multiqc/multiqc_report.html"),
    # expand("dupradar/data/{sample}/.done", sample = SAMPLES),
    # expand("htseq/data/{sample}/.done", sample = SAMPLES),
    expand('feature_counts_readrange/data/{sample}/{generegions}/{readrange}/.done', sample=RIBO_TOTAL_DICT.keys(), generegions=GENEREGIONS+TRNAs, readrange=READRANGES),
    expand("feature_counts/data/{sample}/feature_counts", sample = SAMPLES),
    expand("feature_counts/all_feature_counts"),
    # expand("feature_counts/data/{sample}/{generegions}/{readrange}/.done", sample = SAMPLES, generegions = generegions+TRNAs, readrange = READRANGES),
    # expand("kallisto/data/{sample}/.done", sample = SAMPLES),
    # expand("bigwigs/{group}/{strand}/{istrans}.done",group = SAMPLES,strand=strands,istrans=istransvals),
    # expand("mergedbigwigs/{group}/{strand}/{istrans}.done",group = GROUPS,strand=strands,istrans=istransvals),
    # expand("ribotaper/{sample}/.done",sample=list(RIBO_TOTAL_DICT.keys())),
    # expand('junctioncounts/{sample}/{sample}.junctioncounts.tsv',sample=ribosamples),
    expand("SaTAnn/{sample}/.done", sample = ribosamples+groupnames),
    expand('riboqc/reports/{sample}/riboqcreport.html', sample = ribosamples+groupnames),
    expand('groupedsatan/{group}.fasta', group = groupnames),
    expand('groupedsatan/{group}.variantmod.fasta', group = groupnames),

# expand("ribotapermetaplots/{sample}/.done",sample=list(RIBO_TOTAL_DICT.keys())),
    

MINREADLENGTH=20
MAXREADLENGTH=300
QUALLIM=20
CUTADAPTBIN="~/work/bin/cutadapt"
REMOVE8NBIN="~/work/bin/remove8N_twoarg.pl"


###OOOOOOkay. GFF read does some weird shit including excluding 'transcript' lines.....
rule gffread:
  input: REF,GTF_orig
  output: GTF,CDSGTF,RNAFASTA,CDSFASTA,BED
  run:
    shell(r""" 
      # set -x
      #with filtering output all sequences
      cat {GTF_orig} \
      | sed -r  's/((transcript_id|gene_id|protein_id|ID|Parent|exon_id|havana_gene|havana_transcript)\W+\w+)\.[0-9]+/\1/g' \
      > {GTF}

      cat {GFF_orig} \
      | sed -r  's/((transcript_id|gene_id|protein_id|ID|Parent|exon_id|havana_gene|havana_transcript)\W+\w+)\.[0-9]+/\1/g' \
      > {GFF}


      #needs gff - output exon sequences
      cat {GFF_orig} |  grep -P -e'\texon\t|^\#' | gffread - -F -E -g {REF} -W -w {RNAFASTA} -o /dev/null

      #Note we are now minus the transcript and exon entries for these
      #now make GTF

      #| grep -P -e'\tCDS\t|^\#' 
     #with filtering, output the coding sequences filteirng out the ones that aren't in frame, have a stop codon, are pseudogenes etc.
      
      cat {GFF_orig}  \
        | sed -r  's/((transcript_id|gene_id|protein_id|ID=\w+|Parent)\W+\w+)\.[0-9]+/\1/g' \
        | gffread - -C -V -J --no-pseudo  -F -E -g {REF} \
        -W -w {ANNOBASE}.coding.transcript.fa -x {CDSFASTA} -y {ANNOBASE}.protein.fa -T \
        -o /dev/stdout \
        | awk -v FS="\t" -v OFS="\t" '{{if($3=="CDS"){{$3="exon";print $0}}}}' \
         > {CDSGTF}

      #now make bed
      cat {GTF_orig} | awk '{{print $1,$4,$5,"name",$6,$7}}' > {BED}
      """)

 
rule make_utrs:
  input: GTF=GTF_orig
  output: fputrs='fputrs.gtf',tputrs='tputrs.gtf'
  # script: 'make_utrfiles.R'
  run:
    shell(r"""
      set -ex
      #with filtering output all sequences
      cat {input.GTF}  \
      | awk -v OFS="\t"  '{{if($3=="five_prime_UTR"){{         ;print $0}}}}' \
      | sed -r  's/((transcript_id|gene_id|protein_id|ID=\w+|Parent)\W+\w+)\.[0-9]+/\1/g' \
      > {output.fputrs} 

      cat {input.GTF} \
      | awk -v OFS="\t"  '{{if($3=="three_prime_UTR"){{         ;print $0}}}}' \
      | sed -r  's/((transcript_id|gene_id|protein_id|ID=\w+|Parent)\W+\w+)\.[0-9]+/\1/g' \
      > {output.tputrs}

     
      """) 

rule fastqc:
     input: 'processed_reads/{sample}/.done'
     output: touch('fastqc/data/{sample}/.done')
     threads: 4
     log:'fastqc/reports/{sample}/fastqc.log'
     params:
      reads = lambda wc: [fq.replace('input/','processed_reads/') for fq in SAMPLEFASTQDICT[wc['sample']]],
      outdir = lambda wc: 'fastqc/data/'+wc.sample+'/'
     shell: '''
          OUTDIR=$(dirname {output[0]})
          mkdir -p {params.outdir}
          wait $(for i in {params.reads}; do $( fastqc -o {params.outdir} $i ) & done) 
        '''

rule collect_fastqc:
     input:
          all_results = expand("fastqc/data/{sample}/.done", sample=SAMPLES)
     output:
          result='fastqc/summary/fastqc_summary.tsv',
          log='fastqc/summary/fastqc_summary.log'
     shell:
          r"""
          set -e
          mkdir -p $(dirname {output.result}) 
          {SCRIPTDIR}/collect_fastqc_results.sh -i fastqc/ \
          > {output.result} \
          2> {output.log} 
          """


    
rule star_index:
  input: alljuncfile='junctions/all.tsv',REF=REF
  output: touch('starindex/.done')
  threads: 8
  run:
    shell(r"""
      STAR \
      --runThreadN {threads} \
      --runMode genomeGenerate \
      --genomeDir $(dirname {output}) \
      --sjdbFileChrStartEnd {input.alljuncfile} \
      --genomeFastaFiles {input.REF}
      """)  


def get_fastqops(inputdir,read_pattern,lstring='<( zcat ',rstring=')'):
  import glob as glob
  #get our input files, in either paired end or single end form
  assert '-f' in read_pattern
  fastqops = read_pattern.split('-')[1].replace('f ','',1).strip()
  fastqops = glob.glob(inputdir+'/'+fastqops)
  fastqops.sort()
  assert fastqops
  assert all([os.stat(fastq).st_size!=0 for fastq in fastqops])
  fastqops = ' '.join(fastqops)
  
  fastqops = lstring+fastqops+')'

  if '-q' in read_pattern:
    fastqops2 = read_pattern.split('-')[2].replace('q ','',1).strip()
    fastqops2 = glob.glob(inputdir+'/'+fastqops2)
    assert all([os.stat(fastq).st_size!=0 for fastq in fastqops2])
    fastqops2.sort()
    fastqops2 = ' '.join(fastqops2)
    fastqops2 = lstring+fastqops2+')'
    fastqops += ' ' + fastqops2
  return(fastqops)

rule star:
     input:
          fastqs='processed_reads/{sample}/.done',
          STARINDEX='starindex/.done',
          bowtie_index='bowtie_index/.done',
     output:
          done = touch('star/data/{sample,[^/]+}/.done')
     threads: 8
     run:
          input.STARINDEX=input.STARINDEX.replace('.done','')
          markdup = '' if ASSAY_DICT[wildcards['sample']] == 'ribo' else '-m'
          platform = 'NotSpecified'
          inputdir = os.path.dirname(input['fastqs'])
          outputdir = os.path.dirname(output[0])
          read_pattern = READ_PATTERN_DICT[wildcards['sample']]
          fastqops = get_fastqops(inputdir,read_pattern,lstring='<( zcat ',rstring=')')
          repdir = outputdir.replace('data','reports')
          tophatindex =input['bowtie_index'].replace('.done','')
          
          halfthreads = threads/2
          sortmem = str(int(5000/halfthreads))+'M'

          remap = '1' if ASSAY_DICT[wildcards['sample']] == 'ribo' else ''

          sample = wildcards['sample']
          shell(r"""
            set -x
         MY_TMP_DIR=$(mktemp -d)
        trap "set -x; rm -rf ${{MY_TMP_DIR}}" EXIT KILL TERM INT HUP

         mkdir -p $MY_TMP_DIR
        mkdir -p $MY_TMP_DIR/star
        mkdir -p $MY_TMP_DIR/tophat2

        #--outSAMmultNmax 20 --winAnchorMultimapNmax 50 --outFilterMultimapNmax 20 \

        STAR \
              --genomeDir {input.STARINDEX} \
              --runThreadN {threads} \
              --outSAMunmapped Within \
              --outFilterType BySJout \
              --outMultimapperOrder Random \
              --alignSJoverhangMin 8 \
              --alignSJDBoverhangMin 1 \
              --outFilterMismatchNmax 999 \
              --outFilterMismatchNoverLmax 0.04 \
              --alignIntronMin 20 \
              --alignIntronMax 1000000 \
              --alignMatesGapMax 1000000 \
              --genomeLoad NoSharedMemory \
              --quantMode GeneCounts \
              --outSAMattributes NH HI AS NM MD \
              --outSAMtype BAM  Unsorted\
              --outSAMattrRGline \"ID:{sample}\" \"SM:{sample}\" \"PL:{platform}\" \
              --outFileNamePrefix ${{MY_TMP_DIR}}/star/ \
              --outReadsUnmapped Fastx \
              --readFilesIn {fastqops}
          

          iftophat=
          if [ {remap} -eq 1 ] && [ ${{MY_TMP_DIR}}/star/Unmapped.out.mate* ]; then
            iftophat=true

            tophat2 \
                -p {threads} \
                -z0 \
                -g 100 \
                --output-dir ${{MY_TMP_DIR}}/tophat2 \
                --library-type fr-unstranded \
                --no-coverage-search \
                --transcriptome-index {tophatindex}/transcriptome \
                {tophatindex}/ref \
                ${{MY_TMP_DIR}}/star/Unmapped.out.mate*

            umapped=${{MY_TMP_DIR}}/tophat2/unmapped.bam
            tmapped=
            
            samtools merge \
              -@ {threads}  -f \
             ${{MY_TMP_DIR}}/all.bam \
             ${{MY_TMP_DIR}}/star/Aligned.out.bam \
             ${{MY_TMP_DIR}}/tophat2/*.bam
          else
            cp ${{MY_TMP_DIR}}/star/Aligned.out.bam ${{MY_TMP_DIR}}/all.bam
          fi
          
         samtools sort \
          -@ {halfthreads}\
          -m {sortmem} \
          -T ${{MY_TMP_DIR}} \
          -o {outputdir}/{sample}.bam \
          ${{MY_TMP_DIR}}/all.bam
      
        samtools index {outputdir}/{sample}.bam 

        mkdir -p {repdir}
        samtools stats {outputdir}/{sample}.bam > {repdir}/{sample}.bamstats.txt
        samtools flagstat {outputdir}/{sample}.bam > {repdir}/{sample}.flagstat.log
        samtools idxstats {outputdir}/{sample}.bam > {repdir}/{sample}.idxstats.log
        
        cp  ${{MY_TMP_DIR}}/star/ReadsPerGene.out.tab {outputdir}/ReadsPerGene.out.tab
        cp  ${{MY_TMP_DIR}}/star/SJ.out.tab {outputdir}/
        cp  ${{MY_TMP_DIR}}/star/{{Log.final.out,Log.out}} {repdir}/
        if [ $iftophat ] ;then cp ${{MY_TMP_DIR}}/tophat2/align_summary.txt {repdir} ;fi

          """)
          
rule junctioncounts:
  input: 
    uniquejuncrds = 'junctions/uniqueintrons.rds',
    star = 'star/data/{sample}/.done',
  output:
    'junctioncounts/{sample}/{sample}.junctioncounts.tsv'
  run:
    bamfile = 'star/data/'+wildcards['sample']+'/'+wildcards['sample']+'.bam'  
    shell(r'''
      Rscript ../src/junctioncounts.R {input[0]} {bamfile} {output}
        ''')
         


# TESTSECTION="1:1-1260071"

# rule preprocess_bam:
#     input: get_bam
#     output: "work/{sample}/{sample}.preprocess.bam"
#     shell: r"""
#     samtools view -bh {input} {TESTSECTION} > {output}
#     samtools index {output}
#     """


# rule bwa_mem:
#      input:
#           'processed_reads/{sample}/{sample}.fastq.gz'
#      output:
#           done = touch('star/data/{sample}/.done')
#      threads: 8
#      run:
#           read_pattern = READ_PATTERN_DICT[wildcards['sample']]
#           shell(r"""
#           set -e
#           mkdir -p star/reports/{wildcards.sample}
#           {SCRIPTDIR}/run_my_star.sh -m -@ {threads} -i {input} {read_pattern} -o $(dirname {output}) \
#           -r star/reports/{wildcards.sample} -t {TMPDIR} -l NoSharedMemory -x {STARINDEX} 
#           """)

def get_trackfile(sample,strand,ext):
  return 'bigwigs/'+sample+'/'+strand+'/'+sample+'.'+strand+ext

#TODO strand is a problem here
rule bigwigs:
     input: 
          "star/data/{sample}/.done"
     output:
          done = touch("bigwigs/{sample}/{strand}/{istrans}.done")
     threads: 1
     run:
        bedgraph = get_trackfile(wildcards.sample,wildcards.strand,wildcards.istrans+'.bg')
        bw = get_trackfile(wildcards.sample,wildcards.strand,wildcards.istrans+'.bw')
        bext = BEXTS[wildcards.istrans]
        bam = 'star/data/'+wildcards.sample+'/'+wildcards.sample+bext+'.bam'        
        chromsizes = 'bigwigs/'+wildcards.sample+'/'+wildcards.istrans+'.chrsizes.txt'
        strandsym = STRANDSYMS[wildcards.strand]
        #turns out reversing the y axis breaks bigwig merge....
        # revyaxis = -1 if strandsym is '-' else 1
        revyaxis = 1

        shell(r"""
        set -e
        mkdir -p bigwigs/{wildcards.sample}/{wildcards.strand}
       
        # count="$(cat star/reports/{wildcards.sample}/{wildcards.sample}.bam.bamstats.txt | \
        #  grep -Po "reads mapped:\s+\d+" | \
        #  grep -Po "\d+$" | awk '{{print {revyaxis}*1000000/$0}}' 
        #  )"
        count=1
        samtools view -H {bam} | grep "SQ" | sed -r 's/^@SQ\s+SN://' | sed -r "s/\s+LN:/\t/" > {chromsizes} 

        bedtools genomecov -5 -ibam {bam} -bg -strand {strandsym} -scale $count > {bedgraph}
        bedSort {bedgraph} {bedgraph}
        bedGraphToBigWig {bedgraph} {chromsizes} {bw}
       

        """)
        is_over_size(output[0],1e6)


def getGroupBigWigs(wildcards):
    strand = wildcards['strand']
    samples = GROUP_SAMPLES[wildcards['group']]
    return [get_trackfile(s,strand,wildcards.istrans+'.bw') for s in samples ]


def getGroupBigWigDone(wildcards):
    strand = wildcards['strand']
    samples = GROUP_SAMPLES[wildcards['group']]
    return ['bigwigs/'+s+'/'+strand+'/'+wildcards.istrans+'.done' for s in samples ]


rule mergedbigwigs:
     input: getGroupBigWigDone
     output:
          done = touch('mergedbigwigs/{group}/{strand}/{istrans}.done')
     threads: 1
     run:
        chromsizes = 'bigwigs/'+GROUP_SAMPLES[wildcards.group][0]+'/'+wildcards.istrans+'.chrsizes.txt' 
        bigwigs = getGroupBigWigs(wildcards)
        libnum = len(bigwigs)
        mergedfilebase="mergedbigwigs/"+wildcards.group+"/"+wildcards.strand+"/"+wildcards.group+"."+wildcards.strand+wildcards.istrans


        shell(r"""
        set -xe
        mkdir -p mergedbigwigs/{wildcards.group}/{wildcards.strand}/
        
        bigWigMerge {bigwigs} /dev/stdout | awk 'BEGIN{{OFS="\t"}}{{$4 = $4 /{libnum} ; print $0 }}' > {mergedfilebase}.bg 
        
        bedSort {mergedfilebase}.bg {mergedfilebase}.bg 

        bedGraphToBigWig \
          {mergedfilebase}.bg \
          {chromsizes} \
          {mergedfilebase}.bw

        rm {mergedfilebase}.bg
        """)

rule infer_experiment:
     input:
          'star/data/{sample}/.done'
     output:
          done = touch('infer_experiment/data/{sample}/.done')
     shell:
          (r"""
          set -e
          mkdir -p infer_experiment/data/{wildcards.sample}
          
          infer_experiment.py -r {BED} \
          -i star/data/{wildcards.sample}/{wildcards.sample}.bam \
          -q 255 > infer_experiment/data/{wildcards.sample}/{wildcards.sample}.strand_stat.txt
          """)
     
transcript_gene_map = 'transcript_gene_map.tsv'
gene_transcript_map = 'gene_transcript_map.tsv'

rule make_gene_transcript_map:
  input: GTF
  output: gene_transcript_map,transcript_gene_map
  run:
    shell(r"""
    cat {input} \
      | grep -Pe'\ttranscript\t'  \
      | perl -lane '/transcript_id\W+([\w\.]+)/;$t=$1; $g=/gene_id\W+([\w\.]+)/;$g=$1;print($g,"\t",$t)' \
      | sort | uniq \
      > {gene_transcript_map}

    cat {gene_transcript_map} \
      | awk '{{print $2,$1}}' > {transcript_gene_map}
    """)
    is_nonempty(gene_transcript_map)


rrna_intervals = 'qc/picard_rrna_intervals.txt'
refflat = 'qc/'+ANNOBASE+'.refflat'

rule make_picard_files:
  input: GTF,'star/data/'+SAMPLES[0]+'/.done'
  output: intervals=rrna_intervals,refflat=refflat
  conda: '../envs/picard'
  shell:r"""
         samtools view -H star/data/{SAMPLES[0]}/{SAMPLES[0]}.bam > {output.intervals}
        
         grep -Pe 'gene_type..rRNA.' {input[0]} \
         | awk '$3 =="transcript"' \
         | cut -f 1,4,5,7,9 \
         | perl -lane ' /transcript_id "([^"]+)"/ or die "notranscript_id on $."; print join "\t", (@F[0,1,2,3], $1) ' \
         | sort -k1V -k2n -k3n  - >> {output.intervals}
        
        gtfToGenePred -geneNameAsName2 {GTF} {GTF}.genepred
        cat {GTF}.genepred | awk -vOFS="\t" '{{print $1,$0}}' > {output.refflat}

  """


rule qc:
     input:
          fastqc='fastqc/data/{sample}/.done',
          star='star/data/{sample}/.done',
          refflat = refflat,
          rrna_intervals = rrna_intervals,
     output:
          done=touch('qc/data/{sample}/.done'),
     conda: '../envs/picard'
     resources:
     params:
        singleendflag = lambda wc: ' -singeEnd ' if LIBRARY_DICT[wc['sample']] == 'PAIRED' else '',
        bamfile = lambda wc:'star/data/'+wc['sample']+'/'+wc['sample']+'.bam' 
    
     shell: """
          set -e
          set -xv
          
        OUTDIR=$(dirname {output.done})
        mkdir -p qc/reports/{wildcards.sample}/

        {SCRIPTDIR}/read_statistic_report.sh \
         -l star/reports/{wildcards.sample}/Log.final.out  \
         -g $(dirname {input.fastqc}) \
         -o ${{OUTDIR}}/read_alignment_report.tsv \
         &> qc/reports/{wildcards.sample}/{wildcards.sample}_qc.log 

         picard CollectRnaSeqMetrics -Xms4G \
          I={params.bamfile} \
          O=${{OUTDIR}}/{wildcards.sample}_picard_qc.txt \
          REF_FLAT={refflat} \
          STRAND=FIRST_READ_TRANSCRIPTION_STRAND \
          RIBOSOMAL_INTERVALS={rrna_intervals}
        
        picard CollectAlignmentSummaryMetrics \
          INPUT={params.bamfile} \
          OUTPUT=${{OUTDIR}}/{wildcards.sample}.picard.alignmentmetrics.txt \
          R={REF}

      {SCRIPTDIR}/read_duplication.sh \
        -i {params.bamfile} \
        -o ${{OUTDIR}}/duplication/ \
        &> qc/reports/{wildcards.sample}/{wildcards.sample}_qc.log 

          """
     


multiqcscript = '~/work/Applications/MultiQC/scripts/multiqc'

rule multiqc:
  input:
      expand("fastqc/data/{sample}/.done", sample = SAMPLES),
      expand("star/data/{sample}/.done", sample = SAMPLES),
      expand("qc/data/{sample}/.done", sample = SAMPLES),
      # expand("tophat2/data/{sample}/.done", sample = SAMPLES),
      [f.replace('input','filter_reads') for f in  ribofastqs],
      expand("feature_counts/data/{sample}/feature_counts", sample = SAMPLES),
      'sample_file.txt'
 
  output:
    'multiqc/multiqc_report.html'
  run:
    reportsdirs = list(input)
    reportsdirs=[s.replace('star/data','star/reports') for s in reportsdirs]
    reportsdirs=[s.replace('tophat2/data','tophat2/reports') for s in reportsdirs]
    reportsdirs=[os.path.dirname(s) for s in list(reportsdirs)]
    shell(r"""
      cat sample_file.txt | sed 's/.fastq.gz//g' | sed 's/\t.*\//\t/g' \
      | awk -vOFS='\t' 'BEGIN{{print "fastqname","samplename"}}{{sumsamp[$1] = sumsamp[$1]+1;print $2,$1"_fq"sumsamp[$1]}}' \
      > multiqc/samplenames.txt

      {multiqcscript} {reportsdirs} -fo $(dirname {output[0]}) -c multiqc_config.yaml --sample-names multiqc/samplenames.txt
      """)

#this is going to count reads in each library over 5'UTRS, CDS, and 3' UTRs


rule readlenfilt:
  input: 'star/data/{sample}/.done'
  output:  'readlenfilt/data/{sample}/{readrange}/readlenfilt.bam'
  run:
    minreadlen,maxreadlen = READRANGEDICT[wildcards['readrange']]

    readrangebam = output[0]
    shell(r"""
          set -ex
          samtools view -h $(dirname {input})/{wildcards.sample}.bam \
           | awk '((length($10) >= {minreadlen})&&(length($10) <= {maxreadlen})) || $1 ~ /^@/' \
          | samtools view -S -b - > {readrangebam}
            
      """)


rule feature_counts_readrange:
     input:
          GTF,CDSGTF,'tputrs.gtf','fputrs.gtf',
          readrangefilt='readlenfilt/data/{sample}/{readrange}/readlenfilt.bam'
     output:
          done = touch('feature_counts_readrange/data/{sample,[^/]+}/{generegions}/{readrange}/feature_counts')
     threads: 2
     log: r"""feature_counts_readrange/reports/{sample}/{generegions}/{readrange}/feature_counts.log"""
     run:
          if (wildcards['generegions'] in ['gene','cds']):
            GTF = GTF_DICT[wildcards['sample']]
            groupcol = 'gene_id'
          else:
            GTF = wildcards['generegions']+'.gtf'
            groupcol = 'transcript_id'
          
          protocol = PROTOCOL_DICT[wildcards['sample']]
          if (protocol == 'no'):
               protocol = 0
          elif (protocol == 'yes'):
               protocol = 1
          elif (protocol == 'reverse'):
               protocol = 2
          else:
               sys.exit('Protocol not known!')

          library = LIBRARY_DICT[wildcards['sample']]

          if (library == 'PAIRED'):
               library = '-p'
          else:
               library = ''

          countmultimappers = ' ' 
          
          if (wildcards['generegions']=='tRNAs'):
            featuretype = 'tRNA'
            countmultimappers = '-M --fraction'
          
          elif  (wildcards['generegions']=='tputrs'):
            featuretype = 'exon'
          elif  (wildcards['generegions']=='fputrs'):
            featuretype = 'exon'
          else:
            featuretype = 'exon'


          sample = wildcards['sample']
          generegions = wildcards['generegions']
          readrange = wildcards['readrange']
          rangebam = input['readrangefilt']
          
          shell(r"""
          set -ex
          mkdir -p feature_counts_readrange/data/{sample}/{generegions}/{readrange}/
          mkdir -p feature_counts/reports/{wildcards.sample}/
          featureCounts \
            -T {threads} \
            -t {featuretype} -g {groupcol} \
            -a {GTF} \
            -s {protocol} {library} {countmultimappers} \
            -o feature_counts_readrange/data/{sample}/{generegions}/{readrange}/feature_counts \
            {rangebam} \
             &> feature_counts/reports/{wildcards.sample}/{wildcards.sample}.feature_counts.log

          """)

def get_readrange(wc):
    if wc['sample'] in ribosamples:
          selregion='cds'
          selreadrange='25_31'
    else:
          selregion='cds'
          selreadrange='1_300'
    fcountfile = 'feature_counts_readrange/data/'+wc['sample']+'/'+selregion+'/'+selreadrange+'/feature_counts'
    return   fcountfile

rule feature_counts:
     input:
          get_readrange,
          'star/data/{sample}/.done',GTF,CDSGTF,'tputrs.gtf','fputrs.gtf',
     output:
          'feature_counts/data/{sample,[^/]+}/feature_counts'
     threads: 2
     run:       
        bamfile = 'star/data/'+wildcards['sample']+'/'+wildcards['sample']+'.bam' 

        shell(r"""

          mkdir -p feature_counts/reports/{wildcards.sample}
          mkdir -p feature_counts/data/{wildcards.sample}

          head -n2 {input[0]} \
          | tail -n1 \
          | awk -v FS="\t" -v OFS="\t"  '{{$7= "{bamfile}";print $0}}' \
          >    feature_counts/data/{wildcards.sample}/feature_counts

          tail -n+3 {input[0]} \
          >>   feature_counts/data/{wildcards.sample}/feature_counts
          """)

# #this is going to count reads in each library over 5'UTRS, CDS, and 3' UTRs
rule aggregate_feature_counts:
  input : expand("feature_counts/data/{sample}/feature_counts", sample = SAMPLES),
  output: 'feature_counts/all_feature_counts'
  run:
    fcountfiles = list(input)
    shell(r""" 

       #( (sed '2q;d' {fcountfiles[0]} | cut -f1 && tail -n+3 {fcountfiles[0]}| sort -k1 | cut -f1) > {output})
       tail -n+3 {fcountfiles[0]}| sort -k1 | cut -f1 > {output}
       
       #now for eahc fcount table, join it to the ids
       for fcountfile in $(echo {fcountfiles}); do

          tail -n+3 $fcountfile| sort -k1 | cut -f1,7 | join {output} - | sort -k1 > {output}tmp
          mv {output}tmp {output}
       
       done

      echo "feature_id {SAMPLES}" | cat - {output} > {output}tmp
      mv {output}tmp {output}
    
      """)

rule aggregate_junctioncounts:
  input : expand("junctioncounts/{sample}/{sample}.junctioncounts.tsv", sample = SAMPLES),
  output: 'junctioncounts/all_junctioncounts.tsv'
  run:
    shell(r""" 
          echo {input}  | {{ read -a args; \
            cp ${{args[0]}} {output};  \
            for oth in ${{args[@]:1}} ;do \
              paste  {output}  <( cut -f2 $oth)  > {output}.tmp; 
              mv {output}.tmp {output} ; 
            done ;
          }} ;
      """)

rule segment_periodicity:
  input:
    juncgtf = 'junctions/{junction}.gtf',
    pluspsites = 'riboqc/data/{sample}/_P_sites_plus.bw',
    negpsites = 'riboqc/data/{sample}/_P_sites_plus.bw',
  output:
    'junctionspec/ref_metadata.filt_REF.C3N-02289.filt_L1.tsv'
  run:
    shell(r"""
    Rscript ../src/segment_periodicity.R
""")
twobitfile=REF.replace('.fa','.twobit')

satan_annot_script =  '/fast/work/groups/ag_ohler/dharnet_m/satann_working/Annot_make_bioc_gtf.R'
riboqc_script =  '/fast/work/groups/ag_ohler/dharnet_m/satann_working/analysis_qc_mod_jan2018_12.R'
'my_hg19.twobit'
rule make_riboqc_anno:
  input : GTF
  output: touch('riboqc/annot.done')
  run:
  	outdir = output[0].replace('annot.done','')
  	shell(r"""
  		set -x 
    source activate faToTwoBit
    faToTwoBit my_hg19.fa my_hg19.twobit
    cat {GTF} | perl -lanpe ' s/GL(\d+).1/chrUN_gl\1/g' > matchchrs.{GTF}
 	
    source activate cortexomics
	grep -ve 'GL00'  {GTF} > matchchrs.{GTF}
      mkdir -p riboqc
       R -e ' library(RiboseQC); prepare_annotation_files("{outdir}","{twobitfile}","matchchrs.{GTF}","Homo.sapiens","{ANNOBASE}") '

""")

rule run_riboqc:
  input : 'riboqc/annot.done','star/data/{sample}/.done',REF
  output: touch('riboqc/data/{sample}/.done'),'riboqc/data/{sample}/_for_SaTAnn'
  run:
    annofile = input[0].replace('annot.done',ANNOBASE+'.mainchrs.gtf_Rannot')
    bamfile = 'star/data/'+wildcards['sample']+'/'+wildcards['sample']+'.bam'
    outname = output[0].replace('.done','')
    report_file = 'riboqc/reports/'+wildcards['sample']+'/'+'riboqcreport.html'
    shell(r""" 
        R -e ' library(RiboseQC);RiboseQC::RiboseQC_analysis("{annofile}", bam="{bamfile}",dest_names="{outname}", report_file="{report_file}")' 
    """)


SaTAnn_script = '/fast/work/groups/ag_ohler/dharnet_m/satann_working/Satan_working_hg38_genc25_May8.R'

rule make_chrnamefile:
 input: GTF
 output: 'chrnames.txt'
 shell:r""" cut -d$'\t' -f1 {input} | uniq | sort | uniq | grep -v '#' > {input}"""

# rule run_satann:
#   input : 'riboqc/data/{sample}/.done'
#   output: touch('SaTAnn/{sample}/.done')
#   threads: 10
#   run:
#     shell(r"""
#     mkdir -p $( dirname {output} )

#     if [ ! -s riboqc/data/{wildcards.sample}/_for_SaTAnn ]; then
#       touch $(dirname {output} )/cannot_run_no_psites
#       exit
#     fi

  
#       Rscript {SaTAnn_script} \
#         {REF} \
#         riboqc/{ANNOBASE}.annoout \
#         {threads} \
#         riboqc/data/{wildcards.sample}/_for_SaTAnn \
#         $(dirname {output})/SaTAnn

#       """)
    
rule run_satann_group:
  input :
    lambda wc: 
      ['riboqc/data/'+s+'/.done' for s in samplegroups[wc['groupname']] ]
  output: touch('SaTAnn/{groupname}/.done')
  threads: 20
  run:

    sinputs = [i.replace('.done','_for_SaTAnn') for i in list(input)]
    sinputs = list(filter(os.path.isfile,sinputs))
    nopsites = 'true' if len(sinputs) is 0 else ''
    print(nopsites)
    print(sinputs)

    sinputsjoin = ','.join(sinputs)

    # print('\n\n\n')
    # print( (input[4]))
    # print('\n\n\n')
    # print( sinputs[4])
    # print('\n\n\n')
    shell(r"""
    mkdir -p $( dirname {output} )

     if [ {nopsites} ] ; then
     	touch $(dirname {output})/nopsites
     	exit
     fi

      Rscript {SaTAnn_script} \
        {REF} \
        riboqc/{ANNOBASE}.annoout \
        {threads} \
        {sinputsjoin} \
        $(dirname {output})/SaTAnn

      """)
#&& [[ -s {output} ]]

ribqc_report='/fast/work/groups/ag_ohler/dharnet_m/satann_working/riboseqc.Rmd'
rule knit_riboqc:
  input: 
    riboqcoutput=lambda wc:
      ['riboqc/data/'+s+'/.done' for s in samplegroups[wc['groupname']] ],
    ribqc_report=ribqc_report
  output: 'riboqc/reports/{groupname}/riboqcreport.html'
  run:
    #knitr's file path fuckery means we need to cd and hand it absolute paths, god knows where it's working directory ends up.
    indata = [i.replace('.done','_results_all_toknit') for i in list(input['riboqcoutput'])]
    indata = [os.path.abspath(r) for r in indata]
    indata = ','.join(indata)
    sinputs = ','.join(samplegroups[wildcards['groupname']])
    figpath = os.path.abspath(output[0]).replace('.html','')+'_fig'
    outpath = os.path.abspath(output[0])
    shell(r"""
    cd $(dirname {output})
    R -e 'rmarkdown::render("{input.ribqc_report}",params = list(input_list = "{indata}",input_list_names = "{sinputs}", output_fig_path = "{figpath}"),output_file = "{outpath}")'
""")


VCF_DICT = {'ONVC' : ['../ext_data/vcfs_october_2018/vcfs/0NVC_M07_240517','../ext_data/vcfs_october_2018/vcfs/0NVC_M07_240517/M07_mq_240517.vcf'],
  'OD5P' : ['../ext_data/vcfs_october_2018/vcfs/0D5P_M11_240517/M11_tumor_240517.fasta','../ext_data/vcfs_october_2018/vcfs/0D5P_M11_240517/M11_mq_240517.vcf'],
  'OMM' : ['../ext_data/vcfs_october_2018/vcfs/0MM745_M02_240517/M02_mq_240517.vcf']
}

modfastascript='../src/modify_orf_seqs.R'

rule modify_satannfasta:
  input:
    'groupedsatan/{group}.genomic.gtf',REF
  output:
    'groupedsatan/{group}.variantmod.fasta'
  run:
    vcfs = ','.join(VCF_DICT[wildcards['group']])
    shell(r"""
      Rscript {modfastascript} {input} {vcfs} {output}
       """)


rule make_id_table:
  input: GTF
  output: 'ids.txt'
  shell: r"""R -e 'library(tidyverse,quiet=T); library(rtracklayer,quiet=T);import("{input}") %>%mcols%>%as.data.frame%>% select(gene_id,gene_name)%>%distinct%>%write.table("{output}", col.names=TRUE, row.names=FALSE)'"""

DISPDIFF = 0
rule run_ribodiff:
  input: 'feature_counts/all_feature_counts'
  conda: 'ribodiff'
  output: touch('ribodiff/.done')
  shell: r"""Rscript ../exploration/pipeline/run_ribodiff.R {input} {DISPDIFF} $(dirname {output})"""





#rule run_template:
#  input: templatefinput
#  output: touch('run_template/.done')
#  shell: r"""template {input} $(dirname {output})"""
