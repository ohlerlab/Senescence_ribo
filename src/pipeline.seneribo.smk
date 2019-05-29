import glob
import pandas as pd
from pathlib import Path
import ipdb

def is_nonempty(file):
  assert Path(file).stat().st_size
def is_over_size(file,n):
  assert Path(file).stat().st_size > n
def newfolder(file,newroot):
  file = Path(file)
  assert snakedir in file.parents , str(snakedir) + " doesn't seem to be in parents of " + str(file)
  return str(Path(*(newroot,)+file.relative_to(snakedir).parts[1:]))

shell.executable("bash")
shell.prefix("set -e  pipefail;")
# user set parameter

configfile: "../src/config.yaml"
config['root'] = Path(config['root'])


snakedir = Path(config['root']).resolve() / 'pipeline'
TMPDIR = Path('../tmp')


seqfilesdf = pd.read_csv(config['sample_files']).set_index("sample_id", drop=False)
sampledf = pd.read_csv(config['sample_parameter']).set_index("sample_id", drop=False)

assert sampledf.sample_id.is_unique
assert seqfilesdf.file.is_unique
assert set(seqfilesdf.sample_id) == set(seqfilesdf.sample_id)
for sample in sampledf.sample_id:
  for f in seqfilesdf.loc[[sample],'file']:
    assert snakedir/'input'/sample in Path(f).parents, 'This path should be input/sample/file :'+f


samples = list(sampledf['sample_id'].unique())
fastqs = list(seqfilesdf['file'].unique())
ribosamples = sampledf.sample_id[sampledf['assay']=='ribo']

#for trimming CDS for riboseq
STARTCODTRIM=config['STARTCODTRIM']
STOPCODTRIM=config['STOPCODTRIM']
#local copies of the annotation
REF = snakedir / Path(config['REF_orig']).with_suffix('.fa').name
GTF = snakedir / Path(config['GTF_orig']).with_suffix('.gtf').name
GFF = snakedir / Path(config['GFF_orig']).with_suffix('.gff3').name
CDSGTF = GTF.with_suffix('.cds.gtf')
TRIMCDSGTF = str(CDSGTF).replace('.gtf','.trim_'+str(STARTCODTRIM)+'_'+str(STOPCODTRIM)+'.gtf')
RNAFASTA = GTF.with_suffix('.fa')
CODINGFASTA=GTF.with_suffix('.coding.fa')
PROTEINFASTA=GTF.with_suffix('.protein.fa')
CDSFASTA=GTF.with_suffix('.cds.fa')
BED=GTF.with_suffix('.bed')


rule all:
  input:
    seqfilesdf.file.unique(),
    expand("processed_reads/{sample}/.done", sample = sampledf.sample_id.unique()),
    expand("fastqc/data/{sample}/.done", sample = samples),
    expand("star/data/{sample}/.done", sample = samples),
    expand("qc/data/{sample}/.done", sample = samples),
    ("multiqc/multiqc_report.html"),
    expand("feature_counts/data/{sample}/feature_counts", sample = samples),
    # expand("feature_counts/all_feature_counts"),
    # # expand("bigwigs/{group}/{strand}/{istrans}.done",group = samples,strand=strands,istrans=istransvals),
    # expand("mergedbigwigs/{group}/{strand}/{istrans}.done",group = GROUPS,strand=strands,istrans=istransvals),
    expand('riboqc/reports/{sample}/riboqcreport.html', sample = ribosamples),
    expand('SaTAnn/{sample}/.done', sample = ribosamples),
    expand('cdsmax/{sample}/cdsmax_offsets.tsv', sample = ribosamples),



rule link_in_ref:
  input: config['REF_orig']
  output: REF
  shell:r"""
      ln -fs {input} {REF}
      """

rule link_in_files:
  input: 'input/{sample}/{fastq}'
  output: 'preprocessed_reads/{sample}/{fastq}'
  run:  
    sample = wildcards['sample']
    fastq = wildcards['fastq']
    shell(r"""
      mkdir -p $(dirname {output})
      ln -sf $(readlink -f input/{sample}/{fastq}) {output}
    """)


rule cutadapt_reads:
  input: 'preprocessed_reads/{sample}/{fastq}'
  output: 'cutadapt_reads/{sample}/{fastq}'
  conda: '../envs/cutadapt.yml'
  params: MINREADLENGTH=config['MINREADLENGTH'],MAXREADLENGTH=config['MAXREADLENGTH'],QUALLIM=config['QUALLIM']
  shell: r"""    #   set -evx
      set -e
       mkdir -p cutadapt_reads/{wildcards.sample}/
        zcat {input} \
           | cutadapt \
             -a {params.adaptorseq} \
            --minimum-length {params.MINREADLENGTH} \
            --maximum-length {params.MAXREADLENGTH} \
            -q {params.QUALLIM} - \
        2> cutadapt_reads/{wildcards.sample}/{wildcards.fastq}.cutadaptstats.txt \
        | gzip  > {output}
"""

rule collapse_reads:
    input: 'cutadapt_reads/{sample}/{fastq}'
    output: 'collapse_reads/{sample}/{fastq}'
    params: collapse_reads_script = config['collapse_reads_script']
    run:
        sample = wildcards['sample']
        shell(r"""
       set -evx
     
       mkdir -p collapse_reads/{sample}/
     
       zcat {input}  \
         | {params.collapse_reads_script} {wildcards.sample} \
         2> collapse_reads/{wildcards.sample}/{wildcards.fastq}.collreadstats.txt \
         | cat > {output}
     """)
        is_over_size(output[0],100)

#this finds small filse
 # find collapse_reads/ -name "*.fastq.gz" -size -100M
# find trim_reads/ -name "*.fastq.gz" -size -10M  | xargs ls -latr
# #this finds everything in a certain rule that's less than 10M and then quits
# for i in $(find trim_reads/ -name "*.fastq.gz" -size -10M);do   find . -name $(dirname $i | xargs basename) | grep -v input | grep -v cutadapt; done
rule trim_reads:
    input: 'collapse_reads/{sample}/{fastq}'
    output: 'trim_reads/{sample}/{fastq}'
    params: REMOVE8NBIN = config['REMOVE8NBIN']
    run:
        sample = wildcards['sample']
        shell(r"""
       set -evx
     
       OUTDIR=$(dirname {output})
       mkdir -p  $OUTDIR
     
       {params.REMOVE8NBIN} {input} {output}

       gzip -f {output}
       mv {output}.gz {output}

     """)

rule make_trna_rrna_indices:
  input: GTF,contaminants="../contaminants/contaminants.fa"
  output: touch('tRNA_rRNA_index/tRNA_rRNA_index.done')
  run:
    outprefix = output[0].replace('.done','')
    fafile =outprefix+'.fa'
    
    shell(r"""
      source activate tophat
      #get rRNAs and tRNAs, rename the tRNAs to exons for GenePRed,
      #get the gene type out and stick it front of the transcript id
      #for better names in the fasta
      cp {input.contaminants} {fafile}

       bowtie2-build {fafile} {outprefix} -p {threads}
       
      """)

rule make_bowtie_indices:
  input: REF
  output: touch('bowtie_index/.done')
  threads: 8
  run:
    outprefix = output[0].replace('.done','')
    fafile =outprefix+'.fa'
    shell(r"""
      #get rRNAs and tRNAs, rename the tRNAs to exons for GenePRed,
      #get the gene type out and stick it front of the transcript id
      #for better names in the fasta
       bowtie2-build {REF} bowtie_index/ref 
       cp {REF} bowtie_index/ref.fa
       bowtie2-build {RNAFASTA} bowtie_index/transcriptome
       cp {RNAFASTA} bowtie_index/transcriptome.fa 
       
      """)

collate_idxscript = "../exploration/collate_idx.R"

rule filter_tRNA_rRNA:
    input: 
      'trim_reads/{sample}/{fastq}',
      'tRNA_rRNA_index/tRNA_rRNA_index.done'  
      # filter_index    
    output: 'filter_reads/{sample}/{fastq}',
    threads: 8
    conda: '../envs/bowtie2'
    params: 
      indexname = lambda wc,input: input[1].replace('.done',''),
      outdir = lambda wc,output: os.path.dirname(output[0]) 
    shell:"""
       set -evx

       [ -f {params.outdir} ] && rm -rf {params.outdir}
     
       mkdir -p  {params.outdir}

      bowtie2 \
        -x {params.indexname} \
        -L 20  \
        -p {threads}  \
        -N 0 \
        -U  {input[0]} \
        --un-gz {output[0]} \
        --no-unal \
        2> {output[0]}.alignreport.log > {output[0]}.filtered_reads.sam


        samtools view -bh  {output[0]}.filtered_reads.sam \
        | samtools sort -@ {threads}  > {output[0]}.filtered_reads.bam

      samtools index {output[0]}.filtered_reads.bam
      
      #those which mismatch twice should not be included
      samtools view -hb {params.outdir}/filtered_reads.bam \
      | bamtools filter -tag XM:2-10 -in - -out /dev/stdout \
      | samtools view -H > {output[0]}.mm.sam
      #>> {params.outdir}/unmapped.sam
     
      #group the idx columns stuff is from 
      samtools idxstats {output[0]}.filtered_reads.bam \
      | perl -lanpe 's/^(\S+)_[^_\s]+\t/$1\t/' > {output[0]}.idxtmp

      #Rscript --vanilla {collate_idxscript} {output[0]}.idxtmp {params.indexname}.fa

      samtools stats {output[0]}.filtered_reads.bam > samtools stats {output[0]}.filtered_reads.bam.stats

    """


#takes a file, which has 
def get_processed_files(wc): 
  if wc['sample'] in ribosamples:
    return [ newfolder(fq,'filter_reads') for fq in seqfilesdf.loc[ [wc['sample']],'file' ] ]
  else:
    return seqfilesdf['file'][[wc['sample']]]

#this rule is the 'signal spliter where we go from sample to indiv fastqs
rule link_processed_reads:
  input: get_processed_files
  output: touch('processed_reads/{sample}/.done')
  run:
    shell(r"""
        mkdir -p processed_reads/{wildcards.sample}
        ln -rifs $(readlink -f {input}) processed_reads/{wildcards.sample}/
    """)


rule fastqc:
     input: 'processed_reads/{sample}/.done'
     output: touch('fastqc/data/{sample}/.done')
     threads: 4
     log:'fastqc/reports/{sample}/fastqc.log'
     params:
      reads = lambda wc: [fq.replace('input/','processed_reads/') for fq in seqfilesdf['file'][[wc['sample']]]],
      outdir = lambda wc: 'fastqc/data/'+wc.sample+'/'
     shell: '''
          OUTDIR=$(dirname {output[0]})
          mkdir -p {params.outdir}
          wait $(for i in {params.reads}; do $( fastqc -o {params.outdir} $i ) & done) 
        '''

##################
#process the annotation
##################

rule gffread:
  input: REF=config['REF_orig'],GTF=config['GTF_orig'],GFF=config['GFF_orig']
  output: GTF,CDSGTF,RNAFASTA,CDSFASTA,BED,GFF
  conda: '../envs/gffread.yml'
  shell: r""" 
      ln -s {input.REF} {REF}
      # set -x
      #with filtering output all sequences
      cat {input.GTF} \
      | sed -r  's/((transcript_id|gene_id|protein_id|ID|Parent|exon_id|havana_gene|havana_transcript)\W+\w+)\.[0-9]+/\1/g' \
      > {GTF}

      cat {input.GFF} \
      | sed -r  's/((transcript_id|gene_id|protein_id|ID|Parent|exon_id|havana_gene|havana_transcript)\W+\w+)\.[0-9]+/\1/g' \
      > {GFF}


      #needs gff - output exon sequences
      cat {input.GFF} |  grep -P -e'\texon\t|^\#' | gffread - -F -E -g {REF} -W -w {RNAFASTA} -o /dev/null

      #Note we are now minus the transcript and exon entries for these
      #now make GTF

      #| grep -P -e'\tCDS\t|^\#' 
     #with filtering, output the coding sequences filteirng out the ones that aren't in frame, have a stop codon, are pseudogenes etc.
      
      cat {input.GFF}  \
        | sed -r  's/((transcript_id|gene_id|protein_id|ID=\w+|Parent)\W+\w+)\.[0-9]+/\1/g' \
        | gffread - -C -V -J --no-pseudo  -F -E -g {REF} \
        -W -w {CODINGFASTA} -x {CDSFASTA} -y {PROTEINFASTA} -T \
        -o /dev/stdout \
        | awk -v FS="\t" -v OFS="\t" '{{if($3=="CDS"){{$3="exon";print $0}}}}' \
         > {CDSGTF}

      #now make bed
      cat {input.GTF} | awk '{{print $1,$4,$5,"name",$6,$7}}' > {BED}
    """

# rule trim_CDS:
#   input: GTF,
#   params:
#     STARTCODTRIM=15,
#     STOPCODTRIM=5,
#   output: TRIMCDSGTF
#     outdir= '.'
#   shell:r"""
#     Rscript ../src/trim_coding_sequences.R {input[0]} {params[0:2]} {output[0]}
#   """

rule star_index:
  input: REF=REF,GTF=GTF
  output: touch('starindex/.done')
  threads: 8
  run:
    shell(r"""
      STAR \
      --runThreadN {threads} \
      --runMode genomeGenerate \
      --genomeDir $(dirname {output}) \
      --sjdbGTFfile {input.GTF} \
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
          done = touch('star/data/{sample,[^/]+}/.done'),bam='star/data/{sample}/{sample}.bam'
     threads: 8
     run:
          input.STARINDEX=input.STARINDEX.replace('.done','')
          markdup = '' if sampledf.assay[wildcards['sample']] == 'ribo' else '-m'
          platform = 'NotSpecified'
          inputdir = os.path.dirname(input['fastqs'])
          outputdir = os.path.dirname(output[0])
          read_pattern = sampledf.read_pattern[wildcards['sample']]
          fastqops = get_fastqops(inputdir,read_pattern,lstring='<( zcat ',rstring=')')
          repdir = outputdir.replace('data','reports')
          tophatindex =input['bowtie_index'].replace('.done','')
          
          halfthreads = threads/2
          sortmem = str(int(5000/halfthreads))+'M'

          # remap = '1' if sampledf.assay[wildcards['sample']] == 'ribo' else ''
          remap = '' 
          
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


rrna_intervals = 'qc/picard_rrna_intervals.txt'
refflat = snakedir/ 'qc' / Path(config['GFF_orig']).with_suffix('.refflat').name

rule make_picard_files:
  input: GTF,'star/data/'+samples[0]+'/.done'
  output: intervals=rrna_intervals,refflat=refflat
  conda: '../envs/picard'
  shell:r"""
         samtools view -H star/data/{samples[0]}/{samples[0]}.bam > {output.intervals}
        
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
        singleendflag = lambda wc: ' -e ' if sampledf.loc[wc['sample'],'library_layout'] != 'PAIRED' else '',
        bamfile = lambda wc:'star/data/'+wc['sample']+'/'+wc['sample']+'.bam' ,
        scriptdir = lambda wc:  config['rnaseqpipescriptdir']
     shell: """
          set -exv
          
        OUTDIR=$(dirname {output.done})
        mkdir -p qc/reports/{wildcards.sample}/


        #{params.scriptdir}/run_RNA-SeQC.sh -i {params.bamfile} {params.singleendflag} -t {GTF} -r {REF} -o ${{OUTDIR}}/

        {params.scriptdir}/read_statistic_report.sh \
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

        {params.scriptdir}/read_duplication.sh \
         -i {params.bamfile} \
         -o ${{OUTDIR}}/duplication/

          """
     

def get_multiqc_dirs(wildcards,input):
      reportsdirs = list(input)
      reportsdirs=[s.replace('star/data','star/reports') for s in reportsdirs]
      reportsdirs=[s.replace('tophat2/data','tophat2/reports') for s in reportsdirs]
      reportsdirs=[os.path.dirname(s) for s in list(reportsdirs)]
      assert len(reportsdirs) > 0
      return(reportsdirs)

rule multiqc:
  input:
      expand("fastqc/data/{sample}/.done", sample = samples),
      expand("star/data/{sample}/.done", sample = samples),
      expand("qc/data/{sample}/.done", sample = samples),
      # expand("tophat2/data/{sample}/.done", sample = samples),
      # [f.replace('input','filter_reads') for f in  seqfilesdf.file[ribosamples]],
      expand("feature_counts/data/{sample}/feature_counts", sample = samples),
      # 'sample_file.txt'
  conda: '../envs/multiqc'
  params: 
    multiqcscript = config['multiqcscript'],
    sample_reads_file=config['sample_files'],
    reportsdirs= get_multiqc_dirs,
    sampnames = '--sample-names '+config.get('samplenamesfile') if config.get('samplenamesfile') else ''
  output:
    'multiqc/multiqc_report.html'
  shell:r"""
      cat {params.sample_reads_file} | sed 's/.fastq.gz//g' | sed 's/\t.*\//\t/g' \
      | awk -vOFS='\t' 'BEGIN{{print "fastqname","samplename"}}{{sumsamp[$1] = sumsamp[$1]+1;print $2,$1"_fq"sumsamp[$1]}}' \
      > multiqc/samplenames.txt

      {params.multiqcscript} {params.reportsdirs} -fo $(dirname {output[0]}) {params.sampnames}
      """



rule feature_counts:
     input:
          TRIMCDSGTF=TRIMCDSGTF,
          BAM='star/data/{sample}/{sample}.bam'
     output:
          done = 'feature_counts/data/{sample}/feature_counts'
     threads: 2
     run:
          
          #protocol type
          protocol = sampledf.protocol[wildcards['sample']]
          if (protocol == 'no'):
               protocol = 0
          elif (protocol == 'yes'):
               protocol = 1
          elif (protocol == 'reverse'):
               protocol = 2
          else:
               sys.exit('Protocol not known!')

          library = sampledf.library_layout[wildcards['sample']]

          #get library type
          if (library == 'PAIRED'):
               library = '-p'
          else:
               library = ''

          #multimapper type
          countmultimappers = ' ' 
          
          # if (wildcards['region']=='tRNAs'):
          #   featuretype = 'tRNA'
          #   countmultimappers = '-M --fraction'


          sample = wildcards['sample']
          # region = wildcards['region']
          region = 'CDS'
          # rangebam = input['readrangefilt']

          groupcol = 'gene_id'
          
          shell(r"""
          set -ex
          mkdir -p feature_counts/data/{sample}/
          mkdir -p feature_counts/reports/{sample}/
          featureCounts \
            -Q 50 \
            -T {threads} \
            -t {region} -g {groupcol} \
            -a {input.TRIMCDSGTF} \
            -s {protocol} {library} {countmultimappers} \
            -o $(realpath feature_counts/data/{sample}/feature_counts) \
            $(realpath  {input.BAM}) \
             &> feature_counts/reports/{wildcards.sample}/{wildcards.sample}.feature_counts.log

          """)

rule make_riboqc_anno:
   input : GTF,REF
   output: touch('riboqc/annot.done')
   params:
    speciesname = lambda wc,input: Path(input[1]).name.replace('.fa','') + '.foo',
    annobase = lambda wc,input: input[0].replace('.gtf',''),
    outdir = lambda wc,output: output[0].replace('annot.done',''),
    twobitfile = lambda wc,input: input[1].replace('.fa','.twobit'),
    gtfmatchchrs = lambda wc,input: input[0].replace('.gtf','.matchchrs.gtf')
   conda: '../envs/riboseqc'
   shell:r"""
         set -x
     faToTwoBit {REF} {params.twobitfile}

     samtools faidx {REF}

      awk -vOFS="\t" '{{print $1,0,$2}}'  {REF}.fai | bedtools intersect -b -  -a {GTF} > {params.gtfmatchchrs}

     mkdir -p riboqc
     conda activate R2
     R -e ' library(RiboseQC); prepare_annotation_files("{params.outdir}","{params.twobitfile}","{params.gtfmatchchrs}","{params.speciesname}","{params.annobase}",forge_BS=FALSE,genome_seq=FaFile("{REF}")) '
 """

rule run_riboqc:
   input: 'riboqc/annot.done'
   output: touch('riboqc/data/{sample}/.done'),'riboqc/data/{sample}/_for_SaTAnn','riboqc/reports/{sample}/riboqcreport.html'
   # conda: '../envs/riboseqc'
   params:    
     annofile = lambda wc,input: input[0].replace('annot.done',Path(GTF).name.replace('.gtf','.matchchrs.gtf_Rannot')),
     bamfile = lambda wc,input: 'star/data/'+wc['sample']+'/'+wc['sample']+'.bam',
     outname = lambda wc,output: output[0].replace('.done',''),
     report_file = lambda wc: 'riboqc/reports/'+wc['sample']+'/'+'riboqcreport.html',
   shell:r"""
         set -x
         mkdir -p {params.outname}
         mkdir -p riboqc/reports/{wildcards.sample}
         R -e 'library(RiboseQC);RiboseQC::RiboseQC_analysis("{params.annofile}", bam="{params.bamfile}",rescue_all_rls=TRUE,dest_names="{params.outname}", genome_seq = "{REF}", report_file="{params.report_file}")'
     """

rule run_satann:
  input : 'riboqc/data/{sample}/.done','riboqc/annot.done'
  output: touch('SaTAnn/{sample}/.done')
  # input: lambda wc: 
      # ['riboqc/data/'+s+'/.done' for s in samplegroups[wc['groupname']] ]
  # output: touch('SaTAnn/{groupname}/.done')
  params: 
    for_satannfile = lambda wc,input: 'c("'+('","'.join(['riboqc/data/'+s+'/''_for_SaTAnn' for s in [wc['sample']] ]))+'")',
    annofile = lambda wc,input: 'riboqc/'+Path(GTF).name.replace('.gtf','.matchchrs.gtf_Rannot'),
    outputdir = lambda wc,output: output[0].replace('.done','')
  threads: 10
  shell:r"""
    set -ex
      mkdir -p {params.outputdir}

      R -e 'devtools::load_all("~/work/Applications/SaTann");run_SaTAnn(for_SaTAnn_file = {params.for_satannfile},annotation_file = "{params.annofile}",genome_seq = "{REF}", n_cores = {threads},prefix="{params.outputdir}") '
        
      """



rule cdsmax:
  input: 'star/data/{sample}/{sample}.bam',GTF,REF
  output: 'cdsmax/{sample}/cdsmax_offsets.tsv'
  shell: r"""Rscript --vanilla ../src/offsets_cdsmax.R {input} $(dirname {output[0]})"""

rule cds_offset_count:
  input: 'star/data/{sample}/{sample}.bam',GTF,REF,'cdsmax/{sample}/cdsmax_offsets.tsv'
  output: 'cdsmax/{sample}/gene_id_counts.tsv',
  shell: r"""Rscript --vanilla ../src/cds_offset_count.R {input} $(dirname {output[0]}) {STARTCODTRIM} {STOPCODTRIM}"""




rule aggregate_feature_counts:
  input : expand("feature_counts/data/{sample}/feature_counts", sample = samples),
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

      echo "feature_id {samples}" | cat - {output} > {output}tmp
      mv {output}tmp {output}
    
      """)

rule run_xtail:
  input: 'feature_counts/all_feature_counts'
  output: touch('xtail/.done')
  shell: r"""Rscript ../src/run_xtail.R {input} $(dirname {output})"""



