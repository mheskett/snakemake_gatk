

SAMPLES, = glob_wildcards("/home/exacloud/lustre1/SpellmanLab/havens/elements3/samples/{smp}_R1.fastq") #will make a list of samples based on the files in sample folder which end with this pattern
STARtool="/opt/installed/STAR"
REF="/home/exacloud/lustre1/SpellmanLab/havens/elements3/ref/ucsc.hg19.fasta" #give genome fasta file
REFdir="/home/exacloud/lustre1/SpellmanLab/havens/elements3/ref" #give directory with referance file and index files, should not end in /, just the dir name
ANNOTA="/home/exacloud/lustre1/SpellmanLab/havens/elements3/ref/ucsc.hg.annotations.gtf" #annotation file for the genome reference
MAIL='havensj@ohsu.edu' #please change if using mail reporting function
OUTdir="/home/exacloud/lustre1/SpellmanLab/havens/elements3/align" #gives directory where aligned files will be output, end without /
THREADS = 16

rule all:
    input:
        expand("{OUTdir}/{smp}.Log.final.out", smp=SAMPLES)
#    shell:
#        "mail -s 'complete' {MAIL}"


#Make index files of the referance genome for STARalign, can run by snakemake STARindex.
rule STARindex:
    threads: {THREADS}
    shell:
        "{STARtool} "
        "--runMode genomeGenerate "
        "--runThreadN {threads} "
        "--genomeDir {REFdir} "
        "--genomeFastaFiles {REF}"


#default target, aligns reads to reference genome, requires index files 
rule STARalign:
    input:
        read1="samples/{smp}_R1.fastq",
        read2="samples/{smp}_R2.fastq",
        ref={REFdir}
    output:
        "{OUTdir}/{smp}.Log.final.out"
    threads: {THREADS}
    shell:
        "{STARtool} --runThreadN {threads} "
        "--genomeDir {input.ref} "
        "--readFilesIn {input.read1} "
        "--sjdbGTFfile {ANNOTA} "
        "--outFileNamePrefix {OUTdir}/{wildcards.smp}. "
        "--outSAMtype BAM SortedByCoordinate "
        "--outSAMstrandField intronMotif "
        "--quantMode GeneCounts TranscriptomeSAM "
        "--outSAMattributes NH HI NM MD jM jI "
        

