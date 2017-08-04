configfile: "config.yaml"

SAMPLESDIC = config["samples"]
JAVA="/usr/lib/jvm/jre-1.8.0/bin/java"
PICARD="/opt/installed/picard/picard-tools-1.110"
GATK="/opt/installed/GATK/GenomeAnalysisTK-3.5.jar"
REF="/home/exacloud/lustre1/SpellmanLab/havens/elements3/refs/ucsc.hg19.fasta"  #reference genome fasta file, note need *.fasta.fai and *.dict in same folder
knownSites="/home/exacloud/lustre1/SpellmanLab/heskett/refs/dbsnp_138.hg19.vcf.gz" #note need index file in same folder

#README - information on modification for use:
#when modifing the file please keep note of the spaces between the last charcter and the end "
#to change what rules (and so what processes) are run, modify the targeting function, referance targetingNotes.txt
#to change the which samples are being used for GATK as input adjust the config file
#all sequences with *_R1.fasta will be aligned if STAR is run 
#note when writing rules the length of inputs and outputs must be consistant, uses expand() to make this work out
 

#returns a list of targets (output of run rules) without wildcards 
#change lists in this function to adjust endpoint of script
def targeting():
    SAMPLES = list(SAMPLESDIC.keys())
    #put string of the final target from the rules to be run, using WILD where the sample name would go, assosications are makred in targetNotes.txt
    desiredWildTargets = ["GATK/WILD.filtered.vcf"]
    #if there is not sample name in the output file, put it in this list
    targetList = []
    for name in SAMPLES:
        for tar in desiredWildTargets:
            targetList.append(tar.replace("WILD", name))
    return targetList

#default target of snakemake
rule all:
    input:
        targeting()

onsuccess:
    print("completed steps")

onerror:
    print("error has occured")


#sets up inedx files of reference for STAR alignment
rule STARindex:
    shell:
        "./STARwork/snakemake STARindex"


#run STAR alignment, for each of the reads with *_R1.fasta pattern in samples folder
rule STARalign:
    input:
        "./STARwork/snakemake"


#adds readGroup for GATK pipeline to input bam file
rule replace_rg:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        "mapped/{sample}.rg.sorted.bam"
    log:
        "logs/picard/{sample}.rg.sorted.log"
    #params:
    shell:
        "{JAVA} -Xmx30g -jar {PICARD}/AddOrReplaceReadGroups.jar "
        "I={input} "
        "O={output} "
        "RGLB=lib1 "
        "RGPL=illumina "
        "RGPU=unit1 "
        "RGSM={wildcards.sample} "

#marks duplicates in a bam file from a bam file
rule mark_duplicates:
    input:
        rules.replace_rg.output
    output:
        bam="mapped/{sample}.rg.sorted.markdup.bam",
        metrics="mapped/{sample}.metrics.txt"
    log:
        "logs/picard/{sample}.rg.sorted.markdup.log"
    shell:
        "{JAVA} -Xmx30g -jar {PICARD}/MarkDuplicates.jar "
        "I={input} "
        "O={output.bam} "
        "METRICS_FILE={output.metrics} "
        "VALIDATION_STRINGENCY=LENIENT "
        "CREATE_INDEX=true "

#handels split reads inicated by CIGAR notation-N:skippped region, with correction for GATK in  output in bam file
rule SplitNCigarReads:
    #note: need have {sample}.dict and {sample}.fasta.fai in folder with reads
    input:
        bam=rules.mark_duplicates.output.bam,
        ref={REF}
    output:
        "GATK/{sample}.split.bam"
    log:
        "logs/GATK/{sample}.split.log"
    shell:
        "{JAVA} -Xmx30g -jar {GATK} "
        "-T SplitNCigarReads "
        "-R {input.ref} "
        "-I {input.bam} "
        "-o {output} "
        "-rf ReassignOneMappingQuality "
        "-RMQF 255 "
        "-RMQT 60 "
        "-U ALLOW_N_CIGAR_READS "

#analyze patters of covariation in the sequence in bam, outputs the patterns in table
rule BQSR1:
    input:
        split=rules.SplitNCigarReads.output,
        ref={REF},
        known= {knownSites}
    output:
        "GATK/{sample}.recal_data.table"
    log:
        "logs/GATK/{sample}.recal_data.table.log"
    shell:
        "{JAVA} -Xmx30g -jar {GATK} "
        "-T BaseRecalibrator "
        " -nct 8 "
        " -nt 1 "
        "-R {input.ref} "
        "-I {input.split} "
        "-knownSites {input.known} "
        "-o {output} "
        "--disable_auto_index_creation_and_locking_when_reading_rods"
#uses BQSR1 table to apply recalibration to an input bam and produces corrected bam
rule BQSR2:
    #note: skipped visulization steps
    input:
        split=rules.SplitNCigarReads.output,
        ref={REF},
        tab=rules.BQSR1.output
    output:
        "GATK/{sample}.split.recal_reads.bam"
    log:
        "logs/GATK/{sample}.split.recal_reads.log"
    shell:
        "{JAVA} -Xmx30g -jar {GATK} "
        "-T PrintReads "
        " -nct 4 "
        " -nt 1 "
        "-R {input.ref} "
        "-I {input.split} "
        "-BQSR {input.tab} "
        "-o {output} "


#idenify sites in bam which may be variant written into vcf file     
rule varCalling:
    input:
        reads=rules.BQSR2.output,
        ref={REF},
        known= {knownSites}
    output:
        "GATK/{sample}.raw.vcf"
    log:
        "logs/GATK/{sample}.raw.vcf"
    shell:
        "{JAVA} -Xmx30g -jar {GATK} "
        "-T HaplotypeCaller "
        "-R {input.ref} "
        "-I {input.reads} "
        "--genotyping_mode DISCOVERY "
        "-L {input.known} "
        "--output_mode EMIT_ALL_CONFIDENT_SITES "
        "-o {output} "
        "-stand_call_conf 10.0 "
        "-stand_emit_conf 30.0 "
        "-dontUseSoftClippedBases "


#takes vcf and applies hard-filtering variant calls
rule varFilter:
    input:
        var=rules.varCalling.output,
        ref={REF},
    output:
        "GATK/{sample}.filtered.vcf"
    log:
        "logs/GATK/{sample}.filtered.vcf"
    shell:
        "{JAVA} -Xmx30g -jar {GATK} "
        "-T VariantFiltration "
        "-R {input.ref} "
        "-V {input.var} "
        "-window 35 "
        "-o {output} "
        "-cluster 3 "
        "-filterName FS "
        "-filter 'FS > 30.0' "
        "-filterName QD "
        "-filter 'QD < 2.0' "



    
