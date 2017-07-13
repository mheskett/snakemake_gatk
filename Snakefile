configfile: "config.yaml"

SAMPLESDIC = config["samples"]
JAVA="/usr/lib/jvm/jre-1.8.0/bin/java"
PICARD="/opt/installed/picard/picard-tools-1.110"
GATK="/opt/installed/GATK/GenomeAnalysisTK-3.5.jar"
REF="/home/exacloud/lustre1/SpellmanLab/havens/elements3/data/ucsc.hg19.fasta"  #note need *.fasta.fai and *.dict in same folder
HG19VCF="/home/exacloud/lustre1/SpellmanLab/heskett/refs/dbsnp_138.hg19.vcf.gz" #note need index file in same folder


#README - information on modification for use:
#when modifing the file please keep note of the spaces between the last charcter and the end "
#to change what rules (and so what processes) are run, modify the targeting function, referance targetingNotes.txt
#to change the which samples (that is the portion of the file names which refernce specific samples) are being used as input adjust the config file
#note when writing rules the length of inputs and outputs must be consistant, uses expand() to make this work out
 
def targeting():
    SAMPLES = list(SAMPLESDIC.keys())
    #put string the target from the rules to be run, using WILD where the sample name would go, assosications are makred in targetNotes.txt
    desiredWildTargets = ["GATK/WILD.filtered.vcf"]
    #if there is not sample name in the output file, put it in this list
    targetList = []
    for name in SAMPLES:
        for tar in desiredWildTargets:
            targetList.append(tar.replace("WILD", name))
    return targetList


rule all:
    input:
        targeting()

onsuccess:
    print("completed steps")

onerror:
    print("error has occured")


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

rule BQSR1:
    input:
        split=rules.SplitNCigarReads.output,
        ref={REF},
        known= {HG19VCF}
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
    
rule varCalling:
    input:
        reads=rules.BQSR2.output,
        ref={REF},
        known= {HG19VCF}
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




#rule moveToOut:
#    input:
#        expand("mapped/{sample}.rg.sorted.markdup.bam", sample=config["samples"])
#    output:
#        "outfiles/"
#    shell:
#        "cp {input} {output}"
    


#rule report:
#    input:
#        "outfiles/"
#    output:
#        "outfiles/report.html"
#    run:
#        from snakemake.utils import report
#        report("""
#        attmepts at adding individual moduals
#        """, output[0])

