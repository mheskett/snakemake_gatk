configfile: "config.yaml"

SAMPLESDIC = config["samples"]

def targeting():
    SAMPLES = list(SAMPLESDIC.keys())
    #put string the target from the rules to be run, using WILD where the sample name would go, assosications are makred in targetNotes.txt
    desiredWildTargets = ["mapped/WILD.rg.sorted.bam", "mapped/WILD.rg.sorted.markdup.bam"]
    #if there is not sample name in the output file, put it in this list
    targetList = ["outfiles/report.html"]
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
        "logs/picard/mapped/{sample}.rg.sorted.log"
    params:
        "RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM={sample} VALIDATION_STRINGENCY=SILENT"
    shell:
        "/usr/lib/jvm/jre-1.8.0/bin/java -Xmx30g -jar /opt/installed/picard/picard-tools-1.110/AddOrReplaceReadGroups.jar "
        "I={input} "
        "O={output} "
        "RGLB=lib1 "
        "RGPL=illumina "
        "RGPU=unit1 "
        "RGSM=test "


rule mark_duplicates:
    input:
        #data=lambda wildcards: config["samples"][wildcards.sample],
        data="mapped/{sample}.rg.sorted.bam"
    output:
        bam="mapped/{sample}.rg.sorted.markdup.bam",
        metrics="mapped/{sample}.metrics.txt"
    log:
        "logs/picard/mapped/{sample}.rg.sorted.markdup.log"
    params:
        "REMOVE_DUPLICATES=true"
    shell:
        "/usr/lib/jvm/jre-1.8.0/bin/java -Xmx30g -jar /opt/installed/picard/picard-tools-1.110/MarkDuplicates.jar "
        "I={input.data} "
        "O={output.bam} "
        "METRICS_FILE={output.metrics} "
        "VALIDATION_STRINGENCY=LENIENT "
        "CREATE_INDEX=true "

rule moveToOut:
    input:
        expand("mapped/{sample}.rg.sorted.markdup.bam", sample=config["samples"])
    output:
        "outfiles/"
    shell:
        "cp {input} {output}"
    


rule report:
    input:
        "outfiles/"
    output:
        "outfiles/report.html"
    run:
        from snakemake.utils import report
        report("""
        attmepts at adding individual moduals
        """, output[0])

