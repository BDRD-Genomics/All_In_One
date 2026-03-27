#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
#########################################################################
                            Assemblies
*/

/* metaSPAdes PE trimmed assembly (dependency STAGE 5) */
process MetaSPAdes_Assembly {
    tag { sample_id }
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/spades/meta_pe_trim", mode: 'copy'
    label 'normal'
    conda "${baseDir}/env/spades.yml"
    cpus { cpus_in }
    memory { mem_in }
    //time '36h'
    input:
    //tuple val(sample_id), path(qc_SR_Read1), path(qc_SR_Read2)
    tuple val(sample_id), path(qc_SR_Read1), path(qc_SR_Read2), val(cpus_in), val(mem_in)
    output:
    tuple val(sample_id), file("contigs.fasta"), emit: metaspades_assembly_ch
    tuple val(sample_id), file("*"), emit: extra_metaspades_assembly_ch
    when:
    params.metaspades

    script:
    """
    if [ "${params.longreads}" == "true" ]; then
        echo "skip this process" > stage6_skipped.txt
    elif [ "${params.hybrid}" == "true" ]; then 
        ${params.metaSPAdes} -t ${task.cpus} \
                -1 ${qc_SR_Read1} \
                -2 ${qc_SR_Read2} \
                -m ${task.memory.giga} \
                -o . >> MetaSPAdes_Assembly_hybrid.log 2>&1

        # Filtering
        mv contigs.fasta raw_contigs.fasta 
        bash ${params.scripts}/filter_contigs_minlen.sh raw_contigs.fasta contigs.fasta


    else
        ${params.metaSPAdes} -t ${task.cpus} \
                -1 ${qc_SR_Read1} \
                -2 ${qc_SR_Read2} \
                -m ${task.memory.giga} \
                -o . >> MetaSPAdes_Assembly_shortreads.log 2>&1

        # Filtering
        mv contigs.fasta raw_contigs.fasta 
        bash ${params.scripts}/filter_contigs_minlen.sh raw_contigs.fasta contigs.fasta


    fi
    """
}

process Spades {
    tag { sample_id }
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/spades/meta_pe_trim", mode: 'copy'
    label 'normal'
    conda "${baseDir}/env/spades.yml"
    errorStrategy 'ignore'
    cpus { 80 }
    memory { '300 GB'}
    time '36h'
    input:
    tuple val(sample_id), path(qc_SR_Read1), path(qc_SR_Read2)

    output:
    tuple val(sample_id), file("contigs.fasta"), emit: spades_assembly_ch
    tuple val(sample_id), file("*"), emit: extra_spades_assembly_ch
    when:
    params.spades

    script:
    """
    if [ "${params.longreads}" == "true" ]; then
        echo "skip this process" > SPADES_skipped.txt
    elif [ "${params.hybrid}" == "true" ]; then 
        echo "skip this process" > SPADES_skipped.txt
    else
        ${params.SPAdes} -t ${task.cpus} \
                -1 ${qc_SR_Read1} -2 ${qc_SR_Read2} -o . >> Spades.log 2>&1
    fi
    """
}

process Dragonflye {

    tag { sample_id }

    publishDir "${params.outdir}/${params.project_id}/${sample_id}/dragonflye/", mode: 'copy'
    label 'normal'
    conda "${baseDir}/env/dragonflye.yml"
    errorStrategy 'ignore'
    cpus { 16 }
    memory { '32 GB'}
    time '36h'
    input:
    tuple val(sample_id), val(qc_SR_Read1), val(qc_SR_Read2), val(qc_long_read), val(q2)

    output:
    tuple val(sample_id), file("contigs.fa"), emit: dragonflye_assembly_ch
    tuple val(sample_id), file("*.gfa"), file("*.log"), file("*.txt"), file("*.gz"), emit: extra_dragonflye_assembly_ch

    script:
    """
    mkdir -p status_log
    #eval "\$(command conda 'shell.bash' 'hook' 2> /dev/null)"
    #conda activate dragonflye
    if [[ "${qc_long_read}" == "null" || ! -s "${qc_long_read}" ]]; then
        echo "[dragonflye] ERROR: No long reads provided!" >&2
        echo "FAIL: no long reads" > status_log/dragonflye_assembly.failed
        exit 1
    fi

    if [[ "${params.longreads}" == "true" && "${params.use_gsize}" == "true" ]]; then
        dragonflye --force \\
                   --gsize ${params.gsize} \\
                   --cpus ${task.cpus} \\
                   --tmpdir /dev/shm \\
                   --ram ${task.memory} \\
		   --minreadlen ${q2} \\
		   --minquality ${params.dragonflye_min_quality} \\
	           --keepfiles \\
                   --trim --trimopts "--discard_middle" --depth 100 --opts "--iterations 5" --racon 0 --polypolish 5 --polypolish_careful --nanohq  \\
                   --outdir . \\
                   --reads ${qc_long_read} > dragonflye.log 2>&1

    elif [[ "${params.hybrid}" == "true" && "${params.use_gsize}" == "true" ]]; then
        dragonflye --force \\
                   --gsize ${params.gsize} \\
                   --cpus ${task.cpus} \\
                   --tmpdir /dev/shm \\
                   --ram ${task.memory} \\
                   --minreadlen ${q2} \\
                   --minquality ${params.dragonflye_min_quality} \\
	  	   --keepfiles \\
                   --trim --trimopts "--discard_middle" --depth 100 --opts "--iterations 5" --racon 0 --polypolish 5 --polypolish_careful --nanohq  \\
                   --outdir . \\
                   --reads ${qc_long_read} \\
                   --R1 ${qc_SR_Read1} \\
                   --R2 ${qc_SR_Read2} > dragonflye.log 2>&1

    elif [[ "${params.longreads}" == "true" && "${params.use_gsize}" == "false" ]]; then
        dragonflye --force \\
                   --cpus ${task.cpus} \\
                   --tmpdir /dev/shm \\
                   --ram ${task.memory} \\
                   --minreadlen ${q2} \\
                   --minquality ${params.dragonflye_min_quality} \\
		           --keepfiles \\
                   --trim --trimopts "--discard_middle" --depth 100 --opts "--iterations 5" --racon 0 --polypolish 5 --polypolish_careful --nanohq \\
                   --outdir . \\
                   --reads ${qc_long_read} > dragonflye.log 2>&1

    elif [[ "${params.hybrid}" == "true" && "${params.use_gsize}" == "false" ]]; then
        dragonflye --force --cpus ${task.cpus} \\
                   --tmpdir /dev/shm \\
                   --ram ${task.memory} \\
                   --minreadlen ${q2} \\
                   --minquality ${params.dragonflye_min_quality} \\
		           --keepfiles \\
		           --trim --trimopts "--discard_middle" --depth 100 --opts "--iterations 5" --racon 0 --polypolish 5 --polypolish_careful --nanohq \\
                   --outdir . \\
                   --reads ${qc_long_read} \\
                   --R1 ${qc_SR_Read1} \\
                   --R2 ${qc_SR_Read2} > dragonflye.log 2>&1

    else
        echo "SKIP: dragonflye not enabled or unsupported mode" > dragonflye_skipped.txt
        echo "SKIP" > status_log/dragonflye_assembly.skipped
        exit 0
    fi

    echo "OK" > status_log/dragonflye_assembly.finished
    """
}

process Dragonflye_Medaka {

    tag { sample_id }

    publishDir "${params.outdir}/${params.project_id}/${sample_id}/dragonflye_medaka/", mode: 'copy'
    label 'medaka_gpu'
    conda "${baseDir}/env/dragonflye.yml"
    errorStrategy 'ignore'
    //clusterOptions = [ '--partition="normal"','--gpus=2' ]
    //cpus { 16 }
    //memory { '32 GB'}
    time '36h'

    input:
    tuple val(sample_id), val(qc_SR_Read1), val(qc_SR_Read2), val(qc_long_read), val(q2)

    output:
    tuple val(sample_id), file("contigs.fa"), emit: dragonflye_assembly_medaka_ch
    //tuple val(sample_id), file("*.gfa"), file("*.log"), file("*.txt"), emit: extra_dragonflye_medaka_assembly_ch
    tuple val(sample_id), file("*.gfa"), file("*.log"), file("*.txt"), file("*.gz"), emit: extra_dragonflye_medaka_assembly_ch
    script:
    """
    mkdir -p status_log
    #eval "\$(command conda 'shell.bash' 'hook' 2> /dev/null)"
    #conda activate dragonflye
    if [[ "${qc_long_read}" == "null" || ! -s "${qc_long_read}" ]]; then
        echo "[dragonflye] ERROR: No long reads provided!" >&2
        echo "FAIL: no long reads" > status_log/dragonflye_medaka_assembly.failed
        exit 1
    fi

    if [[ "${params.longreads}" == "true" && "${params.use_gsize}" == "true" && "${params.medaka}" == "true" ]]; then
        dragonflye --force \\
                   --gsize ${params.gsize} \\
                   --cpus ${task.cpus} \\
                   --tmpdir /dev/shm \\
                   --ram ${task.memory} \\
		           --minreadlen ${q2} \\
		           --minquality ${params.dragonflye_min_quality} \\
	               --keepfiles \\
                   --trim --trimopts "--discard_middle" --depth 500 --opts "--iterations 5" --racon 0 --polypolish 5 --polypolish_careful --nanohq --model ${params.medaka_model} \\
                   --outdir . \\
                   --reads ${qc_long_read} > dragonflye.log 2>&1

    elif [[ "${params.hybrid}" == "true" && "${params.use_gsize}" == "true" && "${params.medaka}" == "true" ]]; then
        dragonflye --force \\
                   --gsize ${params.gsize} \\
                   --cpus ${task.cpus} \\
                   --tmpdir /dev/shm \\
                   --ram ${task.memory} \\
                   --minreadlen ${q2} \\
                   --minquality ${params.dragonflye_min_quality} \\
                   --keepfiles \\
                   --trim --trimopts "--discard_middle" --depth 500 --opts "--iterations 5" --racon 0 --polypolish 5 --polypolish_careful --nanohq --model ${params.medaka_model} \\
                   --outdir . \\
                   --reads ${qc_long_read} \\
                   --R1 ${qc_SR_Read1} \\
                   --R2 ${qc_SR_Read2} > dragonflye.log 2>&1

    elif [[ "${params.longreads}" == "true" && "${params.use_gsize}" == "false" && "${params.medaka}" == "true" ]]; then
        dragonflye --force \\
                   --cpus ${task.cpus} \\
                   --tmpdir /dev/shm \\
                   --ram ${task.memory} \\
                   --minreadlen ${q2} \\
                   --minquality ${params.dragonflye_min_quality} \\
                   --keepfiles \\
                   --trim --trimopts "--discard_middle" --depth 500 --opts "--iterations 5" --racon 0 --polypolish 5 --polypolish_careful --nanohq --model ${params.medaka_model} \\
                   --outdir . \\
                   --reads ${qc_long_read} > dragonflye.log 2>&1

    elif [[ "${params.hybrid}" == "true" && "${params.use_gsize}" == "false" && "${params.medaka}" == "true" ]]; then
        dragonflye --force --cpus ${task.cpus} \\
                   --tmpdir /dev/shm \\
                   --ram ${task.memory} \\
                   --minreadlen ${q2} \\
                   --minquality ${params.dragonflye_min_quality} \\
                   --keepfiles \\
                   --trim --trimopts "--discard_middle" --depth 500 --opts "--iterations 5" --racon 0 --polypolish 5 --polypolish_careful --nanohq --model ${params.medaka_model} \\
                   --outdir . \\
                   --reads ${qc_long_read} \\
                   --R1 ${qc_SR_Read1} \\
                   --R2 ${qc_SR_Read2} > dragonflye.log 2>&1

    else
        echo "SKIP: dragonflye not enabled or unsupported mode" > dragonflye_skipped.txt
        echo "SKIP" > status_log/dragonflye_assembly.skipped
        exit 0
    fi

    echo "OK" > status_log/dragonflye_assembly.finished
    """
}

process Unicycler_Assembly {
    tag { sample_id }
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/unicycler/", mode: 'copy'
    label 'normal'
    conda "${baseDir}/env/unicycler.yml"
    cpus { 30 }
    memory { '162 GB'}
    time '36h'
    input:
    tuple val(sample_id), val(qc_SR_Read1), val(qc_SR_Read2), val(qc_long_read), val(q2), val(mode)

    output:
    tuple val(sample_id), file("assembly.fasta"), emit: unicycler_assembly_ch
    tuple val(sample_id), file("*"), emit: extra_unicycler_assembly_ch
    when:
    params.unicycler

    script:
    """
    if [[ "$mode" == "short" ]]; then
    unicycler --spades_options "--threads ${params.threads}" \\
              --short1 ${qc_SR_Read1} \\
              --short2 ${qc_SR_Read2} \\
              --out . \\
              --verbosity 2 \\
              --threads ${task.cpus} \\
              --keep 1 --min_fasta_length 1000 --mode ${params.unicycler_mode}
    elif [[ "$mode" == "hybrid"  ]]; then 
    unicycler --spades_options "--threads ${params.threads}" \\
              --long ${qc_long_read} \\
              --short1 ${qc_SR_Read1} \\
              --short2 ${qc_SR_Read2} \\
              --out . \\
              --verbosity 2 \\
              --threads ${task.cpus} \\
              --keep 1 --min_fasta_length 1000 --mode ${params.unicycler_mode}
    fi
    """
}

process Plasmid_Spades {
    tag { sample_id }
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/spades/plasmid", mode: 'copy'
    label 'normal'
    conda "${baseDir}/env/spades.yml"
    errorStrategy 'ignore'
    input:
    tuple val(sample_id), path(qc_SR_Read1), path(qc_SR_Read2)

    output:
    tuple val(sample_id), file("contigs.fasta"), emit: plasmidspades_assembly_ch
    tuple val(sample_id), file("*"), emit: extra_plasmidspades_assembly_ch
    when:
    params.plasmidspades

    script:
    """
    spades.py --plasmid \
               -1 ${qc_SR_Read1} -2 ${qc_SR_Read2} -t ${params.threads} -o . >> plasmidspades.log 2>&1
    """
}

