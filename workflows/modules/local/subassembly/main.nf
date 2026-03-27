#!/usr/bin/ nextflow

nextflow.enable.dsl=2
/*
#########################################################################
                            Subassembly
*/

process Subassembly {
    tag {sample_id}
    //errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/subassembly/", mode: 'copy'
    label 'normal'
    conda './env/md.yml'

    input:
    tuple val(sample_id), path(qc_SR_Read1), path(qc_SR_Read2)

    output:
    tuple val(sample_id), file("*"), emit: metaspades_assembly_ch

    when:
    params.subassembly

    script:
    """
    if [ "${params.longreads}" == "true" ]; then
        echo "skip this process" > SPADES_skipped.txt
    elif [ "${params.hybrid}" == "true" ]; then 
        echo "skip this process" > SPADES_skipped.txt
    else
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate md   
        seqkit sample -s 5 \
                -1 ${qc_SR_Read1} \
                -n 25000
        seqkit sample -s 5 \
                -1 ${qc_SR_Read2} \
                -n 25000
        touch ${params.outdir}/${params.project_id}/${sample_id}/status_log/subassembly.finished 
    fi    
    """

}
