#!/usr/bin/ nextflow

nextflow.enable.dsl=2
/* 
################################################################
        AMR VF Protein BERT Classification
################################################################ 
*/

process pyrodigal {
    tag {sample_id}
    //errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/protBERT/", mode: 'copy'
    label 'normal'
    conda './env/protBERT.yml'

    input:
    tuple val(sample_id), file("contigs.fasta"),  val(cpus)

    cpus {cpus} // setting slurm allocation dynamically 

    output:
    tuple val(sample_id), file("*"), optional:true, emit: protBERT_ch
    tuple val(sample_id), file("*.faa"), optional:true, emit: prot.faa 

    when:
    params.protBERT

    script:
    """
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate protBERT   
        pyrodigal -p meta -n -m -f gbk -a ${sample_id}.faa -d ${sample_id}.ffn  -s ${sample_id}.coord -o ${sample_id}.gbk -j ${cpus} -i ${sample_id}
        touch ${params.outdir}/${params.project_id}/${sample_id}/status_log/pyrodigal.finished 
    fi    
    """

}
process protBERT{
    tag {sample_id}
    //errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/protBERT/", mode: 'copy'
    label 'normal'
    conda './env/protBERT.yml'

    input:
    tuple val(sample_id), path(prot.faa)

    when:
    params.protBERT

    script:
    """
       	eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
       conda activate protBERT
       ##Command to run protein BERT
    fi
    """

}

