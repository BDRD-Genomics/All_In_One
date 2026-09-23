#!/usr/bin/ nextflow

nextflow.enable.dsl=2

/*
========================================================================================
   Exercise Reporting Workflow
========================================================================================
   Github   : 
   Contact  :     
----------------------------------------------------------------------------------------

*/

process Exercise_Report_Full {
    publishDir "${params.outdir}/${params.project_id}/exercise_report", mode: 'copy'
    label 'lowmem'
    conda "$baseDir/env/pandas_env.yml"
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), val(assembler)
    val dep1
    val	dep2
    val	dep3
    val	dep4

    output:
    path "*.docx", emit: report_out_ch

 
    when:
    params.exercise_report

    script:
    """
    python \\
        ${params.scripts}/exercise_report_Full_v3.py \\
        -s ${sample_id} -i ${params.outdir}/${params.project_id}/ -o . \\
        -p ${params.project_id} -t ${baseDir}/templates/FFBS_OA_Full_Report_Template.docx \\
        -a ${assembler}
    """
}
