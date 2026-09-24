#!/usr/bin/ nextflow

nextflow.enable.dsl=2

/*
========================================================================================
   BlastX
========================================================================================
*/


    process UNMAPPED_LR_TO_FASTA {
        tag { "${sample_id} | ${source_assembler}" }

        publishDir { "${params.outdir}/${params.run_id}/${sample_id}/unmapped_read_analysis/filtered_reads/" },
            mode: 'copy',
            overwrite: true

        errorStrategy 'ignore'

        input:
        tuple val(sample_id),
              val(source_assembler),
              path(unmapped_lr),
              val(min_length),
              val(mode)

        output:
        tuple val(sample_id),
              val(source_assembler),
              path("${sample_id}_${source_assembler}_unmapped_min700.fasta"),
              emit: fasta_ch

        script:
        """
        seqkit fq2fa ${unmapped_lr} \
        | seqkit seq -m ${min_length} \
        > ${sample_id}_${source_assembler}_unmapped_min700.fasta
        """
    }



process BLASTN_NT_contigs {
    tag { "${assembler} | ${sample_id}" }
    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/blast/${assembler}" }, mode: 'copy'
    cpus 64
    memory '128 GB'
    errorStrategy 'ignore'
    input:
    tuple val(sample_id), val(assembler), path(contigs)

    output:
    tuple val(sample_id), val(assembler), file("*_core_nt_blastn.tsv"), emit: blastn_contigs_nt_ch

    when:
    (params.blast_contigs)

    script:
    """
    blastn -query ${contigs} \\
        -db ${params.blast_core_nt} \\
        -outfmt '6 qaccver qlen qstart qend saccver stitle sstart send length slen pident qcovs sstrand gaps evalue bitscore score' \\
        -out ${sample_id}_${assembler}_core_nt_blastn.tsv -evalue 1e-5 -num_threads ${task.cpus} -perc_identity 90 -max_target_seqs 10

    sed -i "1i Query_ID\tQuery_length\tqstart\tqend\tSubject_accession\tTitle\tsstart\tsend\taln_length\tslen\tpercent_id\tTotal_Query_Cov\tsstrand\tgaps\tevalue\tbitscore\tscore" ${sample_id}_${assembler}_core_nt_blastn.tsv
    """
}

process PROKKA {
    tag { "${assembler} | ${sample_id}" }
    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/prokka/" }, mode: 'copy'
    cpus 1
    memory '16 GB'
    errorStrategy 'ignore'
    input:
    tuple val(sample_id), val(assembler), path(contigs)

    output:
    tuple val(sample_id), val(assembler), path("${assembler}"), emit: prokka_all_ch
    tuple val(sample_id), val(assembler), file("${assembler}/*.fna"),file("${assembler}/*.faa"), file("${assembler}/*.gff"),file("${assembler}/*.gbk"), emit: prokka_ch

    when:
    (params.prokka)

    script:
    """
    export SHELL=/bin/bash
    mkdir -p ${assembler} tmp
    export TMPDIR="\$PWD/tmp"
    export TMP="\$TMPDIR"
    export TEMP="\$TMPDIR"
    prokka ${contigs} \\
       --outdir ${assembler} \\
       --prefix ${sample_id}_${assembler} \\
       --cpus ${task.cpus} \\
       --force
    """
}


process BUSCO {
    tag { "${assembler} | ${sample_id}" }
    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/busco/${assembler}" }, mode: 'copy'
    label 'lowmem'
    errorStrategy 'ignore'
    input:
    tuple val(sample_id), val(assembler), path(contigs_fasta)

    output:
    tuple val(sample_id), val(assembler), file("*"), emit: busco_ch

    when:
    (params.busco)

    script:
    """
    busco -i ${contigs_fasta} --mode geno --auto-lineage -f --download_path ${params.busco_db} --out ${sample_id}_${assembler} --offline --cpu ${task.cpus}
    """
}

process AMR_VF_BLAST {
    tag { "${assembler} | ${sample_id}" }
    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/AMR_VF/${assembler}" }, mode: 'copy'
    label 'lowmem'
    errorStrategy 'ignore'
    input:
    tuple val(sample_id), val(assembler), path(contigs_fasta)

    output:
    tuple val(sample_id), val(assembler), file("*"), emit: amr_vf_ch

    when:
    (params.amr_vf)

    script:
    """
    blastn -query ${contigs_fasta} -db ${params.addgeneDB} -outfmt '6 qaccver qlen qstart qend saccver stitle sstart send length slen pident qcovs sstrand gaps evalue bitscore score' -out ./${sample_id}_${assembler}_addgene_plasmids_blastn.tsv -evalue 1e-5 -num_threads 32 -perc_identity 80 -max_target_seqs 10
    blastn -query ${contigs_fasta} -db ${params.snapgeneDB} -outfmt '6 qaccver qlen qstart qend saccver stitle sstart send length slen pident qcovs sstrand gaps evalue bitscore score' -out ./${sample_id}_${assembler}_snapgene_plasmids_blastn.tsv -evalue 1e-5 -num_threads 32 -perc_identity 80 -max_target_seqs 10
    blastn -query ${contigs_fasta} -db ${params.vfdb_nt} -outfmt '6 qaccver qlen qstart qend saccver stitle sstart send length slen pident qcovs sstrand gaps evalue bitscore score' -out ./${sample_id}_${assembler}_vfdb_nt_blastn.tsv -evalue 1e-5 -num_threads 32 -perc_identity 90 -culling_limit 1 -max_target_seqs 1000
    blastn -query ${contigs_fasta} -db ${params.AMR_CDS} -outfmt '6 qaccver qlen qstart qend saccver stitle sstart send length slen pident qcovs sstrand gaps evalue bitscore score' -out ./${sample_id}_${assembler}_AMR_CDS_blastn_80threshold.tsv -evalue 1e-5 -num_threads 32 -perc_identity 80 -max_target_seqs 10
    blastn -query ${contigs_fasta} -db ${params.selectAgents_nuc} -outfmt '6 qaccver qlen qstart qend saccver stitle sstart send length slen pident qcovs sstrand gaps evalue bitscore score' -out ./${sample_id}_${assembler}_select_agents_nuc_blastn.tsv -evalue 1e-5 -num_threads 32 -perc_identity 80 -max_target_seqs 10
    blastp -query ${contigs_fasta} -db ${params.vfdb_prot} -outfmt '6 qaccver qlen qstart qend saccver stitle sstart send length slen pident qcovs sstrand gaps evalue bitscore score' -out ./${sample_id}_${assembler}_vfdb_prot.tsv -evalue 1e-20 -num_threads 32 -culling_limit 1
    blastp -query ${contigs_fasta} -db ${params.AMR_Prot} -outfmt '6 qaccver qlen qstart qend saccver stitle sstart send length slen pident qcovs sstrand gaps evalue bitscore score' -out ./${sample_id}_${assembler}_AMRProt.tsv -evalue 1e-20 -num_threads 32 -culling_limit 1
    blastp -query ${contigs_fasta} -db ${params.selectAgents_prot} -outfmt '6 qaccver qlen qstart qend saccver stitle sstart send length slen pident qcovs sstrand gaps evalue bitscore score' -out ./${sample_id}_${assembler}_select_agents_toxins_prot.tsv -evalue 1e-20 -num_threads 32 -culling_limit 1
    """
}

process AMR_VF_BLASTN_UNMAPPED {
    tag { "${assembler} | ${sample_id}" }
    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/unmapped_read_analysis/AMR_VF/${assembler}" }, mode: 'copy'
    label 'lowmem'
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), val(assembler), path(contigs_fasta)

    output:
    tuple val(sample_id), val(assembler), path("*.tsv"), emit: amr_vf_blastn_ch

    when:
    params.amr_vf

    script:
    """
    blastn -query ${contigs_fasta} -db ${params.addgeneDB} \
        -outfmt '6 qaccver qlen qstart qend saccver stitle sstart send length slen pident qcovs sstrand gaps evalue bitscore score' \
        -out ${sample_id}_${assembler}_addgene_plasmids_blastn.tsv \
        -evalue 1e-5 -num_threads 1 -perc_identity 80 -max_target_seqs 10

    blastn -query ${contigs_fasta} -db ${params.snapgeneDB} \
        -outfmt '6 qaccver qlen qstart qend saccver stitle sstart send length slen pident qcovs sstrand gaps evalue bitscore score' \
        -out ${sample_id}_${assembler}_snapgene_plasmids_blastn.tsv \
        -evalue 1e-5 -num_threads 1 -perc_identity 80 -max_target_seqs 10

    blastn -query ${contigs_fasta} -db ${params.vfdb_nt} \
        -outfmt '6 qaccver qlen qstart qend saccver stitle sstart send length slen pident qcovs sstrand gaps evalue bitscore score' \
        -out ${sample_id}_${assembler}_vfdb_nt_blastn.tsv \
        -evalue 1e-5 -num_threads 1 -perc_identity 90 -culling_limit 1 -max_target_seqs 1000

    blastn -query ${contigs_fasta} -db ${params.AMR_CDS} \
        -outfmt '6 qaccver qlen qstart qend saccver stitle sstart send length slen pident qcovs sstrand gaps evalue bitscore score' \
        -out ${sample_id}_${assembler}_AMR_CDS_blastn_80threshold.tsv \
        -evalue 1e-5 -num_threads 1 -perc_identity 80 -max_target_seqs 10

    blastn -query ${contigs_fasta} -db ${params.selectAgents_nuc} \
        -outfmt '6 qaccver qlen qstart qend saccver stitle sstart send length slen pident qcovs sstrand gaps evalue bitscore score' \
        -out ${sample_id}_${assembler}_select_agents_nuc_blastn.tsv \
        -evalue 1e-5 -num_threads 1 -perc_identity 80 -max_target_seqs 10
    """
}

process MLST {
    tag { "${assembler} | ${sample_id}" }
    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/MLST/${assembler}" }, mode: 'copy'
    label 'lowmem'
    errorStrategy 'ignore'
    input:
    tuple val(sample_id), val(assembler), path(contigs_fasta)

    output:
    tuple val(sample_id), val(assembler), file("${sample_id}_${assembler}_mlst.csv"), emit: mlst_ch

    when:
    (params.mlst)

    script:
    """
    mlst --threads ${task.cpus} --quiet --nopath --full ${contigs_fasta} --csv --outfile ${sample_id}_${assembler}_mlst.csv
    """
}

process GE_SCREEN {
    tag { "${assembler} | ${sample_id}" }
    publishDir { "${params.outdir}/${params.run_id}/GE_screen/${assembler}" }, mode: 'copy'
    label 'lowmem'
    errorStrategy 'ignore'
    input:
    tuple val(sample_id), val(assembler), path(assembly_dir), path(reads_dir)

    output:
    tuple val(sample_id), val(assembler), file("*"), emit: ge_screen_ch

    when:
    (params.ge_screen)

    script:
    """
    mkdir -p assemblies reads

    bash ${params.scripts}/run_ge_screen.sh \\
        --assembly-dir ${assembly_dir} \\
        --reads-dir ${reads_dir} \\
        --outdir ge_screen_results \\
        --amrfinderplus \\
        --amrfinder-db ${params.amrfinderDB} \\
        --amrfinder-annotation-format prodigal \\
        --partition normal \\
        --threads ${task.cpus} \\
        --submit-dashboard \\ 
        --dashboard-wrapper ${params.scripts}/submit_ge_screen_dashboard.sh \\
        --dashboard-dependency afterok \\
        --dashboard-embed-plotly \\
        --dashboard-max-nt-rows 50 \\
        --dashboard-max-univec-rows 50 \\
        --dashboard-max-coverage-rows 100 \\
        --dashboard-max-amrfinder-rows 100
    """
}

process RGI {
    tag { "${assembler} | ${sample_id}" }
    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/rgi/${assembler}" }, mode: 'copy'
    label 'lowmem'
    errorStrategy 'ignore'
    input:
    tuple val(sample_id), val(assembler), path(fna), path(faa)

    output:
    tuple val(sample_id), val(assembler), file("*"), emit: rgi_ch

    when:
    (params.rgi)

    script:
    """
    ${params.rgi_load_cmd}

    rgi main -i ${fna} -o ${sample_id}_${assembler}_rgi_nuc -t contig --clean -n ${task.cpus} --include_loose
    rgi main -i ${faa}  -o ${sample_id}_${assembler}_rgi_prot -t protein --clean -n ${task.cpus}
    """
}

process RGI_UNMAPPED {
    tag { "${assembler} | ${sample_id}" }
    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/rgi/${assembler}" }, mode: 'copy'
    label 'lowmem'
    errorStrategy 'ignore'
    input:
    tuple val(sample_id), val(assembler), path(fna)

    output:
    tuple val(sample_id), val(assembler), file("*"), emit: rgi_ch

    when:
    (params.rgi)

    script:
    """
    ${params.rgi_load_cmd}

    rgi main -i ${fna} -o ${sample_id}_${assembler}_rgi_nuc -t contig --clean -n ${task.cpus} --local
    """
}


process AMRFINDER {
    tag { "${assembler} | ${sample_id}" }
    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/amrfinder/${assembler}" }, mode: 'copy'
    label 'lowmem'
    errorStrategy 'ignore'
    input:
    tuple val(sample_id), val(assembler), path(fna), path(faa), path(gff), path(gbk)

    output:
    tuple val(sample_id), val(assembler), file("*"), emit: amrfinder_ch

    when:
    (params.amrfinder)

    script:
    """
    amrfinder --plus \\
      -p ${faa} \\
      -g ${gff} \\
      -n ${fna} \\
      -a prokka \\
      --name ${sample_id}_${assembler} \\
      -o ${sample_id}_${assembler}_amrfinder.tsv \\
      -d ${params.amrfinderDB} \\
      --threads ${task.cpus}
    """
}

process AMRFINDER_UNMAPPED {
    tag { "${assembler} | ${sample_id}" }
    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/amrfinder/${assembler}" }, mode: 'copy'
    label 'lowmem'
    errorStrategy 'ignore'
    input:
    tuple val(sample_id), val(assembler), path(fna)

    output:
    tuple val(sample_id), val(assembler), file("*"), emit: amrfinder_ch

    when:
    (params.amrfinder)

    script:
    """
    amrfinder --plus \\
      -n ${fna} \\
      --name ${sample_id}_${assembler} \\
      -o ${sample_id}_${assembler}_amrfinder.tsv \\
      -d ${params.amrfinderDB} \\
      --threads 1 
    """
}

process PLASME {
    tag { "${assembler} | ${sample_id}" }
    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/plasme/${assembler}" }, mode: 'copy'
    label 'lowmem'
    errorStrategy 'ignore'
    input:
    tuple val(sample_id), val(assembler), path(contigs_fasta)

    output:
    tuple val(sample_id), val(assembler), file("*"), emit: plasme_detection_ch

    when:
    (params.run_plasme)

    script:
    """
    python ${params.scripts}/PLASMe.py ${contigs_fasta} ${sample_id}_${assembler}_${params.plasme_mode}.out -d ${params.plasme_db} -m ${params.plasme_mode}
    """
}

process PHISPY {
    tag { "${assembler} | ${sample_id}" }
    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/phispy/${assembler}" }, mode: 'copy'
    label 'lowmem'
    errorStrategy 'ignore'
    input:
    tuple val(sample_id), val(assembler), path(gbk)

    output:
    tuple val(sample_id), val(assembler), file("*"), emit: phispy_ch

    when:
    (params.run_phispy)

    script:
    """
    PhiSpy.py --threads ${task.cpus} --output_choice 11 -o ./phispy -p ${sample_id} ${gbk} 
    """
}

process MOBSUITE {
    tag { "${assembler} | ${sample_id}" }
    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/mobsuite/${assembler}" }, mode: 'copy'
    label 'no_fips'
    errorStrategy 'ignore'
    input:
    tuple val(sample_id), val(assembler), path(contigs_fasta)

    output:
    tuple val(sample_id), val(assembler), file("*"), emit: mobsuite_ch

    when:
    (params.run_mobsuite)

    script:
    """
    mob_recon --force --infile ${contigs_fasta} --outdir ./mobsuite -c -s ${sample_id} -p ${sample_id}
    """
}

process AMR_VF_BLAST_plasmids {
    tag { "${assembler} | ${sample_id}" }
    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/AMR_VF/${assembler}" }, mode: 'copy'
    label 'lowmem'
    conda "${baseDir}/env/amrfinderplus.yml"
    errorStrategy 'ignore'
    input:
    tuple val(sample_id), val(assembler), path(contigs_fasta)

    output:
    tuple val(sample_id), val(assembler), file("*"), emit: amr_vf_plasmids_ch

    when:
    (params.amr_vf)

    script:
    """
    blastn -query ${contigs_fasta} -db ${params.addgeneDB} -outfmt '6 qaccver qlen qstart qend saccver stitle sstart send length slen pident qcovs sstrand gaps evalue bitscore score' -out ./${sample_id}_${assembler}_addgene_plasmids_blastn.tsv -evalue 1e-5 -num_threads 32 -perc_identity 80 -max_target_seqs 10
    blastn -query ${contigs_fasta} -db ${params.snapgeneDB} -outfmt '6 qaccver qlen qstart qend saccver stitle sstart send length slen pident qcovs sstrand gaps evalue bitscore score' -out ./${sample_id}_${assembler}_snapgene_plasmids_blastn.tsv -evalue 1e-5 -num_threads 32 -perc_identity 80 -max_target_seqs 10

    sed -i "1i Query_ID\tQuery_length\tqstart\tqend\tSubject_accession\tTitle\tsstart\tsend\taln_length\tslen\tpercent_id\tTotal_Query_Cov\tsstrand\tgaps\tevalue\tbitscore\tscore" ${sample_id}_${assembler}_addgene_plasmids_blastn.tsv
    sed -i "1i Query_ID\tQuery_length\tqstart\tqend\tSubject_accession\tTitle\tsstart\tsend\taln_length\tslen\tpercent_id\tTotal_Query_Cov\tsstrand\tgaps\tevalue\tbitscore\tscore" ${sample_id}_${assembler}_snapgene_plasmids_blastn.tsv
    """
}

process AMR_VF_BLAST_select_agents {
    tag { "${assembler} | ${sample_id}" }
    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/AMR_VF/${assembler}" }, mode: 'copy'
    label 'lowmem'
    conda "${baseDir}/env/amrfinderplus.yml"
    errorStrategy 'ignore'
    input:
    tuple val(sample_id), val(assembler), path(contigs_fasta)

    output:
    tuple val(sample_id), val(assembler), file("*"), emit: amr_vf_select_ch

    when:
    (params.amr_vf)

    script:
    """
    blastn -query ${contigs_fasta} -db ${params.selectAgents_nuc} -outfmt '6 qaccver qlen qstart qend saccver stitle sstart send length slen pident qcovs sstrand gaps evalue bitscore score' -out ./${sample_id}_${assembler}_select_agents_nuc_blastn.tsv -evalue 1e-5 -num_threads 32 -perc_identity 80 -max_target_seqs 10
    blastp -query ${contigs_fasta} -db ${params.selectAgents_prot} -outfmt '6 qaccver qlen qstart qend saccver stitle sstart send length slen pident qcovs sstrand gaps evalue bitscore score' -out ./${sample_id}_${assembler}_select_agents_toxins_prot.tsv -evalue 1e-20 -num_threads 32 -culling_limit 1

    sed -i "1i Query_ID\tQuery_length\tqstart\tqend\tSubject_accession\tTitle\tsstart\tsend\taln_length\tslen\tpercent_id\tTotal_Query_Cov\tsstrand\tgaps\tevalue\tbitscore\tscore" ${sample_id}_${assembler}_select_agents_nuc_blastn.tsv
    sed -i "1i Query_ID\tQuery_length\tqstart\tqend\tSubject_accession\tTitle\tsstart\tsend\taln_length\tslen\tpercent_id\tTotal_Query_Cov\tsstrand\tgaps\tevalue\tbitscore\tscore" ${sample_id}_${assembler}_select_agents_toxins_prot.tsv
    """
}

process AMR_VF_BLAST_AMR {
    tag { "${assembler} | ${sample_id}" }
    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/AMR_VF/${assembler}" }, mode: 'copy'
    label 'lowmem'
    conda "${baseDir}/env/amrfinderplus.yml"
    errorStrategy 'ignore'
    input:
    tuple val(sample_id), val(assembler), path(contigs_fasta)

    output:
    tuple val(sample_id), val(assembler), file("*"), emit: amr_vf_amr_ch

    when:
    (params.amr_vf)

    script:
    """
    blastn -query ${contigs_fasta} -db ${params.AMR_CDS} -outfmt '6 qaccver qlen qstart qend saccver stitle sstart send length slen pident qcovs sstrand gaps evalue bitscore score' -out ./${sample_id}_${assembler}_AMR_CDS_blastn_80threshold.tsv -evalue 1e-5 -num_threads 32 -perc_identity 80 -max_target_seqs 10
    blastp -query ${contigs_fasta} -db ${params.AMR_Prot} -outfmt '6 qaccver qlen qstart qend saccver stitle sstart send length slen pident qcovs sstrand gaps evalue bitscore score' -out ./${sample_id}_${assembler}_AMRProt.tsv -evalue 1e-20 -num_threads 32 -culling_limit 1

    sed -i "1i Query_ID\tQuery_length\tqstart\tqend\tSubject_accession\tTitle\tsstart\tsend\taln_length\tslen\tpercent_id\tTotal_Query_Cov\tsstrand\tgaps\tevalue\tbitscore\tscore" ./${sample_id}_${assembler}_AMR_CDS_blastn_80threshold.tsv
    sed -i "1i Query_ID\tQuery_length\tqstart\tqend\tSubject_accession\tTitle\tsstart\tsend\taln_length\tslen\tpercent_id\tTotal_Query_Cov\tsstrand\tgaps\tevalue\tbitscore\tscore" ./${sample_id}_${assembler}_AMRProt.tsv
    """
}
process AMR_VF_BLAST_VF {
    tag { "${assembler} | ${sample_id}" }
    publishDir { "${params.outdir}/${params.run_id}/${sample_id}/AMR_VF/${assembler}" }, mode: 'copy'
    label 'lowmem'
    conda "${baseDir}/env/amrfinderplus.yml"
    errorStrategy 'ignore'
    input:
    tuple val(sample_id), val(assembler), path(contigs_fasta)

    output:
    tuple val(sample_id), val(assembler), file("*"), emit: amr_vf_vf_ch

    when:
    (params.amr_vf)

    script:
    """
    blastn -query ${contigs_fasta} -db ${params.vfdb_nt} -outfmt '6 qaccver qlen qstart qend saccver stitle sstart send length slen pident qcovs sstrand gaps evalue bitscore score' -out ./${sample_id}_${assembler}_vfdb_nt_blastn.tsv -evalue 1e-5 -num_threads 32 -perc_identity 90 -culling_limit 1 -max_target_seqs 1000
    blastp -query ${contigs_fasta} -db ${params.vfdb_prot} -outfmt '6 qaccver qlen qstart qend saccver stitle sstart send length slen pident qcovs sstrand gaps evalue bitscore score' -out ./${sample_id}_${assembler}_vfdb_prot.tsv -evalue 1e-20 -num_threads 32 -culling_limit 1

    sed -i "1i Query_ID\tQuery_length\tqstart\tqend\tSubject_accession\tTitle\tsstart\tsend\taln_length\tslen\tpercent_id\tTotal_Query_Cov\tsstrand\tgaps\tevalue\tbitscore\tscore" ${sample_id}_${assembler}_vfdb_nt_blastn.tsv
    sed -i "1i Query_ID\tQuery_length\tqstart\tqend\tSubject_accession\tTitle\tsstart\tsend\taln_length\tslen\tpercent_id\tTotal_Query_Cov\tsstrand\tgaps\tevalue\tbitscore\tscore" ${sample_id}_${assembler}_vfdb_prot.tsv
    """
}
