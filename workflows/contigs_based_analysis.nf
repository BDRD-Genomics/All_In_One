#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Imports
include { BLASTN_NT_contigs;PROKKA;BUSCO;AMR_VF_BLAST;MLST;CHIMERIC_DETECTION;RGI;AMRFINDER;PLASME;PHISPY;MOBSUITE } from './modules/local/contig_characterization/main_contigs.nf'
include { AMR_VF_BLAST_plasmids;AMR_VF_BLAST_select_agents;AMR_VF_BLAST_AMR;AMR_VF_BLAST_VF } 			     from './modules/local/contig_characterization/main_contigs.nf'
include { AMRFINDER_UNMAPPED; UNMAPPED_LR_TO_FASTA; RGI_UNMAPPED; AMR_VF_BLASTN_UNMAPPED } 			     from './modules/local/contig_characterization/main_contigs.nf'
include { CheckM_Assemblies;CheckM2_Assemblies;CheckV_Assemblies }                                   	    	     from './modules/local/checkm/main.nf'
include { Map_Reads_2_Contigs }                                                                      		     from './modules/local/map2refseq/main.nf'
//include { Assembly_Workflow as Unmapped_Assembly_Workflow }                                          		     from './assembly.nf'
include {UNMAPPED_CHIMERIC_DETECTION}                                                                                         from './modules/local/contig_characterization/main_contigs.nf'

def contigPickCpus(size_mb) {
    size_mb < 400 ? 12 : size_mb < 800 ? 24 : size_mb < 2000 ? 64 : 128
}

def contigPickMem(size_mb) {
    size_mb < 400 ? '24 GB' : size_mb < 800 ? '48 GB' : size_mb < 2000 ? '128 GB' : '256 GB'
}

workflow CONTIG_based_analysis {

    /*
    TAKE:
    contigs_based_analysis_ch : (sid, fq1, fq2, lr, contigs, mode, assembler)

    EMITS:
    - (from PROKKA) prokka_ch
    - (from BUSCO) busco_ch
    - (from AMR_VF_BLAST) amr_vf_blast_ch
    - (from MLST) mlst_ch
    - (from RGI) rgi_ch
    - (from AMR_FINDER) amrfinder_ch
    - (from CHIMERIC_DETECTION) chimeric_detection_ch
     */

    take:
    contigs_based_analysis_ch

    main:

    def global_mode =
        params.hybrid ? 'hybrid' :
        params.shortreads ? 'short' :
        (params.longreads ? 'long' : null)

    def asm_key = contigs_based_analysis_ch.map { sid, fq1, fq2, lr, contigs, mode2, assembler ->
        tuple(sid, assembler)
    }

    sample_ch = asm_key.map { sid, assembler -> 
		sid
    }
    def sample_var = sample_ch.first()

    def all_tagged = contigs_based_analysis_ch.map { sid, fq1, fq2, lr, contigs, mode2, assembler ->
        double contigs_mb = contigs ? (contigs.size() / 1e6) : 0.0
        def mode = global_mode
        def cpus = contigPickCpus(contigs_mb)
        def mem = contigPickMem(contigs_mb)
        tuple(sid, fq1, fq2, lr, contigs, assembler, mode, cpus, mem)
    }



    def sample_asm_contigs = all_tagged.map { sid, fq1, fq2, lr, contigs, assembler, mode, cpus, mem ->
        tuple(sid, assembler, contigs)
    }

    BLASTN_NT_contigs(sample_asm_contigs)

    PROKKA(sample_asm_contigs)
    def prokka_out_ch = PROKKA.out.prokka_ch
    //prokka_out_ch.view{ "prokka_output_ch: $it" } //sid, assembler, fna, faa, gff
    CHIMERIC_DETECTION(sample_asm_contigs)
    def chimeric_detection_out_ch = CHIMERIC_DETECTION.out.chimeric_detection_ch

    def map2contigs_ch = all_tagged.map { sid, fq1, fq2, lr, contigs, assembler, mode, cpus, mem ->
        tuple(sid, fq1, fq2, lr, contigs, mode, assembler)}
    Map_Reads_2_Contigs(map2contigs_ch)

    //def map2contigs_ch = all_tagged.map {
    //   sid, fq1, fq2, lr, contigs, assembler, mode, cpus, mem ->

    //    tuple(
    //        sid,
    //        fq1 ?: [],
    //        fq2 ?: [],
    //        lr  ?: [],
    //       contigs,
    //        mode,
     //       assembler
     //   )
    //}

    map2contigs_ch.view {
        "Map_Reads_2_Contigs input: ${it}"
    }

    //Map_Reads_2_Contigs(map2contigs_ch)


    def unmapped_prokka_out_ch      = Channel.empty()
    def unmapped_amr_vf_out_ch      = Channel.empty()
    def unmapped_rgi_out_ch         = Channel.empty()
    def unmapped_amrfinder_out_ch   = Channel.empty()
    def unmapped_contigs_out_ch     = Channel.empty()
    def unmapped_chimeric_detection = Channel.empty()

        if (params.test_unmapped_reassembly) {


	def unmapped_long_reads_ch =
    	Map_Reads_2_Contigs.out.unmap2assembly_ch
        	.filter {
            	sid,
            	unmapped_lr,
            	source_assembler,
            	mode ->

            	mode in ['long', 'hybrid']
        	}
        	.map {
            	sid,
            	unmapped_lr,
            	source_assembler,
            	mode ->
                     tuple(
               	          sid,
                          source_assembler,
                          unmapped_lr,
                          700,
                          mode
                )
            }

	unmapped_long_reads_ch.view {
	"UNMAPPED_LR_TO_FASTA input: ${it}"
	}

	UNMAPPED_LR_TO_FASTA(unmapped_long_reads_ch)


        /*
         * Change the label so output directories clearly identify these as
         * unassembled unmapped-read sequences.
         *
         * Result:
         * sid, "unmapped_reads_from_raven", filtered_fasta
         */
        def unmapped_fasta_for_analysis_ch =
            UNMAPPED_LR_TO_FASTA.out.fasta_ch.map {
                sid, source_assembler, unmapped_fasta ->

                tuple(
                    sid,
                    "unmapped_reads_from_${source_assembler}",
                    unmapped_fasta
                )
            }

        /*
         * Both processes consume:
         * sid, label, nucleotide FASTA
         */
        //PROKKA_UNMAPPED(unmapped_fasta_for_analysis_ch)
        AMR_VF_BLASTN_UNMAPPED(unmapped_fasta_for_analysis_ch)

        /*
         * RGI consumes the nucleotide and protein FASTA files produced by
         * Prokka:
         * sid, label, FNA, FAA
         */
        /*
        def unmapped_prokka_for_rgi_ch =
            PROKKA_UNMAPPED.out.prokka_ch.map {
                sid, assembler, fna, faa, gff, gbk ->
                    tuple(sid, assembler, fna, faa)
            }
        */
        //RGI_UNMAPPED(unmapped_prokka_for_rgi_ch)
        RGI_UNMAPPED(unmapped_fasta_for_analysis_ch)
        //AMRFINDER_UNMAPPED(PROKKA_UNMAPPED.out.prokka_ch)
        AMRFINDER_UNMAPPED(unmapped_fasta_for_analysis_ch)
        UNMAPPED_CHIMERIC_DETECTION(unmapped_fasta_for_analysis_ch)
        unmapped_contigs_out_ch =
            UNMAPPED_LR_TO_FASTA.out.fasta_ch.map {
                sid, source_assembler, fasta ->
                    tuple(
                        sid,
                        "unmapped_reads_from_${source_assembler}",
                        fasta
                    )
            }

        //unmapped_prokka_out_ch =
        //    PROKKA_UNMAPPED.out.prokka_ch

        unmapped_amr_vf_out_ch =
            AMR_VF_BLASTN_UNMAPPED.out.amr_vf_blastn_ch

        unmapped_rgi_out_ch =
            RGI_UNMAPPED.out.rgi_ch
        unmapped_amrfinder_out_ch =
            AMRFINDER_UNMAPPED.out.amrfinder_ch
        unmapped_chimeric_detection =
	    UNMAPPED_CHIMERIC_DETECTION.out.unmapped_chimeric_detection_ch
    }
 
    CheckM_Assemblies(sample_asm_contigs)
    CheckM2_Assemblies(sample_asm_contigs)
    CheckV_Assemblies(sample_asm_contigs)


    BUSCO(sample_asm_contigs)
    MLST(sample_asm_contigs)
    //AMR_VF_BLAST(sample_asm_contigs)
    AMR_VF_BLAST_plasmids(sample_asm_contigs)
    AMR_VF_BLAST_select_agents(sample_asm_contigs)
    AMR_VF_BLAST_AMR(sample_asm_contigs)
    AMR_VF_BLAST_VF(sample_asm_contigs)

    rgi_ch = prokka_out_ch.map{ sid, assembler, fna, faa, gff, gbk -> tuple(sid, assembler, fna, faa) }
    rgi_ch.view{ "rgi_ch input: $it" }
    //amrfinder_ch = prokka_out_ch
    phispy_in_ch = prokka_out_ch.map{ sid, assembler, fna, faa, gff, gbk -> tuple(sid, assembler, gbk) } 


    RGI(rgi_ch)
    AMRFINDER(prokka_out_ch)
    PHISPY(phispy_in_ch)
    MOBSUITE(sample_asm_contigs)

    PLASME(sample_asm_contigs)
    //PLACEHOLDER FOR CHIMERIC DETECTION

    emit:
    busco_out_ch = BUSCO.out.busco_ch
    mlst_out_ch = MLST.out.mlst_ch
    //amr_vf_out = AMR_VF_BLAST.out.amr_vf_ch
    amr_vf_plasmids_out = AMR_VF_BLAST_plasmids.out.amr_vf_plasmids_ch
    amr_vf_select_out = AMR_VF_BLAST_select_agents.out.amr_vf_select_ch
    amr_vf_amr_out = AMR_VF_BLAST_AMR.out.amr_vf_amr_ch
    amr_vf_vf_out = AMR_VF_BLAST_VF.out.amr_vf_vf_ch
    rgi_out_ch = RGI.out.rgi_ch
    amrfinder_out_ch = AMRFINDER.out.amrfinder_ch
    blastn_contigs_out_ch = BLASTN_NT_contigs.out.blastn_contigs_nt_ch
    checkm_out_ch = CheckM_Assemblies.out.checkm_all_files_ch
    checkm2_out_ch = CheckM2_Assemblies.out.checkm2_all_files_ch
    checkv_out_ch = CheckV_Assemblies.out.checkv_all_files_ch
    map2assembly_out_ch = Map_Reads_2_Contigs.out.map2assembly_ch
    phispy_out_ch = PHISPY.out.phispy_ch
    mobsuite_out_ch = MOBSUITE.out.mobsuite_ch
    unmap2assembly_ch = Map_Reads_2_Contigs.out.unmap2assembly_ch
    unmapped_contigs_ch = unmapped_contigs_out_ch
    unmapped_prokka_ch = unmapped_prokka_out_ch
    unmapped_amr_vf_blastn_ch = unmapped_amr_vf_out_ch
    unmapped_rgi_ch = unmapped_rgi_out_ch
    unmapped_amrfinder_ch = unmapped_amrfinder_out_ch
    //chimeric_detection_ch = chimeric_detection_ch
    //unmapped_chimeric_detection = unmapped_chimeric_detection
}
