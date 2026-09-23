#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================
		            Readmapping and Variant Calling Workflow
========================================================================================
*/
// Include modules
include { Fasta_preprocess } from './modules/local/readmapping_auto/main.nf'
include { BBMap } from './modules/local/readmapping_auto/main.nf'
include { Minimap } from './modules/local/readmapping_auto/main.nf'
include { Minimap_Hyb } from './modules/local/readmapping_auto/main.nf'
include { Minimap_sr } from './modules/local/readmapping_auto/main.nf'
include { Filter_bams } from './modules/local/readmapping_auto/main.nf'
include { Individual_Covstats } from './modules/local/readmapping_auto/main.nf'
include { Combine_Covstats } from './modules/local/readmapping_auto/main.nf'
include { PileUp } from './modules/local/readmapping_auto/main.nf'
include { IVar } from './modules/local/readmapping_auto/main.nf'



workflow {

workflow CONTIG_based_analysis {

    /*
    TAKE:
    readmapping_input_ch : (sid, fq1, fq2, lr, reference, gff3, interleave)
    EMITS:
    --
     */

    take:
    readmapping_input_ch

    main:



    def sample_info_ch = readmapping_input_ch.map { sid, fq1, fq2, lr, reference, gff3, interleave ->
        def mode = fastq_1 && fastq_2 && long_read ? 'hybrid'
                     : fastq_1 && fastq_2	? 'short'
                     : long_read && interleave_flag ? 'short_interleaved'
                     : long_read               ? 'long'
                     : 'none'
        tuple(sample_id, fastq_1, fastq_2, long_read, reference, gff3, mode)
    }.filter { it[4] != 'none' }

    def outdir=.    
    
    sample_info_ch.view { "Raw sample data: $it" }

    Fasta_preprocess(sample_info_ch)

    short_ch = Fasta_preprocess.out.cleaned_fasta_ch.filter { sample_id, fastq_1, fastq_2, long_read, cleaned_ref, gff3, mode ->
      mode == 'short'
    }
    long_ch = Fasta_preprocess.out.cleaned_fasta_ch.filter { sample_id, fastq_1, fastq_2, long_read, cleaned_ref, gff3, mode ->
      mode == 'long'
    }
    hybrid_ch = Fasta_preprocess.out.cleaned_fasta_ch.filter { sample_id, fastq_1, fastq_2, long_read, cleaned_ref, gff3, mode ->
      mode == 'hybrid'
    }
    short_in_ch = Fasta_preprocess.out.cleaned_fasta_ch.filter { sample_id, fastq_1, fastq_2, long_read, cleaned_ref, gff3, mode ->
      mode == 'short_interleaved'
    }
    

    BBMap(short_ch)
    Minimap(long_ch)
    Minimap_Hyb(hybrid_ch)
    Minimap_sr(short_in_ch)
	
    all_bams = BBMap.out.bbmap_bam_ch.mix(Minimap.out.minimap_bam_ch, Minimap_Hyb.out.minimap_hyb_bam_ch, Minimap_sr.out.minimap_bam_ch)
    

    Filter_bams(all_bams)
    
    Individual_Covstats(all_bams)
    
    collected_covstats = Individual_Covstats.out.covstats_ch.collect(flat: false).map {it.transpose() }
    // collected_covstats.view { "$it" }
    Combine_Covstats(collected_covstats, Fasta_preprocess.out.header_ch)

  
    if (params.run_variant) {
      PileUp(Filter_bams.out.filtered_bam_ch)
      IVar(PileUp.out.pileup_output_ch)
    } 
  
}



