#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Imports
include { Map_Reads_2_Contigs as Map_Reads_2_Contigs_VS } from './modules/local/virusseeker/main.nf'
include { Filter_Putative_Viruses } from './modules/local/virusseeker/main.nf'
include { Generate_VS_Outputs } from './modules/local/virusseeker/main.nf'
include { Merge_Arc } from './modules/local/virusseeker/main.nf'
include { Heatmap_Virusseeker } from './modules/local/virusseeker/main.nf'

workflow VS_Workflow {

    /*
    TAKE:
    vs_calls_input_ch : (sid, fq1, fq2, lr, contigs, mmseqs_parsed, assembler)

    EMITS:
    - (from Generate_VS_Outputs) vs_reports_ch, accurate_read_counts_ch, accurate_read_counts_dir_ch
    - (from Merge_Arc) total_csv, nfam_csv, nrpm_csv
    - (from Heatmap_Virusseeker) heatmap, heatmap_scaled_by_family, heatmap_scaled_by_sample
    - (from mapper) map2assembly_ch_sr, map2assembly_ch_lr
    */
    take:
    vs_calls_input_ch

    main:

    def global_mode =
        params.hybrid ? 'hybrid' :
        params.shortreads ? 'short' :
        (params.longreads ? 'long' : null)

    def pickCpus = { size_gb -> size_gb < 1 ? 2 : size_gb < 5 ? 4 : size_gb < 10 ? 8 : 16 }
    def pickMem = { size_gb -> size_gb < 1 ? '8 GB' : size_gb < 5 ? '16 GB' : size_gb < 10 ? '32 GB' : '64 GB' }

    def asm_key = vs_calls_input_ch.map { sid, fq1, fq2, lr, contigs, parsed, assembler ->
        tuple(sid, assembler) // (sid, assembler)
    }

    def all_map_in_tagged = vs_calls_input_ch.map { sid, fq1, fq2, lr, contigs, parsed, assembler ->
        long s1 = fq1 ? fq1.size() : 0L
        long s2 = fq2 ? fq2.size() : 0L
        double sr_gb = (s1 + s2) / 1e9
        double lr_gb = lr ? (lr.size() / 1e9) : 0.0
    def mode = global_mode ?: (
        ( (fq1 || fq2) && lr ) ? 'hybrid' :
        ( (fq1 || fq2) ) ? 'short' :
        ( lr ) ? 'long' : 'short'
    )
    def cpus = pickCpus(Math.max(sr_gb, lr_gb))
    def mem = pickMem(Math.max(sr_gb, lr_gb))
        tuple(sid, fq1, fq2, lr, contigs, mode, cpus, mem, assembler)
    }

    //def all_map_in_plain = all_map_in_tagged.map { sid, fq1, fq2, lr, contigs, mode, cpus, mem, assembler ->
    //    tuple(sid, fq1, fq2, lr, contigs, mode, cpus, mem)
    //}

    Map_Reads_2_Contigs_VS( all_map_in_tagged )

    def sr_stats_tagged = asm_key
        .join( Map_Reads_2_Contigs_VS.out.map2assembly_ch_sr ) // (sid, assembler, pe_covstats, pe_qc_reads)
        .map { sid, assembler, pe_covstats, pe_qc_reads -> tuple(assembler, sid, pe_covstats, pe_qc_reads) }

    def lr_stats_tagged = asm_key
        .join( Map_Reads_2_Contigs_VS.out.map2assembly_ch_lr )
        .map { sid, assembler, lr_covstats, lr_qc_reads -> tuple(assembler, sid, lr_covstats, lr_qc_reads) }

    def contigs_by_assembler = all_map_in_tagged
        .map { sid, fq1, fq2, lr, contigs, mode, cpus, mem, assembler ->
        tuple(assembler, sid, contigs) // (asm, sid, contigs)
    }

    def mode_by_asm = all_map_in_tagged
        .map { sid, fq1, fq2, lr, contigs, mode, cpus, mem, assembler ->
        tuple(assembler, sid, mode) // (asm, sid, mode)
    }

    // Filter_Putative_Viruses runs once per (sid, assembler) mmseqs_parsed
    //def parsed_plain = vs_calls_input_ch.map { sid, fq1, fq2, lr, contigs, parsed, assembler ->
    //    tuple(sid, parsed)
    //}
    def parsed_tagged = vs_calls_input_ch.map { sid, fq1, fq2, lr, contigs, parsed, assembler ->
        tuple(sid, parsed, assembler ) 
        }
    Filter_Putative_Viruses( parsed_tagged )

    def viral_tagged = Filter_Putative_Viruses.out.viral_blast_tab_ch
        .map { sid, viral_blast_tab, assembler -> tuple(assembler, sid, viral_blast_tab) }
    
    def viral_k = viral_tagged // (asm, sid, viral)
        .map { asm, sid, viral -> tuple([asm, sid], viral) }

    def contigs_k = contigs_by_assembler // (asm, sid, contigs)
        .map { asm, sid, contigs -> tuple([asm, sid], contigs) }

    def mode_k = mode_by_asm // (asm, sid, mode)
        .map { asm, sid, mode -> tuple([asm, sid], mode) }

    // Stats (real)
    def sr_real_k = sr_stats_tagged // (asm, sid, cov, qc)
        .map { asm, sid, cov, qc -> tuple([asm, sid], [cov, qc]) }

    def lr_real_k = lr_stats_tagged
        .map { asm, sid, cov, qc -> tuple([asm, sid], [cov, qc]) }

    def all_keys = viral_tagged
        .map { asm, sid, _ -> tuple([asm, sid], true) }
        .mix( contigs_by_assembler.map { asm, sid, _ -> tuple([asm, sid], true) } )
        .groupTuple()
        .map { key, _ -> key } // emits [asm, sid]

    def sr_default_k = all_keys.map { key -> tuple(key, [null, null]) }
    def lr_default_k = all_keys.map { key -> tuple(key, [null, null]) }

    def sr_any_k = sr_default_k.mix(sr_real_k) // ([asm,sid], [cov,qc] possibly nulls)
        .groupTuple()
        .map { key, pairs ->
    def got = pairs.find { it[0] != null || it[1] != null } ?: pairs[0]
        tuple(key, got)
    }

    def lr_any_k = lr_default_k.mix(lr_real_k)
        .groupTuple()
        .map { key, pairs ->
    def got = pairs.find { it[0] != null || it[1] != null } ?: pairs[0]
        tuple(key, got)
    }

    def gen_vs_all = viral_k
        .join(contigs_k) // ([asm,sid], viral, contigs)
        .join(sr_any_k) // ([asm,sid], viral, contigs, [pe_cov,pe_qc])
        .join(lr_any_k) // ([asm,sid], viral, contigs, [pe_cov,pe_qc], [lr_cov,lr_qc])
        .join(mode_k) // ([asm,sid], viral, contigs, [pe_cov,pe_qc], [lr_cov,lr_qc], mode)
        .map { key, viral_tab, contigs_fa, pe_pair, lr_pair, mode ->
    def (asm, sid) = key
    def (pe_cov, pe_qc) = pe_pair
    def (lr_cov, lr_qc) = lr_pair
        tuple(sid, viral_tab, contigs_fa, pe_cov, pe_qc, lr_cov, lr_qc, mode, asm)
    }
    
    Generate_VS_Outputs( gen_vs_all )
    
    //def merged_arc_input_ch = Generate_VS_Outputs.out.accurate_read_counts_dir_ch.collect()
    //  	.map{ sid, tsv, asm -> tsv } .collect()    

    //merged_arc_input_ch = Generate_VS_Outputs.out.accurate_read_counts_dir_ch.collect()
    //merged_arc_input_ch.view { "merged_arc_input_ch: $it" }
    //def arc_lists_by_asm = Generate_VS_Outputs.out.accurate_read_counts_file_ch.collect()
    //	.groupTuple()
    //def arc_files_list = Generate_VS_Outputs.out.accurate_read_counts_file_ch
    //  	.map { asm, p -> p }
    //	.collect()
    //def arc_files_list = arc_files.collect()
   
    //arc_files_list.view { "merged_arc_input_ch: $it" }
    //def arc_files_by_asm = Generate_VS_Outputs.out.accurate_read_counts_file_ch
    //	.map { sid, tsv, asm -> tuple( asm, tsv ) }
    //def arc_lists_by_asm = arc_files_by_asm.groupTuple()
   
    //def arc_pairs_all = Generate_VS_Outputs.out.accurate_read_counts_file_ch.collect()
    //arc_pairs_all.view { "merged_arc_input_ch: $it" }
    //def arc_files_list = Generate_VS_Outputs.out.accurate_read_counts_file_ch
    //	.map { asm, p -> p }
    //	.collect()
    def arc_lists_by_asm = Generate_VS_Outputs.out.accurate_read_counts_file_ch.groupTuple()
    Merge_Arc( arc_lists_by_asm )
    
    def heatmap_in_ch = Merge_Arc.out.total_csv
        .mix(Merge_Arc.out.nfam_csv)
        .mix(Merge_Arc.out.nrpm_csv)
        .map { f -> tuple(f.baseName, f) }

    Heatmap_Virusseeker( heatmap_in_ch )

    emit:
    map2assembly_ch_sr = Map_Reads_2_Contigs_VS.out.map2assembly_ch_sr
    map2assembly_ch_lr = Map_Reads_2_Contigs_VS.out.map2assembly_ch_lr

    vs_reports_ch = Generate_VS_Outputs.out.vs_reports_ch
    accurate_read_counts_ch = Generate_VS_Outputs.out.accurate_read_counts_ch
    accurate_read_counts_dir_ch= Generate_VS_Outputs.out.accurate_read_counts_file_ch
    total_csv = Merge_Arc.out.total_csv
    nfam_csv = Merge_Arc.out.nfam_csv
    nrpm_csv = Merge_Arc.out.nrpm_csv
    heatmap = Heatmap_Virusseeker.out.heatmap
    heatmap_scaled_by_family = Heatmap_Virusseeker.out.heatmap_scaled_by_family
    heatmap_scaled_by_sample = Heatmap_Virusseeker.out.heatmap_scaled_by_sample
}
