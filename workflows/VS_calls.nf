#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Fixed version generated for short-only, long-only, and hybrid VS mode handling.

// Imports
include { Map_Reads_2_Contigs as Map_Reads_2_Contigs_VS } from './modules/local/virusseeker/main.nf'
include { MMSEQ_AllReads }                                from './mmseqs/mmseqs_reads.nf'
include { MMSEQ_AllReads as MMSEQ_AllReads_LR }           from './mmseqs/mmseqs_reads.nf'
include { Generate_VS_Outputs }                           from './modules/local/virusseeker/main.nf'
include { Merge_Arc }                                     from './modules/local/virusseeker/main.nf'
include { Heatmap_Virusseeker }                           from './modules/local/virusseeker/main.nf'
include { Pull_Unassigned_MMseqs_Blastx }                 from './modules/local/virusseeker/main.nf'
include { Diamondview }                                   from './modules/local/virusseeker/main.nf'
include { Pull_and_Parse_Unmapped_Reads }                 from './modules/local/virusseeker/main.nf'
include { Pull_Unassigned_MMseqs_Reads }                  from './modules/local/virusseeker/main.nf'
include { Fastq_2_Fasta }                                 from './modules/local/virusseeker/main.nf'
include { Parse_MMseqs_Reads }                            from './modules/local/virusseeker/main.nf'
include { Filter_Putative_Viruses_ALL }                   from './modules/local/virusseeker/main.nf'

def vsPickCpus(size_gb) {
    size_gb < 1 ? 12 : size_gb < 5 ? 24 : size_gb < 10 ? 64 : 128
}

def vsPickMem(size_gb) {
    size_gb < 1 ? '24 GB' : size_gb < 5 ? '48 GB' : size_gb < 10 ? '128 GB' : '256 GB'
}

def firstPresent(value) {
    if (value == null) return null
    if (value instanceof List) return value.size() > 0 ? value[0] : null
    value
}

workflow VS_Workflow {



    take:
    vs_calls_input_ch

    main:

    def nt_db_ch = Channel.value(file(params.mmseqs_nt_db))

    def global_mode =
        params.hybrid ? 'hybrid' :
        params.shortreads ? 'short' :
        (params.longreads ? 'long' : null)

    def asm_key = vs_calls_input_ch.map { sid, fq1, fq2, lr, contigs, parsed, assembler, diamond_parsed, diamond_sr, diamond_lr ->
        tuple(sid, assembler)
    }

    def sample_ch = asm_key.map { sid, assembler ->
        sid
    }
    def sample_var = sample_ch.first()

    def all_tagged = vs_calls_input_ch.map { sid, fq1, fq2, lr, contigs, mmseqs_contigs_parsed, assembler, diamond_parsed, diamond_sr, diamond_lr ->
        long s1 = fq1 ? fq1.size() : 0L
        long s2 = fq2 ? fq2.size() : 0L
        double sr_gb = (s1 + s2) / 1e9
        double lr_gb = lr ? (lr.size() / 1e9) : 0.0
        def mode = global_mode ?: (
            ((fq1 || fq2) && lr) ? 'hybrid' :
            ((fq1 || fq2)) ? 'short' :
            (lr) ? 'long' : 'short'
        )
        def cpus = vsPickCpus(Math.max(sr_gb, lr_gb))
        def mem = vsPickMem(Math.max(sr_gb, lr_gb))
        tuple(sid, fq1, fq2, lr, contigs, mmseqs_contigs_parsed, mode, cpus, mem, assembler, diamond_parsed, diamond_sr, diamond_lr)
    }

//    def blastx_daa_ch = all_tagged.map { sid, fq1, fq2, lr, contigs, mmseqs_contigs_parsed, mode, cpus, mem, assembler, diamond_parsed, diamond_sr, diamond_lr ->
//        tuple(sid, diamond_sr, diamond_lr, mode, cpus, mem)
//    }

//    Diamondview(blastx_daa_ch)


    def blastx_daa_ch = all_tagged.flatMap { sid,
                                             fq1,
                                             fq2,
                                             lr,
                                             contigs,
                                             mmseqs_contigs_parsed,
                                             mode,
                                             cpus,
                                             mem,
                                             assembler,
                                             diamond_parsed,
                                             diamond_sr,
                                             diamond_lr ->  

        def rows = []   

        if ((mode == 'short' || mode == 'hybrid') && diamond_sr) {
            rows << tuple(
                sid,
                'short',
                diamond_sr,
                mode,
                cpus,
                mem,
                assembler
            )
        }   

        if ((mode == 'long' || mode == 'hybrid') && diamond_lr) {
            rows << tuple(
                sid,
                'long',
                diamond_lr,
                mode,
                cpus,
                mem,
                assembler
            )
        }   

        return rows
    }   

    blastx_daa_ch.view {
        "DIAMONDVIEW INPUT CLEANED: $it"
    }   

    Diamondview(blastx_daa_ch)  

    def diamondview_out_all =
        Diamondview.out.diamondview_reads_ch
            .map { sid, read_type, assembler, mode, diamondview_tsv ->
                tuple([sid, assembler, mode], tuple(read_type, diamondview_tsv))
            }
            .groupTuple()
            .map { key, items ->    

                def sid       = key[0]
                def assembler = key[1]
                def mode      = key[2]  

                def diamondview_sr = items
                    .findAll { it[0] == 'short' }
                    .collect { it[1] }  

                def diamondview_lr = items
                    .findAll { it[0] == 'long' }
                    .collect { it[1] }  

                tuple(
                    sid,
                    assembler,
                    diamondview_sr,
                    diamondview_lr
                )
            }   

    diamondview_out_all.view {
        "DIAMONDVIEW GROUPED: $it"
    }


    def mmseqs_2_blastx_ch = all_tagged.map { sid, fq1, fq2, lr, contigs, mmseqs_contigs_parsed, mode, cpus, mem, assembler, diamond_parsed, diamond_sr, diamond_lr ->
        tuple(sid, contigs, mmseqs_contigs_parsed, diamond_parsed, mode, cpus, mem, assembler)
    }

    Pull_Unassigned_MMseqs_Blastx(mmseqs_2_blastx_ch)

    def all_map_in_tagged = all_tagged.map { sid, fq1, fq2, lr, contigs, mmseqs_contigs_parsed, mode, cpus, mem, assembler, diamond_parsed, diamond_sr, diamond_lr ->
        tuple(sid, fq1, fq2, lr, contigs, mode, cpus, mem, assembler)
    }

    def params_only_ch = all_tagged.map { sid, fq1, fq2, lr, contigs, mmseqs_contigs_parsed, mode, cpus, mem, assembler, diamond_parsed, diamond_sr, diamond_lr ->
        tuple(sid, mode, cpus, mem, assembler)
    }


    Map_Reads_2_Contigs_VS(all_map_in_tagged)

    Map_Reads_2_Contigs_VS.out.map2assembly_ch_all.view{ "view all map2assembly outputs: $it"}

    unmapped_assembly_ch = Map_Reads_2_Contigs_VS.out.map2assembly_ch_all
        .map { sid, files ->
            def SR_unmapped = files.findAll {it.name.endsWith("_unmapped_SR.fastq.gz")} ?: null
            def LR_unmapped =  files.findAll {it.name.endsWith("_unmapped_LR.fastq.gz")} ?: null
            return [sid, SR_unmapped, LR_unmapped] }.transpose()


    def fq2fa_ch = unmapped_assembly_ch
        .join( params_only_ch )
        .map { sid, unmapped_fq_SR, unmapped_fq_LR, mode, cpus, mem, assembler ->
             tuple(sid, unmapped_fq_SR, unmapped_fq_LR, cpus, mem, mode, assembler)
             }

    fq2fa_ch.view { "testing new fq2fa_ch: $it" }
    Fastq_2_Fasta(fq2fa_ch)

    def mmseq_sr_input_ch = Fastq_2_Fasta.out.fastq2fasta_ch_sr
        .filter { sid, fasta -> fasta != null }
	.map { sid, fasta ->
	    tuple('short', sid, fasta)}

    def mmseq_lr_input_ch = Fastq_2_Fasta.out.fastq2fasta_ch_lr
        .map { sid, fasta ->
            tuple('long', sid, fasta)}

    def mmseqs_input_ch = mmseq_sr_input_ch
        .mix(mmseq_lr_input_ch)

    MMSEQ_AllReads(mmseqs_input_ch, nt_db_ch)

    mmseqs_all_ch   = MMSEQ_AllReads.out.mmseqs_all_reads_ch

    mmseqs_all_ch.view { "RAW MMSEQS_ALL_READS: $it" }


    assembler_by_sample_ch = params_only_ch
	.map { sid, mode, cpus, mem, assembler ->
		tuple(sid, assembler)
	}
    assembler_by_sample_ch.view { row -> "ASSEMBLER BY SAMPLE: ${row}" } 
    parse_mmseqs_ch = MMSEQ_AllReads.out.mmseqs_all_reads_ch
		.filter { row -> 
			row[2] != null &&
			row[2].toString() != 'null' &&
			row[2].toString() != 'No_Read'
		}
		.map {	row ->
		    def read_type = row[0]
		    def sample_id = row[1]
		    def mmseq_out = row[2]

		    tuple(sample_id, read_type, mmseq_out)
		}
    parse_mmseqs_ch.view{ row -> "PARSE READY: ${row}" }
    parse_mmseqs_ch.view{ "*****PARSE MMSEQS INPUT*************: $it" }  


    final_parse_mmseqs_ch = parse_mmseqs_ch
	.join(assembler_by_sample_ch)
	.map { sid, read_type, mmseq_out, assembler -> 
	     tuple( sid, assembler, read_type, mmseq_out)
	}
    final_parse_mmseqs_ch.view{ "*****WHAAAAAAAAAAAAAAAT*************: $it" }
    Parse_MMseqs_Reads( final_parse_mmseqs_ch ) // Output: tuple(sid,parsed)

    Parse_MMseqs_Reads.out.mmseqs_parsed_ch.view { "WHAT WHAT WHAT: $it" } 


    def fasta_read_files = Fastq_2_Fasta.out.fastq2fasta_ch_all
        .map { sid, files ->

            def file_list = files instanceof List ? files : [files]

            def SR_fasta = file_list
                .findAll { it != null && it.name.endsWith("_assembly_unmapped_SR.fasta") }

            def LR_fasta = file_list
                .findAll { it != null && it.name.endsWith("_assembly_unmapped_LR.fasta") }

            tuple(
                sid,
                SR_fasta,
                LR_fasta
            )
        }

    fasta_read_files.view {
        "FASTA READ FILES GROUPED: $it"
    }

//    def parse_mmseqs_out_all = Parse_MMseqs_Reads.out.mmseqs_parsed_ch
//        .map { sid, files ->
//            def file_list1 = files instanceof List ? files : [files]
//            def SR_mmparsed = file_list1.findAll {it.name.endsWith("_sr.mmseqs.parsed")} ?: null
//            def LR_mmparsed = file_list1.findAll {it.name.endsWith("_lr.mmseqs.parsed")} ?: null
//            return [sid, SR_mmparsed, LR_mmparsed] }
    def parse_mmseqs_out_all = Parse_MMseqs_Reads.out.mmseqs_parsed_ch
	.map { sid, assembler, read_type, parsed_file -> 
	    tuple([sid, assembler], tuple(read_type, parsed_file))
	}
	.groupTuple()
	.map { key, parsed_items ->
	    def sid = key[0]
	    def assembler = key[1]

	    def sr_mmseqs_parsed = parsed_items
		.findAll { item -> item[0] == 'short' }
                .collect { item -> item[1] }

            def lr_mmseqs_parsed = parsed_items
                .findAll { item -> item[0] == 'long' }
                .collect { item -> item[1] }

	    tuple(sid, assembler, sr_mmseqs_parsed, lr_mmseqs_parsed)
	}


    parse_mmseqs_out_all.view { row -> "PARSED MMSEQS GROUPED: ${row}" }


    def pull_unassigned_mmseqs_reads_ch =
        fasta_read_files
            .join(parse_mmseqs_out_all)
            .join(diamondview_out_all)
            .join(params_only_ch)
            .flatMap { sid,
                       sr_fasta,
                       lr_fasta,
                       assembler_from_parse,
                       sr_mmseqs_parsed,
                       lr_mmseqs_parsed,
                       assembler_from_diamondview,
                       diamondview_sr,
                       diamondview_lr,
                       mode,
                       cpus,
                       mem,
                       assembler ->

                if (assembler_from_parse != assembler) {
                    throw new IllegalStateException("Assembler mismatch for ${sid}: parse=${assembler_from_parse}, params=${assembler}")
                }

                if (assembler_from_diamondview != assembler) {
                    throw new IllegalStateException("Assembler mismatch for ${sid}: diamondview=${assembler_from_diamondview}, params=${assembler}")
                }

                def rows = []

                def sr_fasta_file   = firstPresent(sr_fasta)
                def lr_fasta_file   = firstPresent(lr_fasta)

                def sr_parsed_file  = firstPresent(sr_mmseqs_parsed)
                def lr_parsed_file  = firstPresent(lr_mmseqs_parsed)

                def sr_diamond_file = firstPresent(diamondview_sr)
                def lr_diamond_file = firstPresent(diamondview_lr)

                if ((mode == 'short' || mode == 'hybrid') &&
                    sr_fasta_file != null &&
                    sr_parsed_file != null &&
                    sr_diamond_file != null) {

                    rows << tuple(
                        sid,
                        'short',
                        sr_fasta_file,
                        sr_parsed_file,
                        sr_diamond_file,
                        mode,
                        cpus,
                        mem,
                        assembler
                    )
                }

                if ((mode == 'long' || mode == 'hybrid') &&
                    lr_fasta_file != null &&
                    lr_parsed_file != null &&
                    lr_diamond_file != null) {

                    rows << tuple(
                        sid,
                        'long',
                        lr_fasta_file,
                        lr_parsed_file,
                        lr_diamond_file,
                        mode,
                        cpus,
                        mem,
                        assembler
                    )
                }

                return rows
            }

    pull_unassigned_mmseqs_reads_ch.view { row -> "PULL INPUT CLEANED: ${row}" } 

    Pull_Unassigned_MMseqs_Reads( pull_unassigned_mmseqs_reads_ch ) // Output: tuple val(sample_id), file("${sample_id}_reads_blastx_unassignedMMseqs.out") unassigned_mmseqs_reads_blastx_ch

    Pull_Unassigned_MMseqs_Reads.out.unassigned_mmseqs_blastx_ch.view {"THIS IS PRE FORCELIST *?*?*?* $it"}

    def pull_and_parse_blastx_reads_ch = 
	Pull_Unassigned_MMseqs_Reads.out.unassigned_mmseqs_blastx_ch
	    .map { sid,
		   read_type,
                   assembler,
		   mode,
		   unassigned_blastx ->
	    
                 tuple(
		      sid,
	              read_type,
		      unassigned_blastx,
                      mode,
                      assembler
	     )
        }
    pull_and_parse_blastx_reads_ch.view { "PULL AND PARSE INPUT: $it " }

    Pull_and_Parse_Unmapped_Reads( pull_and_parse_blastx_reads_ch )


    def parsed_diamond_unassigned_grouped_ch =
	Pull_and_Parse_Unmapped_Reads.out.parsed_diamond_all_ch
	    .map { sid, read_type, assembler, mode, parsed_file ->
		tuple([sid, assembler, mode], tuple(read_type, parsed_file))
	    }
	    .groupTuple()
	    .map { key, items ->
		
		def sid		= key[0]
		def assembler	= key[1]
		def mode 	= key[2]

		def sr_parsed	= items
		    .findAll { it[0] == 'short' }
		    .collect { it[1] }

		def lr_parsed   = items
		    .findAll { it[0] == 'long' }
		    .collect { it[1] }

		tuple(sid, assembler, mode, sr_parsed, lr_parsed)
	    }

    parsed_diamond_unassigned_grouped_ch.view {
	"PARSED DIAMOND UNASSIGNED GROUPED: $it"
    }
  

    def filter_putative_viruses_ch =
        all_tagged
            .join(Pull_Unassigned_MMseqs_Blastx.out.unassigned_mmseqs_parsed_ch)
            .join(parse_mmseqs_out_all)
            .join(parsed_diamond_unassigned_grouped_ch)
            .join(params_only_ch)
            .map { sid,
                   fq1,
                   fq2,
                   lr,
                   contigs_fasta,
                   mmseqs_contigs_parsed,
                   mode_from_all_tagged,
                   cpus_from_all_tagged,
                   mem_from_all_tagged,
                   assembler_from_all_tagged,   

                   diamond_contigs_parsed,
                   diamond_sr_reads,
                   diamond_lr_reads,
                   mmseqs_unassigned_parsed,    

                   assembler_from_mmseqs_reads,
                   mmseqs_sr_parsed,
                   mmseqs_lr_parsed,    

                   assembler_from_diamond_unassigned,
                   mode_from_diamond_unassigned,
                   diamond_sr_parsed,
                   diamond_lr_parsed,   

                   mode,
                   cpus,
                   mem,
                   assembler -> 

                tuple(
                    sid,
                    mmseqs_contigs_parsed,
                    diamond_contigs_parsed,
                    firstPresent(mmseqs_sr_parsed),
                    firstPresent(mmseqs_lr_parsed),
                    firstPresent(diamond_sr_parsed),
                    firstPresent(diamond_lr_parsed),
                    assembler,
                    mode
                )
            }   

    filter_putative_viruses_ch.view {
        "FILTER PUTATIVE CLEAN INPUT: $it"
    }   

    Filter_Putative_Viruses_ALL(filter_putative_viruses_ch)

    Filter_Putative_Viruses_ALL.out.viral_blast_tab_ch.view { "Filter VS: $it"}



    def params_by_sid_ch = params_only_ch
        .map { row ->
            def sid       = row[0]
            def mode      = row[1]
            def cpus      = row[2]
            def mem       = row[3]
            def assembler = row[4]
            tuple(sid, [assembler, mode, cpus, mem])
        }   

    def params_k = params_only_ch
        .map { row ->
            def sid       = row[0]
            def mode      = row[1]
            def cpus      = row[2]
            def mem       = row[3]
            def assembler = row[4]
            tuple([assembler, sid], [mode, cpus, mem])
        }   

    params_k.view { "params_k: $it" }   
    

    def viral_k = Filter_Putative_Viruses_ALL.out.viral_blast_tab_ch
        .map { row ->
            def sid             = row[0]
            def viral_blast_tab = row[1]
            def assembler       = row[2]
            tuple([assembler, sid], viral_blast_tab)
        }   

    viral_k.view { "viral_k: $it" } 
    

    def contigs_k = all_tagged
        .map { row ->
            def sid     = row[0]
            def contigs = row[4]
            tuple(sid, contigs)
        }
        .join(params_by_sid_ch)
        .map { sid, contigs, meta ->
            def assembler = meta[0]
            tuple([assembler, sid], contigs)
        }   

    contigs_k.view { "contigs_k: $it" } 
    
    def reads_sr_real_k = Fastq_2_Fasta.out.fastq2fasta_ch_sr
        .map { row ->
            def sid         = row[0]
            def reads_fasta = row[1]
            tuple(sid, reads_fasta)
        }
        .join(params_by_sid_ch)
        .map { sid, reads_fasta, meta ->
            def assembler = meta[0]
            tuple([assembler, sid], reads_fasta)
        }   

    reads_sr_real_k.view { "reads_sr_real_k: $it" } 
    
    def reads_lr_real_k = Fastq_2_Fasta.out.fastq2fasta_ch_lr
        .map { row ->
            def sid         = row[0]
            def reads_fasta = row[1]
            tuple(sid, reads_fasta)
        }
        .join(params_by_sid_ch)
        .map { sid, reads_fasta, meta ->
            def assembler = meta[0]
            tuple([assembler, sid], reads_fasta)
        }   

    reads_lr_real_k.view { "reads_lr_real_k: $it" } 
    

    def sr_real_k = Map_Reads_2_Contigs_VS.out.map2assembly_ch_all
        .map { sid, outs -> 

            def pe_covstats = null
            def pe_qc_reads = null  

            if (outs instanceof List) {
                pe_covstats = outs.find { f -> f.name.endsWith("_assembly_mapping_SR_covstats.txt") }
                pe_qc_reads = outs.find { f -> f.name == "sr_count.txt" }
            }   

            tuple(sid, [pe_covstats, pe_qc_reads])
        }
        .join(params_by_sid_ch)
        .map { sid, pe_pair, meta ->
            def assembler = meta[0]
            tuple([assembler, sid], pe_pair)
        }   

    sr_real_k.view { "sr_real_k FIXED: $it" }   

    def lr_real_k = Map_Reads_2_Contigs_VS.out.map2assembly_ch_all
        .map { sid, outs -> 

            def lr_covstats = null
            def lr_qc_reads = null  

            if (outs instanceof List) {
                lr_covstats = outs.find { f -> f.name.endsWith("_assembly_mapping_LR_covstats.txt") }
                lr_qc_reads = outs.find { f -> f.name == "lr_count.txt" }
            }   

            tuple(sid, [lr_covstats, lr_qc_reads])
        }
        .join(params_by_sid_ch)
        .map { sid, lr_pair, meta ->
            def assembler = meta[0]
            tuple([assembler, sid], lr_pair)
        }   

    lr_real_k.view { "lr_real_k FIXED: $it" }
 
    def all_keys = viral_k
        .map { row -> tuple(row[0], true) }
        .mix(contigs_k.map { row -> tuple(row[0], true) })
        .mix(params_k.map { row -> tuple(row[0], true) })
        .groupTuple()
        .map { row -> row[0] }  

    all_keys.view { "all_keys: $it" }   
    
    def reads_sr_default_k = all_keys.map { key -> tuple(key, null) }
    def reads_lr_default_k = all_keys.map { key -> tuple(key, null) }   

    def sr_default_k = all_keys.map { key -> tuple(key, [null, null]) }
    def lr_default_k = all_keys.map { key -> tuple(key, [null, null]) } 
    

    def reads_sr_any_k = reads_sr_default_k
        .mix(reads_sr_real_k)
        .groupTuple()
        .map { row ->
            def key    = row[0]
            def values = row[1]
            def got    = values.find { it != null } ?: null
            tuple(key, got)
        }   

    reads_sr_any_k.view { "reads_sr_any_k: $it" }   
    

    def reads_lr_any_k = reads_lr_default_k
        .mix(reads_lr_real_k)
        .groupTuple()
        .map { row ->
            def key    = row[0]
            def values = row[1]
            def got    = values.find { it != null } ?: null
            tuple(key, got)
        }   

    reads_lr_any_k.view { "reads_lr_any_k: $it" }   
    

    def sr_any_k = sr_default_k
        .mix(sr_real_k)
        .groupTuple()
        .map { row ->
            def key   = row[0]
            def pairs = row[1]
            def got   = pairs.find { it[0] != null || it[1] != null } ?: pairs[0]
            tuple(key, got)
        }   

    sr_any_k.view { "sr_any_k: $it" }   
    

    def lr_any_k = lr_default_k
        .mix(lr_real_k)
        .groupTuple()
        .map { row ->
            def key   = row[0]
            def pairs = row[1]
            def got   = pairs.find { it[0] != null || it[1] != null } ?: pairs[0]
            tuple(key, got)
        }   

    lr_any_k.view { "lr_any_k: $it" }   
    

    def generate_vs_ch = viral_k
        .join(contigs_k)
        .join(reads_sr_any_k)
        .join(reads_lr_any_k)
        .join(sr_any_k)
        .join(lr_any_k)
        .join(params_k)
        .map { row ->
            def key                      = row[0]
            def viral_tab                = row[1]
            def contigs_fa               = row[2]
            def contig_short_reads_fasta = row[3]
            def contig_long_reads_fasta  = row[4]
            def pe_pair                  = row[5]
            def lr_pair                  = row[6]
            def meta                     = row[7]   

            def assembler = key[0]
            def sid       = key[1]  

            def pe_covstats = pe_pair[0]
            def pe_qc_reads = pe_pair[1]    

            def lr_covstats = lr_pair[0]
            def lr_qc_reads = lr_pair[1]    

            def mode = meta[0]  

            tuple(
                sid,
                viral_tab,
                contigs_fa,
                contig_short_reads_fasta,
                contig_long_reads_fasta,
                pe_covstats,
                pe_qc_reads,
                lr_covstats,
                lr_qc_reads,
                mode,
                assembler
            )
        }   

    generate_vs_ch.view { "Generate VS Output Channel: $it" }   

    Generate_VS_Outputs(generate_vs_ch)


}
