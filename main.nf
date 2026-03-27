#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ─── Imports ──────────────────────────────────────────────────
include { QC_Workflow }                         from './workflows/qc.nf'
include { MapReads_2_RefSeq }                   from './workflows/MapReads2RefSeq.nf'
include { Assembly_Workflow }                   from './workflows/assembly.nf'
include { BlastX_Reads_Workflow }               from './workflows/BlastX_Reads.nf'
include { Split_BamFiles_Workflow }             from './workflows/dorado_polish.nf'
include { Dorado_Aligner }                      from './workflows/modules/local/dorado_polish/main.nf'
include { Dorado_Polisher }                     from './workflows/modules/local/dorado_polish/main.nf'
include { DAA2INFO_contigs_daa_file }           from './workflows/modules/local/blastx/main.nf'
include { Hecatomb} 		                    from './workflows/modules/local/hecatomb/main.nf'
include { Quast } 	                            from './workflows/modules/local/quast/main.nf'
include { BlastX_contigs }                      from './workflows/modules/local/blastx/main.nf'
include { Meganize_BlastX_Contigs }             from './workflows/modules/local/blastx/main.nf'
include { Parse_BlastX_Contigs }                from './workflows/modules/local/blastx/main.nf'
include { Parse_MMSEQ_Contigs }                 from './workflows/modules/local/mmseq/mmseq_main.nf'
include { Remove_rRNA_Reads_Workflow }          from './workflows/Remove_rRNA_wf.nf'
include { CheckM_Assemblies }                   from './workflows/modules/local/checkm/main.nf'
include { MMSEQ_AllAssemblers }                 from './workflows/mmseqs/mmseq_allassemblers.nf'
//include { Split_Short_Reads }                   from './workflows/modules/local/seqkit_split_reads/main.nf'
//include { Split_Long_Reads }                    from './workflows/modules/local/seqkit_split_reads/main.nf'
include { Megan_ShortReads_WF }                 from './workflows/megan_workflow/main.nf'
include { Megan_LongReads_WF }                  from './workflows/megan_workflow/main.nf'
include { VS_Workflow }				            from './workflows/VS_calls.nf'
include { Read_Distribution }                   from './workflows/modules/local/qc/main.nf'
include { Autocycler_Workflow }                 from './workflows/autocycler_wf.nf'

Channel.value(file(params.mmseqs_nt_db)).set { nt_db_ch }
// ─── Main Workflow ─────────────────────────────────────
workflow {

    def raw_samples_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def fq1 = row.fastq_1?.contains('No_Read') ? null : file(row.fastq_1, checkIfExists: false)
            def fq2 = row.fastq_2?.contains('No_Read') ? null : file(row.fastq_2, checkIfExists: false)
            def lr  = row.long_read?.contains('No_Read') ? null : file(row.long_read, checkIfExists: false)

            def mode = fq1 && fq2 && lr ? 'hybrid'
                     : fq1 && fq2       ? 'short'
                     : lr               ? 'long'
                     : 'none'

            // Estimate size and set resources for QC
            def read_files = [fq1, fq2, lr].findAll()
            def total_size_gb = read_files.sum { it?.size() ?: 0 } / 1e9

            // Dynamically assign cpus and memory
            def cpus = total_size_gb < 1 ? 2 :
                       total_size_gb < 5 ? 8 :
                       total_size_gb < 10 ? 12 : 16
            def mem = total_size_gb < 1 ? '8 GB' :
                       total_size_gb < 5 ? '16 GB' :
                       total_size_gb < 10 ? '32 GB' : '64 GB'
                       
            tuple(row.sample_id, fq1, fq2, lr, mode, cpus, mem)
        }
        .filter { it[4] != 'none' }

    raw_samples_ch.view { "Raw sample data: $it" }

    // QC Workflow
    QC_Workflow(raw_samples_ch)


    // Ensure raw_samples_ch is consumed first (now used by QC_Workflow)

    def pipeline_mode_ch = raw_samples_ch
        .map { it[4] }
        .distinct()
        .collect()
        .map { modes ->
            if (modes.contains('hybrid')) return 'hybrid'
            else if (modes == ['short']) return 'short'
            else if (modes == ['long']) return 'long'
            else throw new IllegalStateException("Unsupported pipeline modes: ${modes}")
        }
    pipeline_mode_ch.view { "Detected pipeline mode: $it" }
    
    def short_ch
    def long_ch
    def short_ch_final
    def long_ch_final
    def sr_ch
    def lr_ch
    def final_long_ch
    def blast_short_ch
    if (params.map2reads) {
        MapReads_2_RefSeq(
            QC_Workflow.out.trimmed_short_ch,
            QC_Workflow.out.trimmed_long_ch
        )   

        MapReads_2_RefSeq.out.host_removed_fq_ch_sr.view { " host_removed_fq_ch_sr: $it" }
        MapReads_2_RefSeq.out.host_removed_fq_ch_lr.view { " host_removed_fq_ch_lr: $it" } 

        if (params.remove_rRNA_reads) {
            def rrna_input_ch = MapReads_2_RefSeq.out.host_removed_fq_ch_sr
                .map { id, sr -> tuple(id, [sr, null]) }
                .mix(
                    MapReads_2_RefSeq.out.host_removed_fq_ch_lr
                        .map { id, lr -> tuple(id, [null, lr]) }
                )
                .groupTuple()
                .map { id, values ->
                    def sr_file = values.collect { it[0] }.find { it != null }
                    def lr_file = values.collect { it[1] }.find { it != null }

                    def files = [sr_file, lr_file].findAll()
                    def total_size_gb = files.sum {it?.size() ?: 0 } / 1e9

                    def cpus = total_size_gb < 1 ? 8 :
                               total_size_gb < 5 ? 16 :
                               total_size_gb < 10 ? 32 :
                               total_size_gb < 20 ? 48 :
                               total_size_gb < 40 ? 64 : 80
                    
                    def mem = total_size_gb < 1 ? '32 GB' :
                               total_size_gb < 5 ? '64 GB' :
                               total_size_gb < 10 ? '128 GB' :
                               total_size_gb < 20 ? '256 GB' :
                               total_size_gb < 40 ? '384 GB' : '512 GB'

                    tuple(id, sr_file, lr_file, cpus, mem)
                }
            rrna_input_ch.view { "RNA input tuple: $it" }    

            Remove_rRNA_Reads_Workflow(rrna_input_ch)   

            short_ch = Remove_rRNA_Reads_Workflow.out.rRNA_host_remove_short
                .filter { sample_id, fq -> fq != null }
                .map { sample_id, fq -> tuple(sample_id, fq, 'short') } 

             blast_short_ch = Remove_rRNA_Reads_Workflow.out.rRNA_host_remove_short
                 .filter { sample_id, fq -> fq != null }
                 .map { sample_id, fq -> tuple(sample_id, fq, 'short') }

            long_ch = Remove_rRNA_Reads_Workflow.out.rRNA_host_remove_long
                .filter { sample_id, fq -> fq != null }
                .map { sample_id, fq -> tuple(sample_id, fq, 'long') }  

            long_ch_final = Remove_rRNA_Reads_Workflow.out.rRNA_host_remove_long
                .map { sample_id, lr -> tuple(sample_id, lr) }  

            //def split_input_ch = Remove_rRNA_Reads_Workflow.out.rRNA_host_remove_short
            //    .filter { sample_id, fq -> fq != null }
            //    .map { sample_id, fq -> tuple(sample_id, fq, 'short') } 

            //Split_Interleave(split_input_ch)    


            short_ch_final = 
                Remove_rRNA_Reads_Workflow.out.rRNA_host_remove_short_r1
                    .join( Remove_rRNA_Reads_Workflow.out.rRNA_host_remove_short_r2, by: 0)
                    .map { sid, r1, r2 -> tuple(sid, r1, r2) }

            //short_ch_final = Split_Interleave.out.split_interleave_ch
            //    .map { sample_id, r1, r2 -> tuple(sample_id, r1, r2) }  

            //short_ch_final = 
            //    MapReads_2_RefSeq.out.host_removed_fq_ch_r1
            //        .join( MapReads_2_RefSeq.out.host_removed_fq_ch_r2, by: 0)
            //        .map { sid, r1, r2 -> tuple(sid, r1, r2) }

	     
            short_ch_final.view { "short_ch_final: $it" } 

        } else {
            short_ch = MapReads_2_RefSeq.out.host_removed_fq_ch_sr
                .filter { sample_id, fq -> fq != null }
                .map { sample_id, fq -> tuple(sample_id, fq, 'short') } 

	    blast_short_ch = MapReads_2_RefSeq.out.host_removed_fq_ch_sr
		.filter { sample_id, fq -> fq != null }
                .map { sample_id, fq -> tuple(sample_id, fq, 'short') }

            long_ch = MapReads_2_RefSeq.out.host_removed_fq_ch_lr
                .filter { sample_id, fq -> fq != null }
                .map { sample_id, fq -> tuple(sample_id, fq, 'long') }  

            long_ch_final = MapReads_2_RefSeq.out.host_removed_fq_ch_lr
                .map { sample_id, lr -> tuple(sample_id, lr) }  

            def split_input_ch = MapReads_2_RefSeq.out.host_removed_fq_ch_sr
                .filter { sample_id, fq -> fq != null }
                .map { sample_id, fq -> tuple(sample_id, fq, 'short') } 

            //Split_Interleave(split_input_ch)    

            //short_ch_final = Split_Interleave.out.split_interleave_ch
            //    .map { sample_id, r1, r2 -> tuple(sample_id, r1, r2) }

            short_ch_final = 
                MapReads_2_RefSeq.out.host_removed_fq_ch_r1
                    .join( MapReads_2_RefSeq.out.host_removed_fq_ch_r2, by: 0)
                    .map { sid, r1, r2 -> tuple(sid, r1, r2) }

            short_ch_final.view { "short_ch_final: $it" } 
        }   

    } else {
            short_ch = QC_Workflow.out.interleaved_trimmed_short_ch
                .filter { sample_id, fq -> fq != null }
                .map { sample_id, fq -> tuple(sample_id, fq, 'short') }     
            
            blast_short_ch = QC_Workflow.out.interleaved_trimmed_short_ch
                .filter { sample_id, fq -> fq != null }
                .map { sample_id, fq -> tuple(sample_id, fq, 'short') }     

            long_ch = QC_Workflow.out.trimmed_long_ch
                .filter { sample_id, fq -> fq != null }
                .map { sample_id, fq -> tuple(sample_id, fq, 'long') }      

            long_ch_final = QC_Workflow.out.trimmed_long_ch
                .filter { sample_id, fq -> fq != null }
                .map { sample_id, fq -> tuple(sample_id, fq) }      

            //def split_input_ch = QC_Workflow.out.interleaved_trimmed_short_ch
            //    .filter { sample_id, fq -> fq != null }
            //    .map { sample_id, fq -> tuple(sample_id, fq, 'short') }     

            //Split_Interleave(split_input_ch)        

            //short_ch_final = Split_Interleave.out.split_interleave_ch
            //    .map { sample_id, r1, r2 -> tuple(sample_id, r1, r2) }

            short_ch_final = 
                QC_Workflow.out.trimmed_short_ch
            short_ch_final.view { "short_ch_final: $it" }             
    }
    // ─── Enforce tuple structure ──────────────────────────────────────────────   

    // Ensure short_ch_final is a 3-tuple (sample_id, R1, R2)
    short_ch_final = short_ch_final.map { it -> tuple(it[0], it[1], it[2]) }    

    // Ensure long_ch_final is a 2-tuple (sample_id, LR)
    long_ch_final = long_ch_final.map { it -> tuple(it[0], it[1]) } 

    // Optional: debug before building hybrid
    short_ch_final.view { "short_ch_final (tuple): $it" }
    long_ch_final.view  { "long_ch_final (tuple): $it" }  
    
    // ─── Build conditional assembly input ─────────────────────────────────────   

    // Short-only input channel: (sample_id, R1, R2, null, null, 'short')
    def short_only_ch = short_ch_final.map { id, r1, r2 ->
        tuple(id, r1, r2, null, null, 'short')
    }   

    // Long-only input channel: (sample_id, null, null, LR, 'long')
    //def long_only_ch = final_long_ch.map { id, lr ->
    //    tuple(id, null, null, lr, q2, 'long')
    //} 

    // ensure sample_id is a plain String in both channels
    short_ch_final = short_ch_final.map { it -> tuple(it[0], it[1], it[2]) }
    long_ch_final  = long_ch_final.map { it -> tuple(it[0], it[1]) }

    // ─── Read Distribution Statistics ─────────────────────────────────────────
    Read_Distribution(long_ch_final)
    q2_value_ch = Read_Distribution.out.read_distribution_ch.map { row ->
	def (sample_id, png_file, summary_file, csv_file) = row
	def lines = csv_file.text.readLines()
	def header = lines[0].trim()
	def q2 = lines[1].trim()
	if (header != 'Q1') {
		throw new RuntimeException("Expected header Q2 but found ${header} in ${csv_file}")
	}
	[ sample_id, q2 ]
    }
    //final_long_ch = long_ch_final.join(q2_value_ch, by: 0 )
    //	.map { sample_id, lr, q2 -> [ sample_id, lr, q2 ] }	
    final_long_ch = long_ch_final.join(q2_value_ch, by: 0) // this is now what get's passed to assembly workflow
    final_long_ch.map { sample_id, lr, q2 -> 
	println "Sample: ${sample_id} LR: ${lr} Q2: ${q2}"
    }
    // Long-only input channel: (sample_id, null, null, LR, 'long')
    def long_only_ch = final_long_ch.map { id, lr, q2 ->
        tuple(id, null, null, lr, q2, 'long')
    }
    short_ch_final.view { " short tuple: $it (${it.getClass()})" }
    final_long_ch.view  { " long tuple:  $it (${it.getClass()})" }
    def hybrid_ch = short_ch_final
        .map { id, r1, r2 -> [ id, ['short', r1, r2] ] }
        .mix(
            final_long_ch.map { id, lr, q2 -> [ id, ['long', lr, q2] ] }
        )
        .groupTuple()
        .filter { id, values ->
            values.find { it[0] == 'short' } && values.find { it[0] == 'long' }
        }
        .map { id, values ->
            def short_vals = values.find { it[0] == 'short' }
            def long_vals  = values.find { it[0] == 'long' }    

            def r1 = short_vals[1]
            def r2 = short_vals[2]
            def lr = long_vals[1]   
            def q2 = long_vals[2]   

            [ id, r1, r2, lr, q2, 'hybrid' ]
        }
        
    hybrid_ch.view { it -> println " hybrid_ch: ${it} (${it.getClass()})" }

    def assembly_input_ch   

    if (params.hybrid) {
        assembly_input_ch = hybrid_ch
    } else if (params.shortreads) {
        assembly_input_ch = short_only_ch
    } else if (params.longreads) {
        assembly_input_ch = long_only_ch
    } else {
        error "You must set one of: --shortreads, --longreads, or --hybrid"
    }
    
    assembly_input_ch.view { " assembly_input_ch: $it" }

    // ─── Call the Autocycler Workflow ────────────────────────────────────────  
    if ( params.run_autocycler ) {
        //def autocycler_in_ch = assembly_input_ch.map { sid, r1, r2, lr, q2, mode ->
        //    tuple(sid, r1, r2, lr, mode) }

        Autocycler_Workflow(
           assembly_input_ch
        )
    }

    // ─── Call the Assembly Workflow ────────────────────────────────────────  

    Assembly_Workflow(
        assembly_input_ch
    )
    short_ch.view { "short_ch: $it" }
    long_ch.view  { "long_ch: $it" }
    if (params.medaka) {
        Assembly_Workflow.out.dragonflye_medaka_assembly_ch.view { " Dragonflye_Medaka output: $it" }
    }

    if ( params.dorado_polish ) {
        log.info " Running Split_BamFiles_Workflow"
        log.info " Using BAM file: ${params.bam_file}"

        def dorado_split_input_ch = Channel.fromPath(params.bam_file, checkIfExists: true)
        dorado_split_input_ch.view { " Input to Split_BamFiles_Workflow: $it" }

        Split_BamFiles_Workflow(dorado_split_input_ch)

        //  Extract sample_id from BAM filename
        def dorado_split_tuples_ch = Split_BamFiles_Workflow.out.dorado_split_output_ch
            .flatten()
            .map { bam_file ->
                def sample_id = bam_file.getBaseName().replaceFirst(/\.bam$/, '')
                tuple(sample_id, bam_file)
            }

        dorado_split_tuples_ch.view { " BAM as tuples: $it" }
        Assembly_Workflow.out.dragonflye_assembly_ch.view { " Dragonflye output: $it" }

        //  Join on sample_id and feed to Dorado_Aligner
        def dorado_aligner_input_ch = dorado_split_tuples_ch
            .join(Assembly_Workflow.out.dragonflye_assembly_ch)
            .map { sample_id, bam_file, contigs_fasta -> tuple(sample_id, bam_file, contigs_fasta) }

        dorado_aligner_input_ch.view { "Joined input to Dorado_Aligner: $it" }

        Dorado_Aligner(dorado_aligner_input_ch)
        Dorado_Aligner.out.dorado_aligner_ch.view { " Dorado aligner output: $it" }
        // Step 1: Get original BAM input (before alignment)
        def dorado_input_bam_ch = dorado_aligner_input_ch
            .map { sample_id, bam_file, contigs_fasta -> tuple(sample_id, bam_file) }       

        // Step 2: Wait for aligner to produce .bam.bai
        def dorado_bai_ch = Dorado_Aligner.out.dorado_aligner_ch
            .map { sample_id, bai_file -> tuple(sample_id, bai_file) }      

        // Step 3: Join BAM + .bai on sample_id
        def dorado_bam_with_index_ch = dorado_input_bam_ch
            .join(dorado_bai_ch)
            .map { sample_id, bam_file, bai_file -> tuple(sample_id, bam_file, bai_file) }      

        // Step 4: Join with contigs.fa from Dragonflye
        def dorado_polisher_input_ch = dorado_bam_with_index_ch
            .join(Assembly_Workflow.out.dragonflye_assembly_ch)
            .map { sample_id, bam_file, bai_file, contigs_fasta -> tuple(sample_id, bam_file, bai_file, contigs_fasta) }        

        dorado_polisher_input_ch.view { "Input to Dorado_Polisher: $it" }     

        Dorado_Polisher(dorado_polisher_input_ch)
        Dorado_Polisher.out.dorado_polisher_ch.view { " Dorado_Polisher output: $it" }
    }


    // CheckM (waits on assembly steps to finish)
    if ( params.dragonflye || params.medaka || params.metaspades || params.unicycler) {
        Channel.empty()
        .mix( params.dragonflye && Assembly_Workflow.out.dragonflye_assembly_ch ? Assembly_Workflow.out.dragonflye_assembly_ch.map { sid, fa -> tuple(sid, 'dragonflye', fa) } : Channel.empty() )
        .mix( params.medaka && Assembly_Workflow.out.dragonflye_medaka_assembly_ch ? Assembly_Workflow.out.dragonflye_medaka_assembly_ch.map { sid, fa -> tuple(sid, 'dragonflye_medaka', fa) } : Channel.empty() )
        .mix( params.unicycler && Assembly_Workflow.out.unicycler_assembly_ch ? Assembly_Workflow.out.unicycler_assembly_ch.map { sid, fa -> tuple(sid, 'unicycler', fa) } : Channel.empty() )
        .mix( params.metaspades && Assembly_Workflow.out.metaspades_assembly_ch ? Assembly_Workflow.out.metaspades_assembly_ch.map { sid, fa -> tuple(sid, 'metaspades', fa) } : Channel.empty() )
        .set { checkm_input_ch }
        CheckM_Assemblies( checkm_input_ch )
    }



    if ( params.nanopore_assembly ) {
        
	// Remove 'mode' field from assembly_input_ch before passing to Hecatomb
	def hecatomb_input_ch = assembly_input_ch.map { id, r1, r2, lr, q2,  mode -> 
	    tuple(id, r1, r2, lr)
	}	

	Hecatomb(hecatomb_input_ch)
        
    // Input to Stage_7a_quast: tuple(sample_id, merged_fasta_path)
	def quast_input_ch = assembly_input_ch
        	.map { sample_id, fq1, fq2, lr, mode -> tuple(sample_id, fq1, fq2, q2, lr) } // remove mode
        	.join(Hecatomb.out.hecatomb_merged_fasta)
        	.map { sample_id, fq1, fq2, lr, merged_fasta ->
            	tuple(sample_id, fq1, fq2, lr, merged_fasta)
      	} 

	Quast(quast_input_ch)
    
    }

    Assembly_Workflow.out.dragonflye_assembly_ch.view { " Dragonflye output: $it" }
    Assembly_Workflow.out.spades_assembly_ch.view { " SPAdes output: $it" }
    Assembly_Workflow.out.metaspades_assembly_ch.view { " MetaSPAdes output: $it" }
    Assembly_Workflow.out.unicycler_assembly_ch.view { " Unicycler output: $it" }
    Assembly_Workflow.out.plasmidspades_assembly_ch.view { " PlasmidSPAdes output: $it" }

if ( params.run_blastx ) {

        Channel.empty()
        // Dragonflye
        .mix( params.dragonflye && Assembly_Workflow.out.dragonflye_assembly_ch ?
        Assembly_Workflow.out.dragonflye_assembly_ch.map { sid, fa -> tuple(sid, 'dragonflye', fa) } : Channel.empty() )
        // Dragonflye + Medaka (only include if you have such a channel and flag)
        .mix( (params.dragonflye && params.medaka && Assembly_Workflow.out.dragonflye_medaka_assembly_ch) ?
        Assembly_Workflow.out.dragonflye_medaka_assembly_ch.map { sid, fa -> tuple(sid, 'dragonflye_medaka', fa) } : Channel.empty() )
        // SPAdes
        .mix( params.spades && Assembly_Workflow.out.spades_assembly_ch ?
        Assembly_Workflow.out.spades_assembly_ch.map { sid, fa -> tuple(sid, 'spades', fa) } : Channel.empty() )
        // metaSPAdes
        .mix( params.metaspades && Assembly_Workflow.out.metaspades_assembly_ch ?
        Assembly_Workflow.out.metaspades_assembly_ch.map { sid, fa -> tuple(sid, 'metaspades', fa) } : Channel.empty() )
        // Unicycler
        .mix( params.unicycler && Assembly_Workflow.out.unicycler_assembly_ch ?
        Assembly_Workflow.out.unicycler_assembly_ch.map { sid, fa -> tuple(sid, 'unicycler', fa) } : Channel.empty() )
        .set { blastx_input_ch }


        BlastX_contigs( blastx_input_ch )

        // (sample_id, assembler, daa) -> (sample_id, assembler, diamondview_file)
        Meganize_BlastX_Contigs( BlastX_contigs.out.blastx_contigs_ch )

        Parse_BlastX_Contigs( Meganize_BlastX_Contigs.out.meganized_blastx_contigs_ch )

        // Debug
        BlastX_contigs.out.blastx_contigs_ch.view { it -> "BLASTX DAA: $it" }
        Parse_BlastX_Contigs.out.parse_blastx_ch.view { it -> "Parsed BLASTX: $it" }
        // Call the workflow with tuple (sid, fq, mode)

        // Normalize SHORT to (sid, path(fq), 'short')
        def final_short_reads_flat =
            (params.shortreads || params.hybrid)
            ? blast_short_ch
                .map     { sid, fq, mode -> tuple(sid, fq, mode) }                    // (sid, fq, 'short')
            : Channel.empty()

        // Normalize LONG to (sid, path(fq), 'long')
        def final_long_reads_flat =
            (params.longreads || params.hybrid)
            ? long_ch_final
                .map     { sid, fq -> tuple(sid, fq, 'long') }                     // (sid, fq, 'long')
            : Channel.empty()
        // optional, but useful for debugging, we can remove this
        final_short_reads_flat.view { "BLASTX short input: $it" }
        final_long_reads_flat.view  { "BLASTX long  input: $it" }
        BlastX_Reads_Workflow(final_short_reads_flat, final_long_reads_flat)    
        // optinal useful for debugging; expected per-part: (sid, mode, daa)
        BlastX_Reads_Workflow.out.short_reads_blastx_channel.view { " Short DAA: $it" }
        BlastX_Reads_Workflow.out.long_reads_blastx_channel.view  { " Long  DAA: $it" }     
        Megan_ShortReads_WF( BlastX_Reads_Workflow.out.short_reads_blastx_channel )
        Megan_LongReads_WF( BlastX_Reads_Workflow.out.long_reads_blastx_channel ) 
    }

    // place holder for VirusSeeker Workflow
    // Build one mixed stream of all enabled assemblers:
    // each must be (asm, sid, contigs.fa)
    def contigs_all = Channel.empty()   

    if (params.dragonflye) {
      contigs_all = contigs_all.mix(
        Assembly_Workflow.out.dragonflye_assembly_ch.map { sid, fa -> tuple('dragonflye', sid, fa) }
      )
    }   

    if (params.medaka) {
      contigs_all = contigs_all.mix(
        Assembly_Workflow.out.dragonflye_medaka_assembly_ch.map { sid, fa -> tuple('dragonflye_medaka', sid, fa) }
      )
    }   

    if (params.metaspades) {
      contigs_all = contigs_all.mix(
        Assembly_Workflow.out.metaspades_assembly_ch.map { sid, fa -> tuple('metaspades', sid, fa) }
      )
    }   

    if (params.unicycler) {
      contigs_all = contigs_all.mix(
        Assembly_Workflow.out.unicycler_assembly_ch.map { sid, fa -> tuple('unicycler', sid, fa) }
      )
    }   

    def mmseqs_dragonflye_contigs_ch = Channel.empty()
    def mmseqs_dragonflye_medaka_contigs_ch = Channel.empty()
    def mmseqs_metaspades_contigs_ch = Channel.empty()
    def mmseqs_unicycler_contigs_ch = Channel.empty()

    if (params.run_mmseqs) {
    MMSEQ_AllAssemblers(contigs_all, nt_db_ch)

    def dragonCh = MMSEQ_AllAssemblers.out.dragon ?: Channel.empty()
    def medakaCh = MMSEQ_AllAssemblers.out.medaka ?: Channel.empty()
    def metaCh = MMSEQ_AllAssemblers.out.meta ?: Channel.empty()
    def uniCh = MMSEQ_AllAssemblers.out.uni ?: Channel.empty()

    def mmseq_parse_in = Channel.empty()
    .mix( dragonCh.map { sid, tsv -> tuple(sid, 'dragonflye', tsv) } )
    .mix( medakaCh.map { sid, tsv -> tuple(sid, 'dragonflye_medaka', tsv) } )
    .mix( metaCh .map { sid, tsv -> tuple(sid, 'metaspades', tsv) } )
    .mix( uniCh .map { sid, tsv -> tuple(sid, 'unicycler', tsv) } )

    Parse_MMSEQ_Contigs( mmseq_parse_in )

    mmseqs_dragonflye_contigs_ch =
    Parse_MMSEQ_Contigs.out.mmseqs_parsed_ch
    .filter { sid, asm, f -> asm == 'dragonflye' }
    .map { sid, asm, f -> tuple(sid, f) }

    mmseqs_dragonflye_medaka_contigs_ch =
    Parse_MMSEQ_Contigs.out.mmseqs_parsed_ch
    .filter { sid, asm, f -> asm == 'dragonflye_medaka' } // note label
    .map { sid, asm, f -> tuple(sid, f) }

    mmseqs_metaspades_contigs_ch =
    Parse_MMSEQ_Contigs.out.mmseqs_parsed_ch
    .filter { sid, asm, f -> asm == 'metaspades' }
    .map { sid, asm, f -> tuple(sid, f) }

    mmseqs_unicycler_contigs_ch =
    Parse_MMSEQ_Contigs.out.mmseqs_parsed_ch
    .filter { sid, asm, f -> asm == 'unicycler' }
    .map { sid, asm, f -> tuple(sid, f) }
}

    if (params.vs) {    
        def reads_min_ch = assembly_input_ch.map { sid, fq1, fq2, lr, q2, mode ->
        tuple(sid, fq1, fq2, lr) // (sid, fq1, fq2, lr)
        }
        def vs_calls_inputs = Channel.empty()
        // DRAGONFLYE
        if (params.dragonflye) {
            def vs_dragon = reads_min_ch
                .join(Assembly_Workflow.out.dragonflye_assembly_ch) // (sid, fq1,fq2,lr, contigs)
                .join(mmseqs_dragonflye_contigs_ch) // (sid, fq1,fq2,lr, contigs, mmseqs_parsed)
                .map { sid, fq1, fq2, lr, contigs, mmseqs_parsed ->
                // add assembler tag as 7th field
                tuple(sid, fq1, fq2, lr, contigs, mmseqs_parsed, 'dragonflye')
            }
            vs_calls_inputs = vs_calls_inputs.mix(vs_dragon)
            }
        // MEDAKA
        if (params.medaka) {
            def vs_medaka = reads_min_ch
                .join(Assembly_Workflow.out.dragonflye_medaka_assembly_ch)
                .join(mmseqs_dragonflye_medaka_contigs_ch)
                .map { sid, fq1, fq2, lr, contigs, mmseqs_parsed ->
                tuple(sid, fq1, fq2, lr, contigs, mmseqs_parsed, 'medaka')
            }
            vs_calls_inputs = vs_calls_inputs.mix(vs_medaka)
            }

        // METASPADES
        if (params.metaspades) {
            def vs_meta = reads_min_ch
                .join(Assembly_Workflow.out.metaspades_assembly_ch)
                .join(mmseqs_metaspades_contigs_ch)
                .map { sid, fq1, fq2, lr, contigs, mmseqs_parsed ->
                tuple(sid, fq1, fq2, lr, contigs, mmseqs_parsed, 'metaspades')
            }
            vs_calls_inputs = vs_calls_inputs.mix(vs_meta)
            }

        // UNICYCLER
        if (params.unicycler) {
            def vs_uni = reads_min_ch
                .join(Assembly_Workflow.out.unicycler_assembly_ch)
                .join(mmseqs_unicycler_contigs_ch)
                .map { sid, fq1, fq2, lr, contigs, mmseqs_parsed ->
                tuple(sid, fq1, fq2, lr, contigs, mmseqs_parsed, 'unicycler')
            }
            vs_calls_inputs = vs_calls_inputs.mix(vs_uni)
        }

        // optional: peek
        vs_calls_inputs.view { "vs_calls_input (all assemblers): $it" }
        VS_Workflow( vs_calls_inputs )
    }

}
