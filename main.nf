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
include { BlastX_contigs }                      from './workflows/modules/local/blastx/main.nf'
include { Meganize_BlastX_Contigs }             from './workflows/modules/local/blastx/main.nf'
include { Parse_BlastX_Contigs }                from './workflows/modules/local/blastx/main.nf'
include { Parse_MMSEQ_Contigs }                 from './workflows/modules/local/mmseq/mmseq_main.nf'
include { Remove_rRNA_Reads_Workflow }          from './workflows/Remove_rRNA_wf.nf'
include { MMSEQ_AllAssemblers }                 from './workflows/mmseqs/mmseq_allassemblers.nf'
include { Megan_ShortReads_WF }                 from './workflows/megan_workflow/main.nf'
include { Megan_LongReads_WF }                  from './workflows/megan_workflow/main.nf'
include { VS_Workflow }				            from './workflows/VS_calls.nf'
include { Read_Distribution }                   from './workflows/modules/local/qc/main.nf'
include { Posttrim_NanoPlot as Posttrim_NanoPlot_After_ReadDistribution } from './workflows/modules/local/qc/main.nf'
include { Autocycler_Workflow }                 from './workflows/autocycler_wf.nf'
include { CONTIG_based_analysis } 		        from './workflows/contigs_based_analysis.nf'
include { LongRead_NanoPlot_QC_Stats }          from './workflows/modules/local/qc/main.nf'
include { Reads_Taxonomic_Classifier_Workflow } from './workflows/reads_taxonomic_classifier.nf'
// ─── Main Workflow ─────────────────────────────────────
workflow {

    def nt_db_ch = Channel.value(file(params.mmseqs_nt_db))

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
            def cpus = total_size_gb < 1 ? 128 :
                       total_size_gb < 5 ? 16 :
                       total_size_gb < 10 ? 24 : 32
            def mem = total_size_gb < 1 ? '128 GB' :
                       total_size_gb < 5 ? '32 GB' :
                       total_size_gb < 10 ? '48 GB' : '64 GB'
                       
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


            short_ch_final = 
                Remove_rRNA_Reads_Workflow.out.rRNA_host_remove_short_r1
                    .join( Remove_rRNA_Reads_Workflow.out.rRNA_host_remove_short_r2, by: 0)
                    .map { sid, r1, r2 -> tuple(sid, r1, r2) }
	     
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

            short_ch_final = 
                MapReads_2_RefSeq.out.host_removed_fq_ch_r1
                    .join( MapReads_2_RefSeq.out.host_removed_fq_ch_r2, by: 0)
                    .map { sid, r1, r2 -> tuple(sid, r1, r2) }

            short_ch_final.view { "short_ch_final: $it" } 
        }   

    } else {

        /*
         * No host-removal path.
         *
         * If --remove_rRNA_reads is true, run rRNA removal directly on
         * the QC outputs. This is the path needed for testing Ribodetector
         * without also enabling --map2reads.
         */
        if (params.remove_rRNA_reads) {

            def rrna_input_ch = QC_Workflow.out.interleaved_trimmed_short_ch
                .map { id, sr -> tuple(id, [sr, null]) }
                .mix(
                    QC_Workflow.out.trimmed_long_ch
                        .map { id, lr -> tuple(id, [null, lr]) }
                )
                .groupTuple()
                .map { id, values ->
                    def sr_file = values.collect { it[0] }.find { it != null }
                    def lr_file = values.collect { it[1] }.find { it != null }

                    def files = [sr_file, lr_file].findAll()
                    def total_size_gb = files.sum { it?.size() ?: 0 } / 1e9

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

            rrna_input_ch.view { "rRNA input tuple from QC outputs: $it" }

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
                .filter { sample_id, fq -> fq != null }
                .map { sample_id, fq -> tuple(sample_id, fq) }

            short_ch_final =
                Remove_rRNA_Reads_Workflow.out.rRNA_host_remove_short_r1
                    .join(Remove_rRNA_Reads_Workflow.out.rRNA_host_remove_short_r2, by: 0)
                    .map { sid, r1, r2 -> tuple(sid, r1, r2) }

            short_ch_final.view { "short_ch_final after rRNA removal: $it" }
        }

        /*
         * Original no-host-removal / no-rRNA-removal path.
         */
        else {
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

            short_ch_final =
                QC_Workflow.out.trimmed_short_ch

            short_ch_final.view { "short_ch_final: $it" }
        }
    }
    // ─── Enforce tuple structure ──────────────────────────────────────────────   

    // Ensure short_ch_final is a 3-tuple (sample_id, R1, R2)
    short_ch_final = short_ch_final.map { it -> tuple(it[0], it[1], it[2]) }    

    // Ensure long_ch_final is a 2-tuple (sample_id, LR)
    long_ch_final = long_ch_final.map { it -> tuple(it[0], it[1]) } 

    // Optional: debug before building hybrid
    short_ch_final.view { "short_ch_final (tuple): $it" }
    long_ch_final.view  { "long_ch_final (tuple): $it" }  
   
    if (params.run_reads_taxonomic_classifier) {

    def taxonomic_short_reads_ch =
        (params.shortreads || params.hybrid) ?
            short_ch_final
                .filter { sid, r1, r2 -> r1 != null && r2 != null }
                .map { sid, r1, r2 ->
                    tuple("${sid}_short", [r1, r2])
                } :
            Channel.empty()

    def taxonomic_long_reads_ch =
        (params.longreads || params.hybrid) ?
            long_ch_final
                .filter { sid, lr -> lr != null }
                .map { sid, lr ->
                    tuple("${sid}_long", [lr])
                } :
            Channel.empty()

    def reads_taxonomic_classifier_input_ch =
        taxonomic_short_reads_ch.mix(taxonomic_long_reads_ch)

    reads_taxonomic_classifier_input_ch.view {
        "reads_taxonomic_classifier_input_ch: ${it}"
    }

    Reads_Taxonomic_Classifier_Workflow(
        reads_taxonomic_classifier_input_ch
    )
    }


 
    // ─── Build conditional assembly input ─────────────────────────────────────   

    // Short-only input channel: (sample_id, R1, R2, null, null, 'short')
    def short_only_ch = short_ch_final.map { id, r1, r2 ->
        tuple(id, r1, r2, null, null, 'short')
    }   

    // ensure sample_id is a plain String in both channels
    short_ch_final = short_ch_final.map { it -> tuple(it[0], it[1], it[2]) }
    long_ch_final  = long_ch_final.map { it -> tuple(it[0], it[1]) }

    // ─── Read Distribution Statistics ─────────────────────────────────────────
    
    Read_Distribution(long_ch_final)
    q2_value_ch = Read_Distribution.out.read_distribution_ch.map { sample_id, png_file, summary_file, csv_file, trimmed_fastq ->
            def lines = csv_file.text.readLines()
            def header = lines[0].trim()
            def q2 = lines[1].trim()

            if (header != 'Q1') {
                throw new RuntimeException("Expected header Q1 but found ${header} in ${csv_file}")
            }

            tuple(sample_id, q2)
        }

    post_distribution_long_ch = Read_Distribution.out.read_distribution_fastq_ch
            .map { sample_id, fq -> tuple(sample_id, fq) }
    if (params.longreads || params.hybrid) {
    	Posttrim_NanoPlot_After_ReadDistribution(post_distribution_long_ch)

    	LongRead_NanoPlot_QC_Stats(
            QC_Workflow.out.nanoplot_pre_ch.collect(),
            Posttrim_NanoPlot_After_ReadDistribution.out.posttrim_nanoplot_ch.collect()
        )
        final_long_ch = post_distribution_long_ch.join(q2_value_ch, by: 0)
        } else {
	
	    final_long_ch = Channel.empty()
     }

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

    // Map2Assembly (waits on assembly steps to finish)
    //if ( params.dragonflye || params.medaka || params.metaspades || params.unicycler || params.raven ) {
    if ( params.characterize_contigs ) {
        def only_reads_ch = assembly_input_ch.map { sid, fq1, fq2, lr, q2, mode -> tuple(sid, fq1, fq2, lr, mode)}
        def characterize_contigs_ch = Channel.empty()
        if (params.dragonflye) {
            def dragonflye_readsAndContigs = only_reads_ch.join(Assembly_Workflow.out.dragonflye_assembly_ch).map {sid, fq1, fq2, lr, mode, contigs -> tuple(sid, fq1, fq2, lr, contigs, mode, 'dragonflye') } 
            characterize_contigs_ch = characterize_contigs_ch.mix(dragonflye_readsAndContigs) }
        if (params.medaka) { 
            def medaka_readsAndContigs =only_reads_ch.join(Assembly_Workflow.out.dragonflye_medaka_assembly_ch).map {sid, fq1, fq2, lr, mode, contigs -> tuple(sid, fq1, fq2, lr, contigs, mode, 'medaka') }
            characterize_contigs_ch = characterize_contigs_ch.mix(medaka_readsAndContigs) }
        if (params.metaspades) { 
            def metaspades_readsAndContigs = only_reads_ch.join(Assembly_Workflow.out.metaspades_assembly_ch).map {sid, fq1, fq2, lr, mode, contigs -> tuple(sid, fq1, fq2, lr, contigs, mode, 'metaspades') }
            characterize_contigs_ch = characterize_contigs_ch.mix(metaspades_readsAndContigs) }
        if (params.unicycler) { 
            def unicycler_readsAndContigs = only_reads_ch.join(Assembly_Workflow.out.unicycler_assembly_ch).map {sid, fq1, fq2, lr, mode, contigs -> tuple(sid, fq1, fq2, lr, contigs, mode, 'unicycler') }
            characterize_contigs_ch = characterize_contigs_ch.mix(unicycler_readsAndContigs) }
        if (params.raven) { 
            def raven_readsAndContigs = only_reads_ch.join(Assembly_Workflow.out.dragonflye_raven_assembly_ch).map {sid, fq1, fq2, lr, mode, contigs -> tuple(sid, fq1, fq2, lr, contigs, mode, 'raven') }
            characterize_contigs_ch = characterize_contigs_ch.mix(raven_readsAndContigs) }

        if (params.myloasm) {
            def myloasm_readsAndContigs = only_reads_ch
                .join(Assembly_Workflow.out.myloasm_assembly_ch)
                .map { sid, fq1, fq2, lr, mode, contigs ->
                    tuple(sid, fq1, fq2, lr, contigs, mode, 'myloasm')
                }

            characterize_contigs_ch = characterize_contigs_ch.mix(myloasm_readsAndContigs)
        }
        characterize_contigs_ch.view { " characterize_contigs_ch input: $it" }

        CONTIG_based_analysis( characterize_contigs_ch )

        }


    Assembly_Workflow.out.dragonflye_assembly_ch.view { " Dragonflye output: $it" }
    Assembly_Workflow.out.spades_assembly_ch.view { " SPAdes output: $it" }
    Assembly_Workflow.out.metaspades_assembly_ch.view { " MetaSPAdes output: $it" }
    Assembly_Workflow.out.unicycler_assembly_ch.view { " Unicycler output: $it" }
    Assembly_Workflow.out.plasmidspades_assembly_ch.view { " PlasmidSPAdes output: $it" }
    Assembly_Workflow.out.dragonflye_raven_assembly_ch.view { " Dragonflye Raven output: $it" }
    Assembly_Workflow.out.myloasm_assembly_ch.view { " Myloasm output: $it" } 

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
	    // Raven
	    .mix( params.raven && Assembly_Workflow.out.dragonflye_raven_assembly_ch ?
        Assembly_Workflow.out.dragonflye_raven_assembly_ch.map { sid, fa -> tuple(sid, 'raven', fa) } : Channel.empty() )
        .mix( params.myloasm && Assembly_Workflow.out.myloasm_assembly_ch ?
        Assembly_Workflow.out.myloasm_assembly_ch.map { sid, fa -> tuple(sid, 'myloasm', fa) } : Channel.empty() )
        .set { blastx_input_ch }
        blastx_input_ch.view { it -> "BLASTX input channel: $it" }

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
    if (params.raven) {
      contigs_all = contigs_all.mix(
        Assembly_Workflow.out.dragonflye_raven_assembly_ch.map { sid, fa -> tuple('raven', sid, fa) }
      )
    }	
    if (params.myloasm) {
      contigs_all = contigs_all.mix(
        Assembly_Workflow.out.myloasm_assembly_ch.map { sid, fa -> tuple('myloasm', sid, fa) }
      )
    }


    def mmseqs_dragonflye_contigs_ch = Channel.empty()
    def mmseqs_dragonflye_medaka_contigs_ch = Channel.empty()
    def mmseqs_metaspades_contigs_ch = Channel.empty()
    def mmseqs_unicycler_contigs_ch = Channel.empty()
    def mmseqs_raven_contigs_ch = Channel.empty()


    if (params.run_mmseqs) {
    MMSEQ_AllAssemblers(contigs_all, nt_db_ch)

    def dragonCh = MMSEQ_AllAssemblers.out.dragon ?: Channel.empty()
    def medakaCh = MMSEQ_AllAssemblers.out.medaka ?: Channel.empty()
    def metaCh = MMSEQ_AllAssemblers.out.meta ?: Channel.empty()
    def uniCh = MMSEQ_AllAssemblers.out.uni ?: Channel.empty()
    def ravenCh = MMSEQ_AllAssemblers.out.raven ?: Channel.empty()

    def mmseq_parse_in = Channel.empty()
    .mix( dragonCh.map { sid, tsv -> tuple(sid, 'dragonflye', tsv) } )
    .mix( medakaCh.map { sid, tsv -> tuple(sid, 'dragonflye_medaka', tsv) } )
    .mix( metaCh .map { sid, tsv -> tuple(sid, 'metaspades', tsv) } )
    .mix( uniCh .map { sid, tsv -> tuple(sid, 'unicycler', tsv) } )
    .mix( ravenCh .map { sid, tsv -> tuple(sid, 'raven', tsv) } )

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

    mmseqs_raven_contigs_ch =
    Parse_MMSEQ_Contigs.out.mmseqs_parsed_ch
    .filter { sid, asm, f -> asm == 'raven' }
    .map { sid, asm, f -> tuple(sid, f) }

}

if (params.vs) {
    def reads_min_ch = assembly_input_ch.map { sid, fq1, fq2, lr, q2, mode ->
        tuple(sid, fq1, fq2, lr)
    }

    def blastx_dragonflye_parsed_ch = Parse_BlastX_Contigs.out.parse_blastx_ch
        .filter { sid, assembler, parsed -> assembler == 'dragonflye' }
        .map    { sid, assembler, parsed -> tuple(sid, parsed) }

    def blastx_medaka_parsed_ch = Parse_BlastX_Contigs.out.parse_blastx_ch
        .filter { sid, assembler, parsed -> assembler == 'dragonflye_medaka' || assembler == 'medaka' }
        .map    { sid, assembler, parsed -> tuple(sid, parsed) }

    def blastx_metaspades_parsed_ch = Parse_BlastX_Contigs.out.parse_blastx_ch
        .filter { sid, assembler, parsed -> assembler == 'metaspades' }
        .map    { sid, assembler, parsed -> tuple(sid, parsed) }

    def blastx_unicycler_parsed_ch = Parse_BlastX_Contigs.out.parse_blastx_ch
        .filter { sid, assembler, parsed -> assembler == 'unicycler' }
        .map    { sid, assembler, parsed -> tuple(sid, parsed) }

    def blastx_raven_parsed_ch = Parse_BlastX_Contigs.out.parse_blastx_ch
        .filter { sid, assembler, parsed -> assembler == 'raven' }
        .map    { sid, assembler, parsed -> tuple(sid, parsed) }


    def sample_keys_for_vs = reads_min_ch
        .map { sid, fq1, fq2, lr ->
            tuple(sid, true)
        }

    def short_daa_default_ch = sample_keys_for_vs
        .map { sid, dummy ->
            tuple(sid, null)
        }

    def long_daa_default_ch = sample_keys_for_vs
        .map { sid, dummy ->
            tuple(sid, null)
        }

    def short_daa_real_ch = Megan_ShortReads_WF.out.meganized_short_daa_ch
        .map { sid, daa ->
            tuple(sid, daa)
        }

    def long_daa_real_ch = Megan_LongReads_WF.out.meganized_long_daa_ch
        .map { sid, daa ->
            tuple(sid, daa)
        }

    def short_daa_any_ch = short_daa_default_ch
        .mix(short_daa_real_ch)
        .groupTuple()
        .map { sid, vals ->
            def got = vals.find { it != null } ?: null
            tuple(sid, got)
        }

    def long_daa_any_ch = long_daa_default_ch
        .mix(long_daa_real_ch)
        .groupTuple()
        .map { sid, vals ->
            def got = vals.find { it != null } ?: null
            tuple(sid, got)
        }

    def reads_daa_pair_ch = short_daa_any_ch
        .join(long_daa_any_ch)
        .map { sid, short_daa, long_daa ->
            tuple(sid, short_daa, long_daa)
        }

    reads_daa_pair_ch.view {
        "VS READS DAA PAIR: $it"
    }


    def vs_calls_inputs = Channel.empty()


    if (params.dragonflye) {
        def vs_dragon = reads_min_ch
            .join(Assembly_Workflow.out.dragonflye_assembly_ch)
            .join(mmseqs_dragonflye_contigs_ch)
            .join(blastx_dragonflye_parsed_ch)
            .join(reads_daa_pair_ch)
            .map { sid,
                   fq1,
                   fq2,
                   lr,
                   contigs,
                   mmseqs_parsed,
                   diamond_parsed,
                   short_daa,
                   long_daa ->

                tuple(
                    sid,
                    fq1,
                    fq2,
                    lr,
                    contigs,
                    mmseqs_parsed,
                    'dragonflye',
                    diamond_parsed,
                    short_daa,
                    long_daa
                )
            }

        vs_calls_inputs = vs_calls_inputs.mix(vs_dragon)
    }


    if (params.medaka) {
        def vs_medaka = reads_min_ch
            .join(Assembly_Workflow.out.dragonflye_medaka_assembly_ch)
            .join(mmseqs_dragonflye_medaka_contigs_ch)
            .join(blastx_medaka_parsed_ch)
            .join(reads_daa_pair_ch)
            .map { sid,
                   fq1,
                   fq2,
                   lr,
                   contigs,
                   mmseqs_parsed,
                   diamond_parsed,
                   short_daa,
                   long_daa ->

                tuple(
                    sid,
                    fq1,
                    fq2,
                    lr,
                    contigs,
                    mmseqs_parsed,
                    'medaka',
                    diamond_parsed,
                    short_daa,
                    long_daa
                )
            }

        vs_calls_inputs = vs_calls_inputs.mix(vs_medaka)
    }


    if (params.metaspades) {
        def vs_meta = reads_min_ch
            .join(Assembly_Workflow.out.metaspades_assembly_ch)
            .join(mmseqs_metaspades_contigs_ch)
            .join(blastx_metaspades_parsed_ch)
            .join(reads_daa_pair_ch)
            .map { sid,
                   fq1,
                   fq2,
                   lr,
                   contigs,
                   mmseqs_parsed,
                   diamond_parsed,
                   short_daa,
                   long_daa ->

                tuple(
                    sid,
                    fq1,
                    fq2,
                    lr,
                    contigs,
                    mmseqs_parsed,
                    'metaspades',
                    diamond_parsed,
                    short_daa,
                    long_daa
                )
            }

        vs_calls_inputs = vs_calls_inputs.mix(vs_meta)
    }


    if (params.unicycler) {
        def vs_uni = reads_min_ch
            .join(Assembly_Workflow.out.unicycler_assembly_ch)
            .join(mmseqs_unicycler_contigs_ch)
            .join(blastx_unicycler_parsed_ch)
            .join(reads_daa_pair_ch)
            .map { sid,
                   fq1,
                   fq2,
                   lr,
                   contigs,
                   mmseqs_parsed,
                   diamond_parsed,
                   short_daa,
                   long_daa ->

                tuple(
                    sid,
                    fq1,
                    fq2,
                    lr,
                    contigs,
                    mmseqs_parsed,
                    'unicycler',
                    diamond_parsed,
                    short_daa,
                    long_daa
                )
            }

        vs_calls_inputs = vs_calls_inputs.mix(vs_uni)
    }


    if (params.raven) {
        def vs_raven = reads_min_ch
            .join(Assembly_Workflow.out.dragonflye_raven_assembly_ch)
            .join(mmseqs_raven_contigs_ch)
            .join(blastx_raven_parsed_ch)
            .join(reads_daa_pair_ch)
            .map { sid,
                   fq1,
                   fq2,
                   lr,
                   contigs,
                   mmseqs_parsed,
                   diamond_parsed,
                   short_daa,
                   long_daa ->

                tuple(
                    sid,
                    fq1,
                    fq2,
                    lr,
                    contigs,
                    mmseqs_parsed,
                    'raven',
                    diamond_parsed,
                    short_daa,
                    long_daa
                )
            }

        vs_calls_inputs = vs_calls_inputs.mix(vs_raven)
    }

    vs_calls_inputs.view {
        "vs_calls_input (all assemblers): $it"
    }

    VS_Workflow(vs_calls_inputs)
}
}