#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {
    Pretrim_fastqc_merged;
    Quality_Control;
    PoreChop;
    Post_trim_fastqc;
    Multiqc_QC_Stats;
    Interleave
} from './modules/local/qc/main.nf'

workflow QC_Workflow {
    take:
    input_samples_ch

    main:
    Interleave(input_samples_ch)
    Pretrim_fastqc_merged(input_samples_ch,Interleave.out.interleave_ch)
    Quality_Control(input_samples_ch)
    PoreChop(input_samples_ch)
    trimmed_short_ch = Quality_Control.out.trimmed_short_ch
    trimmed_long_ch  = PoreChop.out.porechop_ch
    interleaved_trimmed_short_ch  = Quality_Control.out.interleaved_trimmed_short_ch
    // Create channel for fastqc_post_trim
    // interleaved_trimmed_short_ch : (sid, fq)
    // trimmed_long_ch              : (sid, fq)

    def posttrim_input_ch =
        interleaved_trimmed_short_ch
            .map { sid, fq -> tuple(sid, file(fq)) }            // (sid, path)
            .mix(
                trimmed_long_ch
                    .map { sid, fq -> tuple(sid, file(fq)) }    // (sid, path)
            )
            .groupTuple()                                       // (sid, [path, path, ...])
            .map { sid, files -> tuple(sid, files) }            // (sid, List<path>)
    Post_trim_fastqc(posttrim_input_ch)
    Multiqc_QC_Stats(Pretrim_fastqc_merged.out.pretrim_fastqc_ch.collect(),Post_trim_fastqc.out.posttrim_fastqc_ch.collect())
    // Each FastQC emits: (sample_id, file("*"))
    // Strip sample_id and COLLECT to a single emission (List<path>)
    // PRE-trim: (sid, file("*")) where the second element can be a single Path or a List<Path>
    //def fastqc_pre_trim_files_ch = Pretrim_fastqc_merged.out.pretrim_fastqc_ch
    //    .flatMap { sid, f -> (f instanceof List) ? f : [f] }  // emit each Path
    //    .collect()                                            // => one emission: List<Path>

    // POST-trim:
    //def fastqc_post_trim_files_ch = Post_trim_fastqc.out.posttrim_fastqc_ch
    //    .flatMap { sid, f -> (f instanceof List) ? f : [f] }
    //    .collect()

    // Run MultiQC once, in one step
    //Multiqc_QC_Stats(fastqc_pre_trim_files_ch, fastqc_post_trim_files_ch)


    //Read_Distribution(input_samples_ch,interleaved_trimmed_short_ch,trimmed_long_ch)
    emit:
    trimmed_short_ch = trimmed_short_ch
    trimmed_long_ch  = trimmed_long_ch
    interleaved_trimmed_short_ch = interleaved_trimmed_short_ch
}
