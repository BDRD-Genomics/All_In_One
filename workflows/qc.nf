#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {
    Pretrim_fastqc_merged;
    Pretrim_NanoPlot;
    Quality_Control;
    PoreChop;
    Fastp_LongReads;
    Post_trim_fastqc;
    Multiqc_QC_Stats;
    Interleave;
} from './modules/local/qc/main.nf'


workflow QC_Workflow {

    take:
    input_samples_ch

    main:

    Interleave(input_samples_ch)


    Pretrim_fastqc_merged(
        input_samples_ch,
        Interleave.out.interleave_ch
    )


    Pretrim_NanoPlot(input_samples_ch)

    nanoplot_pre_ch = Pretrim_NanoPlot.out.pretrim_nanoplot_ch

    Quality_Control(input_samples_ch)


    PoreChop(input_samples_ch)


    Fastp_LongReads(PoreChop.out.porechop_ch)

    trimmed_short_ch = Quality_Control.out.trimmed_short_ch

    interleaved_trimmed_short_ch = Quality_Control.out.interleaved_trimmed_short_ch

    trimmed_long_ch = Fastp_LongReads.out.fastp_long_ch


    def posttrim_input_ch =
        interleaved_trimmed_short_ch
            .map { sid, fq -> tuple(sid, fq) }
            .mix(
                trimmed_long_ch
                    .map { sid, fq -> tuple(sid, fq) }
            )
            .groupTuple()
            .map { sid, files -> tuple(sid, files) }


    Post_trim_fastqc(posttrim_input_ch)


    Multiqc_QC_Stats(
        Pretrim_fastqc_merged.out.pretrim_fastqc_ch.collect(),
        Post_trim_fastqc.out.posttrim_fastqc_ch.collect()
    )



    emit:

    trimmed_short_ch = trimmed_short_ch

    trimmed_long_ch = trimmed_long_ch

    interleaved_trimmed_short_ch = interleaved_trimmed_short_ch

    nanoplot_pre_ch = nanoplot_pre_ch
}
