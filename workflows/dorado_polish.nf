nextflow.enable.dsl=2

include {
    Dorado_Split_BamFiles;
    Dorado_Aligner;
    Dorado_Polisher
} from './modules/local/dorado_polish/main.nf'  // adjust path if needed

workflow Split_BamFiles_Workflow {

    take:
    dorado_split_input_ch

    main:
    Dorado_Split_BamFiles(dorado_split_input_ch)  //  positional input only — FIXED

    emit:
    dorado_split_output_ch = Dorado_Split_BamFiles.out
}


workflow Dorado_Aligner_Workflow {
    take:
    dorado_split_output_ch
    dragonflye_assembly_ch

    main:
    dorado_aligner_input_ch = dorado_split_output_ch.join(dragonflye_assembly_ch)
    Dorado_Aligner(dorado_aligner_input_ch)

    emit:
    dorado_aligner_ch = Dorado_Aligner.out
}
