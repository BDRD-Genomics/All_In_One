#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {
    Spades;
    MetaSPAdes_Assembly;
    Plasmid_Spades;
    Dragonflye as Dragonflye_LR;
    Dragonflye as Dragonflye_Hybrid;
    Dragonflye_Medaka as Dragonflye_Medaka_LR;
    Dragonflye_Medaka as Dragonflye_Medaka_Hybrid;
    Unicycler_Assembly as Unicycler_Assembly_Short;
    Unicycler_Assembly as Unicycler_Assembly_Hybrid;
    Dragonflye_Raven as Dragonflye_Raven_LR;
    Dragonflye_Raven as Dragonflye_Raven_Hybrid;
    Myloasm
} from './modules/local/assembler/main.nf'

workflow Assembly_Workflow {

    take:
    assembly_input_ch

    main:

    assembly_input_ch.view { "Assembly input (received): $it" }

    short_input_ch  = assembly_input_ch.filter { it[5] in ['short', 'hybrid'] }
    long_input_ch   = assembly_input_ch.filter { it[5] == 'long' }
    hybrid_input_ch = assembly_input_ch.filter { it[5] == 'hybrid' }

    if (params.spades) {
        short_input_ch
            .map { sample_id, fq1, fq2, lr, q2, mode -> tuple(sample_id, fq1, fq2) }
            .set { spades_input_ch }

        Spades(spades_input_ch)
    }

    if (params.plasmidspades) {
        short_input_ch
            .map { sample_id, fq1, fq2, lr, q2, mode -> tuple(sample_id, fq1, fq2) }
            .set { plasmid_spades_input_ch }

        Plasmid_Spades(plasmid_spades_input_ch)
    }

    //if (params.metaspades) {
    //    short_input_ch
    //        .map { sample_id, fq1, fq2, lr, q2, mode -> tuple(sample_id, fq1, fq2) }
    //        .set { metaspades_input_ch }
    //    MetaSPAdes_Assembly(metaspades_input_ch)
    //}

    if (params.metaspades) {
        // short_input_ch is assumed (sid, fq1, fq2, lr, mode)
        def metaspades_input_ch = short_input_ch.map { sample_id, fq1, fq2, lr, q2, mode ->
            long bytes = (fq1 && fq1.exists()) ? fq1.size() : 0L
            double gb  = bytes / (1024.0 * 1024 * 1024)

            int cpus   = (gb > 15.0) ? 254 : params.metaspades_default_cpus
            String mem = (gb > 15.0) ? '996 GB' : params.metaspades_default_mem

            // emit (sid, fq1, fq2, cpus, mem)
            tuple(sample_id, fq1, fq2, cpus, mem)
        }

        // optional: visibility
        metaspades_input_ch.view { sid, fq1, fq2, cpus, mem ->
            def sz = (fq1?.exists() ? String.format('%.2f', fq1.size()/(1024.0*1024*1024)) : '0.00')
            "MetaSPAdes sizing: ${sid} fq1=${sz} GB  ->  cpus=${cpus} mem=${mem}"
        }

	metaspades_input_ch
		.map { sid, fq1, fq2, cpus, mem -> tuple(sid, fq1, fq2, cpus, mem) }
		.set { test_metaspades_ch }		

        MetaSPAdes_Assembly(metaspades_input_ch)
	//MetaSPAdes_Assembly(short_input_ch)

    }


    def dragonflye_ch = Channel.empty()

    if (params.dragonflye) {
        if (long_input_ch) {
            def dragonflye_input_ch_long = 
                   long_input_ch
                     .map { t -> 
		      def ( sample_id, fq1, fq2, lr, q2, mode ) = (t as List )
                      tuple(sample_id, null, null, lr, q2) }
            dragonflye_input_ch_long.view { "Dragonflye Input Channel: $it" }
            Dragonflye_LR(dragonflye_input_ch_long)
            dragonflye_ch = Dragonflye_LR.out.dragonflye_assembly_ch
        }

        if (hybrid_input_ch) {
	     hybrid_input_ch
                .map { sample_id, fq1, fq2, lr, q2, mode -> tuple(sample_id, fq1, fq2, lr, q2) }
                .set { hybrid_dragonflye_input_ch }

            Dragonflye_Hybrid(hybrid_dragonflye_input_ch)
            dragonflye_ch = dragonflye_ch.mix(Dragonflye_Hybrid.out.dragonflye_assembly_ch)
        }
    }

    def dragonflye_medaka_ch = Channel.empty()  

    if (params.medaka) {
        if (long_input_ch) {
            def dragonflye_medaka_input_long =
            long_input_ch
                .map { t -> 
                def ( sample_id, fq1, fq2, lr, q2, mode ) = ( t as List ) 
                tuple(sample_id, null, null, lr, q2) }

            Dragonflye_Medaka_LR(dragonflye_medaka_input_long)
            dragonflye_medaka_ch = Dragonflye_Medaka_LR.out.dragonflye_assembly_medaka_ch
        }   

        if (hybrid_input_ch) {
            hybrid_input_ch
                .map { sample_id, fq1, fq2, lr, q2, mode -> tuple(sample_id, fq1, fq2, lr, q2) }
                .set { dragonflye_medaka_input_hybrid } 

            Dragonflye_Medaka_Hybrid(dragonflye_medaka_input_hybrid)
            dragonflye_medaka_ch = dragonflye_medaka_ch.mix(Dragonflye_Medaka_Hybrid.out.dragonflye_assembly_medaka_ch)
        }
    }

    def dragonflye_raven_ch = Channel.empty()

    if (params.raven) {
        if (long_input_ch) {
            def raven_input_long =
            long_input_ch
                .map { t ->
                def ( sample_id, fq1, fq2, lr, q2, mode ) = ( t as List )
                tuple(sample_id, null, null, lr, q2) }

            Dragonflye_Raven_LR(raven_input_long)
            dragonflye_raven_ch = Dragonflye_Raven_LR.out.dragonflye_assembly_raven_ch
        }

        if (hybrid_input_ch) {
            hybrid_input_ch
                .map { sample_id, fq1, fq2, lr, q2, mode -> tuple(sample_id, fq1, fq2, lr, q2) }
                .set { raven_input_hybrid }

            Dragonflye_Raven_Hybrid(raven_input_hybrid)
            dragonflye_raven_ch = dragonflye_raven_ch.mix(Dragonflye_Raven_Hybrid.out.dragonflye_assembly_raven_ch)
        }
    }
    

    def myloasm_ch = Channel.empty()

    if (params.myloasm) {
        def myloasm_input_ch = assembly_input_ch
            .filter { sample_id, fq1, fq2, lr, q2, mode ->
                mode in ['long', 'hybrid'] && lr != null
            }
            .map { sample_id, fq1, fq2, lr, q2, mode ->
                tuple(sample_id, lr)
            }

        Myloasm(myloasm_input_ch)
        myloasm_ch = Myloasm.out.myloasm_assembly_ch
    }

    def unicycler_ch = Channel.empty()

    if (params.unicycler) {
        if (short_input_ch) {
        short_input_ch
            .map { sample_id, fq1, fq2, lr, q2, mode -> tuple(sample_id, fq1, fq2, lr, q2, mode) }
            .set { unicycler_input_ch }
        Unicycler_Assembly_Short(unicycler_input_ch)
        unicycler_ch = Unicycler_Assembly_Short.out.unicycler_assembly_ch
        }
        
        if (hybrid_input_ch) {
        hybrid_input_ch
            .map { sample_id, fq1, fq2, lr, q2, mode -> tuple(sample_id, fq1, fq2, lr, q2, mode) }
            .set { unicycler_input_ch }

        Unicycler_Assembly_Hybrid(unicycler_input_ch)
        unicycler_ch = Unicycler_Assembly_Hybrid.out.unicycler_assembly_ch
        }
    }

    emit:
    spades_assembly_ch        = params.spades        ? Spades.out.spades_assembly_ch : Channel.empty()
    metaspades_assembly_ch    = params.metaspades    ? MetaSPAdes_Assembly.out.metaspades_assembly_ch : Channel.empty()
    dragonflye_assembly_ch    = dragonflye_ch
    unicycler_assembly_ch     = unicycler_ch
    plasmidspades_assembly_ch = params.plasmidspades ? Plasmid_Spades.out.plasmidspades_assembly_ch : Channel.empty()
    dragonflye_medaka_assembly_ch = params.medaka ? dragonflye_medaka_ch : Channel.empty()
    dragonflye_raven_assembly_ch = params.raven ? dragonflye_raven_ch : Channel.empty()
    myloasm_assembly_ch = params.myloasm ? myloasm_ch : Channel.empty()
}

