// autocycler_wf.nf
nextflow.enable.dsl=2

include { EstimateGenomeSize; SubsampleReads; AssembleFlye; AssembleMyloasm; AssembleMiniasm; CompressAssemblies; ClusterAssemblies; TrimClusters; ResolveClusters; CombineResolved } from './modules/local/autocycler/main.nf'

workflow Autocycler_Workflow {

  /* TAKE: same 6-tuple as Assembly_Workflow
     (sample_id, R1, R2, LR, Q2, mode) */
  take:
    assembly_in_ch

  /*
    ADAPT: build the channels required by the processes.
    - We only proceed for modes that have long reads (long/hybrid) and LR != null
  */
  main:
    def long_only = assembly_in_ch
      .map   { sid, r1, r2, lr, q2, mode -> tuple(sid, lr, mode) }
      .filter{ sid, lr, mode -> lr != null && (mode == 'long' || mode == 'hybrid') }

    // 0) Genome size
    EstimateGenomeSize( long_only )                 // in: (sid, lr, mode)
    // emits: (sid, genome_size.txt)
    def gs_ch = EstimateGenomeSize.out.map { sid, gs_file ->
      def lines = file(gs_file).text.readLines()
      def gs_val = lines[-1].trim()
      tuple(sid, gs_val)
    }
    // 1) Subsample reads using GS
    def gs_joined = long_only
      .join ( EstimateGenomeSize.out, by: 0 )
      .map { sid, lr, mode, gs_val -> 
        tuple(sid, lr, mode, gs_val)
    }    

    SubsampleReads(
      gs_joined
    )
    // emits: (sid, *.fastq*, genome_size)
    SubsampleReads.out.view { "SUB OUT → $it" }         // Expect: [sid, /path.fastq.gz, "4.7m"]

    // 2) Fan-out subsampled sets to assemblers
    //    Expect a separate upstream step to create (sid, subid, subreads, genome_size)
    //    If your SubsampleReads already emits one file per subset, create subids:
    def toAssemblerTuples = SubsampleReads.out.map { sid, subreads, gs ->
        def name = subreads.getName()
        def patterns = [
            /_sub(\d+)\.(?:fastq|fq)(?:\.gz)?$/,     // SampleA_sub03.fastq.gz
            /sub[_-]?(\d+)\.(?:fastq|fq)(?:\.gz)?$/, // SampleA-sub3.fastq.gz
            /[_-](\d{2,})\.(?:fastq|fq)(?:\.gz)?$/,  // SampleA-03.fastq.gz
        ]
        String subid = null
        for (p in patterns) {
            def m = (name =~ p)
            if (m) { subid = m[0][1]; break }
        }
        if (!subid) subid = '01'  // final fallback
        tuple(sid, subid, subreads, gs)
    }
    
    toAssemblerTuples.view { "ASM IN → $it" }           // Expect: [sid, "03", /path.fastq.gz, "4.7m"]

    AssembleFlye(   toAssemblerTuples )
    AssembleMyloasm(toAssemblerTuples )
    AssembleMiniasm(toAssemblerTuples )


    // 3) Collect assembler FASTAs for compression
    def fa_from_flye     = AssembleFlye.out.fa     .map { sid, asm, subid, fa -> tuple(sid, fa) }
    def fa_from_myloasm  = AssembleMyloasm.out.fa  .map { sid, asm, subid, fa -> tuple(sid, fa) }
    def fa_from_miniasm  = AssembleMiniasm.out.fa  .map { sid, asm, subid, fa -> tuple(sid, fa) }

    def all_fa = fa_from_flye
      .mix(fa_from_myloasm)
      .mix(fa_from_miniasm)
      .groupTuple()                                  // (sid, [fa, fa, ...])
      .map { sid, fas -> tuple(sid, fas) }          // (sid, List<path>)

    CompressAssemblies( all_fa )                     // emits: (sid, ac_dir.txt)

    // 4) Cluster → Trim → Resolve (in place, passing ac_dir.txt)
    ClusterAssemblies( CompressAssemblies.out )
    TrimClusters(      ClusterAssemblies.out )
    ResolveClusters(   TrimClusters.out )            // emits: (sid, resolved/)

    // 5) Combine and emit finals
    CombineResolved( ResolveClusters.out )           // emits: (sid, final.gfa?), (sid, final.fasta?)

  emit:
    final_gfa_ch = CombineResolved.out.final_gfa
    final_fa_ch  = CombineResolved.out.final_fa
}
