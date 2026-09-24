nextflow.enable.dsl=2

/*
===============================================================================
 processes_autocycler.nf

===============================================================================
*/

/*  EstimateGenomeSize  */
process EstimateGenomeSize {
  errorStrategy 'ignore'
  tag { sample_id }
  label 'small'
  publishDir { "${params.outdir}/${params.run_id}/${sample_id}/00_genome_size" }, mode: 'copy'

  input:
  tuple val(sample_id), path(lr), val(mode)

  output:
  // Emit keyed TSV; workflow will parse GS and re-key
  tuple val(sample_id), path('genome_size.txt')

  script:
  """
  autocycler helper genome_size --reads ${lr} --threads ${task.cpus} > genome_size.txt 2>&1
  """
}

/*  SubsampleReads  */
process SubsampleReads {
  tag { sample_id }
  label 'small'
  errorStrategy 'ignore'
  publishDir { "${params.outdir}/${params.run_id}/${sample_id}/01_subsamples" }, mode: 'copy'

  input:
  tuple val(sample_id), path(lr), val(mode), val(genome_size)

  output:
  tuple val(sample_id), path('*.fastq*'), val(genome_size)

  script:
  """
  echo "INFO: Trying autocycler subsample (skip_subsample=${params.skip_subsample ?: false}, count=${params.subsample_count ?: 1}, min_depth=${params.min_read_depth ?: 1})" >&2

  if ! ${params.skip_subsample ?: false} && \
     autocycler subsample \
       --reads ${lr} \
       --out_dir . \
       --genome_size ${genome_size} \
       --count ${params.subsample_count ?: 5} \
       --min_read_depth ${params.min_read_depth ?: 25}
  then
    echo "INFO: autocycler subsample succeeded" >&2
  else
    echo "WARN: Subsample skipped or failed; emitting pass-through subset." >&2
    out="${sample_id}_sub01.fastq.gz"
    if [[ "${lr}" == *.gz ]]; then
      ln -s "${lr}" "\$out"
    else
      gzip -c "${lr}" > "\$out"
    fi
  fi
  """
}

/*  AssembleFlye  */
process AssembleFlye {
  tag { "${sample_id}_flye_sub${subid}" }
  //cpus { 16 }
  //memory { '32 GB'}
  //time '36h'
  errorStrategy 'ignore'
  label 'autocycler_mem'
  publishDir { "${params.outdir}/${params.run_id}/${sample_id}/02_assemblies/flye/sub${subid}" }, mode: 'copy'
  input:
  tuple val(sample_id), val(subid), path(subreads), val(genome_size)

  output:
  tuple val(sample_id), val('flye'), val(subid), path("${sample_id}_sub${subid}_flye.gfa"),   optional: true, emit: gfa
  tuple val(sample_id), val('flye'), val(subid), path("${sample_id}_sub${subid}_flye.fasta"), optional: true, emit: fa
  tuple val(sample_id), val('flye'), val(subid), path("${sample_id}_sub${subid}_flye.log"), optional: true, emit: log

  script:
  """
  autocycler helper flye \
    --reads ${subreads} \
    --genome_size ${genome_size} \
    --threads ${task.cpus} \
    --out_prefix ${sample_id}_sub${subid}_flye
  """
}

/*  AssembleMyloasm  */
process AssembleMyloasm {
  tag { "${sample_id}_myloasm_sub${subid}" }
  //cpus { 64 }
  //memory { '128 GB'}
  //time '36h'
  errorStrategy 'ignore'
  label 'autocycler_mem'
  publishDir { "${params.outdir}/${params.run_id}/${sample_id}/02_assemblies/myloasm/sub${subid}" }, mode: 'copy'

  input:
  tuple val(sample_id), val(subid), path(subreads), val(genome_size)

  output:
  tuple val(sample_id), val('myloasm'), val(subid), path("${sample_id}_sub${subid}_myloasm.gfa"),   optional: true, emit: gfa
  tuple val(sample_id), val('myloasm'), val(subid), path("${sample_id}_sub${subid}_myloasm.fasta"), optional: true, emit: fa
  tuple val(sample_id), val('myloasm'), val(subid), path("${sample_id}_sub${subid}_myloasm.log"), optional: true, emit: log

  script:
  """
  autocycler helper myloasm \
    --reads ${subreads} \
    --genome_size ${genome_size} \
    --threads ${task.cpus} \
    --out_prefix ${sample_id}_sub${subid}_myloasm
  """
}

/*  AssembleMiniasm  */
process AssembleMiniasm {
  tag { "${sample_id}_miniasm_sub${subid}" }
  //cpus { 64 }
  //memory { '128 GB'}
  //time '36h'
  errorStrategy 'ignore'
  label 'autocycler_mem'
  publishDir { "${params.outdir}/${params.run_id}/${sample_id}/02_assemblies/miniasm/sub${subid}" }, mode: 'copy'

  input:
  tuple val(sample_id), val(subid), path(subreads), val(genome_size)

  output:
  tuple val(sample_id), val('miniasm'), val(subid), path("${sample_id}_sub${subid}_miniasm.gfa"),   optional: true, emit: gfa
  tuple val(sample_id), val('miniasm'), val(subid), path("${sample_id}_sub${subid}_miniasm.fasta"), optional: true, emit: fa
  tuple val(sample_id), val('miniasm'), val(subid), path("${sample_id}_sub${subid}_miniasm.log"), optional: true, emit: log

  script:
  """
  autocycler helper miniasm \
    --reads ${subreads} \
    --genome_size ${genome_size} \
    --threads ${task.cpus} \
    --out_prefix ${sample_id}_sub${subid}_miniasm
  """
}

/*  CompressAssemblies  */
process CompressAssemblies {
  tag { sample_id }
  label 'medium'
  errorStrategy 'ignore'
  publishDir { "${params.outdir}/${params.run_id}/${sample_id}/03_compress" }, mode: 'copy'

  input:
  tuple val(sample_id), val(fa_list)   // fa_list is list of assembly.fasta paths

  output:
  // Emit a pointer to the base Autocycler working directory (here: the CWD)
  tuple val(sample_id), path('autocycler_dir')

  script:
  def minlen = (params.ac_min_contig_len ?: 1000) as int
  def maxn   = (params.ac_max_contigs_per_assembly ?: 25) as int
  def faLines = fa_list.collect { it.toString() }.join('\n')
  """
  mkdir -p assemblies assemblies_filt
  cat > fa_inputs.txt <<'EOF'
${faLines}
EOF

  copied=0
  while IFS= read -r f; do
    [[ -n "\$f" && -e "\$f" ]] || continue
    cp "\$f" assemblies/\$(basename "\$f")
    copied=\$((copied+1))
  done < fa_inputs.txt

  # Filter to top MAX contigs by length, keep only >= MIN
  for f in assemblies/*.fasta; do
    out="assemblies_filt/\$(basename "\$f")"
    awk -v MIN=${minlen} -v MAX=${maxn} '
      BEGIN{RS=">"; ORS=""; nrec=0}
      NR>1{
        n = split(\$0, lines, "\\n")
        head=lines[1]
        seq=""
        for(i=2;i<=n;i++){ gsub(/[^ACGTNacgtn]/,"",lines[i]); seq=seq lines[i] }
        len=length(seq)
        if(len>=MIN){ nrec++; H[nrec]=head; L[nrec]=len; S[nrec]=seq }
      }
      END{
        # simple sort by length desc
        for(i=1;i<=nrec;i++) for(j=i+1;j<=nrec;j++) if(L[j]>L[i]){
          t=L[i]; L[i]=L[j]; L[j]=t
          t=H[i]; H[i]=H[j]; H[j]=t
          t=S[i]; S[i]=S[j]; S[j]=t
        }
        limit=(nrec<MAX)?nrec:MAX
        for(i=1;i<=limit;i++) print ">"H[i]"\\n"S[i]"\\n"
      }
    ' "\$f" > "\$out"
  done

  # Run in *current directory* as the autocycler dir
  autocycler compress -i assemblies_filt -a autocycler_dir || true
  """
}


/*  ClusterAssemblies  */
process ClusterAssemblies {
  tag { sample_id }
  label 'medium'
  errorStrategy 'ignore'
  publishDir { "${params.outdir}/${params.run_id}/${sample_id}/04_cluster" }, mode: 'copy'

  input:
  tuple val(sample_id), path(ac_dir)

  output:
  // Re-emit a local pointer to the same base Autocycler directory
  tuple val(sample_id), path('autocycler_dir/clustering')

  script:
  """
  autocycler cluster --autocycler_dir ${ac_dir}
  """
}


/*  TrimClusters  */
process TrimClusters {
  tag { sample_id }
  label 'small'
  errorStrategy 'ignore'
  publishDir { "${params.outdir}/${params.run_id}/${sample_id}/05_trim" }, mode: 'copy'

  input:
  tuple val(sample_id), path(ac_dir)   // base Autocycler dir from previous step

  output:
  // Re-emit a *fresh* handle created by THIS task (not the staged input)
  tuple val(sample_id), path('clustering')

  script:
  """
  for c in ${ac_dir}/qc_pass/cluster_*; do
    autocycler trim -c "\$c"
  done

  """
}


/*  ResolveClusters  */
process ResolveClusters {
  tag { sample_id }
  label 'small'
  errorStrategy 'ignore'
  publishDir { "${params.outdir}/${params.run_id}/${sample_id}/06_resolve" }, mode: 'copy'

  input:
  tuple val(sample_id), path(ac_dir)

  output:
  tuple val(sample_id), path('clustering')

  script:
  """
  for c in ${ac_dir}/qc_pass/cluster_*; do
    autocycler resolve -c "\$c"
  done
  """
}

/*  CombineResolved  */
process CombineResolved {
  tag { sample_id }
  label 'tiny'
  errorStrategy 'ignore'
  publishDir { "${params.outdir}/${params.run_id}/${sample_id}/07_final" }, mode: 'copy'

  input:
  tuple val(sample_id), path(resolved_dirs)

  output:
  tuple val(sample_id), path('autocycler_out'), optional: true, emit: final_dir
  tuple val(sample_id), file('consensus_assembly.fasta'), optional: true, emit: final_fa
  tuple val(sample_id), file('consensus_assembly.gfa'), optional: true, emit: final_gfa
  tuple val(sample_id), file('consensus_assembly.yaml'), optional: true, emit: final_yaml

  script:
  """
  autocycler combine -a autocycler_out -i ${resolved_dirs}/qc_pass/cluster_*/5_final.gfa
  """
}
