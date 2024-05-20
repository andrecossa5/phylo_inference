// FILTER module

nextflow.enable.dsl = 2

// Collapse output
process collapse_output {

    tag "${sample}_${type}"
    publishDir "${params.outdir}/${sample}", mode: 'copy'

    input:
        tuple val(sample), val(type), path(files)

    output:
        tuple val(sample), val(type), path("${sample}_${type}.csv"), emit: csv 

    script:
    """
    outfile="${sample}_${type}.csv"
    files=(${files})
    cat "\${files[0]}" > \$outfile
    for f in "\${files[@]:1}"; do
        tail -n +2 "\$f" >> \$outfile
    done
    """

}

