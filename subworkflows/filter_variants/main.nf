// FILTER_VARIANTS

// Include here
nextflow.enable.dsl = 2
include { FILTER_AFM } from "./modules/filter_afm.nf"

// 
 
//----------------------------------------------------------------------------//
// FILTER_VARIANTS subworkflow
//----------------------------------------------------------------------------//


// Collapse
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


//


// FILTER_VARIANTS subworkflow
workflow FILTER_VARIANTS {
    
    take:
        ch_samples  

    main:

        // Run each AFM filtering job
        def counter = 1 
        ch_jobs = ch_samples
            .combine(params.min_site_cov)
            .combine(params.min_var_quality)
            .combine(params.min_frac_negative)
            .combine(params.min_n_positive)
            .combine(params.af_confident_detection)
            .combine(params.min_n_confidently_detected)
            .combine(params.min_median_af)    
            .map { it -> 
                def result = [counter++, *it] 
                result 
            }
        FILTER_AFM(ch_jobs)

        // Collapse outputs
        ch_grouped = FILTER_AFM.out.stats
            .groupTuple(by: 0)
            .flatMap { sample, jobPaths, datasetPaths, varsPaths ->
                [
                    [sample, 'job', jobPaths],
                    [sample, 'dataset', datasetPaths],
                    [sample, 'vars', varsPaths]
                ]
            }
        collapse_output(ch_grouped)
            
    emit:
        results = collapse_output.out.csv

}