// FILTER_VARIANTS

// Include here
nextflow.enable.dsl = 2
include { FILTER } from "./modules/filter_with_metrics.nf"

// 
 
//----------------------------------------------------------------------------//
// FILTER_VARIANTS subworkflow
//----------------------------------------------------------------------------//

// FILTER_VARIANTS subworkflow
workflow FILTER_VARIANTS {
    
    take:
        ch_samples  

    main:

        // Prep each job input
        def counter = 1 
        ch_jobs = ch_samples
            .combine(params.min_site_cov)
            .combine(params.min_var_quality)
            .combine(params.min_frac_negative)
            .combine(params.min_n_positive)
            .combine(params.low_confidence_af)
            .combine(params.high_confidence_af)
            .combine(params.min_prevalence_low_confidence_af)
            .combine(params.min_cells_high_confidence_af)
            .filter { it -> it[6] > 10 * it[5] }      
            .map { it -> 
                def result = [counter++, *it] 
                result 
            }
  
        // Run each job 
        FILTER(ch_jobs)
            
    emit:
        // results = FILTER.out.stats
        results = ch_jobs

}