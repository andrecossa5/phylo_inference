// mito_unsupervised

// Include here
nextflow.enable.dsl = 2
include { CASSIOPEIA } from "./modules/cassiopeia.nf"

// 
 
//----------------------------------------------------------------------------//
// CASSIOPEIA subworkflow
//----------------------------------------------------------------------------//

workflow cassiopeia_workflow {
    
    take:
        ch_samples 

    main:

        // Create options
        options = ch_samples
        .combine(params.filtering)
        .combine(params.metric)
        .combine(params.solver)
        .combine(params.boot_strategy)

        CASSIOPEIA(options)

    emit:
        results = CASSIOPEIA.out.results

}
