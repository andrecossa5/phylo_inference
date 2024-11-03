// preprocess

// Include here
nextflow.enable.dsl = 2
include { MAESTER } from "./modules/maester.nf"
include { DISTANCES } from "./modules/distances.nf"
include { DISTANCE_METRICS } from "./modules/distance_metrics.nf"

//
 
//----------------------------------------------------------------------------//
// preprocess subworkflow
//----------------------------------------------------------------------------//

// preprocess subworkflow
workflow preprocess {
    
    take:
        ch_jobs 

    main:
     
        // Process input channel based on tracer system
        if (params.tracing_system == "MAESTER") {
            MAESTER(ch_jobs)   
            ch_afm = MAESTER.out.afm
        } 
        
        // ...
        
        else {
            println('Only valid option is MAESTER so far...')
            // println('Provide valid tracing system option! (e.g., MAESTER, RedeeM, scmtATAC, Cas9, scWGS)')
        }

        // Calculate distances
        replicates = Channel.of( 1..(params.n_boot_distances-1) ).concat(Channel.of( "observed" )) 
        DISTANCES(ch_afm.combine(replicates))
        DISTANCE_METRICS(DISTANCES.out.distances.groupTuple(by: [0,1]))

    emit:

        input = DISTANCES.out.distances
        
} 
