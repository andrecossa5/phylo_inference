// preprocess

// Include here
nextflow.enable.dsl = 2
include { MITO } from "./modules/mito.nf"
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
     
        // Process input AFMs, based on input scLT_system
        if (params.scLT_system == "MAESTER" || params.scLT_system == "RedeeM") {
            MITO(ch_jobs)   
            ch_afm = MITO.out.afm
        }
        else {
            println('Only valid option are MAESTER and RedeeM so far...')
            // println('Provide valid tracing system option! (e.g., MAESTER, RedeeM, scmtATAC, Cas9, scWGS)')
        } 
        
        // Calculate distances
        replicates = Channel.of( 1..(params.n_boot_distances-1) ).concat(Channel.of( "observed")) 
        DISTANCES(ch_afm.combine(replicates))
        DISTANCE_METRICS(DISTANCES.out.distances.groupTuple(by: [0,1]))

    emit:

        input = DISTANCES.out.distances
        
} 
