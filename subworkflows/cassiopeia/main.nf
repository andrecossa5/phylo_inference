// CASSIOPEIA

// Include here
nextflow.enable.dsl = 2
include { BOOTSTRAP } from "./modules/bootstrap.nf"
include { BUILD_CASSIOPEIA } from "./modules/build_cassiopeia.nf"
include { SUPPORT } from "./modules/support.nf"

// 
 
//----------------------------------------------------------------------------//
// CASSIOPEIA subworkflow
//----------------------------------------------------------------------------//

//

// CASSIOPEIA subworkflow
workflow CASSIOPEIA {
    
    take:
        ch_input   

    main: 

        ch_bootstrap = ch_input
                        .combine(params.boot_strategy)
                        .combine(
                            Channel.of( 1..params.n_boot )
                            .concat(Channel.of( "observed" ))
                        )
        BOOTSTRAP(ch_bootstrap)
        BUILD_CASSIOPEIA(BOOTSTRAP.out.bootstrapped_input.combine(params.solver))
        SUPPORT(BUILD_CASSIOPEIA.out.tree.groupTuple(by: [0,1,3,5]))

    emit:

        results = SUPPORT.out.results

}


