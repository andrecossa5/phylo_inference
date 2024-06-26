// IQTREE

// Include here
nextflow.enable.dsl = 2
include { BUILD_IQTREE } from "./modules/build_iqtree.nf"

// 
 
//----------------------------------------------------------------------------//
// IQTREE subworkflow
//----------------------------------------------------------------------------//

//

// IQTREE subworkflow
workflow IQTREE {
    
    take:
        ch_input   

    main: 
 
        BUILD_IQTREE(ch_input)

    emit:
        
        tree = BUILD_IQTREE.out.tree

}


  