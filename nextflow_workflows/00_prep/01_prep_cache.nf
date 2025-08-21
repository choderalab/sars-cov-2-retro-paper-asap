#!/usr/bin/env nextflow
include {
    DEDUPLICATE_LIGANDS
    PREP_CACHE_FOR_DOCKING
    GENERATE_COMBINED_LIGAND_FILES
    GENERATE_SPLIT_LIGAND_FILES
} from "./modules.nf"

workflow {
    cache_ch = Channel.fromPath("${params.fragalysisCache}", type: 'dir', checkIfExists: true)
    frag_ch = Channel.fromPath("${params.curatedFragalysis}", type: 'dir', checkIfExists: true)
    // check to see what params.scripts is being used
    println "Using scripts path: ${params.scripts}"
    DEDUPLICATE_LIGANDS(frag_ch, cache_ch)
    PREP_CACHE_FOR_DOCKING(DEDUPLICATE_LIGANDS.out.fixed_cache)
    GENERATE_COMBINED_LIGAND_FILES(PREP_CACHE_FOR_DOCKING.out.fixed_cache)
    GENERATE_SPLIT_LIGAND_FILES(GENERATE_COMBINED_LIGAND_FILES.out.ligandFile3d, GENERATE_COMBINED_LIGAND_FILES.out.ligandFile2d)
}
