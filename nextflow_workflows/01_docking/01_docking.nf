#!/usr/bin/env nextflow
include {
    CROSS_DOCK_BY_LIGAND
    CROSS_DOCK_BY_LIGAND_MULTIPOSE
} from "./modules.nf"
params.take = -1
params.test_cache = "/data1/choderaj/paynea/asap-datasets/full_cross_dock_v2/mpro_fragalysis-04-01-24_curated_cache_fixed"
params.test_dir = "/data1/choderaj/paynea/asap-datasets/full_cross_dock_v2/mpro_fragalysis-04-01-24_curated"

workflow RUN_LIGAND_DOCKING {
    take:
        num_poses
        posit_method
        pairwise_selector

    main:
        // load in input structure dir
        // input_dir = Channel.fromPath("${params.curatedFragalysis}", type: 'dir')
        input_dir = Channel.fromPath("${params.test_dir}", type: 'dir')
    //     cache_dir = Channel.fromPath("${params.dataPath}/${params.fixedFragalysisCache}", type: 'dir')
        cache_dir = Channel.fromPath("${params.test_cache}", type: 'dir')

        // Create a channel for each ligand file and flatten it
        ligand_files = Channel
            .fromPath("${params.ligandFiles}/${params.split2dligandFiles}/*.sdf")
            .map { file ->
                def id = file.name.toString().find(/([a-zA-Z0-9_-]+)\.sdf/) { match, code -> code }
                return tuple(id, file)
            }

        // Count ligand_files
        ligand_files.count().view { count -> "Total ligand files found: $count" }

        log.info "Taking ${params.take ?: 'all'} ligand files"

        // Create combinations for parallel processing
        docking_combinations = input_dir
            .combine(cache_dir)
            .combine(ligand_files)
            .take(params.take)

    //     Run FRED and POSIT docking in parallel for each combination
        CROSS_DOCK_BY_LIGAND(
            docking_combinations,
            posit_method,
            pairwise_selector,
            num_poses
    )
}
workflow RUN_LIGAND_DOCKING_MULTIPOSE {
    take:
        num_poses
        posit_method
        pairwise_selector

    main:
        // load in input structure dir
        // input_dir = Channel.fromPath("${params.curatedFragalysis}", type: 'dir')
        input_dir = Channel.fromPath("${params.test_dir}", type: 'dir')
    //     cache_dir = Channel.fromPath("${params.dataPath}/${params.fixedFragalysisCache}", type: 'dir')
        cache_dir = Channel.fromPath("${params.test_cache}", type: 'dir')

        // Create a channel for each ligand file and flatten it
        ligand_files = Channel
            .fromPath("${params.ligandFiles}/${params.split2dligandFiles}/*.sdf")
            .map { file ->
                def id = file.name.toString().find(/([a-zA-Z0-9_-]+)\.sdf/) { match, code -> code }
                return tuple(id, file)
            }

        // Count ligand_files
        ligand_files.count().view { count -> "Total ligand files found: $count" }

        log.info "Taking ${params.take ?: 'all'} ligand files"

        // Create combinations for parallel processing
        docking_combinations = input_dir
            .combine(cache_dir)
            .combine(ligand_files)
            .take(params.take)

    //     Run FRED and POSIT docking in parallel for each combination
        CROSS_DOCK_BY_LIGAND_MULTIPOSE(
            docking_combinations,
            posit_method,
            pairwise_selector,
            num_poses
    )
}

workflow POSIT_MULTIPOSE {
    RUN_LIGAND_DOCKING_MULTIPOSE(
        50,
        'ALL',
        'PairwiseSelector'
    )
}
workflow FRED_MULTIPOSE {
    RUN_LIGAND_DOCKING_MULTIPOSE(
        50,
        'FRED',
        'PairwiseSelector'
    )
}

workflow POSIT_SINGLE_POSE {
    RUN_LIGAND_DOCKING(
        1,
        'ALL',
        'PairwiseSelector'
    )
}

workflow FRED_SINGLE_POSE {
    RUN_LIGAND_DOCKING(
        1,
        'FRED',
        'PairwiseSelector'
    )
}

workflow {
    // Run workflows
    POSIT_MULTIPOSE()
    POSIT_SINGLE_POSE()
    FRED_SINGLE_POSE()
}