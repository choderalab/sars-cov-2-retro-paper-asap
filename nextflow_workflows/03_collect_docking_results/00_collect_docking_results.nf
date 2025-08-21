#!/usr/bin/env nextflow

include {
    CALCULATE_RMSD
    COMBINE_AND_PROCESS_RESULTS
    CONVERT_TO_DOCKING_DATA_MODEL
} from "./modules.nf"

workflow CALCULATE_RMSD_WORKFLOW {
    take:
        name

    main:
        docked_dirs = Channel
            .fromPath("${params.dockedFiles}/${name}/*docked", type: 'dir')
            .map { results_dir ->
                def id = results_dir.name.toString().find(/([a-zA-Z0-9_-]+)\_docked/) { match, code -> code }
                return tuple(name, id, results_dir)
            }
        ligand_file_3d = Channel
            .fromPath("${params.ligandFiles}/${params.ligandFile3d}", type: 'file')

        input_pairs = docked_dirs.combine(ligand_file_3d)
        CALCULATE_RMSD(input_pairs)

    emit:
        rmsd_csvs = CALCULATE_RMSD.out.rmsd_csv
}

workflow PROCESS_RESULTS_WORKFLOW {
    take:
        name
        rmsd_csvs

    main:
        input_csvs = rmsd_csvs.collect()
        name_ch = Channel.value(name)
        COMBINE_AND_PROCESS_RESULTS(
            input_csvs,
            name_ch
        )
}

// Individual dataset workflows for RMSD calculation
workflow CALCULATE_FRED_RMSD {
    CALCULATE_RMSD_WORKFLOW('FRED_1_poses')
}

workflow CALCULATE_ALL_MULTIPOSE_RMSD {
    CALCULATE_RMSD_WORKFLOW('ALL_50_poses')
}

workflow CALCULATE_ALL_SINGLE_POSE_RMSD {
    CALCULATE_RMSD_WORKFLOW('ALL_1_poses')
}

// Individual dataset workflows for processing results
workflow PROCESS_FRED_RESULTS {
    PROCESS_RESULTS_WORKFLOW('FRED_1_poses', Channel.fromPath("${params.dockedLigandRMSDs}/FRED_1_poses/*.csv"))
}

workflow PROCESS_ALL_SINGLE_POSE_RESULTS {
    PROCESS_RESULTS_WORKFLOW('ALL_1_poses', Channel.fromPath("${params.dockedLigandRMSDs}/ALL_1_poses/*.csv"))
}

workflow PROCESS_ALL_MULTIPOSE_RESULTS {
    PROCESS_RESULTS_WORKFLOW('ALL_50_poses', Channel.fromPath("${params.dockedLigandRMSDs}/ALL_50_poses/*.csv"))
}

workflow CALCULATE_FRED_MULTIPOSE_RMSD {
    CALCULATE_RMSD_WORKFLOW('FRED_50_poses')
}

workflow PROCESS_FRED_MULTIPOSE_RESULTS {
    PROCESS_RESULTS_WORKFLOW('FRED_50_poses', Channel.fromPath("${params.dockedLigandRMSDs}/FRED_50_poses/*.csv"))
}

// Example workflow entries
workflow calculate_rmsd {
    CALCULATE_FRED_RMSD()
    CALCULATE_ALL_SINGLE_POSE_RMSD()
}

workflow process_results {
    PROCESS_FRED_RESULTS()
    PROCESS_ALL_SINGLE_POSE_RESULTS()
    PROCESS_ALL_MULTIPOSE_RESULTS()
}