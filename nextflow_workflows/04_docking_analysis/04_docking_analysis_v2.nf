#!/usr/bin/env nextflow
include {
    CREATE_EVALUATOR_FACTORY_SETTINGS
    CREATE_EVALUATORS
    RUN_EVALUATORS
    RUN_EVALUATORS_LIGHTWEIGHT
    COMBINE_EVALUATIONS
    CREATE_MULTIPOSE_EVALUATORS
    CREATE_EVALUATORS_MODULAR
} from "./modules.nf"
params.K = 10

workflow RUN_DOCKING_ANALYSIS {
    take:
        name
        docking_results_parquet
        docking_results_json
        evaluator_settings

    main:
        CREATE_EVALUATORS(
            name,
            evaluator_settings,
            docking_results_parquet,
            docking_results_json
        )

        // Create channel from JSON files only after evaluator creation
        eval_inputs_ch = CREATE_EVALUATORS.output.evaluator_json_directory
            .flatMap { dir -> file("${dir}/*.json") }
            .buffer(size: params.K)

        RUN_EVALUATORS(
            name,
            docking_results_parquet,
            docking_results_json,
            eval_inputs_ch,
        )

        // Collect all evaluator results before combining
        all_results = RUN_EVALUATORS.output.evaluator_results
            .flatten()
            .collect()

        COMBINE_EVALUATIONS(
            name,
            all_results
        )
}
dataset_names = [
"posit_single_pose": "ALL_1_poses",
"fred_single_pose": "FRED_1_poses",
"posit_multipose": "ALL_50_poses",
"fred_multipose": "FRED_50_poses",
]

def results = [:]
dataset_names.each { label, name ->
    def parquet = "${params.combinedDockingResultsPath}/${name}.parquet"
    def json = "${params.combinedDockingResultsPath}/${name}.json"

    results[label] = [  // Store directly in map with label as key
        name: name,
        docking_results_parquet: parquet,
        docking_results_json: json
    ]
}

// Print dataset definitions
println "\nDataset Definitions:"
println "===================="
results.each { label, data ->
    println """
    Name: ${label}
    Label: ${data.name}
    Parquet: ${data.docking_results_parquet}
    JSON: ${data.docking_results_json}
    -------------------"""
}

settings_map = [
"datesplit": "reference_split_comparison.yaml",
"x_to_x": "x_to_x_scaffold_split.yaml",
"x_to_x_5": "x_to_x_scaffold_split_5_refs.yaml",
"x_to_y": "x_to_y_scaffold_split.yaml",
"x_to_y_5": "x_to_y_scaffold_split_5_refs.yaml",
"x_to_not_x": "x_to_not_x_scaffold_split.yaml",
"not_x_to_x": "not_x_to_x_scaffold_split.yaml",
"not_x_to_x_5": "not_x_to_x_scaffold_split_5_refs.yaml",
"ecfp4": "increasing_similarity_ecfp4.yaml",
"mcs": "increasing_similarity_mcs.yaml",
"tc": "increasing_similarity_tanimoto_combo_aligned.yaml",]

def settings = [:]
settings_map.each { label, filename ->
    settings[label] = [label: label, filename: "${params.evaluator_configs}/${filename}"]
}

workflow CREATE_EVALUATOR_FACTORY_SETTINGS_WORKFLOW {
    CREATE_EVALUATOR_FACTORY_SETTINGS()
}

workflow RUN_ANALYSIS {
    take:
        result
        setting
    main:
        RUN_DOCKING_ANALYSIS(
            "${result.name}_${setting.label}",
            result.docking_results_parquet,
            result.docking_results_json,
            setting.filename,
        )
}

workflow DATESPLIT_POSIT {RUN_ANALYSIS(results.posit_single_pose, settings.datesplit)}
workflow DATESPLIT_FRED {RUN_ANALYSIS(results.fred_single_pose, settings.datesplit)}
workflow X_TO_X_POSIT {RUN_ANALYSIS(results.posit_single_pose, settings.x_to_x)}
workflow X_TO_X_POSIT_5_REFS {RUN_ANALYSIS(results.posit_single_pose, settings.x_to_x_5)}
workflow X_TO_Y_POSIT {RUN_ANALYSIS(results.posit_single_pose, settings.x_to_y)}
workflow X_TO_Y_POSIT_5_REFS {RUN_ANALYSIS(results.posit_single_pose, settings.x_to_y_5)}
workflow NOT_X_TO_X_POSIT {RUN_ANALYSIS(results.posit_single_pose, settings.not_x_to_x)}
workflow NOT_X_TO_X_POSIT_5_REFS {RUN_ANALYSIS(results.posit_single_pose, settings.not_x_to_x_5)}
workflow X_TO_NOT_X_POSIT {RUN_ANALYSIS(results.posit_single_pose, settings.x_to_not_x)}
workflow INCREASING_SIMILARITY_TC_ALIGNED_POSIT {
    RUN_ANALYSIS(results.posit_single_pose, settings.tc)
}
workflow INCREASING_SIMILARITY_TC_ALIGNED_FRED {
    RUN_ANALYSIS(results.fred_single_pose, settings.tc)
}
workflow INCREASING_SIMILARITY_MCS_POSIT {
    RUN_ANALYSIS(results.posit_single_pose, settings.mcs)
}
workflow INCREASING_SIMILARITY_MCS_FRED {
    RUN_ANALYSIS(results.fred_single_pose, settings.mcs)
}
workflow INCREASING_SIMILARITY_ECFP4_POSIT {
    RUN_ANALYSIS(results.posit_single_pose, settings.ecfp4)
}
workflow INCREASING_SIMILARITY_ECFP4_FRED {
    RUN_ANALYSIS(results.fred_single_pose, settings.ecfp4)
}
workflow INCREASING_SIMILARITY_ECFP4_POSIT_MULTIPOSE {
    RUN_ANALYSIS(results.posit_single_pose, settings.ecfp4)
}
workflow posit_scaffold_splits {
    X_TO_X_POSIT()
    X_TO_X_POSIT_5_REFS()
    X_TO_Y_POSIT()
    X_TO_Y_POSIT_5_REFS()
    NOT_X_TO_X_POSIT()
    NOT_X_TO_X_POSIT_5_REFS()
    X_TO_NOT_X_POSIT()
}
workflow analyze_posit {
    DATESPLIT_POSIT()
    X_TO_X_POSIT()
    X_TO_X_POSIT_5_REFS()
    X_TO_Y_POSIT()
    X_TO_Y_POSIT_5_REFS()
    NOT_X_TO_X_POSIT()
    NOT_X_TO_X_POSIT_5_REFS()
    X_TO_NOT_X_POSIT()
    INCREASING_SIMILARITY_TC_ALIGNED_POSIT()
    INCREASING_SIMILARITY_MCS_POSIT()
    INCREASING_SIMILARITY_ECFP4_POSIT()
    INCREASING_SIMILARITY_ECFP4_POSIT_MULTIPOSE()
}
workflow analyze_fred {
    DATESPLIT_FRED()
    INCREASING_SIMILARITY_TC_ALIGNED_FRED()
    INCREASING_SIMILARITY_MCS_FRED()
    INCREASING_SIMILARITY_ECFP4_FRED()
}
workflow {
    analyze_posit()
    analyze_fred()
}

workflow POSIT_MULTIPOSE_ANALYSIS {
    name = "posit_multipose_analysis"
    CREATE_MULTIPOSE_EVALUATORS(name)

    // Create channel from JSON files only after evaluator creation
    eval_inputs_ch = CREATE_MULTIPOSE_EVALUATORS.output.evaluator_json_directory
        .flatMap { dir -> file("${dir}/*.json") }
        .buffer(size: 1)

    RUN_EVALUATORS(
        name,
        results.posit_multipose.docking_results_parquet,
        results.posit_multipose.docking_results_json,
        eval_inputs_ch,
    )

    // Collect all evaluator results before combining
    all_results = RUN_EVALUATORS.output.evaluator_results
        .flatten()
        .collect()

    COMBINE_EVALUATIONS(
        name,
        all_results
    )
}
workflow FRED_MULTIPOSE_ANALYSIS {
    name = "fred_multipose_analysis"
    CREATE_MULTIPOSE_EVALUATORS(name)

    // Create channel from JSON files only after evaluator creation
    eval_inputs_ch = CREATE_MULTIPOSE_EVALUATORS.output.evaluator_json_directory
        .flatMap { dir -> file("${dir}/*.json") }
        .buffer(size: 1)

    RUN_EVALUATORS(
        name,
        results.fred_multipose.docking_results_parquet,
        results.fred_multipose.docking_results_json,
        eval_inputs_ch,
    )

    // Collect all evaluator results before combining
    all_results = RUN_EVALUATORS.output.evaluator_results
        .flatten()
        .collect()

    COMBINE_EVALUATIONS(
        name,
        all_results
    )
}
workflow SCAFFOLD_DATE_SPLIT {
    name = "posit_scaffold_date_split"
    script_path = "${params.scripts}/create_evaluators_scaffold_datesplit.py"

    CREATE_EVALUATORS_MODULAR(
        name,
        script_path,
        results.posit_single_pose.docking_results_parquet,
        results.posit_single_pose.docking_results_json,
    )

    // Create channel from JSON files only after evaluator creation
    eval_inputs_ch = CREATE_EVALUATORS_MODULAR.output.evaluator_json_directory
        .flatMap { dir -> file("${dir}/*.json") }
        .buffer(size: 1)

    RUN_EVALUATORS_LIGHTWEIGHT(
        name,
        results.posit_single_pose.docking_results_parquet,
        results.posit_single_pose.docking_results_json,
        eval_inputs_ch,
    )

    // Collect all evaluator results before combining
    all_results = RUN_EVALUATORS_LIGHTWEIGHT.output.evaluator_results
        .flatten()
        .collect()

    COMBINE_EVALUATIONS(
        name,
        all_results
    )
}
workflow POSIT_REVERSE_SIMILARITY_SPLIT {
    name = "posit_reverse_similarity_split"
    script_path = "${params.scripts}/create_reverse_similarity_split_evaluators.py"

    CREATE_EVALUATORS_MODULAR(
        name,
        script_path,
        results.posit_single_pose.docking_results_parquet,
        results.posit_single_pose.docking_results_json,
    )

    // Create channel from JSON files only after evaluator creation
    eval_inputs_ch = CREATE_EVALUATORS_MODULAR.output.evaluator_json_directory
        .flatMap { dir -> file("${dir}/*.json") }
        .buffer(size: 1)

    RUN_EVALUATORS_LIGHTWEIGHT(
        name,
        results.posit_single_pose.docking_results_parquet,
        results.posit_single_pose.docking_results_json,
        eval_inputs_ch,
    )

    // Collect all evaluator results before combining
    all_results = RUN_EVALUATORS_LIGHTWEIGHT.output.evaluator_results
        .flatten()
        .collect()

    COMBINE_EVALUATIONS(
        name,
        all_results
    )
}