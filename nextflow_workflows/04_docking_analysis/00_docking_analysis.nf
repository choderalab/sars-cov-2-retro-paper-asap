#!/usr/bin/env nextflow
include {
    CREATE_EVALUATOR_FACTORY_SETTINGS
    CREATE_EVALUATORS
    RUN_EVALUATORS
    COMBINE_EVALUATIONS
} from "./modules.nf"

import groovy.yaml.YamlSlurper

params.K = 10
params.analysis_config = "docking_analysis_config.yaml"

// Load configuration using YamlSlurper
def config = new YamlSlurper().parse(file(params.analysis_config))

// Helper function to generate workflow names
def getWorkflowName(analysis, dataset, variant = null) {
    def parts = [analysis, dataset]
    if (variant) parts.add(variant)
    return parts.join('_')
}

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

        eval_inputs_ch = CREATE_EVALUATORS.output.evaluator_json_directory
            .flatMap { dir -> file("${dir}/*.json") }
            .buffer(size: params.K)

        RUN_EVALUATORS(
            name,
            docking_results_parquet,
            docking_results_json,
            eval_inputs_ch,
        )

        all_results = RUN_EVALUATORS.output.evaluator_results
            .flatten()
            .collect()

        COMBINE_EVALUATIONS(
            name,
            all_results
        )
}

// Create workflow definitions first
def workflow_definitions = []
config.analyses.each { analysis_name, analysis_config ->
    analysis_config.each { variant_name, variant_config ->
        variant_config.enabled_datasets.each { dataset_id ->
            def dataset = config.datasets[dataset_id]

            if (!dataset) {
                println "\nWarning: Dataset '${dataset_id}' not found in config, skipping workflow"
                return  // skip this iteration
            }

            def dataset_name = dataset.name
            def workflow_name = getWorkflowName(analysis_name, dataset_name, variant_name)
            def dataset_parquet = "${params.combinedDockingResultsPath}/${dataset_name}.parquet"
            def dataset_json = "${params.combinedDockingResultsPath}/${dataset_name}.json"
            def settings_file = "${params.evaluator_configs}/${variant_config.settings}"

            workflow_definitions << [
                name: workflow_name,
                parquet: dataset_parquet,
                json: dataset_json,
                settings: settings_file
            ]
        }
    }
}

// Print workflow definitions
println "\nWorkflow Definitions:"
println "===================="
workflow_definitions.each { def workflow_def ->
    println """
    Name: ${workflow_def.name}
    Parquet: ${workflow_def.parquet}
    JSON: ${workflow_def.json}
    Settings: ${workflow_def.settings}
    -------------------"""
}
println "====================\n"

// Convert workflow definitions to a channel
def workflow_ch = Channel.fromList(workflow_definitions)

workflow RUN_ANALYSIS {
    take:
        setup_analysis_check

    main:
        RUN_DOCKING_ANALYSIS(
            workflow_ch.map { it.name },
            workflow_ch.map { it.parquet },
            workflow_ch.map { it.json },
            workflow_ch.map { it.settings }
        )
}

workflow CREATE_EVALUATOR_FACTORY_SETTINGS_WORKFLOW {
    CREATE_EVALUATOR_FACTORY_SETTINGS()
}

workflow {
    RUN_ANALYSIS(CREATE_EVALUATOR_FACTORY_SETTINGS_WORKFLOW().out.collect())
}