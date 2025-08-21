process CREATE_EVALUATOR_FACTORY_SETTINGS{
    publishDir "${params.evaluator_configs}"
    conda "${params.harbor}"
    tag "create-evaluator-factory-settings"
    cache false
    memory { 4.GB }
    time { 10.m }
    label 'cpushort'

    output:
    path("*.yaml"), emit: evaluator_configs

    script:
    """
    python3 "${params.scripts}"/create_evaluator_factory_settings.py
    """
}

process CREATE_EVALUATORS {
    publishDir "${params.evaluationResults}", mode: 'copy', overwrite: false
    conda "${params.harbor}"
    tag "create-evaluators ${name}"
    memory { 32.GB }
    time { 20.m }
    label 'cpushort'
    cache 'lenient'

    input:
    val(name)
    path(settings)
    path(docking_results_parquet)
    path(docking_results_json)

    output:
    path("${name}/*"), emit: evaluator_json_directory

    script:
    """
    python3 "${params.scripts}"/create_evaluators.py \
    --input-parquet "${docking_results_parquet}" \
    --settings "${settings}" \
    --output "${name}" \
    """
}
process CREATE_MULTIPOSE_EVALUATORS {
    publishDir "${params.evaluationResults}", mode: 'copy', overwrite: false
    conda "${params.harbor}"
    tag "create-evaluators ${name}"
    memory { 32.GB }
    time { 20.m }
    label 'cpushort'
    cache 'lenient'

    input:
    val(name)

    output:
    path("${name}/*"), emit: evaluator_json_directory

    script:
    """
    python3 "${params.scripts}"/create_multipose_evaluators.py \
    --output "${name}" \
    """
}
process CREATE_EVALUATORS_MODULAR {
    publishDir "${params.evaluationResults}", mode: 'copy', overwrite: false
    conda "${params.harbor}"
    tag "create-evaluators ${name}"
    label 'local'

    input:
    val(name)
    path(script_path)
    path(docking_results_parquet)
    path(docking_results_json)

    output:
    path("${name}/*"), emit: evaluator_json_directory

    script:
    """
    python3 "${script_path}" \
    --output "${name}" \
    --input-parquet "${docking_results_parquet}" \
    """
}
process RUN_EVALUATORS {
    conda "${params.harbor}"
    tag "run-evaluators ${name}"
    errorStrategy = { task.exitStatus in [137,140,143,247] ? 'retry' : 'terminate' }
    maxRetries 3
    // Dynamic memory allocation
    memory { task.attempt > 1 ? (2 ** (task.attempt - 1)) * 256.GB : 256.GB }
    // Dynamic time allocation
    time { task.attempt > 1 ? (2 ** (task.attempt - 1)) * 1.h : 1.h }
    // set n cpus to request
    cpus 32
    'lenient'

    input:
    val(name)
    path(docking_results_parquet)
    path(docking_results_json)
    path("evaluator_jsons_*")


    output:
    path("*.csv"), emit: evaluator_results

    script:
    """
    python3 "${params.scripts}"/run_evaluators.py \
    evaluator_jsons_* \
    --input-parquet "${docking_results_parquet}" \
    --n-cpus 32
    """
}
process RUN_EVALUATORS_LIGHTWEIGHT {
    conda "${params.harbor}"
    tag "run-evaluators ${name}"
    errorStrategy = { task.exitStatus in [137,140,143,247] ? 'retry' : 'terminate' }
    maxRetries 3
    // Dynamic memory allocation
    memory { task.attempt > 1 ? (2 ** (task.attempt - 1)) * 32.GB : 32.GB }
    // Dynamic time allocation
    time { task.attempt > 1 ? (2 ** (task.attempt - 1)) * 1.h : 1.h }
    // set n cpus to request
    cpus 8
    'lenient'

    input:
    val(name)
    path(docking_results_parquet)
    path(docking_results_json)
    path("evaluator_jsons_*")


    output:
    path("*.csv"), emit: evaluator_results

    script:
    """
    python3 "${params.scripts}"/run_evaluators.py \
    evaluator_jsons_* \
    --input-parquet "${docking_results_parquet}" \
    --n-cpus 8
    """
}
process COMBINE_EVALUATIONS {
    publishDir "${params.evaluationResults}", mode: 'copy', overwrite: true
    conda "${params.harbor}"
    tag "combine-evaluations ${name}"

    input:
    val(name)
    path("evaluator_results_*")

    output:
    path("${name}_combined_results.csv"), emit: combined_results

    script:
    """
    python "${params.scripts}/combine_evaluation_results.py" \
    evaluator_results_* \
    "${name}_combined_results.csv"
    """
}