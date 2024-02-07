process apply_gates {
    container "${params.container}"
    publishDir "${params.output_directory}", mode: 'copy', overwrite: true

    input:
    path input_fcs
    path input_wsp

    output:
    path "*"

    script:
    template "apply_gates.R"
}

workflow {
    Channel
        .fromPath(
            "${params.input_fcs}".split(',').toList(),
            checkIfExists: true
        )
        .ifEmpty { error "No files found matching '${params.input_fcs}'" }
        .toSortedList()
        .set { input_fcs }

    input_wsp = file(params.input_wsp, checkIfExists: true)

    apply_gates(
        input_fcs,
        input_wsp
    )
}