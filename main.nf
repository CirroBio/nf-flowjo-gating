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

process filter_beads {
    container "${params.container}"
    publishDir "${params.output_directory}", mode: 'copy', overwrite: true

    input:
    path input_fcs

    output:
    path "*"

    script:
    template "filter_beads.R"
}

process make_anndata {
    container "${params.container}"
    publishDir "${params.output_directory}", mode: 'copy', overwrite: true

    input:
    path "*"

    output:
    path "output.h5ad"

    script:
    template "make_anndata.py"
}

process make_plots {
    container "${params.container}"
    publishDir "${params.output_directory}", mode: 'copy', overwrite: true

    input:
    path "input.?.h5ad"

    output:
    path "*"

    script:
    template "make_plots.py"
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

    if(params.input_wsp) {
        input_wsp = file(params.input_wsp, checkIfExists: true)

        apply_gates(
            input_fcs,
            input_wsp
        )
        output = apply_gates.out
    }else{
        filter_beads(input_fcs)
        output = filter_beads.out
    }
    
    // Make an AnnData from the gated data
    make_anndata(output)

    // If there are any additional summary files provided,
    // include them in the input for the next process
    if (params.input_anndata) {
        make_plots(
            make_anndata.out.merge(
                Channel
                    .fromPath(
                        "${params.input_anndata}".split(',').toList(),
                        checkIfExists: true
                    )
                    .toSortedList()
            )
        )
    } else {
        make_plots(make_anndata.out)
    }
}