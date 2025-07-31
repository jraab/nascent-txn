/*
========================================================================================
    Main workflow functions
========================================================================================
*/

class WorkflowMain {

    //
    // Citation string for pipeline
    //
    public static String citation(workflow) {
        return "If you use ${workflow.manifest.name} for your analysis please cite:\\n\\n" +
            "* The pipeline\\n" +
            "  https://doi.org/10.5281/zenodo.XXXXXX\\n\\n" +
            "* The Nextflow framework\\n" +
            "  https://doi.org/10.1038/nbt.3820\\n\\n" +
            "* Software dependencies\\n" +
            "  https://github.com/${workflow.manifest.name}/blob/master/CITATIONS.md"
    }

    //
    // Print help to screen if required
    //
    public static String helpMessage(workflow, params, log) {
        def command = "nextflow run ${workflow.manifest.name} --sample_sheet samplesheet.csv --genome hg38 -profile docker"
        def help_string = ''
        help_string += "\\n" + workflow.manifest.name + " v${workflow.manifest.version}\\n"
        help_string += "Usage: ${command}\\n\\n"
        help_string += "Mandatory arguments:\\n"
        help_string += "  --sample_sheet [file]        Path to sample sheet CSV file\\n"
        help_string += "  --genome [str]                Reference genome (hg38, mm10)\\n\\n"
        help_string += "Optional arguments:\\n"
        help_string += "  --outdir [path]               Output directory (default: results)\\n"
        help_string += "  --spike_in_genome [str]       Spike-in genome for normalization\\n"
        help_string += "  --help                        Show this help message\\n\\n"
        help_string += citation(workflow) + '\\n'
        return help_string
    }

    //
    // Print parameter summary log to screen
    //
    public static String paramsSummaryLog(workflow, params, log) {
        def summary_log = ''
        summary_log += "\\n" + workflow.manifest.name + " v${workflow.manifest.version}\\n"
        summary_log += "===================================\\n"
        summary_log += "Sample sheet    : ${params.sample_sheet}\\n"
        summary_log += "Genome          : ${params.genome}\\n"
        summary_log += "Output dir      : ${params.outdir}\\n"
        if (params.spike_in_genome) {
            summary_log += "Spike-in genome : ${params.spike_in_genome}\\n"
        }
        summary_log += "===================================\\n"
        summary_log += citation(workflow) + '\\n'
        return summary_log
    }

    //
    // Validate parameters and print summary to screen
    //
    public static void initialise(workflow, params, log) {
        // Print help to screen if required
        if (params.help) {
            log.info helpMessage(workflow, params, log)
            System.exit(0)
        }

        // Print parameter summary log to screen
        log.info paramsSummaryLog(workflow, params, log)
    }

    //
    // Check resource limits
    //
    public static def check_max(obj, type, params) {
        if (type == 'memory') {
            try {
                if (obj.compareTo(params.max_memory as nextflow.util.MemoryUnit) == 1)
                    return params.max_memory as nextflow.util.MemoryUnit
                else
                    return obj
            } catch (all) {
                println "   ### ERROR ###   Max memory '${params.max_memory}' is not valid! Using default value: $obj"
                return obj
            }
        } else if (type == 'time') {
            try {
                if (obj.compareTo(params.max_time as nextflow.util.Duration) == 1)
                    return params.max_time as nextflow.util.Duration
                else
                    return obj
            } catch (all) {
                println "   ### ERROR ###   Max time '${params.max_time}' is not valid! Using default value: $obj"
                return obj
            }
        } else if (type == 'cpus') {
            try {
                return Math.min( obj, params.max_cpus as int )
            } catch (all) {
                println "   ### ERROR ###   Max cpus '${params.max_cpus}' is not valid! Using default value: $obj"
                return obj
            }
        }
    }
}