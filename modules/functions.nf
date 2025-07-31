/*
========================================================================================
    Utility functions for TT-seq pipeline
========================================================================================
*/

import groovy.json.JsonBuilder

//
// Function to initialise default values and to generate a Groovy Map of available tools and their versions
//
public static Map toolVersions() {
    def tool_versions = [:]
    tool_versions['nextflow'] = "$workflow.nextflow.version"
    return tool_versions
}

//
// Function to generate file path for publishing
//
public static String getPathFromPublishDir(Map args) {
    def publish_dir = args.publish_dir
    def filename = args.filename
    if (publish_dir instanceof Map) {
        def path = publish_dir.path instanceof Closure ? publish_dir.path.call(args) : publish_dir.path
        def mode = publish_dir.mode ?: 'copy'
        def pattern = publish_dir.pattern ?: '*'
        return "${path}/${filename}"
    } else {
        return "${publish_dir}/${filename}"
    }
}

//
// Function to save/publish module results
//
public static String getProcessName(String name) {
    return name.split(':')[-1]
}

//
// Function to validate sample sheet
//
def validateSampleSheet(csvFile) {
    def summary_stats = [:]
    csvFile.eachLine { line ->
        def cols = line.split(',')
        if (cols.size() < 3) {
            exit 1, "ERROR: Sample sheet must have at least 3 columns: sample,fastq_1,fastq_2"  
        }
        def sample = cols[0]
        def fastq1 = cols[1]
        def fastq2 = cols.size() > 2 ? cols[2] : ''
        
        if (!fastq1 || !file(fastq1).exists()) {
            exit 1, "ERROR: FastQ file 1 does not exist: ${fastq1}"
        }
        if (fastq2 && !file(fastq2).exists()) {
            exit 1, "ERROR: FastQ file 2 does not exist: ${fastq2}"  
        }
        
        summary_stats[sample] = [
            'fastq_1': fastq1,
            'fastq_2': fastq2,
            'single_end': fastq2 ? false : true
        ]
    }
    return summary_stats
}