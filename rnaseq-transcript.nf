#!/usr/bin/env nextflow

// Copyright (C) 2026 IARC/WHO
// This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
// See the GNU General Public License for more details <http://www.gnu.org/licenses/>.

nextflow.enable.dsl = 2

// --------------------------------------------------
// PARAMETERS
// --------------------------------------------------

params.input_folder = '.'
params.input_file = null
params.output_folder = "."
params.mem  = 2
params.cpu  = 2
params.gtf  = null
params.prepDE_input = null
params.readlength = 75
params.twopass  = null
params.annot_organism = "Homo sapiens"
params.annot_genome   = "Unknown"
params.annot_provider = "Unspecified"
params.annot_version  = "Unspecified"
params.ref = "Unspecified"
params.help = null

//Header for the IARC tools - logo generated using the following page : http://patorjk.com/software/taag  (ANSI logo generator)
def IARC_Header (){
     return  """
#################################################################################
# ██╗ █████╗ ██████╗  ██████╗██████╗ ██╗ ██████╗ ██╗███╗   ██╗███████╗ ██████╗  #
# ██║██╔══██╗██╔══██╗██╔════╝██╔══██╗██║██╔═══██╗██║████╗  ██║██╔════╝██╔═══██╗ #
# ██║███████║██████╔╝██║     ██████╔╝██║██║   ██║██║██╔██╗ ██║█████╗  ██║   ██║ #
# ██║██╔══██║██╔══██╗██║     ██╔══██╗██║██║   ██║██║██║╚██╗██║██╔══╝  ██║   ██║ #
# ██║██║  ██║██║  ██║╚██████╗██████╔╝██║╚██████╔╝██║██║ ╚████║██║     ╚██████╔╝ #
# ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═════╝ ╚═╝ ╚═════╝ ╚═╝╚═╝  ╚═══╝╚═╝      ╚═════╝  #
# Nextflow pipelines for cancer genomics.########################################
"""
}

// --------------------------------------------------
// FILE DEFINITION
// --------------------------------------------------
gtf = params.gtf ? file(params.gtf) : null
if (!params.gtf) error "GTF file is required"
prepDE_input = Channel.value(params.prepDE_input)

// --------------------------------------------------
// INPUT CHANNELS
// --------------------------------------------------

if (params.input_file) {

    bam_files = Channel
        .fromPath(params.input_file)
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.ID, row.readlength as Integer, file(row.bam)) }

} else {
	 bam_files = Channel
    	.fromPath("${params.input_folder}/*.bam")
    	.ifEmpty { error "No BAM files found" }
    	.map { bam -> tuple(bam.baseName, params.readlength, bam) }
    } 

bam_1pass = bam_files
bam_2pass = bam_files

// --------------------------------------------------
// PROCESSES
// --------------------------------------------------

process STRINGTIE_1STPASS {
    cpus params.cpu
    memory "${params.mem}G"
    tag { sample_id }

    input:
    tuple val(sample_id), val(readlength), path(bam)
    path gtf

    output:
    tuple path("${sample_id}"), val(readlength), emit: st1
    path "*.log", emit: logs

    publishDir params.output_folder, mode: 'copy',
        saveAs: { f ->
    def fname = f instanceof Path ? f.name : f.toString()
    fname.endsWith('.log')
        ? "logs/${fname}"
        : params.twopass
            ? "intermediate_files/sample_folders/ST1of2passes/${fname}"
            : "intermediate_files/sample_folders/ST1pass/${fname}"
}

    script:

	def opts
	def log

	if (params.twopass) {
    	opts = "-o ${sample_id}/${sample_id}_1of2passes_ST.gtf -e -B -A ${sample_id}/${sample_id}_pass1_gene_abund.tab"
    	log  = "${sample_id}_1of2passes.log"
	} else {
    	opts = "-o ${sample_id}/${sample_id}_1pass_ST.gtf -e -B -A ${sample_id}/${sample_id}_pass1_gene_abund.tab"
    	log  = "${sample_id}_1pass.log"
	}

    """
    stringtie ${opts} -p ${params.cpu} -G ${gtf} -l ${sample_id} ${bam}
    cp .command.log ${log}
    """
}

process MERGE_GTF {

    cpus params.cpu
    memory "${params.mem}G"

    input:
    path st_dirs
    path gtf

    output:
    path "stringtie_merged.gtf", emit: merged_gtf
    path "gffcmp_merged*", emit: gffcmp

    publishDir "${params.output_folder}/gtf", mode: 'copy'

    script:
    """
    ls */*_ST.gtf > mergelist.txt
    stringtie --merge -p ${params.cpu} -G ${gtf} -o stringtie_merged.gtf mergelist.txt
    gffcompare -r ${gtf} -G -o gffcmp_merged stringtie_merged.gtf
    """
}

process STRINGTIE_2NDPASS {

    cpus params.cpu
    memory "${params.mem}G"
    tag { sample_id }

    input:
    tuple val(sample_id), val(readlength), path(bam)
    path merged_gtf
    path gtf

    output:
    tuple path(sample_id), val(readlength), emit: st2
    path "*.log", emit: logs

    publishDir params.output_folder, mode: 'copy',
        saveAs: { f ->
            f.name.endsWith('.log')
                ? "logs/${f.name}"
                : "intermediate_files/sample_folders/ST2pass/${f.name}"
        }

    script:
    """
    stringtie -o ${sample_id}/${sample_id}_2pass_ST.gtf -e -B -A ${sample_id}/${sample_id}_pass2_gene_abund.tab \
        -p ${params.cpu} -G ${merged_gtf} -l ${sample_id} ${bam}

cp .command.log ${sample_id}_2pass.log
    """
}

process PREPDE {

    cpus params.cpu
    memory "${params.mem}G"
    tag { readlength }

    input:
    tuple path(st_dirs), val(readlength)
    val prep_input

    output:
    path "*count_matrix*.csv", emit: count_matrices

    publishDir "${params.output_folder}/intermediate_files/expr_matrices", mode: 'copy'

    script:

	def input_dir = (prep_input && prep_input != "null") ? prep_input : "."
    def suffix    = params.twopass ? "_2pass" : ""

    """
	prepDE.py -i ${input_dir} -l ${readlength} \
        -g gene_count_matrix${suffix}_l${readlength}.csv \
        -t transcript_count_matrix${suffix}_l${readlength}.csv
    """
}

process BALLGOWN {

    cpus params.cpu
    memory "${params.mem}G"

    input:
    path st_dirs

    output:
    path "*_matrix*.csv", emit: norm_matrices
    path "*.rda", emit: rdata

   publishDir params.output_folder, mode: 'copy',
        saveAs: { f ->
            def fname = f.toString().tokenize('/').last()
            fname.endsWith('.csv') ?
                "expr_matrices/${fname}" :
                "Robjects/${fname}"
        }

    script:
	def suffix = params.twopass ? "_2pass" : ""
    """
    Rscript ${baseDir}/bin/create_matrices.R ${suffix}
    """
}

process SUMMARIZEDEXPERIMENT {

    cpus params.cpu
    memory "${params.mem}G"

    input:
    path st_dirs
    path norm_matrices
    path count_matrices
    path gtf
    path merged_gtf

    output:
    path "*.rda"
    path "*_matrix*.csv"

    publishDir params.output_folder, mode: 'copy',
        saveAs: { f ->
            def fname = f.toString().tokenize('/').last()
            fname.endsWith('.csv')
                ? "expr_matrices/${fname}"
                : "Robjects/${fname}"
        }

    script:
	def suffix = params.twopass ? "_2pass" : "_1pass"
    def gtf2use
    def gtf_path

 	if (merged_gtf) {
        gtf2use  = gtf
        gtf_path = params.gtf
    } else {
        gtf2use  = merged_gtf
        gtf_path = "${params.output_folder}/gtf/${merged_gtf.getName()}"
    }

    """
    Rscript ${baseDir}/bin/create_summarizedExperiment.R \
        ${gtf2use} \
        "${params.annot_provider}" \
        "${params.annot_version}" \
        "${params.annot_genome}" \
        "${params.annot_organism}" \
        ${gtf_path} \
        ${suffix} \
        ${params.ref}
    """
}

// --------------------------------------------------
// WORKFLOW
// --------------------------------------------------

workflow {

  		log.info IARC_Header()
// --------------------------------------------------
// INFO / HELP
// --------------------------------------------------

log.info ""
log.info "----------------------------------------------------------------------------------------------------------------------"
log.info "RNAseq-transcript-nf 2.2: gene- and transcript-level expression quantification from RNA sequencing data with StringTie"
log.info "----------------------------------------------------------------------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it under certain conditions; see LICENSE for details."
log.info "----------------------------------------------------------------------------------------------------------------------"
log.info ""

// --------------------------------------------------
// INFO / HELP
// --------------------------------------------------

if (params.help) {
    log.info "--------------------------------------------------------"
    log.info "  USAGE                                                 "
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "nextflow run iarcbioinfo/rnaseq-transcript-nf [-profile <docker/singularity/conda>] [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments:"
    log.info '    --input_folder   FOLDER                  Folder containing RNA-seq BAM files whose expression to quantify.'
    log.info '    --gtf            FILE                    Annotation file.'
    log.info ""
    log.info "Optional arguments:"
	log.info '    --input_file     FILE                    File in TSV format containing columns "ID" (sample ID), "bam" '
	log.info '											   (path to RNA-seq BAM file), and "readlength" (sample read length)'
    log.info '    --output_folder  STRING                  Output folder (default: .).'
	log.info '    --readlength     STRING                  Mean read length for count computation (default: 75).'
	log.info '    --prepDE_input   FILE					   File given to script prepDE from StringTie (default: none).'
	log.info '    --annot_organism, '
	log.info '    --annot_genome, '
	log.info '    --annot_provider,'
	log.info '    --annot_version,'
	log.info '    --ref 		   STRINGS				   metainformation stored in SummarizedExperiment R object'
	log.info '					    					   (default: "Homo sapiens","hg38", Unknown, Unknown, Unknown)'
    log.info '    --cpu            INTEGER                 Number of cpu used by bwa mem and sambamba (default: 2).'
    log.info '    --mem            INTEGER                 Size of memory used for mapping (in GB) (default: 2).' 
    log.info ""
    log.info "Flags:"
    log.info "--twopass                                    Enable StringTie 2pass mode"
    log.info ""
    exit 0
} else {
/* Software information */
   log.info "input_folder	= ${params.input_folder}"
   log.info "input_file		= ${params.input_file}"
   log.info "cpu			= ${params.cpu}"
   log.info "mem			= ${params.mem}"
   log.info "readlength		= ${params.readlength}"
   log.info "output_folder	= ${params.output_folder}"
   log.info "gtf			= ${params.gtf}"
   log.info "twopass		= ${params.twopass}"
   log.info "params.prepDE_input = ${params.prepDE_input}"
   log.info "annot_organism	= ${params.annot_organism}"
   log.info "annot_genome	= ${params.annot_genome}"
   log.info "annot_provider	= ${params.annot_provider}"
   log.info "annot_version	= ${params.annot_version}"
   log.info "ref			= ${params.ref}"
   log.info "help:            ${params.help}"
}

// RUN PROCESSES

    STRINGTIE_1STPASS(bam_1pass, gtf)

    if (params.twopass) {

        MERGE_GTF(
    		STRINGTIE_1STPASS.out.st1
        	.map { it[0] }
        	.collect(),
			gtf)

        STRINGTIE_2NDPASS(
            bam_2pass,
            MERGE_GTF.out.merged_gtf,
            gtf)

        st_final_ch = STRINGTIE_2NDPASS.out.st2
		merged_gtf4se_ch = MERGE_GTF.out.merged_gtf

    } else {
        st_final_ch = STRINGTIE_1STPASS.out.st1
		merged_gtf4se_ch = Channel.value(null)
    }

    grouped = st_final_ch
    .groupTuple(by: 1)
    .map { readlength, dirs -> tuple(dirs, readlength) }

    PREPDE(grouped, Channel.value(params.prepDE_input))

   BALLGOWN(
    st_final_ch
        .map { it[0] }
        .collect()
)
	
	SUMMARIZEDEXPERIMENT(
        st_final_ch.map { it[0] }.collect(),
        BALLGOWN.out.norm_matrices,
        PREPDE.out.count_matrices.collect(),
        gtf,
        merged_gtf4se_ch
    )
}
