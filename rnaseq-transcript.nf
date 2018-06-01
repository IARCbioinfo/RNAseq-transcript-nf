#! /usr/bin/env nextflow

// Copyright (C) 2017 IARC/WHO

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


params.input_folder = '.'
params.output_folder= "."
params.mem  = 2
params.cpu  = 2
params.gtf  = null

params.twopass  = null

params.help = null

log.info ""
log.info "--------------------------------------------------------"
log.info "  <PROGRAM_NAME> <VERSION>: <SHORT DESCRIPTION>         "
log.info "--------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "--------------------------------------------------------"
log.info ""

if (params.help) {
    log.info "--------------------------------------------------------"
    log.info "  USAGE                                                 "
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "nextflow run iarcbioinfo/rnaseq-transcript-nf [-with-docker] [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments:"
    log.info '    --input_folder   FOLDER                  Folder containing BAM or fastq files to be aligned.'
    log.info '    --gtf          FILE                    Annotation file.'
    log.info ""
    log.info "Optional arguments:"
    log.info '    --output_folder     STRING                Output folder (default: results_alignment).'
    log.info '    --cpu          INTEGER                 Number of cpu used by bwa mem and sambamba (default: 2).'
    log.info '    --mem          INTEGER                 Size of memory used for mapping (in GB) (default: 2).' 
    log.info ""
    log.info "Flags:"
    log.info "--<FLAG>                                                    <DESCRIPTION>"
    log.info ""
    exit 0
} else {
/* Software information */
   log.info "input_folder = ${params.input_folder}"
   log.info "cpu          = ${params.cpu}"
   log.info "mem          = ${params.mem}"
   log.info "output_folder= ${params.output_folder}"
   log.info "help:                               ${params.help}"
}

if (file(params.input_folder).listFiles().findAll { it.name ==~ /.*bam/ }.size() > 0){
       println "BAM files found, proceed with transcript quantification"; mode ='bam'
       bam_files = Channel.fromPath( params.input_folder+'/*.bam'
       //bai_files = Channel.fromPath( params.input_folder+'/*.bai'
       )
}else{
       println "ERROR: input folder contains no fastq nor BAM files"; System.exit(1)
}

gtf = file(params.gtf)
bam_files.into { bam_files_41stpass; bam_files_42ndpass }

// 1st pass identifies new transcripts for each BAM file
process StringTie1stpass {
	cpus params.cpu
	memory params.mem+'G'
	tag { file_tag }

	input:
	file bam from bam_files_41stpass
	file gtf
	
	output:
	file("${file_tag}_ST.gtf") into STgtfs
	if( params.twopass == null ){
      	    publishDir "${params.output_folder}/${file_tag}", mode: 'copy'
	}
	shell:
	if(params.twopass==null){
	  STopts=''
	}else{
	  STopts="-e -B -A !{file_tag}_pass1_gene_abund.tab "
	}
	file_tag=bam.baseName
    	'''
    	stringtie !{STopts} -o !{file_tag}_ST.gtf -p !{params.cpu} -G !{gtf} -l !{file_tag} !{file_tag}.bam
    	'''
}

if(params.twopass){
// Merges the list of transcripts of each BAM file
process mergeGTF {
	cpus params.cpu
	memory params.mem+'G'
	tag { "merge" }

	input:
	file gtfs from STgtfs.collect()
	file gtf

	output: 
	file("stringtie_merged.gtf") into merged_gtf
	file("gffcmp_merged*") into gffcmp_output
	publishDir "${params.output_folder}", mode: 'copy', pattern: '{gffcmp_merged*}' 

	shell:
	'''
	ls *_ST.gtf > mergelist.txt
	stringtie --merge -p !{params.cpu} -G !{gtf} -o stringtie_merged.gtf mergelist.txt
	gffcompare -r !{gtf} -G -o gffcmp_merged stringtie_merged.gtf
	'''
}

// Quantifies transcripts identified in 1st pass in each sample
process StringTie2ndpass {
	cpus params.cpu
	memory params.mem+'G'
	tag { file_tag }

	input:
	file bam from bam_files_42ndpass
	file merged_gtf
	file gtf

	output: 
	file("${file_tag}_gene_abund.tab") into tabs
	file("*.ctab") into ctabs
	file("${file_tag}_merged.gtf") into gtf_merged
	publishDir "${params.output_folder}/${file_tag}", mode: 'copy'

	shell:
	file_tag=bam.baseName
	'''
	stringtie -e -B -p !{params.cpu} -G !{merged_gtf} -o !{file_tag}_merged.gtf -A !{file_tag}_gene_abund.tab !{file_tag}.bam
	'''
}
}