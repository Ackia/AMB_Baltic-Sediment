#!/usr/bin/env nextflow
/*
 * pipeline input parameters
 */
params.reads = ""
params.outdir = ""
params.cpus = "2"
params.mem = "8GB"


// requires --reads for Assembly
if (params.reads == '') {
    exit 1, '--reads is a required paramater for Meta-Assembly pipeline'
}

// requires --outdir for Assembly
if (params.outdir == '') {
    exit 1, '--outdir is a required paramater for Meta-Assembly pipeline'
}

println """\
         Hybrid Assembly- N F   P I P E L I N E
         ===================================
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()

         reads_atropos_pe = Channel
                      .fromFilePairs(params.reads + '*_{1,2}.fastq.gz', size: 2, flat: true)
process trimming_pe {
                          publishDir params.outdir, mode: 'copy'

                          input:
                              set val(id), file(read1), file(read2) from reads_atropos_pe

                          output:
                              set val(id), file("${id}_R1.fastq"), file("${id}_R2.fastq") into trimmed_reads_pe

                          script:
                              """
                              mkdir trimmed
                              atropos -a TGGAATTCTCGGGTGCCAAGG -B AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT \
                                  -T 4 -m 50 --max-n 0 -q 20,20 -pe1 $read1 -pe2 $read2 \
                                  -o ${id}_R1.fastq -p ${id}_R2.fastq
                              """
                      }
process fastqc {
                          publishDir params.outdir, mode: 'copy'

                          input:
                              file reads from trimmed_reads_pe.collect()

                          output:
                              file "*_fastqc.{zip,html}" into fastqc_results

                          script:
                              """
                              fastqc $reads
                              """
                      }

process multiqc {
                          publishDir params.outdir, mode: 'copy'

                          input:
                              file 'fastqc/*' from fastqc_results.collect()

                          output:
                              file 'multiqc_report.html'

                          script:
                              """
                              multiqc .
                              """
}
