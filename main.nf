#!/usr/bin/env nextflow
/*
 * pipeline input parameters
 */
params.reads = ""
params.outdir = ""
params.cpus = "12"
params.mem = "2"


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
         outdir       : params.outdir
         """
         .stripIndent()

         reads_atropos_pe = Channel
                      .fromFilePairs(params.reads + '*_{1,2}.fastq.gz', size: 2, flat: true)
process trimming_pe {
                          publishDir params.outdir/trimmed, mode: 'copy'

                          input:
                              set val(id), file(read1), file(read2) from reads_atropos_pe

                          output:
                              set val(id), file("${id}_R1_trimmed.fastq"), file("${id}_R2_trimmed.fastq") into trimmed_reads_pe

                          script:
                              """
                              atropos -a TGGAATTCTCGGGTGCCAAGG -B AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT \
                                  -T 4 -m 50 --max-n 0 -q 20,20 -pe1 $read1 -pe2 $read2 \
                                  -o ${id}_R1_trimmed.fastq -p ${id}_R2_trimmed.fastq
                              """
                      }

trimmed_reads_pe.into {reads_for_fastq; reads_for_megahit; reads_for_spades; reads_for_metabin_1; reads_for_metabat; reads_for_checkm}

process fastqc {
                          publishDir params.outdir/fastq, mode: 'copy'

                          input:
                              file reads from reads_for_fastq.collect()

                          output:
                              file "*_fastqc.{zip,html}" into fastqc_results

                          script:
                              """
                              fastqc $reads
                              """
                      }

process multiqc {
                          publishDir params.outdir/multiqc, mode: 'copy'

                          input:
                              file 'fastqc/*' from fastqc_results.collect()

                          output:
                              file 'multiqc_report.html'

                          script:
                              """
                              multiqc .
                              """
}
process megahit {
                            publishDir params.outdir/megahit, mode: 'copy'

                            input:
                                set val(id), file(read1), file(read2) from reads_for_megahit

                            output:
                                file("${id}_megahit/${id}.contigs.fa") into megahit_result

                            script:
                                """
                                megahit  -t $params.cpus -o ${id}_megahit --out-prefix ${id} -1 $read1 -2 $read2
                                """
}
process metabat {
                            publishDir params.outdir/metabat, mode: 'copy'

                            input:
                                file'megahitassembly' from megahit_result
                                set val(id), file(read1), file(read2) from reads_for_metabat

                            output:
                                file"megahitassembly.metabat-bins1500" into metabat_results


                            script:
                                """
                                bowtie2-build $megahitassembly ${id}.contigs
                                bowtie2 -x ${id}.contigs -1 $read1 -2 $read2 | \
                                samtools view -bS -o megahit_to_sort.bam
                                samtools sort megahit_to_sort.bam -o ${id}_megahit.bam
                                samtools index ${id}_megahit.bam
                                runMetaBat.sh -m 1500 $megahitassembly ${id}_megahit.bam
                                """
}
process checkm {
                            publishDir params.outdir/checkm, mode: 'copy'

                            input:
                            file'metabatresult' from metabat_results
                            set val(id), file(read1), file(read2) from reads_for_checkm
                            output:
                            file"${id}_checkM" into checkm_results


                            script:
                            """
                            checkm lineage_wf -x fa -t $params.cpus $metabatresult ${id}_checkM
                            """
}
