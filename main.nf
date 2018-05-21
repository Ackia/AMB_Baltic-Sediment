#!/usr/bin/env nextflow
/*
 * pipeline input parameters
 */
params.reads = ""
params.outdir = ""
params.cpus = "2"
params.mem = "12"


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
                              set val(id), file("${id}_R1_trimmed.fastq"), file("${id}_R2_trimmed.fastq") into trimmed_reads_pe

                          script:
                              """
                              atropos -a TGGAATTCTCGGGTGCCAAGG -B AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT \
                                  -T 4 -m 50 --max-n 0 -q 20,20 -pe1 $read1 -pe2 $read2 \
                                  -o ${id}_R1_trimmed.fastq -p ${id}_R2_trimmed.fastq
                              """
                      }

trimmed_reads_pe.into {reads_for_fastq; reads_for_megahit; reads_for_spades; reads_for_metabin_1; reads_for_metabin_2}

process fastqc {
                          publishDir params.outdir, mode: 'copy'

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
process megahit {
                            publishDir params.outdir, mode: 'copy'

                            input:
                                set val(id), file(read1), file(read2) from reads_for_megahit

                            output:
                                set val(id), file"${id}.contigs.fa" into megahit_result

                            script:
                                """
                                megahit  -t $params.cpus -o ${id}_megahit --out-prefix ${id} -1 $read1 -2 $read2
                                """
}
process metaspades {
                            publishDir params.outdir, mode: 'copy'

                            input:
                                set val(id), file(read1), file(read2) from reads_for_spades

                            output:
                                set val(id), file'assembly.fasta' into spades_result

                            script:
                                """
                                spades.py -o ${id}_spades --meta -1 $read1 -2 $read2 -t $params.cpus
                                """

}
process metabat {
                            publishDir params.outdir, mode: 'copy'

                            input:
                                set val(id), file(megahitassembly), from megahit_result


                            output:
                                file'assembly.fasta' into spades_result
                                file(read1), file(read2) from reads_for_metabin

                            script:
                                """
                                bowtie2-build megahitassembly final.contigs
                                bowtie2 -x final.contigs -1 tara_reads_R1.fastq.gz -2 tara_reads_R2.fastq.gz | \
                                samtools view -bS -o tara_to_sort.bam
                                samtools sort tara_to_sort.bam -o tara.bam
                                samtools index tara.bam
                                runMetaBat.sh -m 1500 final.contigs.fa tara.bam
                                """
}
