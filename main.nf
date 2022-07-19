
nextflow.enable.dsl = 2

include { FASTP } from './modules/fastp.nf'
include { CHECK_STRANDNESS } from './modules/check_strandness.nf'
include { GENERATE_DECOY_TRANSCIPTROME ; SALMON_INDEX_REFERENCE ; SALMON_ALIGN_QUANT } from './modules/salmon.nf'
include { MERGE_QUANT } from './modules/merge_quant.nf'

log.info """\
         RNAseq analysis using NextFlow 
         =============================
         genome: ${params.reference_genome}
         annot : ${params.reference_annotation}
         reads : ${params.reads}
         outdir: ${params.outdir}
         """
         .stripIndent()
 
params.outdir = 'results'

process split_fastq {
    
    input: 
    tuple val(name), path(fastq)

    output:
    tuple val(name), path("*${name}*R1*.f*q"), path("*${name}*R2*.f*q")

    script:
    """ 
    ${params.baseDir}/bin/splitFastq -i ${fastq[0]} -n ${params.split} -o ${fastq[0].getBaseName()} 
    ${params.baseDir}/bin/splitFastq -i ${fastq[1]} -n ${params.split} -o ${fastq[1].getBaseName()} 

    """
}

process split_fastq_unzipped {

    input:
    tuple val(name), path(fastq)

    output:
    tuple val(name), path("${name}_1-${/[0-9]/*params.suffix_length}.fq"), path("${name}_2-${/[0-9]/*params.suffix_length}.fq")

    script:
    """
    cat ${fastq[0]} | split \\
        -a ${params.suffix_length} \\
        -d \\
        -l ${params.num_lines} \\
        - \\
        ${fastq[0].getBaseName()}- \\
        --additional-suffix=".fq"

    cat ${fastq[1]} | split \\
        -a ${params.suffix_length} \\
        -d \\
        -l ${params.num_lines} \\
        - \\
        ${fastq[1].getBaseName()}- \\
        --additional-suffix=".fq" 
    """
}

workflow {

    Channel.fromFilePairs(params.reads, checkIfExists: true).set{ read_pairs_unsplit_ch }

    CHECK_STRANDNESS( read_pairs_unsplit_ch, params.reference_cdna, params.reference_annotation_ensembl )
    FASTP( read_pairs_unsplit_ch )
    
    split_fastq(FASTP.out.sample_trimmed) \
	    | map { name, fastq, fastq1 -> tuple( groupKey(name, fastq.size()), fastq, fastq1 ) } \
        | transpose() \
        | view()
        | set{ read_pairs_ch }
    
    GENERATE_DECOY_TRANSCIPTROME( params.reference_genome, params.reference_cdna )
    SALMON_INDEX_REFERENCE( GENERATE_DECOY_TRANSCIPTROME.out.decoy, GENERATE_DECOY_TRANSCIPTROME.out.gentrome )
    SALMON_ALIGN_QUANT( CHECK_STRANDNESS.out.first(), read_pairs_ch, SALMON_INDEX_REFERENCE.out, params.reference_annotation )
 MERGE_QUANT( SALMON_ALIGN_QUANT.out.quant_file.collect() )
}

