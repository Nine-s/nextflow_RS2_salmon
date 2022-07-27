process MERGE_QUANT {
    label 'python'
    publishDir params.outdir

    input:
    file out_bam
    
    output:
    path("*"), emit: gathered_quant
    
    script:
    """
    echo ${out_bam}
    ls > files.txt
    python ${params.baseDir}/bin/merge_quant.py 
    """
}
