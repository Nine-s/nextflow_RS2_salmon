process MERGE_QUANT {
    label 'python'
    publishDir params.outdir

    input:
    file out_bam
    
    output:
    path("merged_quant.sf"), emit: gathered_quant
    
    script:
    """
    python3 ${params.baseDir}/bin/merge_quant.py 
    """
}
