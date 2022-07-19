process GENERATE_DECOY_TRANSCIPTROME {
    label 'bash'
    publishDir params.outdir
    
    input:
    path(reference)
    path(transcriptome)

    output:
    path("decoy.txt"), emit:decoy
    path("gentrome.fa"), emit:gentrome
    
    //https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
    //grep "^>" <(gunzip -c GRCm38.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > decoys.txt
    //    cat gencode.vM23.transcripts.fa.gz GRCm38.primary_assembly.genome.fa.gz > gentrome.fa.gz
    script:
    """
    grep "^>" <${reference} | cut -d " " -f 1 > decoy.txt
    sed -i -e 's/>//g' decoy.txt
    cat ${transcriptome} ${reference} > gentrome.fa
    """
}

process SALMON_INDEX_REFERENCE {
    label 'salmon'
    publishDir params.outdir
    
    input:
    path(decoy)
    path(gentrome)

    output:
    path("salmon_index/")

    //https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
    script:
    """
    salmon index -t ${gentrome} -d ${decoy} -p 12 -i salmon_index
    """
    //NOTE: --gencode flag is for removing extra metdata in the target header separated by | from the gencode reference. You can skip it if using other references.
}


process SALMON_ALIGN_QUANT {
    label 'salmon'
    publishDir params.outdir
 
    input:
    env STRANDNESS
    tuple val(sample_name), path(reads_1), path(reads_2)
    path(index)
    path(annotation)

    output:
    path("*-transcripts_quant"), emit: quantification
    path("*-transcripts_quant/*-quant.sf"), emit: quant_file

    shell:
    '''
    filename="!{reads_1}"
    prefix=${filename%-*} 
    newFolderName="${prefix}-transcripts_quant"
    newFileName="${prefix}-quant.sf"
    
    if [[ ($STRANDNESS == "firststrand") ]]; then 
	salmon quant --index !{index} -l ISR -1 !{reads_1} -2 !{reads_2} --validateMappings -o transcripts_quant
        mv transcripts_quant/quant.sf transcripts_quant/${newFileName}
        mv transcripts_quant ${newFolderName}
       elif [[ ($STRANDNESS == "secondstrand") ]]; then 
        salmon quant --index !{index} -l ISF -1 !{reads_1} -2 !{reads_2} --validateMappings -o transcripts_quant
        mv transcripts_quant/quant.sf transcripts_quant/${newFileName}
        mv transcripts_quant ${newFolderName}
	
       elif [[ $STRANDNESS == "unstranded" ]]; then
		salmon quant --index !{index} -l IU -1 !{reads_1} -2 !{reads_2} --validateMappings -o transcripts_quant
                mv transcripts_quant/quant.sf transcripts_quant/${newFileName}
        mv transcripts_quant ${newFolderName}

	else  
		echo $STRANDNESS > error_strandness.txt
		echo "strandness cannot be determined" >> error_strandness.txt
	fi

   '''
 
}
