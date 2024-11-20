process ALIGN {
	tag "ALIGN on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/${sample.run}/mapped/", mode:'copy'
	
	input:
	tuple val(name), val(sample), path(reads)

	output:
    tuple val(name), val(sample), path("${name}.rmdup.sorted.bam"), path("${name}.rmdup.sorted.bai")

	script:
	rg = "\"@RG\\tID:${name}\\tSM:${name}\\tLB:${name}\\tPL:ILLUMINA\""
	"""
	echo ALIGN $name
	source activate bwa
	bwa mem -R ${rg} -t $task.cpus ${params.refindex} $reads | samtools view -Sb - | sambamba sort /dev/stdin -o ${name}.rmdup.sorted.bam
        samtools index ${name}.rmdup.sorted.bam ${name}.rmdup.sorted.bai
	"""
}

process VARDICT {
       tag "VARDICT on $name" 
       publishDir "${params.outDirectory}/${sample.run}/varianty/", mode:'copy'

        input:
        tuple val(name), val(sample), path(bam), path(bai)

        output:
        tuple val(name), val(sample), path("${name}.vcf")

        script:
        """
        echo VARDICT $name
        source activate vardict
        vardict -G ${params.refindex} -f 0 -b $bam -N ${name} \
                -c 1 -S 2 -E 3 -g 4 \
                ${params.varbed} | Rscript --vanilla ${params.teststrandbias} | perl ${params.var2vcf_valid} -f 0 -N ${name} -A > ${name}.vcf
        """
}
process NORMALIZACE {
        tag "NORMALIZACE on $name" 
        publishDir "${params.outDirectory}/${sample.run}/varianty/", mode:'copy'
        container "staphb/bcftools:1.10.2"

        input: 
        tuple val(name), val(sample), path(vardict)

        output:
        tuple val(name), val(sample), path("${name}.bcf2.vcf")

        script:
        """
        echo NORMALIZACE $name
        bgzip $vardict 
        tabix ${name}.vcf.gz
        bcftools norm -m-both -o ${name}.bcf1.vcf.gz ${name}.vcf.gz
        bcftools norm -f ${params.refindex} -o ${name}.bcf2.vcf ${name}.bcf1.vcf.gz
        """
}

 process ANOTACE {
       tag "ANOTACE on $name"
       publishDir "${params.outDirectory}/${sample.run}/varianty/", mode:'copy'
       container "staphb/bcftools:1.10.2"
 
        input:
        tuple val(name), val(sample), path(normalizovany)
  
        output:
        tuple val(name), val(sample), path("${name}.bcf2.vcf.hg19_multianno.vcf.gz"), path("${name}.final.txt")

        script:
        """
        echo ANOTACE $name

        ${params.annovar} -vcfinput $normalizovany /mnt/hdd2/reference/humandb/ -buildver hg19 -protocol refGeneWithVer,ensGene,1000g2015aug_all,1000g2015aug_eur,avsnp150,gnomad_exome,cadd13,clinvar_20220320,popfreq_max_20150413,dbnsfp33a -operation gx,g,f,f,f,f,f,f,f,f -nastring . -otherinfo -polish -xreffile ${params.gene_fullxref.txt} -arg '-splicing 5 -exonicsplicing',,,,,,,,, --remove
        bgzip ${name}.bcf2.vcf.hg19_multianno.vcf
        tabix ${name}.bcf2.vcf.hg19_multianno.vcf.gz
  
        ${params.gatk} --java-options "-Xmx4g" VariantsToTable -R ${params.refindex}  --show-filtered  -V ${name}.bcf2.vcf.hg19_multianno.vcf.gz -F CHROM -F POS -F REF -F ALT -GF GT -GF DP -GF VD -GF ALD -GF AF -F Func.refGeneWithVer -F Gene.refGeneWithVer -F GeneDetail.refGeneWithVer -F ExonicFunc.refGeneWithVer -F AAChange.refGeneWithVer -F 1000g2015aug_all -F 1000g2015aug_eur -F avsnp150 -F CADD13_RawScore -F CADD13_PHRED -F CLNSIG -F PopFreqMax -F gnomAD_exome_ALL -F gnomAD_exome_NFE -F SIFT_score -F SIFT_converted_rankscore -F SIFT_pred -F Polyphen2_HDIV_score -F Polyphen2_HDIV_rankscore -F Polyphen2_HDIV_pred -F Polyphen2_HVAR_score -F Polyphen2_HVAR_rankscore -F Polyphen2_HVAR_pred -F Expression-egenetics.refGene --output ${name}.final.txt
        """
}

process COVERAGE {
          tag "COVERAGE on $name"
       publishDir "${params.outDirectory}/${sample.run}/mapped/", mode:'copy'
         container "staphb/samtools:1.20"

        input:
        tuple val(name), val(sample), path(bam), path(bai)

        output:
        tuple val(name), val(sample), path("${name}.coveragefin.txt")

        script:
        """
        echo ANOTACE $name
        source /home/kristina/miniconda3/bin/activate kamila
        samtools bedcov ${params.varbed} $bam -d 20 > ${name}.COV
        awk '{print \$5/(\$3-\$2)}'  ${name}.COV >  ${name}.COV-mean
        awk '{print (\$6/(\$3-\$2))*100"%"}' ${name}.COV > ${name}-procento-nad-20
        paste ${name}.COV-mean ${name}-procento-nad-20 > vysledek
        echo "chr" "start" "stop" "name" ${name}.COV-mean ${name}-procento-nad-20 > hlavicka
        sed -i 's/ /\t/'g hlavicka
        paste ${params.varbed} vysledek > coverage
        cat hlavicka coverage > ${name}.coveragefin.txt
        sed -i -e "s/\r//g" ${name}.coveragefin.txt
        """
}
