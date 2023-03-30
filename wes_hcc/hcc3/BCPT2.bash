ref=/storage/zhangyanxiaoLab/zhangliwen/genome/galtk/Homo_sapiens_assembly38.fasta
snp=/storage/zhangyanxiaoLab/zhangliwen/genome/galtk/dbsnp_146.hg38.vcf
indel=/storage/zhangyanxiaoLab/zhangliwen/genome/galtk/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
	
for sample in BCPT2

do
	echo $sample


	if [ ! -f ${sample}_marked_fixed.bam ]
	then
		gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" FixMateInformation \
			-I ${sample}_marked.bam \
			-O ${sample}_marked_fixed.bam \
			-SO coordinate \
			1>${sample}_log.fix 2>&1
		samtools index ${sample}_marked_fixed.bam
	fi

	if [ ! -f ${sample}_recal.table ]
	then
		echo ${sample} recal
		gatk --java-options "-Xmx20G -Djava.io.tmpdir=./"  BaseRecalibrator \
			-R $ref  \
			-I ${sample}_marked_fixed.bam  \
			--known-sites $snp \
			--known-sites $indel \
			-O ${sample}_recal.table \
			1>${sample}_log.recal 2>&1
	fi



	if [ ! -f ${sample}_bqsr.bam ]
	then
		echo ${sample} bqsr
		gatk --java-options "-Xmx20G -Djava.io.tmpdir=./"   ApplyBQSR \
			-R $ref  \
			-I ${sample}_marked_fixed.bam  \
			-bqsr ${sample}_recal.table \
			-O ${sample}_bqsr.bam \
			1>${sample}_log.ApplyBQSR  2>&1
	fi


done

