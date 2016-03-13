TUMOR_BAM=`echo $1 | sed "s:\.bam::g"`

export SAMTOOLS='/hematology/tools/samtools/samtools'
export PINDEL='/hematology/tools/pindel/pindel'

$SAMTOOLS view -b ${TUMOR_BAM}.bam chr13:28606093-28610689 1> ${TUMOR_BAM}_pindel.bam
$SAMTOOLS index ${TUMOR_BAM}_pindel.bam

#-----------------------------------------------------------------
# Run Pindel
#-----------------------------------------------------------------

echo "${TUMOR_BAM}_pindel.bam	600	${TUMOR_BAM}" > ${TUMOR_BAM}_conf.txt
$PINDEL -i ${TUMOR_BAM}_conf.txt -f $REFERENCE -o ${TUMOR_BAM}_pindel -c 'chr13:28565921-28702250' -r false --report_interchromosomal_events false
rm ${TUMOR_BAM}_conf.txt

$PINDEL2VCF -p ${TUMOR_BAM}_pindel_TD -r $REFERENCE -R RefGene -d 20140901 -v ${TUMOR_BAM}_pindel_TD.vcf

$PINDEL2VCF -p ${TUMOR_BAM}_pindel_SI -r $REFERENCE -R RefGene -d 20140901 -v ${TUMOR_BAM}_pindel_SI.vcf

exit 0