# SAMPLE_NAME="8435A3"
# FASTQ1="/data/rawdata/BGI/iGT_drSNP_250324/20250321-T7S00934/T7S00934-RDMT02-QN68242041JX-25226012-A175V1-25R08434A1_R1.fastq.gz"
# FASTQ2="/data/rawdata/BGI/iGT_drSNP_250324/20250321-T7S00934/T7S00934-RDMT02-QN68242041JX-25226012-A175V1-25R08434A1_R2.fastq.gz"
# OUTDIR="/data/mengxf/Project/KML250324_drSNP_iGT/work/250324/result"
SAMPLE_NAME=8434B3
FASTQ1=/data/rawdata/BGI/iGT_drSNP_250324/20250321-T7S00934/T7S00934-RDMT02-QN68242041JX-25226012-A175V1-25R08434B3_R1.fastq.gz
FASTQ2=/data/rawdata/BGI/iGT_drSNP_250324/20250321-T7S00934/T7S00934-RDMT02-QN68242041JX-25226012-A175V1-25R08434B3_R2.fastq.gz
OUTDIR=/data/mengxf/Project/KML250324_drSNP_iGT/work/250324/result

THREADS=16
REFERENCE="/data/mengxf/Database/reference/hg19/hg19.fa"
BED="/data/mengxf/Project/KML250324_drSNP_iGT/work/250324/target/A175loci-hg19.bed"
RSID_ANNO="/data/mengxf/Project/KML250324_drSNP_iGT/work/250324/target/A175loci-hg19.freq-anno_with_ref_alt.tsv"

# * 进入Conda环境
source /home/mengxf/miniforge3/bin/activate basic

# # 质控
# mkdir -p $OUTDIR/qc/fastqc/$SAMPLE_NAME
# fastqc -t $THREADS -f fastq -o $OUTDIR/qc/fastqc/$SAMPLE_NAME $FASTQ1 $FASTQ2

# mkdir -p $OUTDIR/qc/trimmed/$SAMPLE_NAME/
# fastp -q 15 -u 40 -l 25 --thread $THREADS \
#     --cut_right --cut_window_size 20 --cut_mean_quality 30 --correction \
#     -i $FASTQ1 \
#     -I $FASTQ2 \
#     -o $OUTDIR/qc/trimmed/$SAMPLE_NAME/$SAMPLE_NAME.clean.R1.fq \
#     -O $OUTDIR/qc/trimmed/$SAMPLE_NAME/$SAMPLE_NAME.clean.R2.fq \
#     -j $OUTDIR/qc/trimmed/$SAMPLE_NAME/$SAMPLE_NAME.json \
#     -h $OUTDIR/qc/trimmed/$SAMPLE_NAME/$SAMPLE_NAME.html

# # 比对
# mkdir -p $OUTDIR/align/$SAMPLE_NAME
# # bowtie2
# # ! 提高 properly paired, MAPQ
# bowtie2 -p $THREADS \
#     --local --no-unal --fr --no-discordant --no-mixed --score-min L,0,0.6 --ma 2 --mp 4 \
#     -x ${REFERENCE} \
#     -1 $OUTDIR/qc/trimmed/$SAMPLE_NAME/$SAMPLE_NAME.clean.R1.fq \
#     -2 $OUTDIR/qc/trimmed/$SAMPLE_NAME/$SAMPLE_NAME.clean.R2.fq |
#     samtools view -@ $THREADS -hbS - |
#     samtools sort -@ $THREADS -o $OUTDIR/align/$SAMPLE_NAME/$SAMPLE_NAME.aligned.bam -
# # GATK4 Best Practices Workflows: Data pre-processing for variant discovery
# # AddOrReplaceReadGroups
# gatk AddOrReplaceReadGroups \
#     I=$OUTDIR/align/$SAMPLE_NAME/$SAMPLE_NAME.aligned.bam \
#     O=$OUTDIR/align/$SAMPLE_NAME/$SAMPLE_NAME.aligned.addrg.bam \
#     RGID=$SAMPLE_NAME RGLB=$SAMPLE_NAME RGPL=ILLUMINA RGPU=unit1 RGSM=$SAMPLE_NAME
# # SortAndFixTags
# gatk SortSam \
#     --INPUT $OUTDIR/align/$SAMPLE_NAME/$SAMPLE_NAME.aligned.addrg.bam \
#     --OUTPUT /dev/stdout \
#     --SORT_ORDER coordinate --CREATE_INDEX false --CREATE_MD5_FILE false |
#     gatk SetNmMdAndUqTags \
#         --INPUT /dev/stdin \
#         --OUTPUT $OUTDIR/align/$SAMPLE_NAME/$SAMPLE_NAME.aligned.sorted.bam \
#         --CREATE_INDEX true --CREATE_MD5_FILE true \
#         --REFERENCE_SEQUENCE $REFERENCE
# # BaseRecalibrator
# gatk BaseRecalibrator \
#     -R $REFERENCE -L $BED \
#     -I $OUTDIR/align/$SAMPLE_NAME/$SAMPLE_NAME.aligned.sorted.bam \
#     --use-original-qualities \
#     -O $OUTDIR/align/$SAMPLE_NAME/$SAMPLE_NAME.recal_data.csv \
#     --known-sites /data/mengxf/Database/GATK4/hg19/1000G_phase1.snps.high_confidence.hg19.sites.bgzip.vcf.gz \
#     --known-sites /data/mengxf/Database/GATK4/hg19/dbsnp_138.hg19.bgzip.vcf.gz \
#     --known-sites /data/mengxf/Database/GATK4/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.bgzip.vcf.gz
# # ApplyBQSR
# gatk ApplyBQSR \
#     -R $REFERENCE -L $BED \
#     -I $OUTDIR/align/$SAMPLE_NAME/$SAMPLE_NAME.aligned.sorted.bam \
#     -O $OUTDIR/align/$SAMPLE_NAME/$SAMPLE_NAME.aligned.recalibrated.bam \
#     -bqsr $OUTDIR/align/$SAMPLE_NAME/$SAMPLE_NAME.recal_data.csv \
#     --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
#     --add-output-sam-program-record \
#     --create-output-bam-md5 \
#     --use-original-qualities
# # 靶区域统计
# samtools depth -b ${BED} \
#     $OUTDIR/align/$SAMPLE_NAME/$SAMPLE_NAME.aligned.recalibrated.bam \
#     >$OUTDIR/align/$SAMPLE_NAME/$SAMPLE_NAME.aligned.recalibrated.bam.depth

# # 调用变异
# mkdir -p $OUTDIR/vcf/$SAMPLE_NAME
# # freebayes 敏感性高，深度与IGV一致. 不报告indel, mnp, complex位点
# freebayes --throw-away-indel-obs --throw-away-mnps-obs --throw-away-complex-obs \
#     --use-duplicate-reads --report-monomorphic \
#     --min-alternate-fraction 0.0001 --min-alternate-count 3 \
#     -t $BED -f ${REFERENCE} \
#     $OUTDIR/align/$SAMPLE_NAME/$SAMPLE_NAME.aligned.recalibrated.bam \
#     >$OUTDIR/vcf/$SAMPLE_NAME/$SAMPLE_NAME.freebayes.vcf

# SNP注释
# 1. chr:start - rs_id 对应表, 注释上rs号
# 2. vcf
# 3. depth 深度统计文件, 阴性时补充深度
python /data/mengxf/Project/KML250324_drSNP_iGT/work/250324/jupyter/scripts/freebayes_snp_anno.py \
    ${RSID_ANNO} \
    $OUTDIR/vcf/$SAMPLE_NAME/$SAMPLE_NAME.freebayes.vcf \
    $OUTDIR/vcf/$SAMPLE_NAME/$SAMPLE_NAME.snp_anno.tsv
