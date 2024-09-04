dir="/home/wangmeiyue/tools/cellranger-7.1.0"
for i in CS_singlecell_rep1 CS_singlecell_rep2;do
{
	rm -rf ${i}
	${dir}/cellranger count --id ${i} --fastqs=/data/wangmeiyue/wheat_singlecell_pgDQ/01CS_scRNA-seq/02cellranger/02_cellranger/rawdata --sample=${i} --transcriptome=/data/wangmeiyue/genome/iwgsc_v1.0/cellranger/wheat6_with_chrC_and_chrM_cellranger --nosecondary --localcores=40 --expect-cells=10000
	}&
done
wait
