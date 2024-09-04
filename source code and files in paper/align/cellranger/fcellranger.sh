dir="cellranger-7.1.0"
for i in CS_singlecell_rep1 CS_singlecell_rep2;do
{
	rm -rf ${i}
	${dir}/cellranger count --id ${i} --fastqs=rawdata --sample=${i} --transcriptome=genome/iwgsc_v1.0/cellranger/wheat6_with_chrC_and_chrM_cellranger --nosecondary --localcores=40 --expect-cells=10000
	}&
done
wait
