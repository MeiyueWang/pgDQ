for i in CS_singlecell_rep1 CS_singlecell_rep2;do
{
	samtools view -@ 20 -h ${i}.bam | grep -E "^\@|NH:i:1" | awk 'BEGIN{FS="\t"} $5==255' |samtools view -@ 20 -t /data/wangmeiyue/genome/iwgsc_v1.0/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta.fai -o ${i}_uniq.bam -
	samtools view -@ 20 -h ${i}.bam | grep -E -v "^\@|NH:i:1" | awk 'BEGIN{FS="\t"} $5!=255' |samtools view -@ 20 -t /data/wangmeiyue/genome/iwgsc_v1.0/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta.fai -o ${i}_multi.bam -
	samtools stats ${i}_uniq.bam >${i}_uniq.bam.stats 2>&1
	samtools stats ${i}_multi.bam >${i}_multi.bam.stats 2>&1
	}&
done
wait
