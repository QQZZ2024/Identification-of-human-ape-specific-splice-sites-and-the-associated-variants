###Quality control of RNA-seq data

num=$(cat /home/id/samplesid.txt) #including five tissues-across humans, chimpanzees, rhesus macaques, and mice
for ID in $num
do
fastp -i /home/fastq/${ID}.fastq -o /home/fastq/${ID}_filter.fastq -h /home/fastq/${ID}_report.html -j /home/fastq/${ID}_report.json -q 30
done

#####################################################

###RNA-seq reads alignment
#human
hisat2-build /home/genome/hg38.fa /home/hisat/index/human/human_index
num=$(cat /home/id/humanid.txt)
for ID in $num
do
hisat2 -p 4 -x /home/hisat/index/human/human_index -U /home/fastq/${ID}_filter.fastq -S /home/hisat/${ID}_hisat2_human.sam
samtools sort -@ 2 -O bam -o /home/hisat/${ID}_hisat2_human.bam /home/hisat/${ID}_hisat2_human.sam
done

#chimpanzee
hisat2-build /home/genome/panTro6.fa /home/hisat/index/chimp/chimp_index
num=$(cat /home/id/chimpid.txt)
for ID in $num
do
hisat2 -p 4 -x /home/hisat/index/chimp/chimp_index -U /home/fastq/${ID}_filter.fastq -S /home/hisat/${ID}_hisat2_chimp.sam
samtools sort -@ 2 -O bam -o /home/hisat/${ID}_hisat2_chimp.bam /home/hisat/${ID}_hisat2_chimp.sam
done

#rhesus macaque
hisat2-build /home/genome/rheMac10.fa /home/hisat/index/rhesus/rhesus_index
num=$(cat /home/id/rhesusid.txt)
for ID in $num
do
hisat2 -p 4 -x /home/hisat/index/rhesus/rhesus_index -U /home/fastq/${ID}_filter.fastq -S /home/hisat/${ID}_hisat2_rhesus.sam
samtools sort -@ 2 -O bam -o /home/hisat/${ID}_hisat2_rhesus.bam /home/hisat/${ID}_hisat2_rhesus.sam
done

#mouse
hisat2-build /home/genome/mm39.fa /home/kyoku/data3/hisat/index/mouse/mouse_index
num=$(cat /home/kyoku/data3/run/id/mouseid.txt)
for ID in $num
do
hisat2 -p 4 -x /home/hisat/index/mouse/mouse_index -U /home/fastq/${ID}_filter.fastq -S /home/hisat/${ID}_hisat2_mouse.sam
samtools sort -@ 2 -O bam -o /home/hisat/${ID}_hisat2_mouse.bam /home/hisat/${ID}_hisat2_mouse.sam
done

#####################################################

###identify candidate splice site by extracting jucntion reads
#human
num=$(cat /home/id/humanid.txt)
for ID in $num
do
samtools view /home/hisat/${ID}_hisat2_human.bam | awk -F "\t" '{if($6~"^[0-9]*M[0-9]*N[0-9]*M$") print}' | awk -F "\t" '{split ($6,T,"[M-N]");$6=T[1] OFS T[2] OFS T[3]}1' OFS="\t" | awk -F "\t" 'BEGIN{OFS="\t"}{$1=$1"_"FNR;$6=$4+$6-1;$7=$6+$7+1;$8=$7+$8-1}{print $0}' | awk -F "\t" '{print$1"\t"$3"\t"$4"\t"$6"\t"$7"\t"$8}' > /home/junction/${ID}_human.bed
done

#chimpanzee
num=$(cat /home/kyoku/data3/run/id/chimpid.txt)

for ID in $num
do

samtools view /home/kyoku/data3/run/bam/${ID}_hisat2_chimp.bam | awk -F "\t" '{if($6~"^[0-9]*M[0-9]*N[0-9]*M$") print}' | awk -F "\t" '{split ($6,T,"[M-N]");$6=T[1] OFS T[2] OFS T[3]}1' OFS="\t" | awk -F "\t" 'BEGIN{OFS="\t"}{$1=$1"_"FNR;$6=$4+$6-1;$7=$6+$7+1;$8=$7+$8-1}{print $0}' | awk -F "\t" '{print$1"\t"$3"\t"$4"\t"$6"\t"$7"\t"$8}' > /home/kyoku/data3/run/mnm/${ID}_chimp.bed

done


