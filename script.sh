###Quality control of RNA-seq data

num=$(cat /home/id/samplesid.txt) #including five tissues-across humans, chimpanzees, rhesus macaques, and mice
for ID in $num
do
fastp -i /home/fastq/${ID}.fastq -o /home/fastq/${ID}_filter.fastq -h /home/fastq/${ID}_report.html -j /home/fastq/${ID}_report.json -q 30
done

#####################################################

###RNA-seq reads alignment
#build index
hisat2-build /home/genome/hg38.fa /home/hisat/index/human/human_index
hisat2-build /home/genome/panTro6.fa /home/hisat/index/chimp/chimp_index
hisat2-build /home/genome/rheMac10.fa /home/hisat/index/rhesus/rhesus_index
hisat2-build /home/genome/mm39.fa /home/kyoku/data3/hisat/index/mouse/mouse_index

#such as in human
num=$(cat /home/id/humanid.txt)
for ID in $num
do
hisat2 -p 4 -x /home/hisat/index/human/human_index -U /home/fastq/${ID}_filter.fastq -S /home/hisat/${ID}_hisat2_human.sam
samtools sort -@ 2 -O bam -o /home/hisat/${ID}_hisat2_human.bam /home/hisat/${ID}_hisat2_human.sam
done

#The scripts used for other species were similar

#####################################################

###identify candidate splice site by extracting jucntion reads
#such as in human
num=$(cat /home/id/humanid.txt)
for ID in $num
do
samtools view /home/hisat/${ID}_hisat2_human.bam | awk -F "\t" '{if($6~"^[0-9]*M[0-9]*N[0-9]*M$") print}' | awk -F "\t" '{split ($6,T,"[M-N]");$6=T[1] OFS T[2] OFS T[3]}1' OFS="\t" | awk -F "\t" 'BEGIN{OFS="\t"}{$1=$1"_"FNR;$6=$4+$6-1;$7=$6+$7+1;$8=$7+$8-1}{print $0}' | awk -F "\t" '{print$1"\t"$3"\t"$4"\t"$6"\t"$7"\t"$8}' > /home/junction/${ID}_human.bed
done

#The scripts used for other species were similar

#####################################################

###Cross-species genomic coordinate comparisons using LiftOver
#such as in human

#Convert splice site positions from one species to their orthologous positions in another
num=$(cat /home/id/humanid.txt)
for ID in $num
do
cat /home/junction/${ID}_human.bed | awk -F "\t" '{print $2"\t"$4-5"\t"$4+5}' | sort | uniq > /home/junction/${ID}_human_donor.bed
cat /home/junction/${ID}_human.bed | awk -F "\t" '{print $2"\t"$5-5"\t"$5+5}' | sort | uniq > /home/junction/${ID}_human_acceptor.bed
bedtools intersect -a /home/junction/${ID}_human_donor.bed -b /home/genome/hg38_trans.gtf -wa -wb | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$7}' | sort | uniq | awk -F "\t" 'BEGIN{OFS="\t"}{$4=$4"\t""id_"FNR}{print $0}' | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$5"\t""100""\t"$4}' > /home/junction/${ID}_human_donor_gene.bed
bedtools intersect -a /home/junction/${ID}_human_acceptor.bed -b /home/genome/hg38_trans.gtf -wa -wb | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$7}' | sort | uniq | awk -F "\t" 'BEGIN{OFS="\t"}{$4=$4"\t""id_"FNR}{print $0}' | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$5"\t""100""\t"$4}' > /home/junction/${ID}_human_acceptor_gene.bed
liftOver /home/junction/${ID}_human_donor_gene.bed /home/liftover/hg38ToPanTro6.over.chain.gz /home/liftover/${ID}_human_chimp_donor.bed unMapped
liftOver /home/junction/${ID}_human_donor_gene.bed /home/liftover/hg38ToRheMac10.over.chain.gz /home/liftover/${ID}_human_rhesus_donor.bed unMapped
liftOver /home/junction/${ID}_human_donor_gene.bed /home/liftover/hg38ToMm39.over.chain.gz /home/liftover/${ID}_human_mouse_donor.bed unMapped
liftOver /home/junction/${ID}_human_acceptor_gene.bed /home/liftover/hg38ToPanTro6.over.chain.gz /home/liftover/${ID}_human_chimp_acceptor.bed unMapped
liftOver /home/junction/${ID}_human_acceptor_gene.bed /home/liftover/hg38ToRheMac10.over.chain.gz /home/liftover/${ID}_human_rhesus_acceptor.bed unMapped
liftOver /home/junction/${ID}_human_acceptor_gene.bed /home/liftover/hg38ToMm39.over.chain.gz /home/liftover/${ID}_human_mouse_acceptor.bed unMapped
done

#Convert positions back to the original genome
num=$(cat /home/id/humanid.txt)
for ID in $num
do
liftOver /home/liftover/${ID}_human_chimp_donor.bed /home/liftover/panTro6ToHg38.over.chain.gz /home/liftover/back/${ID}_human_chimp_donor.bed unMapped
liftOver /home/liftover/${ID}_human_rhesus_donor.bed /home/liftover/rheMac10ToHg38.over.chain.gz /home/liftover/back/${ID}_human_rhesus_donor.bed unMapped
liftOver /home/liftover/${ID}_human_mouse_donor.bed /home/liftover/mm39ToHg38.over.chain.gz /home/liftover/back/${ID}_human_mouse_donor.bed unMapped
liftOver /home/liftover/${ID}_human_chimp_acceptor.bed /home/liftover/panTro6ToHg38.over.chain.gz /home/liftover/back/${ID}_human_chimp_acceptor.bed unMapped
liftOver /home/liftover/${ID}_human_rhesus_acceptor.bed /home/liftover/rheMac10ToHg38.over.chain.gz /home/liftover/back/${ID}_human_rhesus_acceptor.bed unMapped
liftOver /home/liftover/${ID}_human_mouse_acceptor.bed /home/liftover/mm39ToHg38.over.chain.gz /home/liftover/back/${ID}_human_mouse_acceptor.bed unMapped
done

#Retain only the splice sites with consistent positions in both forward and reverse mappings
num=$(cat /home/kyoku/data3/run/id/humanid.txt)
for ID in $num
do
awk -F "\t" 'NR==FNR{a[$4]=$0;next}{$7=a[$4];print}' /home/liftover/${ID}_human_chimp_donor.bed /home/junction/${ID}_human_donor_gene.bed |  tr ' ' '\t' | awk '$7!=""' | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$7"\t"$8"\t"$9"\t"$12}' > /home/liftover/retain/${ID}_tmp1.bed
awk -F "\t" 'NR==FNR{a[$4]=$0;next}{$10=a[$4];print}' /home/liftover/back/${ID}_human_chimp_donor.bed /home/liftover/retain/${ID}_tmp1.bed | tr ' ' '\t' | awk -F "\t" '{if($1==$10 && $2==$11 && $3==$12) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' > /home/liftover/retain/${ID}_tmp2.bed
awk -F "\t" 'NR==FNR{a[$4]=$0;next}{$10=a[$4];print}' /home/liftover/${ID}_human_rhesus_donor.bed /home/liftover/retain/${ID}_tmp2.bed | tr ' ' '\t' | awk '$10!=""' | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$15}' > /home/liftover/retain/${ID}_tmp3.bed
awk -F "\t" 'NR==FNR{a[$4]=$0;next}{$14=a[$4];print}' /home/liftover/back/${ID}_human_rhesus_donor.bed /home/liftover/retain/${ID}_tmp3.bed | tr ' ' '\t' | awk -F "\t" '{if($1==$14 && $2==$15 && $3==$16) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}' > /home/liftover/retain/${ID}_tmp4.bed
awk -F "\t" 'NR==FNR{a[$4]=$0;next}{$14=a[$4];print}' /home/liftover/${ID}_human_mouse_donor.bed /home/liftover/retain/${ID}_tmp4.bed | tr ' ' '\t' | awk '$14!=""' | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$19}' > /home/liftover/retain/${ID}_tmp5.bed
awk -F "\t" 'NR==FNR{a[$4]=$0;next}{$18=a[$4];print}' /home/liftover/back/${ID}_human_mouse_donor.bed /home/liftover/retain/${ID}_tmp5.bed |  tr ' ' '\t' | awk -F "\t" '{if($1==$18 && $2==$19 && $3==$20) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17}' > /home/liftover/retain/${ID}_human_donor.bed
awk -F "\t" 'NR==FNR{a[$4]=$0;next}{$7=a[$4];print}' /home/liftover/${ID}_human_chimp_acceptor.bed /home/junction/${ID}_human_acceptor_gene.bed |  tr ' ' '\t' | awk '$7!=""' | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$7"\t"$8"\t"$9"\t"$12}' > /home/liftover/retain/${ID}_tmp1.bed
awk -F "\t" 'NR==FNR{a[$4]=$0;next}{$10=a[$4];print}' /home/liftover/back/${ID}_human_chimp_acceptor.bed /home/liftover/retain/${ID}_tmp1.bed | tr ' ' '\t' | awk -F "\t" '{if($1==$10 && $2==$11 && $3==$12) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' > /home/liftover/retain/${ID}_tmp2.bed
awk -F "\t" 'NR==FNR{a[$4]=$0;next}{$10=a[$4];print}' /home/liftover/${ID}_human_rhesus_acceptor.bed /home/liftover/retain/${ID}_tmp2.bed | tr ' ' '\t' | awk '$10!=""' | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$15}' > /home/liftover/retain/${ID}_tmp3.bed
awk -F "\t" 'NR==FNR{a[$4]=$0;next}{$14=a[$4];print}' /home/liftover/back/${ID}_human_rhesus_acceptor.bed /home/liftover/retain/${ID}_tmp3.bed | tr ' ' '\t' | awk -F "\t" '{if($1==$14 && $2==$15 && $3==$16) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}' > /home/liftover/retain/${ID}_tmp4.bed
awk -F "\t" 'NR==FNR{a[$4]=$0;next}{$14=a[$4];print}' /home/liftover/${ID}_human_mouse_acceptor.bed /home/liftover/retain/${ID}_tmp4.bed | tr ' ' '\t' | awk '$14!=""' | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$19}' > /home/liftover/retain/${ID}_tmp5.bed 
awk -F "\t" 'NR==FNR{a[$4]=$0;next}{$18=a[$4];print}' /home/liftover/back/${ID}_human_mouse_acceptor.bed /home/liftover/retain/${ID}_tmp5.bed |  tr ' ' '\t' | awk -F "\t" '{if($1==$18 && $2==$19 && $3==$20) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17}' > /home/liftover/retain/${ID}_human_acceptor.bed
done

#The scripts used for other species were similar

#####################################################

###identify candidate species-specific splice sites
#such as in human
cat /home/liftover/retain/*_human_donor.bed| sort | uniq > /home/splice/human_donor.bed
cat /home/liftover/retain/*_human_acceptor.bed | sort | uniq > /home/splice/human_acceptor.bed
cat /home/junction/*_human.bed | awk -F "\t" '{print$2"\t"$4}' | sort | uniq > /home/junction/allhuman_donor.bed
cat /home/junction/*_human.bed | awk -F "\t" '{print$2"\t"$5}' | sort | uniq > /home/junction/allhuman_acceptor.bed

#The scripts used for other species were similar

library(dplyr)
library(tidyr)
human1 <- read.table(file = "/home/splice/human_donor.bed")
human2 <- read.table(file = "/home/splice/human_acceptor.bed")
human<-rbind(human1,human2)
human<-human[,c(1,2,4,5,7,8,10,11)]
index<-duplicated(human)
human<-human[!index,]
human$human<-paste(human[,1],human[,2],sep="_")
human$chimp<-paste(human[,3],human[,4],sep="_")
human$rhesus<-paste(human[,5],human[,6],sep="_")
human$mouse<-paste(human[,7],human[,8],sep="_")
all1_1<-read.table(file = "/home/junction/allhuman_donor.bed")
all1_2<-read.table(file = "/home/junction/allhuman_acceptor.bed")
all1<-rbind(all1_1,all1_2)
all2_1<-read.table(file = "/home/junction/allchimp_donor.bed")
all2_2<-read.table(file = "/home/junction/allchimp_acceptor.bed")
all2<-rbind(all2_1,all2_2)
all3_1<-read.table(file = "/home/junction/allrhesus_donor.bed")
all3_2<-read.table(file = "/home/junction/allrhesus_acceptor.bed")
all3<-rbind(all3_1,all3_2)
all4_1<-read.table(file = "/home/junction/allmouse_donor.bed")
all4_2<-read.table(file = "/home/junction/allmouse_acceptor.bed")
all4<-rbind(all4_1,all4_2)
all1$new<-paste(all1[,1],all1[,2],sep="_")
all2$new<-paste(all2[,1],all2[,2],sep="_")
all3$new<-paste(all3[,1],all3[,2],sep="_")
all4$new<-paste(all4[,1],all4[,2],sep="_")
hc<-intersect(human$chimp,all2$new)
res_hc<-human[which(human$chimp %in% hc),]
res_hc2<-human[which(!(human$chimp %in% hc)),]
res_hc$lab1<-"chimp"
res_hc2$lab1<-"nochimp"
human<-rbind(res_hc,res_hc2)
hc<-intersect(human$rhesus,all3$new)
res_hc<-human[which(human$rhesus %in% hc),]
res_hc2<-human[which(!(human$rhesus %in% hc)),]
res_hc$lab2<-"rhesus"
res_hc2$lab2<-"norhesus"
human<-rbind(res_hc,res_hc2)
hc<-intersect(human$mouse,all4$new)
res_hc<-human[which(human$mouse %in% hc),]
res_hc2<-human[which(!(human$mouse %in% hc)),]
res_hc$lab3<-"mouse"
res_hc2$lab3<-"nomouse"
human<-rbind(res_hc,res_hc2)
write.table(human,"/home/specific/human.bed",sep="\t",quote=F,row.names=F,col.names=F)

#The scripts used for other species were similar

human<-read.table("/home/specific/human.bed",header = F,sep="\t",quote = "",fill = T)
human$lab<-paste("human",human[,13],human[,14],human[,15],sep="_")
chimp<-read.table("/home/specific/chimp.bed",header = F,sep="\t",quote = "",fill = T)
chimp$lab<-paste(chimp[,13],"chimp",chimp[,14],chimp[,15],sep="_")
rhesus<-read.table("/home/specific/rhesus.bed",header = F,sep="\t",quote = "",fill = T)
rhesus$lab<-paste(rhesus[,13],rhesus[,14],"rhesus",rhesus[,15],sep="_")
mouse<-read.table("/home/specific/mouse.bed",header = F,sep="\t",quote = "",fill = T)
mouse$lab<-paste(mouse[,13],mouse[,15],mouse[,14],"mouse",sep="_")
result<-rbind(human[,c(9:12,16)],chimp[,c(9:12,16)],rhesus[,c(9:12,16)],mouse[,c(9:12,16)])
index<-duplicated(result)
result<-result[!index,]
a<-as.data.frame(table(result[,1]))
id<-as.data.frame(a[which(a[,2]==1),1])
colnames(id)<-"id"
pin<-left_join(id,result,by=c("id"="V9"))
a<-as.data.frame(table(pin[,2]))
id<-as.data.frame(a[which(a[,2]==1),1])
colnames(id)<-"id"
pin2<-left_join(id,pin,by=c("id"="V10"))
a<-as.data.frame(table(pin2[,3]))
id<-as.data.frame(a[which(a[,2]==1),1])
colnames(id)<-"id"
pin3<-left_join(id,pin2,by=c("id"="V11"))
a<-as.data.frame(table(pin3[,4]))
id<-as.data.frame(a[which(a[,2]==1),1])
colnames(id)<-"id"
pin4<-left_join(id,pin3,by=c("id"="V12"))
h<-pin4[which(pin4$lab=="human_nochimp_norhesus_nomouse"),1:4]
inter<-pin4[which(pin4$lab=="human_chimp_norhesus_nomouse"),1:4]
write.table(h,"/home/specific/human_specific.txt",sep="\t",quote=F,row.names=F,col.names=F)
write.table(inter,"/home/specific/ape_specific.txt",sep="\t",quote=F,row.names=F,col.names=F)

#####################################################

###Pipeline for filtering species-specific splice sites
num=$(cat /home/id/humanid.txt)
for ID in $num
do
awk -F "\t" 'NR==FNR{a[$1,$2]=$0;next}{$13=a[$1,$2];print}' /home/specific/human_specific.txt /home/junction/${ID}_human_donor_gene.bed | tr ' ' '\t' | awk -F "\t" '{if($13!=""){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}}' | sort | uniq > /home/specific/${ID}_human_donor.bed
awk -F "\t" 'NR==FNR{a[$1,$2]=$0;next}{$13=a[$1,$2];print}'/home/specific/human_specific.txt /home/junction/${ID}_human_acceptor_gene.bed | tr ' ' '\t' | awk -F "\t" '{if($13!=""){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}}' | sort | uniq > /home/specific/${ID}_human_acceptor.bed
done
#Read coverage (The number of junction reads ≥3)
