###Quality control of RNA-seq data

num=$(cat /home/id/samplesid.txt) #including five tissues-across humans, chimpanzees, rhesus macaques, and mice
for ID in $num
do
fastp -i /home/fastq/${ID}.fastq -o /home/fastq/${ID}_filter.fastq -h /home/fastq/${ID}_report.html -j /home/fastq/${ID}_report.json -q 30
done

###RNA-seq reads alignment
#Build index
hisat2-build /home/genome/hg38.fa /home/hisat/index/human/human_index
hisat2-build /home/genome/panTro6.fa /home/hisat/index/chimp/chimp_index
hisat2-build /home/genome/rheMac10.fa /home/hisat/index/rhesus/rhesus_index
hisat2-build /home/genome/mm39.fa /home/kyoku/data3/hisat/index/mouse/mouse_index
#Alignment, such as in human
num=$(cat /home/id/humanid.txt)
for ID in $num
do
hisat2 -p 4 -x /home/hisat/index/human/human_index -U /home/fastq/${ID}_filter.fastq -S /home/hisat/${ID}_hisat2_human.sam
samtools sort -@ 2 -O bam -o /home/hisat/${ID}_hisat2_human.bam /home/hisat/${ID}_hisat2_human.sam
done

###Identify candidate splice site by extracting jucntion reads, such as in human
num=$(cat /home/id/humanid.txt)
for ID in $num
do
samtools view /home/hisat/${ID}_hisat2_human.bam | awk -F "\t" '{if($6~"^[0-9]*M[0-9]*N[0-9]*M$") print}' | awk -F "\t" '{split ($6,T,"[M-N]");$6=T[1] OFS T[2] OFS T[3]}1' OFS="\t" | awk -F "\t" 'BEGIN{OFS="\t"}{$1=$1"_"FNR;$6=$4+$6-1;$7=$6+$7+1;$8=$7+$8-1}{print $0}' | awk -F "\t" '{print$1"\t"$3"\t"$4"\t"$6"\t"$7"\t"$8}' > /home/junction/${ID}_human.bed
done

###Cross-species genomic coordinate comparisons using LiftOver
#Convert splice site positions from one species to their orthologous positions in another, such as human to chimpanzee
num=$(cat /home/id/humanid.txt)
for ID in $num
do
cat /home/junction/${ID}_human.bed | awk -F "\t" '{print $2"\t"$4-5"\t"$4+5}' | sort | uniq > /home/junction/${ID}_human_donor.bed
cat /home/junction/${ID}_human.bed | awk -F "\t" '{print $2"\t"$5-5"\t"$5+5}' | sort | uniq > /home/junction/${ID}_human_acceptor.bed
bedtools intersect -a /home/junction/${ID}_human_donor.bed -b /home/genome/hg38_trans.gtf -wa -wb | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$7}' | sort | uniq | awk -F "\t" 'BEGIN{OFS="\t"}{$4=$4"\t""id_"FNR}{print $0}' | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$5"\t""100""\t"$4}' > /home/junction/${ID}_human_donor_gene.bed
bedtools intersect -a /home/junction/${ID}_human_acceptor.bed -b /home/genome/hg38_trans.gtf -wa -wb | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$7}' | sort | uniq | awk -F "\t" 'BEGIN{OFS="\t"}{$4=$4"\t""id_"FNR}{print $0}' | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$5"\t""100""\t"$4}' > /home/junction/${ID}_human_acceptor_gene.bed
liftOver /home/junction/${ID}_human_donor_gene.bed /home/liftover/hg38ToPanTro6.over.chain.gz /home/liftover/${ID}_human_chimp_donor.bed unMapped
liftOver /home/junction/${ID}_human_acceptor_gene.bed /home/liftover/hg38ToPanTro6.over.chain.gz /home/liftover/${ID}_human_chimp_acceptor.bed unMapped
done
#Convert positions back to the original genome, such as chimpanzee back to human
num=$(cat /home/id/humanid.txt)
for ID in $num
do
liftOver /home/liftover/${ID}_human_chimp_donor.bed /home/liftover/panTro6ToHg38.over.chain.gz /home/liftover/back/${ID}_human_chimp_donor.bed unMapped
liftOver /home/liftover/${ID}_human_chimp_acceptor.bed /home/liftover/panTro6ToHg38.over.chain.gz /home/liftover/back/${ID}_human_chimp_acceptor.bed unMapped
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

###Identify candidate species-specific splice sites
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

###Pipeline for filtering species-specific splice sites
#One side of the jucntion reads is annotated
num=$(cat /home/id/humanid.txt)
for ID in $num
do
awk -F "\t" 'NR==FNR{a[$1,$2]=$0;next}{$7=a[$2,$4];print}' /home/specific/human_specific.txt /home/junction/${ID}_human.bed | tr ' ' '\t' | awk -F "\t" '{if($7!=""){print $2"\t"$4"\t"$5}}' | sort | uniq > /home/filter/tmp/${ID}_human_donor.bed
awk -F "\t" 'NR==FNR{a[$1,$2]=$0;next}{$7=a[$2,$5];print}' /home/specific/human_specific.txt /home/junction/${ID}_human.bed | tr ' ' '\t' | awk -F "\t" '{if($7!=""){print $2"\t"$4"\t"$5}}' | sort | uniq > /home/filter/tmp/${ID}_human_acceptor.bed
done

cat /home/filter/tmp/*_human_donor.bed | sort | uniq > /home/filter/tmp/allhuman_donor.bed
cat /home/filter/tmp/*_human_acceptor.bed | sort | uniq > /home/filter/tmp/allhuman_acceptor.bed

library(dplyr)
library(tidyr)
gtf<-read.table("/home/genome/hg38.gtf",header = F,sep="\t",quote = "",fill = T)
qian<-read.table("/home/filter/tmp/allhuman_donor.bed",header = F,sep="\t",quote = "",fill = T)
pin<-left_join(qian,gtf,by=c("V1"="V1","V3"="V3"))
pin<-pin[,c(1,2,3,6)]
index<-duplicated(pin)
pin<-pin[!index,]
qian<-na.omit(pin)
pin<-NULL
for(i in 1:nrow(qian)){
res1<-gtf[which(qian[i,1] == gtf[,1] & qian[i,2] >= gtf[,3] & qian[i,2] <= gtf[,4]),6]
res2<-gtf[which(qian[i,1] == gtf[,1] & qian[i,3] >= gtf[,3] & qian[i,3] <= gtf[,4]),6]
res<-intersect(res1,res2)
if(length(res)>0){
pin<-rbind(qian[i,],pin)
}
print(i)
}
qian<-pin
qian$id<-paste(qian[,1],qian[,2],sep="_")
id<-paste(qian[,1],qian[,2],sep="_")
index<-duplicated(id)
id<-id[!index]
for(i in 1:length(id)){
res<-qian[which(qian$id==id[i]),4]
index<-duplicated(res)
res<-res[!index]
if(length(res)>1){
print(i)
}
}
hou<-read.table("/home/filter/tmp/allhuman_acceptor.bed",header = F,sep="\t",quote = "",fill = T)
pin<-left_join(hou,gtf,by=c("V1"="V1","V2"="V4"))
pin<-pin[,c(1,2,3,6)]
index<-duplicated(pin)
pin<-pin[!index,]
hou<-na.omit(pin)
pin<-NULL
for(i in 1:nrow(hou)){
res1<-gtf[which(hou[i,1] == gtf[,1] & hou[i,2] >= gtf[,3] & hou[i,2] <= gtf[,4]),6]
res2<-gtf[which(hou[i,1] == gtf[,1] & hou[i,3] >= gtf[,3] & hou[i,3] <= gtf[,4]),6]
res<-intersect(res1,res2)
if(length(res)>0){
pin<-rbind(hou[i,],pin)
}
print(i)
}
hou<-pin
hou$id<-paste(hou[,1],hou[,3],sep="_")
id<-paste(hou[,1],hou[,3],sep="_")
index<-duplicated(id)
id<-id[!index]
for(i in 1:length(id)){
res<-hou[which(hou$id==id[i]),4]
index<-duplicated(res)
res<-res[!index]
if(length(res)>1){
print(i)
}
}
hou<-hou[,c(1,3,4)]
write.table(qian,"/home/filter/tmp/allhuman_donor2.bed",sep="\t",quote=F,row.names=F,col.names=F)
write.table(hou,"/home/filter/tmp/allhuman_acceptor2.bed",sep="\t",quote=F,row.names=F,col.names=F)

num=$(cat /home/id/humanid.txt)
for ID in $num
do
awk -F "\t" 'NR==FNR{a[$1,$2,$3]=$0;next}{$13=a[$1,$2,$3];print}' /home/filter/tmp/allhuman_donor2.bed /home/filter/tmp/${ID}_human_donor.bed | tr ' ' '\t' | awk '$13!=""'  | awk -F "\t" '{print $13"\t"$14"\t"$15"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' > /home/filter/first/${ID}_human_donor.bed
awk -F "\t" 'NR==FNR{a[$1,$2,$3]=$0;next}{$13=a[$1,$2,$3];print}' /home/filter/tmp/allhuman_acceptor2.bed /home/filter/tmp/${ID}_human_acceptor.bed | tr ' ' '\t' | awk '$13!=""'  | awk -F "\t" '{print $13"\t"$14"\t"$15"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' > /home/filter/first/${ID}_human_acceptor.bed
done

#Read coverage (The number of junction reads ≥3) & Read quality (Overhang ≥5 bp)
num=$(cat /home/id/humanid.txt)
for ID in $num
do
awk -F "\t" 'NR==FNR{a[$1,$2]=$0;next}{$7=a[$2,$4];print}' /home/filter/first/${ID}_human_donor.bed /home/junction/${ID}_human.bed | tr ' ' '\t' | awk -F "\t" '{if($7!=""){print $0"\t"$4-$3+1"\t"$6-$5+1}}' | awk -F "\t" '{if($9>4 && $10>4) print $1"\t"$2"\t"$4}' | awk -F "\t" '{count[$2"\t"$3]++}{print$0"\t"count[$2"\t"$3]}' | awk -F "\t" '{if($4>2){print$2"\t"$3}}' | sort | uniq > /home/filter/tmp/${ID}_human_donor.bed
awk -F "\t" 'NR==FNR{a[$1,$2]=$0;next}{$7=a[$2,$5];print}' /home/filter/first/${ID}_human_acceptor.bed /home/junction/${ID}_human.bed | tr ' ' '\t' | awk -F "\t" '{if($7!=""){print $0"\t"$4-$3+1"\t"$6-$5+1}}' | awk -F "\t" '{if($9>4 && $10>4) print $1"\t"$2"\t"$5}' | awk -F "\t" '{count[$2"\t"$3]++}{print$0"\t"count[$2"\t"$3]}' | awk -F "\t" '{if($4>2){print$2"\t"$3}}' | sort | uniq > /home/filter/tmp/${ID}_human_acceptor.bed
awk -F "\t" 'NR==FNR{a[$1,$2]=$0;next}{$13=a[$1,$2];print}' /home/filter/tmp/${ID}_human_donor.bed /home/filter/first/${ID}_human_donor.bed | tr ' ' '\t' | awk -F "\t" '{if($13!=""){print$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}}' > /home/filter/second/${ID}_human_donor.bed
awk -F "\t" 'NR==FNR{a[$1,$2]=$0;next}{$13=a[$1,$2];print}' /home/filter/tmp/${ID}_human_acceptor.bed /home/filter/first/${ID}_human_acceptor.bed | tr ' ' '\t' | awk -F "\t" '{if($13!=""){print$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}}' > /home/filter/second/${ID}_human_acceptor.bed
done

#Observed in over half of the samples in a tissue
num=$(cat /home/id/humanid.txt)
for ID in $num
do
cat <(cat /home/junction/${ID}_human.bed | awk -F "\t" '{print$2"\t"$4}') <(cat /home/junction/${ID}_human.bed | awk -F "\t" '{print$2"\t"$5}') | sort | uniq | awk -F "\t" '{print$0"\t"$1"_"$2"_"$3}' > /home/junction/${ID}_human_splice.bed
done
cat /home/filter/second/${ID}_human_donor.bed | awk -F "\t" '{print$1"_"$2"_"$3}' | sort | uniq > /home/filter/tmp/human.txt
cat /home/filter/tmp/human.txt  | while read line1
do
cc_h=`awk -F "\t" '{print$1"_human_splice.bed"}' /home/id/cerebellar_cortex_id.bed | xargs -I {} grep -l ${line1} /home/junction/{} | sort | uniq | wc -l`
k_h=`awk -F "\t" '{print$1"_human_splice.bed"}' /home//id/kidney_id.bed | xargs -I {} grep -l ${line1} /home/junction/{} | sort | uniq | wc -l`
m_h=`awk -F "\t" '{print$1"_human_splice.bed"}' /home/id/muscle_id.bed | xargs -I {} grep -l ${line1} /home/junction/{} | sort | uniq | wc -l`
pc_h=`awk -F "\t" '{print$1"_human_splice.bed"}' /home/id/prefrontal_cortex_id.bed | xargs -I {} grep -l ${line1} /home/junction/{} | sort | uniq | wc -l`
pvc_h=`awk -F "\t" '{print$1"_human_splice.bed"}' /home/id/primary_visual_cortex_id.bed | xargs -I {} grep -l ${line1} /home/junction/{} | sort | uniq | wc -l`
echo -e ${cc_h}"\t"${k_h}"\t"${m_h}"\t"${pc_h}"\t"${pvc_h} >> /home/filter/third/human_count.txt
done
library(dplyr)
library(tidyr)
h<-read.table("/home/filter/third/human_count.txt",header = F,sep="\t",quote = "",fill = T)
chu<-h[which(h[,2]>3 | h[,3]>3 | h[,4]>3 | h[,5]>3 | h[,6]>3),]
write.table(chu,"/home/filter/third/human_count_overhalf.txt",sep="\t",quote=F,row.names=F,col.names=F)

#Focused on dinucleotides
#Based on "human_count_overhalf.txt" file, add the position and strand information in each species to "/home/filter/fourth/human_tmp.txt"
cat /home/filter/fourth/human_tmp.txt | awk -F "\t" '{if($3=="+") print}' | awk -F "\t" '{print $1"\t"$2-3"\t"$2+6"\t"$3}' > /home/filter/fourth/human_donor_plus.bed
cat /home/filter/fourth/human_tmp.txt | awk -F "\t" '{if($3=="-") print}' | awk -F "\t" '{print $1"\t"$2-3"\t"$2+18"\t"$3}' > /home/filter/fourth/human_acceptor_minus.bed
cat /home/filter/fourth/human_tmp.txt | awk -F "\t" '{if($3=="+") print}' | awk -F "\t" '{print $1"\t"$2-19"\t"$2+2"\t"$3}' > /home/filter/fourth/human_acceptor_plus.bed
cat /home/filter/fourth/human_tmp.txt | awk -F "\t" '{if($3=="-") print}' | awk -F "\t" '{print $1"\t"$2-7"\t"$2+2"\t"$3}' > /home/filter/fourth/human_donor_minus.bed
seqtk subseq /home/genome/hg38.fa /home/filter/fourth/human_donor_plus.bed > /home/filter/fourth/human_donor_plus.fa
seqtk subseq /home/genome/hg38.fa /home/filter/fourth/human_donor_minus.bed > /home/filter/fourth/human_donor_minus.fa
seqtk subseq /home/genome/hg38.fa /home/filter/fourth/human_acceptor_plus.bed > /home/filter/fourth/human_acceptor_plus.fa
seqtk subseq /home/genome/hg38.fa /home/filter/fourth/human_acceptor_minus.bed > /home/filter/fourth/human_acceptor_minus.fa
seqtk seq -r /home/filter/fourth/human_donor_minus.fa > /home/filter/fourth/human_donor_plus2.fa
seqtk seq -r /home/filter/fourth/human_acceptor_minus.fa > /home/filter/fourth/human_acceptor_plus2.fa
cat /home/filter/fourth/human_donor_plus.fa /home/filter/fourth/human_donor_plus2.fa > /home/filter/fourth/human_donor.fa
cat /home/filter/fourth/human_acceptor_plus.fa /home/filter/fourth/human_acceptor_plus2.fa > /home/filter/fourth/human_acceptor.fa
#The scripts used for other species were similar

library(tidyr)
library(dplyr)
hd<-read.table("/home/filter/fourth/human_donor.fa",header = F,sep="\t",quote = "",fill = T)
ha<-read.table("/home/filter/fourth/human_acceptor.fa",header = F,sep="\t",quote = "",fill = T)
cd<-read.table("/home/filter/fourth/chimp_donor.fa",header = F,sep="\t",quote = "",fill = T)
ca<-read.table("/home/filter/fourth/chimp_acceptor.fa",header = F,sep="\t",quote = "",fill = T)
rd<-read.table("/home/filter/fourth/rhesus_donor.fa",header = F,sep="\t",quote = "",fill = T)
ra<-read.table("/home/filter/fourth/rhesus_acceptor.fa",header = F,sep="\t",quote = "",fill = T)
md<-read.table("/home/filter/fourth/mouse_donor.fa",header = F,sep="\t",quote = "",fill = T)
ma<-read.table("/home/filter/fourth/mouse_acceptor.fa",header = F,sep="\t",quote = "",fill = T)
all<-read.table("/home/filter/fourth/human_tmp.txt",header = F,sep="\t",quote = "",fill = T)
id<-as.data.frame(hd[seq(1,nrow(hd),2),])
seq<-as.data.frame(hd[seq(0,nrow(hd),2),])
hd1<-cbind(id,seq)
hd1[,1]<-gsub(">","",hd1[,1])
colnames(hd1)<-c("a","b")
hd1<-hd1 %>% separate(a, into = c("chr","pos"), sep = ":")
hd1<-hd1 %>% separate(pos, into = c("pos1","pos2"), sep = "-")
hd1$new1<-as.numeric(hd1$pos1) + 2
hd1$new2<-as.numeric(hd1$pos2) - 2
id<-as.data.frame(ha[seq(1,nrow(ha),2),])
seq<-as.data.frame(ha[seq(0,nrow(ha),2),])
ha1<-cbind(id,seq)
ha1[,1]<-gsub(">","",ha1[,1])
colnames(ha1)<-c("a","b")
ha1<-ha1 %>% separate(a, into = c("chr","pos"), sep = ":")
ha1<-ha1 %>% separate(pos, into = c("pos1","pos2"), sep = "-")
ha1$new1<-as.numeric(ha1$pos1) + 2
ha1$new2<-as.numeric(ha1$pos2) - 2
#The scripts used for other species were similar
h<-NULL
for(i in 1:nrow(hd1)){
  res1<-all[which(all[,2]==hd1[i,5] & all[,1]==hd1[i,1]),]
  res2<-all[which(all[,2]==hd1[i,6] & all[,1]==hd1[i,1]),]
  res<-rbind(res1,res2)
  if(nrow(res)==1){
    res3<-cbind(res,hd1[i,1:4])
  }
  else{
    if(res1[1,3]!=res2[1,3]){
      print(i)
    }
    else{
      if(res1[1,3]=="+"){
        res<-res1
        res3<-cbind(res,hd1[i,1:4])
      }
      else{
        res<-res2
        res3<-cbind(res,hd1[i,1:4])
      }
    }
  }
  h<-rbind(h,res3)
}
for(i in 1:nrow(ha1)){
  res1<-all[which(all[,2]==ha1[i,5] & all[,1]==ha1[i,1]),]
  res2<-all[which(all[,2]==ha1[i,6] & all[,1]==ha1[i,1]),]
  res<-rbind(res1,res2)
  res3<-cbind(res,ha1[i,1:4])
  h<-rbind(h,res3)
}
hc<-NULL
for(i in 1:nrow(cd1)){
  res1<-h[which(h[,5]==cd1[i,5] & h[,4]==cd1[i,1]),]
  res2<-h[which(h[,5]==cd1[i,6] & h[,4]==cd1[i,1]),]
  res<-rbind(res1,res2)
  if(nrow(res)==1){
    res3<-cbind(res,cd1[i,1:4])
  }
  else{
    print(i)
    if(res1[1,6]!=res2[1,6]){
      print("no")
    }
    else{
      if(res1[1,6]=="+"){
        res<-res1
        res3<-cbind(res,cd1[i,1:4])
      }
      else{
        res<-res2
        res3<-cbind(res,cd1[i,1:4])
      }
    }
  }

  hc<-rbind(hc,res3)
}
for(i in 1:nrow(ca1)){
  res1<-h[which(h[,5]==ca1[i,5] & h[,4]==ca1[i,1]),]
  res2<-h[which(h[,5]==ca1[i,6] & h[,4]==ca1[i,1]),]
  res<-rbind(res1,res2)
  res3<-cbind(res,ca1[i,1:4])
  hc<-rbind(hc,res3)
}
hcr<-NULL
for(i in 1:nrow(rd1)){
  res1<-hc[which(hc[,8]==rd1[i,5] & hc[,7]==rd1[i,1]),]
  res2<-hc[which(hc[,8]==rd1[i,6] & hc[,7]==rd1[i,1]),]
  res<-rbind(res1,res2)
if(nrow(res)==1){
  res3<-cbind(res,rd1[i,1:4])
}
else{
  print(i)
  if(res1[1,9]!=res2[1,9]){
    print("no")
  }
  else{
    if(res1[1,9]=="+"){
      res<-res1
      res3<-cbind(res,rd1[i,1:4])
    }
    else{
      res<-res2
      res3<-cbind(res,rd1[i,1:4])
    }
  }
}
hcr<-rbind(hcr,res3)
}
for(i in 1:nrow(ra1)){
  res1<-hc[which(hc[,8]==ra1[i,5] & hc[,7]==ra1[i,1]),]
  res2<-hc[which(hc[,8]==ra1[i,6] & hc[,7]==ra1[i,1]),]
  res<-rbind(res1,res2)
  res3<-cbind(res,ra1[i,1:4])
  hcr<-rbind(hcr,res3)
}
result<-NULL
for(i in 1:nrow(md1)){
  res1<-hcr[which(hcr[,11]==md1[i,5] & hcr[,10]==md1[i,1]),]
  res2<-hcr[which(hcr[,11]==md1[i,6] & hcr[,10]==md1[i,1]),]
  res<-rbind(res1,res2)
if(nrow(res)==1){
  res3<-cbind(res,md1[i,1:4])
}
else{
  print(i)
  if(res1[1,12]!=res2[1,12]){
    print(i)
  }
  else{
    if(res1[1,12]=="+"){
      res<-res1
      res3<-cbind(res,md1[i,1:4])
    }
    else{
      res<-res2
      res3<-cbind(res,md1[i,1:4])
    }
  }
}
result<-rbind(result,res3)
}
for(i in 1:nrow(ma1)){
  res1<-hcr[which(hcr[,11]==ma1[i,5] & hcr[,10]==ma1[i,1]),]
  res2<-hcr[which(hcr[,11]==ma1[i,6] & hcr[,10]==ma1[i,1]),]
  res<-rbind(res1,res2)
  res3<-cbind(res,ma1[i,1:4])
  result<-rbind(result,res3)
}
colnames(result)<-c(letters,"a1","b1")
library(stringdist)
scorec<-NULL
scorer<-NULL
scorem<-NULL
for(i in 1:nrow(result)){
  seq1<-stringdist::stringsim(tolower(result[i,16]),tolower(result[i,20]))
  seq2<-stringdist::stringsim(tolower(result[i,16]),tolower(result[i,24]))
  seq3<-stringdist::stringsim(tolower(result[i,16]),tolower(result[i,28]))
  scorec<-rbind(scorec,seq1)
  scorer<-rbind(scorer,seq2)
  scorem<-rbind(scorem,seq3)
  print(i)
}
final<-cbind(result,scorec,scorer,scorem)
write.table(final,"/home/filter/fourth/human_identity.bed",sep="\t",quote = F,row.names = F,col.names = F)

#Retain the variants with AF >0.1 as candidate variants associated with the identified splice-site usage
bedtools intersect -a /home/filter/fourth/human_identity.bed -b /home/gnomad/gnomad.v4.1_af.vcf -wa -wb > /home/filter/fourth/human_af.bed

all1<-read.table("/home/filter/fourth/human_identity.bed",header = F,sep="\t",quote = "",fill = T)
all<-read.table("/home/filter/fourth/human_af.bed",header = F,sep="\t",quote = "",fill = T)
donor<-all[which(nchar(all[,4])==9),]
donor[,18]<-gsub("AF=","",donor[,18])
donor[,18]<-as.numeric(donor[,18])
donor<-donor[which(donor[,18]>0.1),]
donor<-donor[which(nchar(donor[,14])==1),]
donor<-donor[which(nchar(donor[,15])==1),]
all2<-donor
pin1<-NULL
for(i in 1:nrow(all2)){
if(all2[i,10]=="+"){
res<-all2[i,]
res2<-all2[i,]
ins<-all2[i,12]-all2[i,2]+1
substr(res[1,4],ins,ins)<-res[1,15]
res$lab<-"mut"
res2$lab<-"ref"
res3<-rbind(res,res2)
} else{
res<-all2[i,]
res2<-all2[i,]
ins<-all2[i,12]-all2[i,2]
ins<-9-ins
t<-res[1,15]
if(t=="C"){
t<-"G"
} else if(t=="G"){
t<-"C"
} else if(t=="A"){
t<-"T"
} else if(t=="T"){
t<-"A"
}
substr(res[1,4],ins,ins)<-t
res$lab<-"mut"
res2$lab<-"ref"
res3<-rbind(res,res2)
}
pin1<-rbind(pin1,res3)
}
d<-pin1
donor<-all[which(nchar(all[,4])==21),]
donor[,18]<-gsub("AF=","",donor[,18])
donor[,18]<-as.numeric(donor[,18])
donor<-donor[which(donor[,18]>0.1),]
donor<-donor[which(nchar(donor[,14])==1),]
donor<-donor[which(nchar(donor[,15])==1),]
all2<-donor
pin1<-NULL
for(i in 1:nrow(all2)){
if(all2[i,10]=="+"){
res<-all2[i,]
res2<-all2[i,]
ins<-all2[i,12]-all2[i,2]+1
substr(res[1,4],ins,ins)<-res[1,15]
res$lab<-"mut"
res2$lab<-"ref"
res3<-rbind(res,res2)
} else{
res<-all2[i,]
res2<-all2[i,]
ins<-all2[i,12]-all2[i,2]
ins<-21-ins
t<-res[1,15]
if(t=="C"){
t<-"G"
} else if(t=="G"){
t<-"C"
} else if(t=="A"){
t<-"T"
} else if(t=="T"){
t<-"A"
}
substr(res[1,4],ins,ins)<-t
res$lab<-"mut"
res2$lab<-"ref"
res3<-rbind(res,res2)
}
pin1<-rbind(pin1,res3)
}
a<-pin1
pin1<-rbind(d,a)
pin<-full_join(all1,pin1,by = c("V13"="V1","V14"="V2","V15"="V3"))
pin<-pin[,c(1:35,39,40,41,42,43,46,47)]
chimp<-NULL
for(i in 1:nrow(pin)){
if(is.na(pin[i,32])){
  if(pin[i,20]==pin[i,24]){
  h<-toupper(unlist(strsplit(pin[i,16],split = "")))
  c<-toupper(unlist(strsplit(pin[i,20],split = "")))
  pos1<-which(h != c)
  if(length(pos1)==1){
  diff1<-data.frame(Position = pos1, string1 = h[pos1], string2 = c[pos1])
  diff<-cbind(pin[i,],diff1)
  chimp<-rbind(chimp,diff)
}
}
}
 else{
  if(pin[i,20]==pin[i,24]){
  h<-toupper(unlist(strsplit(pin[i,32],split = "")))
  c<-toupper(unlist(strsplit(pin[i,33],split = "")))
  pos1<-which(h != c)
  if(length(pos1)==1){
  diff1<-data.frame(Position = pos1, string1 = h[pos1], string2 = c[pos1])
  diff<-cbind(pin[i,],diff1)
  chimp<-rbind(chimp,diff)}
}
}
}
you<-chimp[!is.na(chimp[,32]),]
you<-you[,c(1:12,13:15,17:19,21:23,25:27,32:45)]
no<-chimp[is.na(chimp[,32]),]
no<-no[,c(1:12,13:15,17:19,21:23,25:27,16,20,24,28,43:45)]
you<-you[,c(1:28,36:38)]
colnames(you)<-c(letters,"a1","b1","c1","d1","e1")
colnames(no)<-c(letters,"a1","b1","c1","d1","e1")
result<-rbind(you,no)
acc1<-result[which(nchar(result[,25])==21 & result[,29]==18),]
acc2<-result[which(nchar(result[,25])==21 & result[,29]==17),]
do1<-result[which(nchar(result[,25])==9 & result[,29]==4),]
do2<-result[which(nchar(result[,25])==9 & result[,29]==5),]
final<-rbind(acc1,acc2,do1,do2)
pin<-full_join(final,all,by = c("a"="V8","b"="V9","c"="V10"))
pin$id<-paste(pin[,1],pin[,2],sep="_")
id<-pin$id
index<-duplicated(id)
id<-id[!index]
shuchu<-NULL
for(i in 1:length(id)){
res<-pin[which(pin$id==id[i]),]
res<-na.omit(res)
if(nrow(res)>0){
res[,46]<-gsub("AF=","",res[,46])
res[,46]<-as.numeric(res[,46])
if(nchar(res[1,25])==21){
if(res[1,3]=="+"){
res$ins<-res[,40]-res[,14]+1
res<-res[which(res$c1==res$ins),]
if(nrow(res)>0){
af1<-sum(res[,46])
af<-1-af1
} else{
af<-1
}
} else{
res$ins<-res[,40]-res[,14]
res$ins<-21-res$ins
res<-res[which(res$c1==res$ins),]
if(nrow(res)>0){
af1<-sum(res[,46])
af<-1-af1
} else{
af<-1
}
}
} else if(nchar(res[1,25])==9){
if(res[1,3]=="+"){
res$ins<-res[,40]-res[,14]+1
res<-res[which(res$c1==res$ins),]
if(nrow(res)>0){
af1<-sum(res[,46])
af<-1-af1
} else{
af<-1
}
} else{
res$ins<-res[,40]-res[,14]
res$ins<-9-res$ins
res<-res[which(res$c1==res$ins),]
if(nrow(res)>0){
af1<-sum(res[,46])
af<-1-af1
} else{
af<-1
}
}
}
ti<-pin[which(pin$id==id[i]),c(1:31)]
ti$af<-af
} else{
ti<-pin[which(pin$id==id[i]),c(1:31)]
ti$af<-1
}
shuchu<-rbind(shuchu,ti)
}
index<-duplicated(shuchu)
shuchu<-shuchu[!index,]
shuchu<-na.omit(shuchu)
write.table(shuchu,"/home/filter/fourth/human_onediff.txt",sep="\t",quote=F,row.names=F,col.names=F)
