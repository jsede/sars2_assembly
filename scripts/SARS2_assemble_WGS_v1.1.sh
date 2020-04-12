#!/bin/bash

## 	Author: John-Sebastian eden
##	Contact: js.eden@sydney.edu.au
## 	Date: April 2020
##	Version: v1.1 with WGS primer clipping
##	Comments: The WGS primer clipping protocol is only recommended when coverage is even

#set working directory
wrk=~/Desktop/SARS2_assembly

#set project folder
prj=$wrk/$1

#set database folder
db=$wrk/database

#set virus for filtering
virus=SARS2
blast=NC_045512

#set compute resources where m=GB t=CPU
m=16
t=12

#change into project folder
cd $prj

#quality trim reads with bbduk/38.79
while read i; do
bbduk.sh -Xmx"$m"g threads="$t" in=$prj/1_reads/"$i"_1.fastq.gz in2=$prj/1_reads/"$i"_2.fastq.gz out=stdout.fq ref=phix k=31 hdist=1 | bbduk.sh -Xmx"$m"g threads="$t" interleaved=true in=stdin.fq out=stdout.fq ref=adapters ktrim=r k=17 mink=3 | bbduk.sh -Xmx"$m"g threads="$t" interleaved=true in=stdin.fq out=$prj/1_reads/"$i"_trim_1.fastq.gz out2=$prj/1_reads/"$i"_trim_2.fastq.gz minlen=50 qtrim=rl trimq=23 entropy=0.7
done < $prj/IDs.list

#normalise read coverage for denovo assembly
while read i; do
bbnorm.sh -Xmx"$m"g threads="$t" in=$prj/1_reads/"$i"_trim_1.fastq.gz in2=$prj/1_reads/"$i"_trim_2.fastq.gz out=$prj/1_reads/"$i"_trim_norm_1.fastq.gz out2=$prj/1_reads/"$i"_trim_norm_2.fastq.gz target=100;
done < $prj/IDs.list

#clip WGS primers from ends of reads with cutadapt/1.8.3 and python/2.7.13
while read i; do
cutadapt -a file:$db/WGS.primers.clips.$2.1.fa -A file:$db/WGS.primers.clips.$2.1.fa -o $prj/1_reads/"$i"_trim_temp_1.fastq.gz -p $prj/1_reads/"$i"_trim_temp_2.fastq.gz $prj/1_reads/"$i"_trim_1.fastq.gz $prj/1_reads/"$i"_trim_2.fastq.gz
cutadapt -g file:$db/WGS.primers.clips.$2.2.fa -G file:$db/WGS.primers.clips.$2.2.fa -o $prj/1_reads/"$i"_trim_clip_1.fastq.gz -p $prj/1_reads/"$i"_trim_clip_2.fastq.gz $prj/1_reads/"$i"_trim_temp_1.fastq.gz $prj/1_reads/"$i"_trim_temp_2.fastq.gz
rm $prj/1_reads/"$i"_trim_temp*
done < $prj/IDs.list

#create folder for de novo assembly
mkdir $prj/2_contigs

#denovo assemble trimmed reads with megahit/1.2.9 and rename output contigs with library name
while read i; do
megahit --k-list 21,27,35,45,55 -t "$t" -1 $prj/1_reads/"$i"_trim_norm_1.fastq.gz -2 $prj/1_reads/"$i"_trim_norm_2.fastq.gz -o $prj/2_contigs/megahit_"$i" --out-prefix "$i" --min-contig-len 5000;
sed "s/>k/>"$i"_k/g" $prj/2_contigs/megahit_"$i"/"$i".contigs.fa > $prj/2_contigs/megahit_"$i"/"$i".contigs.megahit.fa;
done < $prj/IDs.list

#export denovo assembly contig stats
while read i; do
grep ">" $prj/2_contigs/megahit_"$i"/"$i".contigs.megahit.fa | cut -d ">" -f2 | tr ' ' '\t' | sed 's/flag\=//g' | sed 's/multi\=//g' | sed 's/len\=//g' > $prj/2_contigs/megahit_"$i"/"$i".contigs.megahit.stats
done < $prj/IDs.list

#create folder for blasting
mkdir $prj/3_filter

#combine denovo contigs into single file
cat $prj/2_contigs/megahit_*/*.contigs.megahit.fa > $prj/3_filter/all_denovo.contigs.megahit.fasta

#blast to local SARS2 database
export BLASTDB=$db
blastn -query $prj/3_filter/all_denovo.contigs.megahit.fasta -db "$blast".ref -evalue 1E-10 -num_threads "$t" -out $prj/3_filter/all_denovo.contigs.megahit.blastn_"$virus".txt -outfmt "6 qseqid qlen stitle sstart send pident length evalue sstrand"

#get top virus blast results for each contig
awk -F$'\t' '!seen[$1]++' $prj/3_filter/all_denovo.contigs.megahit.blastn_"$virus".txt > $prj/3_filter/all_denovo.contigs.megahit.blastn_"$virus".top.txt

#take column containing contig names
cut -f1 $prj/3_filter/all_denovo.contigs.megahit.blastn_"$virus".top.txt > $prj/3_filter/"$virus"_draft_genomes."$1".list 

#retrieve sequences
seqtk subseq $prj/3_filter/all_denovo.contigs.megahit.fasta $prj/3_filter/"$virus"_draft_genomes."$1".list > $prj/3_filter/"$virus"_draft_genomes."$1".fa

#align to reference
mafft --thread $t --reorder --adjustdirection --maxiterate 10 --add $prj/3_filter/"$virus"_draft_genomes."$1".fa $db/"$blast"_WGS.fasta > $prj/3_filter/"$virus"_draft_genomes."$1".ref_aligned.fa

#unalign draft genome alignment
sed '/^>/! s/\-//g' $prj/3_filter/"$virus"_draft_genomes."$1".ref_aligned.fa | sed 's/_R_//g' | awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' > $prj/3_filter/"$virus"_draft_genomes."$1".ref_unaligned.fa

#retrieve draft denovo genomes
seqtk subseq $prj/3_filter/"$virus"_draft_genomes."$1".ref_unaligned.fa $prj/3_filter/"$virus"_draft_genomes."$1".list > $prj/3_filter/"$virus"_draft_genomes."$1".ref_unaligned.reoriented.fa

#create folder for blasting
mkdir $prj/4_remap

#retrieve individual draft sequences for mapping
while read i; do
grep -A1 ">${i}" $prj/3_filter/"$virus"_draft_genomes."$1".ref_unaligned.reoriented.fa > $prj/4_remap/"$i".draft.fa
done < $prj/IDs.list

#map trimmed reads to draft genome
while read i; do
bbmap.sh -Xmx"$m"g threads="$t" in=$prj/1_reads/"$i"_trim_clip_1.fastq.gz in2=$prj/1_reads/"$i"_trim_clip_2.fastq.gz outm=$prj/4_remap/"$i".remapped.bam ref=$prj/4_remap/"$i".draft.fa;
done < $prj/IDs.list

#sort bam file
while read i; do
samtools sort -o $prj/4_remap/"$i".remapped.sorted.bam $prj/4_remap/"$i".remapped.bam && samtools index $prj/4_remap/"$i".remapped.sorted.bam
done < $prj/IDs.list

#call final consensus
while read i; do
samtools mpileup -d 1000 -A -Q 0 $prj/4_remap/"$i".remapped.sorted.bam | ivar consensus -p $prj/4_remap/"$i".final -q 20 -t 0
done < $prj/IDs.list
