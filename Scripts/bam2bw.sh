#!/bin/bash
#by sm934
# convertes sam file into sorted bamfiles as well as bw files (chrsize needed)
#./sam2bw.sh file.bam genome.chrsizes sort
#output: file.bam, file.bam.bai, file.bw

#echo $1 ${1%.bam} ' -> echo $1'
#$1 = .bam
#$2 = .chrsizes
#$3 = sorting flag
PATH=/applications/samtools/samtools-1.2/:/applications/bedtools/bedtools2/bin/:$PATH
if [ "$3" = "sort" ]; then
  echo "sorting bam files (tmp file deleted afterward)"
  samtools sort $1 ${1%.bam}_sort
  #mv ${1%.bam}_sort.bam $1
  #samtools index ${1%.bam}_sort.bam
else
  echo "no sorting of bam file"
fi
genomeCoverageBed -bg -ibam ${1%.bam}_sort.bam -g $2 | sort -k1,1 -k2,2n  > ${1%.bam}.bedgraph
bedGraphToBigWig ${1%.bam}.bedgraph $2 ${1%.bam}.bw
rm ${1%.bam}.bedgraph
rm ${1%.bam}_sort.bam

