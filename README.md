# Helferlein
Small scripts manly for NGS data analysis

# Description of the scripts

## getPolyAsites.pl

* Purpose: finds and accumulates none-mappable overhangs in RNA-seq data.

* Input: 
  * single-end sequencing data
  * locally mapped sam file (e.g., STAR --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 40 --genomeDir ./STARIDX --readFilesIn reads.fa --alignEndsType Local --alignSoftClipAtReferenceEnds Yes)
  * use only uniquely mapped reads
  * removed PCR duplicates

* Options:

* Output:
