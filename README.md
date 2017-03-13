# Helferlein
Small scripts manly for NGS data analysis

# Description of the scripts

## getMappingOverhang.pl

* Purpose: finds and accumulates none-mappable overhangs in RNA-seq data.

Extracts soft clipped bases from an input sam file, constructs, characterizes and evaluates consensus sequence of the overhang.

* Synapsis 
getMappingOverhang.pl [-l INT -c INT -cc INT] [--help|?] [--man] < FILE.sam

* Input: 
  * single-end sequencing data
  * locally mapped sam file (e.g., STAR --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 40 --genomeDir ./STARIDX --readFilesIn reads.fa --alignEndsType Local --alignSoftClipAtReferenceEnds Yes)
  * use only uniquely mapped reads
  * removed PCR duplicates

* Options:

* -l

The minimal length of softclipped overhang to be considered

* -c

The minimal coverage at the overhang anchor site to be considered

* -cc

The minimal coverage within the overhang to be used for consensus sequence quality evaluation

* Output

Output is a table, holding the following infos

* chr chromosome name

* position  chromosome position

* strand	strand [+-]

* length  length of the overhang

* consensusSeq  consensus sequence of the overhang

* abs_conflictPosition  number of positions in the overhang where at least one reads did not agree with the deduced consensus sequence

* rel_conflictPosition  ratio of positions in the overhang where at least one reads did not agree with the deduced consensus sequence

* rel_BaseCountA	ratio of bases in all overhang sequences which are A

* rel_BaseCountConflicting  ratio of bases in all overhang sequences which do not agree with the deduced consensus sequence (only positions with more than -cc INT reads are considered for Numerator and Denominator)

* Usage
./getMappingOverhang.pl -l 4 -c 3 -cc 2 < FILE.sam > OUT.tsv
samtools view FILE.bam | ./getMappingOverhang.pl -l 4 -c 3 -cc 2 > OUT.tsv

* Author
Fabian Amman, fabian@tbi.univie.ac.at
