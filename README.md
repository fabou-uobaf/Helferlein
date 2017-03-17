# Helferlein
Small scripts manly for NGS data analysis

# Description of the scripts

## getMappingOverhang.pl

### Purpose: 

Finds and accumulates none-mappable overhangs in RNA-seq data.

Extracts soft clipped bases from an input sam file, constructs, characterizes and evaluates consensus sequence of the overhang.

### Synapsis 
getMappingOverhang.pl [-l INT -c INT -cc INT -s [-1,+1]] [--help|?] [--man] < FILE.sam

### Input: 
  * single-end sequencing data
  * locally mapped sam file (e.g., STAR --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 40 --genomeDir ./STARIDX --readFilesIn reads.fa --alignEndsType Local --alignSoftClipAtReferenceEnds Yes)
  * use only uniquely mapped reads
  * removed PCR duplicates

### Options:

* -l <INT>

The minimal length of softclipped overhang to be considered

* -c <INT>

The minimal coverage at the overhang anchor site to be considered

* -cc <INT>

The minimal coverage within the overhang to be used for consensus sequence quality evaluation

* -s <STRING>

Defines if the used sequencing library is of a "++,--" or a "+-,-+" set up. Use +1 for the former and -1 for the latter. 

* -e <STRING>

Consider only 3' [3], 5' [5[, or both [53] overhangs; (Default:3)

### Output

Output is a table, holding the following infos

* chr 

chromosome name

* position  

chromosome position

* strand	

strand [+-]

* length  

length of the overhang

* absolut read counts

number of reads showing overhang at position in focus

* RPM read counts

reads per million (see absolute read counts) supproting the overhang. 

* consensusSeq  

consensus sequence of the overhang

* abs_conflictPosition  

number of positions in the overhang where at least one reads did not agree with the deduced consensus sequence

* rel_conflictPosition  

ratio of positions in the overhang where at least one reads did not agree with the deduced consensus sequence

* rel_BaseCountA	

ratio of bases in all overhang sequences which are A

* rel_BaseCountConflicting  

ratio of bases in all overhang sequences which do not agree with the deduced consensus sequence (only positions with more than -cc INT reads are considered for Numerator and Denominator)

### Usage
./getMappingOverhang.pl -l 4 -c 3 -cc 2 < FILE.sam > OUT.tsv
samtools view FILE.bam | ./getMappingOverhang.pl -l 4 -c 3 -cc 2 > OUT.tsv

### Author
Fabian Amman, fabian@tbi.univie.ac.at

## deduplicationSeqMap.pl

### Purpose

Checks a sorted input file for consequtive reads with the same mapping characteristics (chromosome, position, CIGAR string, strand) and the same read sequence. If all above characterisitcs are the same as the preceeding read, the read will not be reported. If one of the features differ, the read is reported. On STDERR short stats are reported (number of total reads, kept reads, removed reads)

### Synapsis 

deduplicationSeqMap.pl <FILE.sam> 

### Input: 

Sorted single end SAM file

### Options:

NONE

### Output

SAM file.

### Usage

samtools view -h FILE.bam | ./deduplicationSeqMap.pl 2> DeDup.stats | samtools view -bS - | samtools sort - FILE.dedup

### Author
Fabian Amman, fabian@tbi.univie.ac.at
