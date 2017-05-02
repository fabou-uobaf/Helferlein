# Helferlein
Small scripts manly for NGS data analysis

# Description of the scripts

## getMappingOverhang.pl

### Purpose: 

Finds and accumulates none-mappable overhangs in RNA-seq data.

Extracts soft clipped bases from an input sam file, constructs, characterizes and evaluates consensus sequence of the overhang.

### Synapsis:
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

### Output:

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

### Usage:
./getMappingOverhang.pl -l 4 -c 3 -cc 2 < FILE.sam > OUT.tsv
samtools view FILE.bam | ./getMappingOverhang.pl -l 4 -c 3 -cc 2 > OUT.tsv

### Author:
Fabian Amman, fabian@tbi.univie.ac.at

## deduplicationSeqMap.pl

### Purpose:

Checks a sorted input file for consequtive reads with the same mapping characteristics (chromosome, position, CIGAR string, strand) and the same read sequence. If all above characterisitcs are the same as the preceeding read, the read will not be reported. If one of the features differ, the read is reported. On STDERR short stats are reported (number of total reads, kept reads, removed reads)

### Synapsis:

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

## pairedEndStrandPhasing.pl

### Purpose: 

Switches the strand of alignments in a bam file. Can treat pair/mates and single-end reads differently. 

Useful if sequencing libraries are of type "1+-,1-+,2++,2--" (Lexogene library) or "1++,1--,2+-,2-+" (default Illumina library) and follow up tools consider pair/mates as single reads with their respective strand information.

### Synapsis:

pairedEndStrandPhasing.pl --bams <BAM1>,<BAM2>,...,<BAMn> [--switchfirst] [--switchsecond] [--switchSE] [--samtools <PATH>] [--man] [--help]

### Input:

	* bam files

### Options:

* --bams <FILE>

Comma-separated list of bam files to be processed.

* --switchfirst 

Flag if alignments classified as first in pair should be strand switched (Default: no action)

* --switchsecond

Flag if alignments classified as second in pair should be strand switched (Default: no action)

* --switchSE

Flag if alignments classified as single-end reads (no read pair) should be strand switched (Default: no action)

* --samtools <STRING>

Path to samtools which are used to read and write bam files. Tested with samtools version 1.3 and 1.4. (Default: samtools specified in the \$PATH environment.

* --[help|?]

Print the short version of the man-page.

* --[man]

Print man-page.

### Output:

	* bam files with the same name as the input but with the *.strandPhased.bam suffix

### Usage:

perl ./pairedEndStrandPhasing.pl --bams Ecoli_rep1.bam,Ecoli_rep2.bam --switchfirst --switchsecond --switchSE --samtools ${HOME}/.local/bin/samtools

### Author:

Fabian Amman, fabian@tbi.univie.ac.at
