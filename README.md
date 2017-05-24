# Helferlein
Small scripts mainly for NGS data analysis

# Description of the scripts

## getMappingOverhang.pl

### Purpose: 

Finds and accumulates none-mappable overhangs in RNA-seq reads.

Extracts soft clipped bases from an input sam file, constructs, characterizes and evaluates consensus sequence of the overhang.

### Synapsis:
getMappingOverhang.pl [-l INT -c INT -cc INT -s [-1,+1]] [--help|?] [--man] < FILE.sam

### Input: 
* single-end sequencing data (if used with paired end, a work around is to extract /1 or /2 read first)
* locally mapped sam file (e.g., STAR --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 40 --genomeDir ./STARIDX --readFilesIn reads.fa --alignEndsType Local --alignSoftClipAtReferenceEnds Yes, or bwa -L 4)
* use only uniquely mapped reads (advised, not mandatory, depending on the analysis)
* removed PCR duplicates (advised, not mandatory, depending on the analysis)

### Options:

* -l <INT>

The minimal length of softclipped overhang to be considered

* -c <INT>

The minimal coverage at the overhang anchor site to be reported

* -cc <INT>

The minimal coverage within the overhang to be used for consensus sequence quality evaluation

* -s <STRING>

Defines if the used sequencing library is of a "++,--" or a "+-,-+" set up. Use +1 for the former and -1 for the latter. 

* -e <STRING>

Consider only 3' [3], 5' [5], or both [53] overhangs; (Default:3)

### Output:

Output is a table, holding the following infos in respective column

* chr 

chromosome name

* position  

chromosome position

* strand	

strand [+-]

* length  

length of the overhang

* absolut read counts

number of reads with an considered (see argument -l) overhang at the overhang anchor point

* RPM read counts

reads per million ([absolute read counts]/[total mapped reads / 10^6]) supporting the overhang. 

* consensusSeq  

consensus sequence of the overhang, determined by majority vote in each column of an naive alignment (without any gaps)

* abs_conflictPosition  

number of positions in the overhang where at least one reads did not agree with the deduced consensus sequence

* rel_conflictPosition  

ratio of positions in the overhang where at least one reads did not agree with the deduced consensus sequence

* rel_BaseCountA	

ratio of bases in all overhang sequences which are A

* rel_BaseCountConflicting  

ratio of bases in all overhang sequences (in the reads) which do not agree with the deduced consensus sequence (only positions with more than -cc INT reads are considered for Numerator and Denominator)

### Usage:
./getMappingOverhang.pl -l 4 -c 3 -cc 2 < FILE.sam > OUT.tsv
samtools view FILE.bam | ./getMappingOverhang.pl -l 4 -c 3 -cc 2 > OUT.tsv

### Author:
Fabian Amman, fabian@tbi.univie.ac.at

## deduplicationSeqMap.pl

### Purpose:

Checks a sorted input file for consecutive reads with the same mapping characteristics (chromosome, position, CIGAR string, strand) and the identical read sequence. If all above characterisitcs are the same as the preceeding read, the read will not be reported. If one of the features differ, the read is reported. On STDERR short stats are reported (number of total reads, kept reads, removed reads)

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

Toggle if alignments classified as first in pair should be strand switched (Default: no action)

* --switchsecond

Toggle if alignments classified as second in pair should be strand switched (Default: no action)

* --switchSE

Toggle if alignments classified as single-end reads (no read pair) should be strand switched (Default: no action)

* --samtools <STRING>

Path to samtools which are used to read and write bam files. Tested with samtools version 1.3 and 1.4. (Default: samtools specified in the \$PATH environment.

* --[help|?]

Print the short version of the man-page.

* --[man]

Print man-page.

### Output:

* bam files with the same name as the input but with *.strandPhased.bam suffix; existing files with the same name will probably be overwritten.

### Usage:

perl ./pairedEndStrandPhasing.pl --bams Ecoli_rep1.bam,Ecoli_rep2.bam --switchfirst --switchsecond --switchSE --samtools ${HOME}/.local/bin/samtools

### Author:

Fabian Amman, fabian@tbi.univie.ac.at


## FisherExactEnrichment.R

### Purpose: 

Generally makes a two.sided Fisher's Exact test for each line of a table, where each line holds all values for a 2x2 contingency table. Eventually, the obtained p-values are corrected for multiple testing by the method of Benjamin and Hochberg. 
Particularly, useful to test enrichment or depletion of certain functions in gene sets.

### Synapsis:

R CMD BATCH "--args <FILE>" ../bin/	.R FisherExactEnrichment.Rlog

### Input:

Tab-separated table with the following header:

* functionCode		e.g. arGOC or GO term id
* functionName		e.g. arGOC function name or GO term name
* totalGenes		number of genes in the genome/background (total population to be tested - can also be only expressed genes, genes on the chip, etc.)
* positivGenes		number of genes in the gene subset to be tested for enrichment (e.g., differentially regulated genes)
* totalGenesWithFunction	number of genes in the genome/background associated with given function
* positivGenesWithFunction	number of genes in the gene subset associated with given function

### Options:

none

### Output:

Input table with additional columns:

* enrichmentFactor	oberved number of genes with with given function in sub set devided by expected number. Greater 1 for enrichment, less than 1 for depletion.
* enrichmentSign	verbose classification ('enriched', 'depleted', 'same') for enrichment factors (>3/2, <2/3, <3/2 && >2/3), respectively
* p.value		Fisher's exact two-sided p-value
* q.value		P-value adjusted for multiple testing by Benjamini and Hochberg

### Usage:

R CMD BATCH "--args functionCountTable.csv" ../bin/	.R FisherExactEnrichment.Rlog

### Author:

Fabian Amman, fabian@tbi.univie.ac.at

## NCBI_featureTable2rntptt.pl

### Purpose: 

Usefull to produce old style *.rnt *.ptt files from newer file format feature table as currently provided by NCBI.

### Synapsis:

NCBI_featureTable2rntptt.pl -ft <FILE> [-man] [-help]

### Input: 

NCBI provided feature tables for replicon annotation. Must be gziped.

### Options:

* -ft <FILE>

NCBI provided feature tables for replicon annotation. Must be gziped.

### Output:

Two files with the same basename as the input, with the extensions .rnt (for ncRNA)  and .ptt (for protein coding genes).

### Usage:

wget ftp://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Escherichia_coli/reference/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_feature_table.txt.gz
perl NCBI_featureTable2rntptt.pl -ft GCF_000005845.2_ASM584v2_feature_table.txt.gz

### Author:
Fabian Amman, fabian@tbi.univie.ac.at

## Pfropfen

### Purpose:

Calls positions from RNA-seq data with error patterns indicating RNA modifications. Base substitutions (default), indels and reverse transcription abortion events can be considered. For each considered event the background error model is defined as the percentage of observed errors given the observed reference base. Each event is tested if it occurs more often than expected given the background model based on the binomial distribution. Each single deduced p-value is merged into a final p-value for each site using the Fisher's method. If replicas are provided this calculation is done over all events and samples for the particular site. Significant sites after multiple testing correction are reported in the VCF file format.

### Synapsis:

pfropfen --bams <BAM FILE 1>::<BAM FILE 2>::...::<BAM FILE N> --fasta <REFERENCE FASTA FILE> [--delta <FLOAT> --cov <INT> --qual <FLAOT> --noterm|term --noindel|indel --winsor <INT> --pval <FLOAT> --data <FILE NAME> --out <VCF FILE NAME> -verbose]

### Input:

One or more bam files and the associated referenec genome n fasta format.

### Options:

* -bams <FILE::FILE::...>

List of input bam files to be analysed. If several files are provided, separate the file names with '::'

* -fasta <FILE>

Reference genome sequence in fasta format.

* -delta <FLOAT>

Only sites with less substitutions then the specifies value are considered for the background error model. (Default = 0.5)

* -qual <FLOAT>

Only bases with quality score above the specified value are considered. (Default = 20)

* -cov <INT>

Only sites with a coverage higher than the specified value are considered for modification calling. (Default = 4)

* -term|-noterm

Toggle if premature read termination are also considered. (Default: no)

* -indel|-noindel

Toggle if observed indels are also considered as modification trace. (Default: no)

* -winsor <INT>

The final p-value is calculated with Fisher's Method from the single p-values for each single modification trace (i.e., each base substitution, indels, RT abbortion events). Fisher's Method is thereby only applied to the winsorized set of p-values, meaning the N highest and lowest p-values are remove beforehand. If L replica bam files are provided, L*N are removed. (Default: N = 1)

* -pval <FLOAT>

Only sites with a multiple testing corrected p-value below this value are reported. (Default = 0.01)

* -data <FILE>

File name were observed raw counts are verbosely written to. (Default = <Pfropfen_dataTable.csv>)

* -out <FILE>

File name there the VCF file of the detected, significant modification sites are reported. (Default = <Pfropfen.vcf>)

* -verbose

Toggle report more detailed diagnosics via STDERR. (Default off)

* --[help|?]

Print the short version of the man-page.

* --[man]

Print man-page.

### Ouput:

List of called modification sites in VCF file format.

### Usage:

./pfropfen --bams test1.bam::test2.bam --fasta test.fa --delta 0.5 --cov 4 --qual 20 --noterm --noindel --winsor 1 --pval 0.01 --data Pfropfen_dataTable.csv --out Pfropfen.vcf -verbose 1

### Author:
Fabian Amman, fabian@tbi.univie.ac.at
