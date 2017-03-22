# parse comand line argument for basename
args     <- commandArgs(TRUE)
basename <- args[length(args)]
print(basename)

# read in data table
filename <- paste(basename, sep="")
data <- read.table(file=filename, sep="\t", header=T)

# define functions
FisherVector <- function(X) {fisher.test(matrix(c(as.numeric(X[3])-as.numeric(X[4]), as.numeric(X[4]), as.numeric(X[5])-as.numeric(X[6]), as.numeric(X[6])), ncol=2, byrow=F), alternative = "two.sided")$p.value}
Enrichment   <- function(X) {(as.numeric(X[6])/as.numeric(X[5]))/(as.numeric(X[4])/as.numeric(X[3]))}

# calculate scores
data$enrichmentFactor <- apply(data, 1, Enrichment)
data$enrichmentSign <- ifelse(data$enrichmentFactor > 2/3 & data$enrichmentFactor < 3/2, 'same', ifelse(data$enrichmentFactor>1, 'enriched', 'depleted'))
data$p.value <- apply(data, 1, FisherVector)
data$q.value <- p.adjust(data$p.value, method ="fdr")

# print output
resname <- paste(basename, ".Rres", sep="")
write.table(as.data.frame(data),file=resname, sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

## input table with 6 columns
## function abbrev.
### function name
### number of genes in genome
### number of genes in selectio set
### number of genes with function in genome
### number of genes with function in selection set
