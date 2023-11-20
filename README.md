# DADA2

Load required libraries
+ library(dada2)
+ library(phyloseq)
+ library(Biostrings)
+ library(ggplot2)
+ library(randomcoloR)
library(dplyr)

Set working directory 
setwd("/Users/hanlanmcdougall/Desktop/Assignment3")
path <- "/Users/hanlanmcdougall/Desktop/Assignment3"
list.files(path)

Forward fastq filenames have the format: SAMPLENAME_R1_001.fastq.gz and reverse fastq files have the format: SAMPLENAME_R2_001.fastq.gz
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names=TRUE))
fnRs <- sort(list.files(path,pattern="_R2_001.fastq.gz", full.names=TRUE))

Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

Visualize quality of forward reads
plotQualityProfile(fnFs[1:3])

Visualize quality of the reverse reads
plotQualityProfile(fnRs[1:3])

Place the filtered files in the filtered subdirectory and name them 
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

Now want to filter and trim the reads based on the above calls 
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,200),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

Want to see the number of reads that were processed and how many passed 
head(out)

Look at error rates and plot 
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

For each of the samples, want to see how many unique sequences there are 
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

Merge the paired reads 
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

Inspect the merger data.frame from the first sample
head(mergers[[1]])

Construct a sequence table (seqtab)
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

Inspect distribution of sequence lengths 
table(nchar(getSequences(seqtab)))

Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

Not a necessary step, but helpful to track reads through the pipeline😀
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

Final look to see how many sequences have been processed and have made it to the final step 
head(track)

Assign taxonomy to 16S database 
taxa <- assignTaxonomy(seqtab.nochim, "/Users/hanlanmcdougall/Desktop/Assignment3/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

Save csv of the taxa 
write.csv (taxa, file= "/Users/hanlanmcdougall/Desktop/Assignment3/Taxa.csv")
write.csv (seqtab.nochim, file= "/Users/hanlanmcdougall/Desktop/Assignment3/Nochim.csv")

Make objects of the files for easier access 
taxa <- read.csv(file="/Users/hanlanmcdougall/Desktop/Assignment3/Taxa.csv", sep=',', row.names=1)
seqtab.nochim <- read.csv(file="/Users/hanlanmcdougall/Desktop/Assignment3/Nochim.csv", sep=',', row.names=1)

Convert to matrices so that Phyloseq can be used 
seqtab.nochim<-as.matrix(seqtab.nochim)
taxa<-as.matrix(taxa)

Transpose no.chim file so all data is facing the same way so ASV abundance and taxa info can be combined
flipped_seqtab.nochim<- as.data.frame(t(seqtab.nochim))
dim(flipped_seqtab.nochim)

Merge the files 
Combine flipped seqtab.nochim with the taxa file and name as object ASVabund using the cbind command 
ASVabund <- cbind(flipped_seqtab.nochim, taxa)

use the write.csv command to save the file 
write.csv(ASVabund, file="/Users/hanlanmcdougall/Desktop/Assignment3/ASVabund.csv")

samples.out <- rownames(seqtab.nochim) #naming the samples as the row names
samdf <- data.frame(samples.out) #making the sample names a data frame
rownames(samdf) <- samples.out #taking samples out and making it a function with names

Hand off data to phyloseq 
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(samdf), tax_table(taxa))

use biostrings function to use the word ASV and # instead of seeing the full sequence - makes it easier computing wise
dna <-Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <-merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0 ("ASV", seq(ntaxa(ps)))

To do a graph with only ggplot, have to take the phyloseq data into a dataframe using psmelt function
ps.table<-psmelt(ps)

now want to make a factor of the phyla column so that we can graph it 
ps.table$Phylum <- factor(ps.table$Phylum)

**Graphing**

Want colours to make it pretty! 
palette.phyla <- distinctColorPalette(32) #32 phyla
palette.order <- distinctColorPalette(127) #127 orders

Make a graph of relative abundance of phyla in the two samples and negative control 
Barplot.phyla <- ggplot (data=ps.table, mapping=aes (x=Sample, y=Abundance)) + 
  geom_bar(aes(fill=Phylum), stat="identity", position="fill") +
  scale_fill_manual(values=c(palette.phyla)) + labs(y="Relative abundance", title="Relative Abundance of Phyla in Canadian Arctic permafrost samples") +  scale_x_discrete(labels=c("Hummock", "Negative Control", "Trough")) + theme(legend.key.size = unit(0.4, 'cm'), legend.text = element_text(size=5)) + theme(legend.text = element_text (face="italic")) 

#Now want to make a graph of relative abundance for order 
Barplot.Order <- ggplot (data=ps.table, mapping=aes (x=Sample, y=Abundance)) + 
  geom_bar(aes(fill=Order), stat="identity", position="fill") + scale_fill_manual(values=c(palette.order)) + labs(y="Relative Abundance", title="Relative Abundance of Order in Canadian Arctic permafrost samples") + scale_x_discrete(labels=c("Hummock", "Negative Control", "Trough")) + theme(legend.key.size = unit(0.4, 'cm'), legend.text = element_text(size=4)) + theme(legend.text = element_text (face="italic"))
