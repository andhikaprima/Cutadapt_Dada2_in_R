# https://neof-workshops.github.io/Metabarcoding_6xxzqz/Course/01-Metabarcoding.html

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("dada2", version = "3.16")

library(dada2)
library(Biostrings)
library(ShortRead)

# Folder management
input.path <- file.path("1_demux")
if(!dir.exists(input.path)) dir.create(input.path)

output.path <- file.path("2_output")
if(!dir.exists(output.path)) dir.create(output.path)

### Loading data into R and DADA2
## Saving the paths for our input and output directories
# We are now nearly ready to input our data but first we need to set the main paths we will be using which contain the input files (the raw sequencing data files) and where we want to save the output files we generate (the ‘Metabarcoding’ directory we made above).
raw.path <- "0_rawsequences"
input.path <- "1_demux"
output.path <- "2_output"


## Demultiplexing with cutadapter
# Specify the path to the cutadapt software so R knows where to find and run it.
cutadapt <- "/Users/andhikaprima/Library/Python/3.7/bin/cutadapt"


# We can use the system2 command to run a shell command from within R.
# Specify the path to the output files where we want to put the 1_demux output files.
system2(cutadapt, args = "--help")

system2(cutadapt, args = c("-e 0.15", "--no-indels",
                           "-g file:0_rawsequences/Lib1_fwd.fasta",
                           "-G file:0_rawsequences/Lib1_rev.fasta",
                           "-o 1_demux/{name1}-{name2}_R1.fastq",
                           "-p 1_demux/{name1}-{name2}_R2.fastq",
                           "0_rawsequences/Lib1_R1.fastq",
                           "0_rawsequences/Lib1_R2.fastq"))

# get all the names of .fastq files
files <- list.files(input.path, pattern = "\\.fastq", full.names = TRUE)
# extract the basename
filenms <- basename(files)
#create a logical vector where it matches the prefix, suffix from `-`
i1 <- sub("(\\w+)-.*", "\\1",filenms) == 
  sub(".*-(\\w+)_.*", "\\1", filenms)
# subset the files to keep
files_to_keep <- filenms[i1]
# remove all the rest of the files
file.remove(files[!i1])
new_files <- file.path(input.path, sub(".*-", "", files_to_keep))
# rename those files
file.rename(from = files[i1], to = new_files)

# Remove unknown sequences
Unknown <- list.files(input.path, pattern = "\\unknown_", full.names = TRUE)
file.remove(Unknown)

# You should see a list of paired fastq sequencing files listed. Forward and reverse fastq filenames have format: SAMPLENAME_L001_R1_001.fastq and SAMPLENAME_L001_R2_001.fastq. This shows that we have set our input path correctly. There is also a directory called ‘taxonomy’ and a csv file of sample information, both of which we will use later.

## Inputting the forward and reverse reads
# We will assign the input path and specify all files which end in _L001_R1_001.fastq as our forward reads (stored as the variable fnFs) and _L001_R2_001.fastq as our reverse reads (stored as the variable fnRs).
fnFs <- sort(list.files(input.path, pattern = "_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(input.path, pattern = "_R2.fastq", full.names = TRUE))

### Identifying and removing primers
# Now you will carry out primer removal. DADA2 requires the primers you used to be trimmed off the forward and reverse of your reads. We will use the software cutadapt for this from within R.

# Read in your forward and reverse primer sequences for our MiFish-U primerset. The MiFish-U primers do not contain any degenerate bases but cutadapt is able to handle any degenerate bases you may have in primer sequences for your own projects.
FWD <- "CCNGAYATRGCNTTYCCNCG"
REV <- "TANACYTCNGGRTGNCCRAARAAYCA"

## Primer orientation checking
# It might be useful to check the orientation of the primers on our sequences so we know the best way to remove them. We know how our libraries were created and sequenced so we are expecting that the forward primer should be located in a forward orientation in the forward read and the reverse primer should be located in the forward orientation in the reverse read.

# For your own data you might be less certain as to which orientation your primers are in your reads.

# We use the following code to make a list of all orientations of our primer sequences (forward, complement, reverse and reverse completment).
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString 
                                  #objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

FWD.orients <- allOrients(FWD)
FWD.orients # print all orientations of the forward primer to the console
REV.orients <- allOrients(REV)
REV.orients # print all orientations of the reverse primer to the console

# It is a good idea to pre-filter your data before searching for the primer sequences to remove Ns (ambiguously called bases) from the dataset and improve mapping of the primers. To do this we can use the filterAndTrim function. For more information on filterAndTrim and the parameters you can specify type into the console.
?filterAndTrim
# (Press ‘q’ to exit the help and return to the R console).

# For now we are only filtering Ns but we will come back to this function later after removing our primers. Specify the path to a new output directory to contain the N filtered files. We will call this directory filtN, and a copy of each of the forward and reverse fastq files with any ‘N’ containing sequences removed will be written here.
fnFs.filtN <- file.path(output.path, "1_filtN", basename(fnFs)) 
fnRs.filtN <- file.path(output.path, "1_filtN", basename(fnRs))

# Then run filterAndTrim with maxN = 0 (sequences with more than 0 Ns will be discarded) and multithread=TRUE (all input files are filtered in parallel, which is faster).
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

# For this function we have specified:
# - the path to our raw forward read fatsqs fnFs
# - the path to write our N-trimmed forward fastqs fnFs.filtN
# - the path to our raw reverse read fatsqs fnRs
# - the path to write our N-trimmed reverse fastqs fnRs.filtN
# - and the function-specific options for maxN and multithread

# This will take a few minutes to run. When it is finished the cursor will stop flashing and you will see the > prompt and can carry on with the next command.

# The following function will use all possible primer combinations to count the number of times a primer is found in the forward and reverse read in each orientation.
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

## Cutadapt
# Now we know the orientation of our primers we can use the software cutadapt to remove the primers from our samples.

# See https://cutadapt.readthedocs.io/en/stable/ for more information on how to install and run cutadapt.

# Specify the path to the cutadapt software so R knows where to find and run it.
cutadapt <- "/Users/andhikaprima/Library/Python/3.7/bin/cutadapt"

# We can use the system2 command to run a shell command from within R.
# Specify the path to the output files where we want to put the cutadapt output files.
path.cut <- file.path(output.path, "2_PrimerRemoval")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

# Here we specify the options needed by cutadapt in order to trim the forward orientation of the forward and reverse primer off the forward and reverse read and the reverse complement off the end of the forward and reverse reads.
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 

# What do the -g -a -G and -A parameters mean? (Hint: use the following command to check the helpfile for cutadapt)
system2(cutadapt, args = "--help")

# In addition to trimming the primers from the reads we will also specify a couple of extra useful parameters.
# --discard-untrimmed this tells cutadapt to discard any read where the primers haven’t been trimmed off. This is especially important for our data as our files at the moment contain sequences amplified using both MiFish-U and 12S-V5 primer sets. We only want to keep sequences matching the MiFish-U primer set for this analysis.
# --minimum-length 60 discard reads shorter than 60bp. This will remove unexpected short reads and help speed up further analysis.

# Run cutadapt. (Note. The cutadapt step is time intensive so this might take a little while to run, probably about 15 minutes, and a lot of output will be written to the screen while it is running).
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, 
                             # -n 2 required to remove FWD and REV
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i], # input files
                             "--discard-untrimmed",
                             "--minimum-length 250"))
}

# We can now check whether all the primers have been removed using the primerHits function we specified earlier.
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# We now have no primers remaining in our file.

### Checking the quality of your data
# Those of you who attended our Introduction to sequencing data and quality control course would have used FastQC to check the quality of the data. DADA2 has its own quality control option, which plots a similar read length by quality figure.

# To run first import the cutadapt files and extract the sample names and then run the plotQualityProfile function

# Specify the paths and file names the forward and reverse primer cleaned files 
cutFs <- sort(list.files(path.cut, pattern = "_R1.fastq", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2.fastq", full.names = TRUE))

# Extract sample names
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

# check the quality for the first file
plotQualityProfile(cutFs[1:1])

# To interpret this plot, the gray-scale heatmap shows the the frequency of each quality score along the forward read lengths. The green line is the median quality score and the orange lines are the quartiles. The red line at the bottom of the plot represents the proportion of reads of that particular length.
# The quality is very good for our forward reads. You can also see that the majority of the forward reads are ~130bp long after having the primers removed with cutadapt.

# Now check the quality of the reverse file in the same way
plotQualityProfile(cutRs[1:1])

# The reverse reads also look like they are good quality.
# To check the quality of the second and third fastq files we would type.
plotQualityProfile(cutFs[2:3])
plotQualityProfile(cutRs[2:3])

### Cleaning your data
# We will now filter our data to remove any poor quality reads.
# First set the path to a directory to store the filtered output files called filtered.

filtFs <- file.path(path.cut, "../3_filtered", basename(cutFs))
filtRs <- file.path(path.cut, "../3_filtered", basename(cutRs))

# Now run filterAndTrim. This time we use the standard filtering parameters:
# maxN=0 After truncation, sequences with more than 0 Ns will be discarded. (DADA2 requires sequences contain no Ns)
# truncQ = 2 Truncate reads at the first instance of a quality score less than or equal to 2
# rm.phix = TRUE Discard reads that match against the phiX genome
# maxEE=c(2, 2) After truncation, reads with higher than 2 “expected errors” will be discarded
# minLen = 60 Remove reads with length less than 60 (note these should have already been removed by cutadapt)
# multithread = TRUE input files are filtered in parallel

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), 
                     truncQ = 2, minLen = 250, rm.phix = TRUE, compress = TRUE, 
                     multithread = TRUE)
out

# Some samples have very low read numbers after this filtering step. These could be poor quality samples but we also have negatives controls in this dataset so we would expect these to contain zero or very few reads.

### Identification of ASVs
## Generate an error model

# First we need to model the error rates of our dataset using both the forward and reverse reads. Each dataset will have a specific error-signature with errors introduced by PCR amplification and sequencing.
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

# We can use the plotErrors function to check the estimated error rates.
plotErrors(errF, nominalQ = TRUE)

# You will see the Warning messages: 1: Transformation introduced infinite values in continuous y-axis 2: Transformation introduced infinite values in continuous y-axis
# Do not worry about these warning messages. This is a message from the plotting function to let you know that there were some zero values in the data plotted (which turn into infinities on the log-scale). This is expected, it results from the fact that not every combination of error type (e.g. A->C) and quality score (e.g. 33) is observed in your data which is normal.

# Interpreting the plots.
# The error rates for each possible transition (e.g. A→C, A→G) are shown
# Red line - expected based on the quality score. (These are plotted when nominalQ = TRUE is included in the plot command)
# Black line - estimate
# Black dots - observed

# What we are expecting to see here is that the observed match up with estimates. Here we can see that the black dots track well with the black line. We can also see that the error rates drop with increasing quality score as we would expect. Sanity checks complete, we can proceed with the analysis.
# If you are worried about what your error plots look like when you fit your own dataset, one possible way to to improve the fit is to try increasing the number of bases the function is using (the default is 100 million)

## Dereplication
# The next step is to dereplicate identical reads. This is a common step in many workflows used for processing amplicons. This saves time and processing power as identical reads are collapsed together. For example instead of processing 100 identical sequences, only one is processed but the original number (i.e. 100) is associated with it.

exists <- file.exists(filtFs) 
# check that all the samples are still present after filtering
derepFs <- derepFastq(filtFs[exists], verbose=TRUE)
derepRs <- derepFastq(filtRs[exists], verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names[exists]
names(derepRs) <- sample.names[exists]

## Inferrence of ASVs
# Now we are ready to infer the ASVs in our dataset. To do this DADA2 uses the error models created above to infer the true sample composition (follow this link for more details).
# Here we will run the inference algorithm on single samples to save time but it is also possible to run samples together in pseudo-pools to increase the ability to identify ASVs of low abundance. Low abundance ASVs in a sample may be filtered out when run separately but if they are found in higher numbers in another sample then the chance of that ASV being real increases. Pseudo-pooling increases the chances of these “real” low abundance ASVs being retained within samples. Whether or not to use the pseudo-pooling option will depend on your dataset and experimental design (see https://benjjneb.github.io/dada2/pseudo.html#Pseudo-pooling for more information).

dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

# It is important to note that you want to run the error and inference steps on datasets generated from a single Illumina run. Data generated from different runs can have different error structures.

## Merging paired end reads
# Up until now we have carried out all filtering, error correction and inference on the forward and reverse reads separately. It is now time to merge the two files. By default the minimum overlap allowed between the two samples is 12 bp.

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

## Making our ASV matrix
# Now it is time to make the counts table. Each column represents a single ASV and each row is an individual sample.

seqtab <- makeSequenceTable(mergers)
dim(seqtab) 

# Q. How many ASVs are in our matrix?

## Chimera detection and removal
# The last step in generating our ASV matrix is to detect and remove any chimeric sequences.
# Chimeric sequences are formed when two or more biological sequences are joined together. This is fairly common in amplicon sequencing. DADA2 uses a method whereby it combines the left and right segments of abundant reads and compares these with lower abundant sequences. Any low abundant sequences that match are removed.

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", 
                                    multithread=TRUE, verbose=TRUE)

dim(seqtab.nochim)

# create file
write.csv(seqtab.nochim, "2_output/1_seqtab.nochim.csv", row.names=TRUE)

# Q. How many ASVs remain after filtering out chimeras?
# Although this is a large proportion of sequence variants it should be a smaller proportion of the total sequences.
sum(seqtab.nochim)/sum(seqtab)

# In these data ~ 17% of the merged sequence reads were identified as chimeric.
# Check the range of ASV lengths:
table(nchar(getSequences(seqtab.nochim)))

# Sequences that are much longer or shorter than the expected amplicon length could be the result of non-specific priming, and can be removed from your sequence table (e.g. seqtab.len <- seqtab[,nchar(colnames(seqtab)) %in% 150:180]). This is analogous to “cutting a band” in-silico to get amplicons of the targeted length, e.g. the above command would keep the amplicons between 150 bp and 180 bp and remove all other lengths of sequences. This is dependent on the amplicon and you should only remove these if you do not expect any length variation, or after further investigation of what these short or long sequences are. For this example we will proceed without any amplicon length filtering.

## Sequence tracking sanity check
# The last thing to do in this section is to track the number of sequences through the pipeline to check whether everything has run as expected and whether there are any steps where we loose a disproportionate number of sequences. If we end up with too few reads to run further analysis we can use this table to identify any step which might require further investigation and optimisation.

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), 
               sapply(dadaRs, getN), 
               sapply(mergers, getN), 
               rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", 
                     "denoisedF", "denoisedR", 
                     "merged", "nonchim")
rownames(track) <- sample.names
track

### Assigning taxonomy
# To assign taxonomy we will use a custom reference database containing fish sequences available for the MiFish-U amplicon region for species found in Lake Tanganyika and its broader catchment area.

# For more information on how to generate your own custom reference database please see the supplementary materials on day 2.

taxa <- assignTaxonomy(seqtab.nochim, 
                       "DB_simple/references_CCCPbarcodes_ap.fasta", 
                       multithread=T, verbose = T)

taxa.print <- taxa
#rownames(taxa.print) <- NULL
head(taxa.print)
taxa.print

# create file
write.csv(taxa.print, "2_output/2_taxa.print.csv", row.names=TRUE)

# This displays the taxonomic assignment of the first six ASVs at taxonomic levels from kingdom to species. The second ASV has not been assigned to any taxon in the reference database. Of the other five, two have been assigned at the species level, two at the genus level and one at the family level.


# ## Merge dataset
# library(data.table)
# library(tibble)
# 
# # transpose reads data
# seqtab.nochim_df <- data.frame(seqtab.nochim, check.names = F)
# seqtab.nochim_df <- rownames_to_column(seqtab.nochim_df, var = "sequence")
# sequence <- colnames(seqtab.nochim_df)
# seqtab.nochim_df = rbind(sequence,seqtab.nochim_df)
# seqtab.nochim_t <- transpose(seqtab.nochim_df) 
# names(seqtab.nochim_t) <- lapply(seqtab.nochim_t[1, ], as.character)
# seqtab.nochim_t <- seqtab.nochim_t[-1,] 
# 
# # preparing taxa data
# taxa.print_df <- data.frame(taxa.print, check.names = T)
# taxa.print_df <- rownames_to_column(taxa.print_df, var = "sequence")
# 
# # Merger
# FullTable <- merge(taxa.print_df, seqtab.nochim_t, by = "sequence", all = T)
# rownames(FullTable) <- sprintf("OTU%06d", 1:nrow(FullTable))
# FullTable <- rownames_to_column(FullTable, var = "OTU")
# 
# # create file
# write.csv(FullTable, "2_output/3_FullTable.csv", row.names=F)



### Exporting data
# It is useful to be able to export our ASV sequences and a matrix table of counts from DADA2/R as we might want to visualise or analyse these using other software.

# Our ASV sequences and counts per sample are stored in the object seqtab.nochim. The ASVs are not named, so first let’s name them (ASV_1, ASV_2, etc.).

# The column names of seqtab.nochim are actually the ASV sequences, 
# so extract these and assign them to `mifish_seqs`
ASV_seqs <- colnames(seqtab.nochim)

# Make a new variable for ASV names, `mifish_headers`, 
#with length equal to the number of ASVs
ASV_headers <- vector(dim(seqtab.nochim)[2], mode="character")

# Fill the vector with names formatted for a fasta header (>ASV_1, >ASV_2, etc.)
for (i in 1:dim(seqtab.nochim)[2]) {
  ASV_headers[i] <- paste(">ASV", i, sep="_")
}

## Fasta file
# Now we have our sequences and names as variables we can join them and make a fasta file.

ASV_fasta <- c(rbind(ASV_headers, ASV_seqs))
write(ASV_fasta, "2_output/4_ASVs.fa")

# You should now have this fasta file in your working directory on the server.

## Sequence count matrix
# Next make a table of sequence counts for each sample and ASV.
# First transpose the `seqtab.nochim` and assign this to the variable `ASVtab`
ASV_tab <- t(seqtab.nochim)

# Name each row with the ASV name, omitting the '>' used in the fasta file
row.names(ASV_tab) <- sub(">", "", ASV_headers)
write.table(ASV_tab, "2_output/5_ASV_counts.tsv", sep="\t", quote=F, col.names=NA)

## Table of taxon names
# Lastly, if we’ve used dada2 to assign taxonomy we can make a table of taxon names for each ASV.

# Replace the row names in `taxa` with the ASV names, 
# omitting the '>' used for the fasta file.
rownames(taxa) <- gsub(pattern=">", replacement="", x=ASV_headers)
write.table(taxa, "2_output/6_ASV_taxonomy.tsv", sep = "\t", quote=F, col.names=NA)

# You should now have a tsv file of taxonomic assignments for each ASV in your working directory.

#merger
library(data.table)
library(tibble)

OTU_table <-cbind(taxa, ASV_tab, ASV_seqs)
OTU_table <- data.frame(OTU_table)
OTU_table <- tibble::rownames_to_column(OTU_table, "ASV")
write.table(OTU_table, "2_output/7_ASV_table.tsv", sep = "\t", quote=F, col.names=NA)


### Further analysis
# We will use three more handy R packages in order to further explore our data, phyloseq and vegan to explore the diverity in our data, plus the visualisation package ggplot2 for graphics. Load these libraries before proceeding with the steps below.

library(phyloseq)
library(vegan)
library(ggplot2)

# load metadata
meta <- read.table("0_rawsequences/sample_details_sheet.csv", row.names = 1)

## Rarefaction curves
# The more deeply we sequence a sample the more species we discover, until this accumulation levels off as we have found all or as many ASVs/species as we are going to find. This becomes a problem as samples are often sequenced at different depths so will have reached different points in the curve. This is a common difficulty as pooling and sequencing equal amounts of each sample can be tricky.

# One way to visualise this is to plot a rarefaction curve for each sample.
rarecurve(seqtab.nochim, step=100, col= meta$Col, lwd=2, ylab="ASVs", label=F)
# add a vertical line to represent the fewest sequences in any sample
abline(v=(min(rowSums(seqtab.nochim))))

## Alpha diversity
# Measures of alpha diversity are used to describe diversity within a sample.
# 
# We will use the R package phyloseq to plot alpha diversity. For this example we will proceed with unnormalised data. It is generally recommended not to normalise count data before calculating alpha diversity measures in the phyloseq FAQ.
# 
# First make a phyloseq object. To do this we first read in our ASV, taxonomy and metadata tables before making the plyloseq object phylo.

seqtab.nochim2<-t(as.data.frame(seqtab.nochim))
phylo_asv <- otu_table(seqtab.nochim2, taxa_are_rows=TRUE)
phylo_tax <- tax_table(taxa)
phylo_samples <- sample_data(meta)

phylo <- phyloseq(phylo_asv, phylo_tax, phylo_samples)

sample_names(phylo)
rank_names(phylo)
sample_variables(phylo) 

# Plot the two methods of calculating alpha diversity.
plot_richness(phylo, 
              measures=c("Shannon", "Simpson"), 
              color = "Position")
plot_richness(phylo, x="Position", measures=c("Shannon", "Simpson"), 
              color = "Position") + geom_boxplot()

## Beta diversity
# Beta diversity compares the difference in diversity between two sites, or to put it another way it calculates the number of species that are not the same in the two sites.

# We will normalise the data before running the beta diversity calculation. We will transform the data into proportions to be used for Bray-Curtis distances.

ps.prop <- transform_sample_counts(phylo, function(otu) otu/sum(otu))

# Ploting
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color="Position", title="Bray NMDS")

# We will then calculate the Bray–Curtis distances using the distance function and perform a PERMANOVA (permutational multivariate analysis of variance) using the adonis function from Vegan to check whether the separation of samples by site is significantly different.
bray.dist<-distance(ps.prop, method="bray")
sampledf <- data.frame(sample_data(phylo))
adonis(bray.dist ~ Position, data = sampledf)

# The PERMANOVA results suggest that there is a statistical difference in communities between sites.
# Lastly, let’s plot the proportion of ASV sequences within each sample that belong to different taxonomic families.
plot_bar(ps.prop, fill = "Family")+
  geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack")+
  facet_grid(~Position, scales = "free", space = "free")

save.image("2_output/Neof_dada2.RData")