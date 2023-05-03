Final\_class\_project
================
Christopher Ebuka Ojukwu
4/23/2023


## Consensus Peaks Creation
``` r
# filepath to import peaks
broadpeakfilepath <- "/scratch/Shares/rinnclass/CLASS_2023/data/data/peaks"

#import peaks
#peak_list <- import_peaks(consensus_file_path = broadpeakfilepath)

# Creating unique DBP object for create_consensus_peaks_from_reduced
dbp <- unique(sapply(names(peak_list), function(x) {
  unlist(strsplit(x, "_"))[1]
}))

# now run our function consensus_from_reduced
consensus_list <- lapply(dbp, consensus_from_reduced, peak_list)

# adding names to the GRange list
names(consensus_list) <- dbp
```

## Finding out the number of peaks
``` r
# creating list of num_peaks per dbp
num_peaks <- sapply(consensus_list, length)

# plotting
hist(num_peaks, breaks = 1000)
hist(num_peaks, breaks = 1000, xlim = c(0,3000))

ggsave("/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/final_project/analysis/01_create_consensus_peaks/figures/hist_num_peaks.pdf")
```
<img src="figures/numpeaks_1.png" width="1888" />

## Looking at just the first 3000
<img src="figures/numpeaks_2.png" width="1888" />

##Filtering the consensus list 
``` r
# Filter out any chip data less 1,000 peaks == filtered consensus peaks
filtered_consensus_list <- consensus_list[sapply(consensus_list, length) > 1000]

# saving 
save(filtered_consensus_list, file = "results/filtered_consensus_list.RData")

# keeping track of DBPs lost
lost_dbps <- names(consensus_list[sapply(consensus_list, length) < 1000]) %>% as.data.frame()

# saving 
write.table(lost_dbps, "results/lost_dbps.csv")
```

``` r
# How does peak number and genome coverage compare #what does this mean?
for(i in 1:length(filtered_consensus_list)) {
  rtracklayer::export(filtered_consensus_list[[i]], 
                      paste0("results/filtered_consensus_peaks/", 
                            names(filtered_consensus_list)[i], 
                             "_filtered_consensus_peaks.bed"))
}

for(i in 1:length(consensus_list)) {
rtracklayer::export(consensus_list[[i]], paste0("results/consensus_peaks/", names(consensus_list)[i],"_consensus_peaks.bed")) }
```

##Exporting the consensus list and filtered consensus list
``` r
consensus_file_list <- list.files("/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/final_project/analysis/01_create_consensus_peaks/results/consensus_peaks", full.names = T, pattern = ".bed")

# lapply (for loop) across consensus file list to add colnames
# The actual col names for .broadPeak are: chr, start, end, name, score, strand
peaks <- lapply(consensus_file_list, read.table, col.names = c("chr", "start", "end", "name", "score", "strand"))
names(peaks) <- dbps

# double check order by looking at consensus_file_list is same order as dbps
# make chromosomes of interest object
canonical_chr <- c(paste0("chr", 1:22), "chrM", "chrX", "chrY")

# let's use lapply with filter funciton to cannonical_chr
peaks <- lapply(peaks, function(x) x %>% filter(chr %in% canonical_chr))

# now that these are all nice and clean let's export:
new_filenames <- paste0("consensus/", names(peaks), "_consensus.bed")

for(i in 1:length(peaks)) {
  write.table(peaks[[i]], new_filenames[[i]],
              sep = "\t", col.names = FALSE, row.names = FALSE,
              quote = FALSE, append = TRUE)
}
```

##Loading and exploring genome features
``` r
# What is the distribution of promoter overlaps versus gene-bodies (hint hist)
gencode_gr <- rtracklayer::import("/scratch/Shares/rinnclass/CLASS_2023/data/da#ta/genomes/gencode.v32.annotation.gtf")

# gencode genes
gencode_genes <- gencode_gr[gencode_gr$type == "gene"] 

# mrna_genes
mrna_genes <- gencode_genes[gencode_genes$gene_type %in% "protein_coding"]

# lncrna_genes
lncrna_genes <- gencode_genes[gencode_genes$gene_type %in% "lncRNA"] 

# mrna_lncrna_genes
mrna_lncrna_genes <- gencode_genes[gencode_genes$gene_type %in% c("protein_coding","lncRNA")]

# lncrna_mrna_promoters
lncrna_mrna_promoters <- promoters(mrna_lncrna_genes, upstream = 1000, downstream = 1000)

# lncrna_gene_ids
lncrna_gene_ids <- mrna_lncrna_genes$gene_id[mrna_lncrna_genes$gene_type == "lncRNA"]

# mrna_gene_ids
mrna_gene_ids <-mrna_lncrna_genes$gene_id[mrna_lncrna_genes$gene_type == "protein_coding"]
```

## Tracking filtered consensus peaks
``` r
#numpeaks to track peak properties
num_peaks_df <- data.frame("dbp" = names(filtered_consensus_list),
                           #"num_peaks" = sapply(filtered_consensus_list, #length))

# total genome covered by peaks
num_peaks_df$total_peak_length <- sapply(filtered_consensus_list, function(x) sum(width(x)))

# creating number of promoter overlaps entry
promoter_peak_counts <- count_peaks_per_feature(lncrna_mrna_promoters, filtered_consensus_list, type = "counts")

# creating promoter peak_occurence for clustering - Metaplots later.
promoter_peak_matrix <- count_peaks_per_feature(lncrna_mrna_promoters, filtered_consensus_list, type = "occurrence")

# saving
write.table(promoter_peak_matrix, "results/promoter_peak_occurrence_matrix.tsv")

# read back in
promoter_peak_occurrence_matrix <- read.table("results/promoter_peak_occurrence_matrix.tsv")

#summing rows to get total number of promoter overlaps
num_peaks_df$peaks_overlapping_promoters <- rowSums(promoter_peak_counts)

# lncrna promoter overlaps 
num_peaks_df$peaks_overlapping_lncrna_promoters <- rowSums(promoter_peak_counts[,lncrna_gene_ids])

# mrna promoter overlaps
num_peaks_df$peaks_overlapping_mrna_promoters <- rowSums(promoter_peak_counts[,mrna_gene_ids])

# Finding overlaps with gene_bodies (will take a few minutes again)
# Note this takes several minutes
genebody_peak_counts <- count_peaks_per_feature(mrna_lncrna_genes, 
                                                filtered_consensus_list, 
                                                type = "counts")

# All gene bodies overlaps
num_peaks_df$peaks_overlapping_genebody <- rowSums(genebody_peak_counts)

# lncRNA gene bodies 
num_peaks_df$peaks_overlapping_lncrna_genebody <- rowSums(genebody_peak_counts[,lncrna_gene_ids])

# mRNA gene bodies
num_peaks_df$peaks_overlapping_mrna_genebody <- rowSums(genebody_peak_counts[,mrna_gene_ids])

#plotting
hist(promoter_peak_counts, genebody_peak_counts)
plot(genebody_peak_counts, main="Promoter Peak Counts vs. Genebody Peak Counts",
    xlab="Promoter Peak Counts", ylab="Genebody Peak Counts", pch=20, col="blue")
# Add a grid
grid()
# Add a legend
legend("topright", legend="Data", pch=20, col="blue")
```

## Finding out proteins that are TFs and most recurring dbd
``` r
# download TF annotations to results
url <- "https://www.cell.com/cms/10.1016/j.cell.2018.01.029/attachment/ede37821#-fd6f-41b7-9a0e-9d5410855ae6/mmc2.xlsx"
destination_for_url <- "results/TF_annotations.xlsx"

# to download we can use download.file
download.file(url, destination_for_url)

# reading in TF annotations 
human_tfs <- readxl::read_excel("/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/final_project/analysis/01_create_consensus_peaks/results/TF_annotations.xlsx",
                                sheet = 2, skip = 1)

# let's rename the 4th column to indicate if it is a TF.
names(human_tfs)[4] <- "is_tf"
names(human_tfs)[3] <- "dbp"

# now let's intersect gene names that are in our ChIP data and has TF identity.
length(which(tolower(num_peaks_df$dbp) %in% tolower(human_tfs$dbp)))
 #407 of the 430 have matching gene_names - not bad

#human_tfs <- human_tfs[tolower(human_tfs$Name) %in% tolower(num_peaks_df$dbp), 1:4]
# adding new column names
names(human_tfs) <- c("ensembl_id",
                      #"dbp",
                      #"dbd",
                      #"tf")

# merging into num_peaks_df
num_peaks_df <- merge(num_peaks_df, human_tfs, all.x = T)
num_yes <- sum(num_peaks_df$tf == "Yes", na.rm = TRUE)
# Print the result
cat("Number of 'yes' occurrences:", num_yes) #344

# Extract the column as a character vector
column <- as.character(num_peaks_df$dbd)

# Use table() function to count the occurrences of each word in the column
dbd_counts <- table(column)

# Find the most recurring word
most_recurring_dbd <- names(dbd_counts [dbd_counts  == max(dbd_counts )])

# Print the result
cat("Most recurring dbd:", most_recurring_dbd) #C2H2 ZF
```

    ## Most recurring dbd: C2H2 ZF
    
 ## Creating promoter peak occurrence matrix
``` r
# running count_peaks_per_feature
promoter_peak_occurence <- count_peaks_per_feature(lncrna_mrna_promoters, filtered_consensus_list, 
                                               type = "occurrence")

# Let's double check that all lncrna & mrna genes are accounted for:
stopifnot(all(colnames(promoter_peak_occurence) == lncrna_mrna_promoters$gene_id))

# saving
write.table(promoter_peak_occurence, "/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/final_project/analysis/01_create_consensus_peaks/results/lncrna_mrna_promoter_peak_occurence_matrix.tsv")

# Now let's use the 'data.frame()' fucntion. Set up a bunch of colnames and populate them.
peak_occurence_df <- data.frame("gene_id" = colnames(promoter_peak_occurence),
                                "gene_name" = lncrna_mrna_promoters$gene_name,
                                "gene_type" = lncrna_mrna_promoters$gene_type,
                                "chr" = lncrna_mrna_promoters@seqnames,   
                                "1kb_up_tss_start" = #lncrna_mrna_promoters@ranges@start,
                                "strand" = lncrna_mrna_promoters@strand,
                                "number_of_dbp" = #colSums(promoter_peak_occurence))

# saving
write_csv(peak_occurence_df, "/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/final_project/analysis/01_create_consensus_peaks/results/peak_occurence_dataframe.csv")
```

## Finding out which dbps are bound
``` r
# dbps on promoters object
DBPs_on_promoter <- lncrna_mrna_promoters %>%
                    as.data.frame() %>%
  dplyr::select(gene_id, gene_name)

# creating promoter dbps by pivot longer of promoter_peak_occurence_matrix
promoter_dbps <- promoter_peak_occurence %>%
  as.data.frame() %>%
  rownames_to_column("dbp") %>%
pivot_longer(2:ncol(.), names_to = "gene_id", values_to = "occurrence") %>%
  filter(occurrence == 1) %>%
  dplyr::select(-occurrence) %>%
  left_join(DBPs_on_promoter)
  
# checking Firre promoter
firre_promoter <- promoter_dbps %>%
  filter(gene_name == "FIRRE")

# XIST promoter (should be off since male cells)
XIST_promoter <- promoter_dbps %>%
  filter(gene_name == "XIST")

# GAPDH
GAPDH_promoter <- promoter_dbps %>%
  filter(gene_name == "GAPDH")

# saving
#promoter_dbps_df <- promoter_dbps %>% as.data.frame()

#write.csv(promoter_dbps, "/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/final_project/analysis/01_create_consensus_peaks/results/promoter_dbps.csv")

#to look at peaks per dbp
ggplot(num_peaks_df, aes(x = num_peaks)) + 
 geom_histogram(bins = 70)
```

## Gene similarities via clustering
``` r
#let's look at the promoter peak occurence matrix (takes a min to load)
# read.table will make this a data.frame but will keep row names.
promoter_peak_occurence_matrix <- read.table("/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/final_project/analysis/01_create_consensus_peaks/results/lncrna_mrna_promoter_peak_occurence_matrix.tsv")

# Converting to a matrix format for correlation analysis
promoter_peak_occurence_matrix <- as.matrix(promoter_peak_occurence_matrix)

# Now let's create distance matrix for each pairwise vector comparison for the whole matrix now.
peak_occurence_dist <- dist(promoter_peak_occurence_matrix, method = "binary")
peak_occurence_dist
```

    ##                       ADNP       AFF4      AHDC1        AHR      AKAP8
    ## AFF4            0.84962956                                            
    ## AHDC1           0.86369545 0.67243894                                 
    ## AHR             0.86263072 0.28007320 0.72729204                      
    ## AKAP8           0.86919299 0.33174052 0.74377765 0.19031814           
    ## ARID3A          0.86975556 0.70724817 0.75399604 0.74109197 0.76502909
    ## ARID4A          0.86744920 0.30931887 0.74102826 0.18317529 0.17068923
    .
    ##                     ARID3A     ARID4A     ARID4B     ARID5B      ARNTL                                                       
    ## ARID4A          0.76001870                                            
    ## ARID4B          0.74706889 0.14889819                                 
    ## ARID5B          0.75483716 0.19233234 0.10701280                      
    ## ARNTL           0.73026170 0.17870848 0.18059777 0.22155653           
    ## ASH2L           0.72578717 0.62847294 0.65121780 0.67130634 0.59765684
    ## ATAD3A          0.78115249 0.24654088 0.27286223 0.29528944 0.26983506
    .
    ##                      ASH2L     ATAD3A       ATF2       ATF3       ATF5                                                            
    ## ATAD3A          0.63936668                                            
    ## ATF2            0.60611563 0.27392823                                 
    ## ATF3            0.79785331 0.82735442 0.81148137                      
    ## ATF5            0.56501295 0.38154367 0.33297445 0.79345766           
    ## ATF6            0.55245958 0.32083897 0.20037279 0.78743710 0.29623581
    ## ATF7            0.60854530 0.28887140 0.21384074 0.81083821 0.31512061
    .
    ##                       ATF6       ATF7       BCL3       BCL6    BHLHE40                                                         
    ## ATF7            0.26012314                                            
    ## BCL3            0.28745170 0.23947516                                 
    ## BCL6            0.28155587 0.25443896 0.13345249                      
    ## BHLHE40         0.95470048 0.95711132 0.96330568 0.96337651           
    ## BRCA1           0.89052373 0.90264442 0.91816036 0.91748038 0.96924767
    ## CAMTA2          0.81085298 0.82343247 0.85089007 0.85156848 0.94408369
    .
    ##                      BRCA1     CAMTA2       CBFB       CBX5      CEBPA                                                     
    ## CAMTA2          0.88997257                                            
    ## CBFB            0.90524609 0.83029412                                 
    ## CBX5            0.91787777 0.85189649 0.18279303                      
    ## CEBPA           0.90764115 0.84323235 0.27474656 0.30487805           
    ## CEBPB           0.95157299 0.93615894 0.89751013 0.90539607 0.87324265
    ## CEBPD           0.91026446 0.83861752 0.14397923 0.17511073 0.24208304
    .
    ##                      CEBPB      CEBPD      CEBPG      CERS6       CHD2          
    ## CEBPD           0.88665321                                            
    ## CEBPG           0.88037692 0.18043464                                 
    ## CERS6           0.95763331 0.91800419 0.91782690                      
    ## CHD2            0.89957585 0.52344932 0.52650004 0.89397558           
    ## CREB3           0.91532329 0.69270942 0.69511415 0.77269758 0.67002865
    ## CREM            0.92106021 0.85622602 0.85231764 0.87620690 0.80212932
    .
    ##                      CREB3       CREM       CTCF        DBP      DDIT3
    .            
    ## CREM            0.81888909                                            
    ## CTCF            0.88632686 0.93014371                                 
    ## DBP             0.67776129 0.84461078 0.93544087                      
    ## DDIT3           0.77379830 0.81544202 0.91739050 0.77639553           
    ## DLX6            0.56186745 0.81821156 0.83833996 0.80267938 0.79105646
    ## DMAP1           0.67877362 0.84785340 0.82555932 0.87144122 0.84904266
    .
    ##                       DLX6      DMAP1      DMTF1      DNMT1       DPF2
    ## DMAP1           0.38029807                                            
    ## DMTF1           0.62574278 0.69434407                                 
    ## DNMT1           0.50954882 0.29300643 0.74825132                      
    ## DPF2            0.59769103 0.62319677 0.65915648 0.69515994           
    ## DRAP1           0.38184845 0.11885104 0.69640108 0.30751760 0.62365360
    ## DZIP1           0.46419049 0.45232944 0.63838384 0.52326221 0.64139942
    .
    .
    .

``` r
# Ready to cluster
bin_hier <- hclust(peak_occurence_dist, method = "complete")

# we now can plot the information in bin_hier
plot(hist(bin_hier$height), main="Peak Occurence Distribution",
     xlab="Peak Occurence", ylab="Frequency", pch=20,col="blue")
pdf("figures/dbp_hclust_dendro.pdf", height = 24, width = 50)
plot(hist(bin_hier$height), main="Peak Occurence Distribution",
     xlab="Peak Occurence", ylab="Frequency", pch=20,col="blue")
dev.off()
```
<img src="figures/dbp_hclust_dendro-2.jpg" width="1888" />

```r

#Now let's use GGDENDRO (dendrogram) package that will plot the branch lengths #that indicate how similar two samples are 
#ggdendro::ggdendrogram(bin_hier, rotate = FALSE,  size = 3, 
                       #theme_dendro = TRUE) +
   # 90 degree rotation to right
   #coord_flip() +
   #scale_y_continuous() +
   # adds label
   #scale_x_continuous(position = "top") +
   # subsection to labels in order of clustering
   # ? seq_along
   #scale_x_continuous(breaks = seq_along(bin_hier$labels[bin_hier$order]),
                      
                      # adding labels that are in the order 'column'
             #labels = bin_hier$labels[bin_hier$order], position = "top",
             #expand = c(0,0)) +
   #theme(axis.text.x = element_text(angle = 90, hjust  = 1)) + 
   #theme(axis.text.y = element_text(angle = 0,hjust = 1)) +
   #scale_y_reverse(expand = c(0.01, 0)) +
   #theme(
     #plot.background = element_blank(),
     #panel.grid.major = element_blank(),
   #panel.grid.minor = element_blank(),
     #panel.border = element_blank()
   #)
# Nice let's save this using 'ggsave': since we used ggplot
#ggsave("/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/final_project/analysis/03_clustering/figures/ggdendro_plot.pdf", height = 100, width = 22, limitsize = F)
```
<img src="figures/ggdendro_plot.jpg" width="1888" />

## clustering lncrna and mrna separatel
``` r
# loading lncRNA promoters
#lncrna_mrna_promoters <- rtracklayer::import("/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/my_results/my_gene_annotation/lncrna_mrna_promoters.gtf")

# now we can split into lncRNA and mRNA
#lncrna_promoters <- lncrna_mrna_promoters[lncrna_mrna_promoters$gene_type == "lncRNA"]

# mrna promoters
#mrna_promoters <- lncrna_mrna_promoters[lncrna_mrna_promoters$gene_type == "protein_coding"]

# Now we will use indexing to separate peak_occurrence_matrix in lncRNA and mRNA.
# note we are indexing and using indexing in indexing to get what we want.
#lncrna_peak_occurence <- promoter_peak_occurence_matrix[,lncrna_promoters$gene_id]

# we do the clutering the same as above or we can combine the dist and hclust:
#bin_hier_lncrna <- hclust(dist(lncrna_peak_occurence, method = "binary"))

# Now plot with ggdendro
#ggdendro::ggdendrogram(bin_hier_lncrna, rotate = T,  size = 3)
# Now let's save this figure

#ggsave("figures/lncrna_hclust_binary_dist.pdf", height = 49, width = 6)

#getting the gencode_gr again
gencode_gr <- rtracklayer::import("/scratch/Shares/rinnclass/CLASS_2023/data/data/genomes/gencode.v32.annotation.gtf")

# now just type == gene annotations
#gencode_genes <- gencode_gr[gencode_gr$type == "gene"] 
#table(gencode_gr$type)

# exporting all genes file (we will save all the .Rdata too at the end)
#rtracklayer::export(gencode_genes, "/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/final_project/analysis/03_clustering/results/gene_annotation/gencode_genes.gtf")

# mRNA genes (called "protein_coding") in this version of gencode changes sometimes !
#mrna_genes <- gencode_genes[gencode_genes$gene_type %in% "protein_coding"] 
#rtracklayer::export(mrna_genes, "/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/final_project/analysis/03_clustering/results/gene_annotation/mrna_genes.gtf")
#table(gencode_genes$gene_type)

# now doing a second index for lncRNA:
#lncrna_genes <- gencode_genes[gencode_genes$gene_type %in% "lncRNA"] 
#rtracklayer::export(lncrna_genes, "/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/final_project/analysis/03_clustering/results/gene_annotation/lncrna_genes.gtf")

# both mRNA and lncRNA annotations together.
#mrna_lncrna_genes <- gencode_genes[gencode_genes$gene_type %in% c("protein_coding","lncRNA")]
#rtracklayer::export(mrna_lncrna_genes, "/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/final_project/analysis/03_clustering/results/gene_annotation/mrna_lncrna_genes.gtf")

# starting annotation file that we will use moving forward.
#lncrna_mrna_genes <- rtracklayer::import("/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/final_project/analysis/03_clustering/results/gene_annotation/mrna_lncrna_genes.gtf")

# Nice that loaded so much faster -- lets see as DF
#lncrna_mrna_genes_df <- lncrna_mrna_genes %>% as.data.frame()

# creating and exporting promoter annotations
#lncrna_mrna_promoters <- promoters(lncrna_mrna_genes, upstream = 1000, downstream = 1000)

# check right size
#width(lncrna_mrna_promoters)
#rtracklayer::export(lncrna_mrna_promoters, "/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/final_project/analysis/03_clustering/results/gene_annotation/lncrna_mrna_promoters.gtf")

# Clustering of lncRNA and mRNA seperately
#Now let's compare how lncRNAs and mRNAs cluster differnently?
#Strategy: load annotations then index peak_occurrence for specific annotations
# lncRNA promoters ggdendro
#lncrna_mrna_promoters <- rtracklayer::import("/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/final_project/analysis/01_create_consensus_peaks/results/gene_annotations/lncrna_mrna_promoters.gtf")

# now we can split into lncRNA and mRNA
#lncrna_promoters <- lncrna_mrna_promoters[lncrna_mrna_promoters$gene_type == "lncRNA"]

# mrna promoters
#mrna_promoters <- lncrna_mrna_promoters[lncrna_mrna_promoters$gene_type == "protein_coding"]

# Now we will use indexing to separate peak_occurrence_matrix in lncRNA and mRNA.
# note we are indexing and using indexing in indexing to get what we want.
#lncrna_peak_occurence <- promoter_peak_occurence_matrix[,lncrna_promoters$gene_id]

# we do the clutering the same as above or we can combine the dist and hclust:
#bin_hier_lncrna <- hclust(dist(lncrna_peak_occurence, method = "binary"))
# Now plot with ggdendro

#ggdendro::ggdendrogram(bin_hier_lncrna, rotate = T,  size = 3)
# Now let's save this figure
#ggsave("/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/final_project/analysis/03_clustering/figures/lncrna_hclust_binary_dist.pdf", height = 49, width = 6)

# mRNA promter ggdendro
#Now for mRNA
# same strategy used very often:
# indexing into larger file and filtering to a file of stuff you want.
#mrna_peak_occurence <- promoter_peak_occurence_matrix[,mrna_promoters$gene_id]
# getting the distance matrix for only mRNA promoters  
#bin_hier_mrna <- hclust(dist(mrna_peak_occurence, method = "binary"))
# plotting with ggdendro
#ggdendro::ggdendrogram(bin_hier, rotate = TRUE,  size = 3)
#they look alike because the combination of the overlaps are the same for mRNA and lncRNA i.e. similar binding profiles.
# saving
#ggsave("/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/final_project/analysis/03_clustering/figures/mrna_hclust_binary_dist.pdf", height = 44, width = 6)

#save(bin_hier_mrna, mrna_peak_occurence, bin_hier_lncrna, lncrna_peak_occurence, lncrna_promoters, mrna_promoters, lncrna_mrna_promoters, lncrna_mrna_genes, lncrna_mrna_genes_df, mrna_lncrna_genes, lncrna_genes, mrna_genes, gencode_genes, file = "results/clustering_features.RData")
```

## Creating metaplots initial objects
``` r
# loading in needed files from 01_peak_features
#load("/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/final_project/analysis/01_create_consensus_peaks/results/peak_features.RData", verbose = T)

#load("/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/final_project/analysis/03_clustering/results/clustering_features.RData", verbose = T)


consensusPeakPath <- "/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/final_project/analysis/01_create_consensus_peaks/results/consensus"
consensus_peaks_files <- list.files(consensusPeakPath, 
                                             pattern = "*.bed",
                                             full.names = TRUE)

# lapply with import function to make a list of GRanges
consensus_peaks <- lapply(consensus_peaks_files, rtracklayer::import)

#file.path(basepath, peak_path)
# adding DBP name to each list of GRanges from the file name
# note gsub requires no space around or argument ( | )

names(consensus_peaks) <- gsub("/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/final_project/analysis/01_create_consensus_peaks/results/consensus/|_consensus.bed","", consensus_peaks_files)
```

## Metaplot for all dps on lncRNA and mRNA promoters
``` r
# making promoters (again, to be safe)
#lncrna_promoters <- lncrna_mrna_promoters[lncrna_mrna_promoters$gene_type == "lncRNA"]
#mrna_promoters <- lncrna_mrna_promoters[lncrna_mrna_promoters$gene_type == "protein_coding"]
```

## Metaplot for lncRNA
``` r
#setting up lncrna DF.
#lncrna_metaplot_df <- data.frame(x = integer(), dens = numeric(), dbp = character())

# for loop to populate DF with overlap density in lncrna promoters
#for(i in 1:length(consensus_peaks)) {
  #print(names(consensus_peaks)[[i]])
  #tmp_df <- profile_tss(consensus_peaks[[i]], lncrna_mrna_promoters = lncrna_promoters)
  #tmp_df$dbp <- names(consensus_peaks)[[i]]
  #lncrna_metaplot_df <- bind_rows(lncrna_metaplot_df, tmp_df)
#}

#lncrna_metaplot_df$gene_type <- "lncRNA"

#pdf("metaplots/lncrna_promoter_metaplot.pdf", height = 200, width = 32)
#par(cex=2.5)

#ggplot(lncrna_metaplot_df, 
       #aes(x = x, y = dens, color = gene_type )) +
  #geom_vline(xintercept = 0, lty = 2) + 
  #geom_line(size = 2.5) + 
  #facet_wrap(dbp ~ ., scales = "free_y") +
  
  
  #ggtitle("lncrna metaplot promoter") + 
  #scale_x_continuous(breaks = c(-1000, 0, 1000),
                     #labels = c("-1kb", "TSS", "+1kb"),
                     #name = "") + 
  #ylab("Peak frequency") +
  #scale_color_manual(values = c("#424242","#a8404c"))

#dev.off()
```

## Metaplot for mRNA
``` r
#mrna_metaplot_df <- data.frame(x = integer(), dens = numeric(), dbp = character())

#for(i in 1:length(consensus_peaks)) {
  #print(names(consensus_peaks)[[i]])
  #tmp_df <- profile_tss(consensus_peaks[[i]], lncrna_mrna_promoters = mrna_promoters)
  #tmp_df$dbp <- names(consensus_peaks)[[i]]
  #mrna_metaplot_df <- bind_rows(mrna_metaplot_df, tmp_df)
#}

#mrna_metaplot_df$gene_type <- "mRNA"

#pdf("metaplots/mrna_promoter_metaplot.pdf", height = 150, width = 32)
#par(cex=2.3)

#ggplot(mrna_metaplot_df, 
       #aes(x = x, y = dens, color = gene_type )) +
  #geom_vline(xintercept = 0, lty = 2) + 
  #geom_line(size = 1.5) + 
  #facet_wrap(dbp ~ ., scales = "free_y") +
  
  
  #ggtitle("mrna metaplot promoter") + 
  #scale_x_continuous(breaks = c(-1000, 0, 1000),
                     #labels = c("-1kb", "TSS", "+1kb"),
                     #name = "") + 
  #ylab("Peak frequency") +
  #scale_color_manual(values = c("#424242","#a8404c"))

#dev.off()
```

##  Metaplot for both mRNA and lncRNA
``` r
#combined_metaplot_profile <- bind_rows(mrna_metaplot_df, lncrna_metaplot_df)

#pdf("figures/combined_promoter_metaplot.pdf", height = 150, width = 32)
#par(cex=2.3)

# plotting - NOTE facet wrap by dbp !
#ggplot(combined_metaplot_profile, 
       #aes(x = x, y = dens, color = gene_type )) +
  #geom_vline(xintercept = 0, lty = 2) + 
  #geom_line(size = 1.5) + 
  #facet_wrap(dbp ~ ., scales = "free_y") +
  
  
  #ggtitle("Combined Promoter Metaplot") + 
  #scale_x_continuous(breaks = c(-1000, 0, 1000),
                     #labels = c("-1kb", "TSS", "+1kb"),
                     #name = "") + 
  #ylab("Peak frequency") +
  #scale_color_manual(values = c("#424242","#a8404c"))

#dev.off()
```

## plotting the density of DBP localization events
``` r
#peak_occurence_df <- read.csv("/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/final_project/analysis/01_create_consensus_peaks/results/peak_occurence_dataframe.csv")

#ggplot(peak_occurence_df, aes(x = number_of_dbp)) +
#geom_density(alpha = 0.2, color = "#424242", fill = "#424242") +
  #theme_paperwhite() +
  #xlab(expression("Number of DBPs")) +
  #ylab(expression("Density")) +
  #ggtitle("Promoter binding events",
          #subtitle = "mRNA and lncRNA genes")

# saving
#ggsave("/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/final_project/analysis/metaplots/num_binding_events_per_promoter.pdf")
```

## Beginning of RNASeq expression
``` r
# First let's read in the sample sheet to know what is what
#samplesheet <- read_rds("/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/05_R_analyses/05_RNAseq/01_differential_expression/results/final_samplesheet.rds")

# reading in salmon tpm
salmon_tpm <- read.csv("/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/05_R_analyses/05_RNAseq/00_RNAseq_download_NF_core_pipeline/00_NF_CORE_RNAseq_Pipeline_run/results/salmon/salmon_merged_gene_tpm.csv")

#colnames(salmon_tpm)[1] <- "gene_id"

#tpm <- salmon_tpm %>% 
  #pivot_longer(cols = 2:ncol(.), names_to = "sample_id", values_to = "tpm") %>%
  #merge(samplesheet) %>%
  #group_by(gene_id, condition) %>%
  #summarize(tpm = mean(tpm, na.rm = T)) %>%
  #pivot_wider(names_from = condition, values_from = tpm, names_prefix = "tpm_")
  #merge(samplesheet) %>%
    
#View(peak_occurence_df)
    
# Now we can merge in the TPM data to this data.frame
# This will merge on the on teh gene_id column
#promoter_features_df <- merge(peak_occurence_df, tpm)

# saving this file
```

## Number of dbps bound on density plot
``` r
# let's plot the number of DBPs bound as a densuty plot.
#ggplot(promoter_features_df, aes(x = number_of_dbp)) +
  #geom_density() 

# Let's save
#ggsave("/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/final_project/analysis/04_binding_vs_expression/figures/DBP_binding_density_plot.pdf")
```

## Heatmap of genes variable across samples
``` r
# First we need to the tpm DF into a matrix
#tpm_matrix <- tpm %>% 
  #column_to_rownames("gene_id") %>%
  #as.matrix()

# And z-scale each row.
#tpm_scaled <- t(scale(t(tpm_matrix)))

# We want to be sure that we're not removing rows that have NA in only one fraction.
# And right now we're not sure.
#summary(is.na(tpm_matrix))

# checking TPM scaled
#summary(is.na(tpm_scaled))

# clearing all possible NAs
#tpm_scaled <- tpm_scaled[complete.cases(tpm_scaled),]

# Create heatmap using pheatmap function
#heatmap <- pheatmap(tpm_scaled, show_rownames = FALSE)

# Save the heatmap in a PDF file
#pdf("/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/final_project/analysis/metaplots/tpm_heatmap.pdf", height = 20, width = 22)
#heatmap
#dev.off()
```

## Plotting binding vs expression
``` r
#promoter_features_df <- read.csv("/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/final_project/analysis/04_binding_vs_expression/results/promoter_feature_df_tpm.csv")

# Plotting binding versus expression
#ggplot(promoter_features_df, 
       #aes(y = log2(tpm_homo_sapiens_hepg2 + 0.001), x = number_of_dbp, color = gene_type)) + 
  #geom_point(data = promoter_features_df %>% filter(tpm_homo_sapiens_hepg2 > 0.001),
             #shape = 17, alpha = 0.7) +
  
  # Adding a generative additive model 
  #geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cs")) +
  # this adds the statistics from the gam to the figure
  #stat_cor() +
  #geom_smooth(method = "lm") +
  # this is just to make it look nice.
  #scale_x_continuous(expand = c(0,0)) +
  # adding colors manually
  #scale_color_manual(values = c("#a8404c", "#424242"), name = "Gene type") + 
  # title
  #ggtitle("Expression vs. promoter binding events") + 
  # labeling axes
  #xlab(expression('Number of TFs')) +
  #ylab(expression(log[2](TPM))) 
#save plot
#ggsave("/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/final_project/analysis/04_binding_vs_expression/figures/binding_vs_expression.pdf")
```

## Heatmap of nuclear versus cytoplasmic expression
``` r
# Now let's make a similar plot for nuclear RNA abundance versus #DBPs bound to their promoter
#ggplot(promoter_features_df, 
       #aes(y = log2(tpm_homo_sapiens_nuclear_fraction + 0.001), x = number_of_dbp, color = gene_type)) + 
  #geom_point(data = promoter_features_df %>% filter(tpm_homo_sapiens_nuclear_fraction > 0.001),
             #shape = 17, alpha = 0.7) +
  #geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cs")) +
  #stat_cor() +
  #scale_x_continuous(expand = c(0,0)) +
  #scale_color_manual(values = c("#a8404c", "#424242"), name = "Gene type") + 
  #ggtitle("Nuclear Expression vs. promoter binding events") + 
  #xlab(expression('Number of TFs')) +
  #ylab(expression(log[2](TPM))) 
# saving figure
#ggsave("/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/final_project/analysis/04_binding_vs_expression/figures/nuclear_expression-vs-promoter_binding.pdf")
```

## Is there a difference between cyto RNAs versus DBPs on promoter?
``` r
#ggplot(promoter_features_df, 
      # aes(y = log2(tpm_homo_sapiens_cytosolic_fraction + 0.001), x = number_of_dbp, color = gene_type)) + 
  # note if we take off the log2 everything is very lowly expressed
  #geom_point(data = promoter_features_df %>% filter(tpm_homo_sapiens_cytosolic_fraction > 0.001),
             #shape = 17, alpha = 0.7) +
  #geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cs")) +
  #stat_cor() +
  #scale_x_continuous(expand = c(0,0)) +
  #scale_color_manual(values = c("#a8404c", "#424242"), name = "Gene type") + 
  #ggtitle("Cytoplasmic Expression vs. promoter binding events") + 
  #xlab(expression('Number of TFs')) +
  #ylab(expression(log[2](TPM))) 
# saving figure
#ggsave("/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/final_project/analysis/04_binding_vs_expression/figures/cytoplasmic_expression-vs-promoter_binding.pdf")
```

## This is another way to look at it
``` r
# the main function here is color by gene_type in promoter_features_df
#ggplot(promoter_features_df, aes(x = log2(tpm_homo_sapiens_hepg2 + 0.01), color = gene_type))+
 # geom_density()
# saving figure
#ggsave("/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/final_project/analysis/04_binding_vs_expression/figures/mrna_lncrna_tpm_total_rna.pdf")
# let's also do the same for nuclear since lncRNAs are typically more nuclear
#ggplot(promoter_features_df, aes(x = log2(tpm_homo_sapiens_nuclear_fraction + 0.01), color = gene_type))+
  #geom_density()
# saving figure
#ggsave("/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/final_project/analysis/04_binding_vs_expression/figures/mrna_lncrna_tpm_nuclear.pdf")
```

## Extracting proteins without expressions from high-binding proteins
``` r
# filter to super binders
promoter_occurence_df<- peak_occurence_df %>%
  mutate(superbinders = peak_occurence_df$number_of_dbp > 200)
# setting col of superbinders
promoter_occurence_df<- peak_occurence_df %>%
  mutate(superbinder_or_not = ifelse(promoter_occurence_df$superbinders == TRUE, "Superbinder", "notsuperbinder"))

# superbinder promoters
superbinder_promoters <- subset(promoter_occurence_df, superbinder_or_not == "Superbinder")
superbinder_promoters <- dplyr::select(superbinder_promoters, "gene_id")


# non super binder proms
non_superbinder_promoters<- subset(promoter_occurence_df, superbinder_or_not == "notsuperbinder")
non_superbinder_promoters <- dplyr::select(non_superbinder_promoters, "gene_id")


# subet mRNA and lncRNA promoters by super binders
super_lncrna_mrna_gr <- lncrna_mrna_promoters[lncrna_mrna_promoters$gene_id %in% superbinder_promoters$gene_id]
non_super_lncrna_mrna_gr <- lncrna_mrna_promoters[lncrna_mrna_promoters$gene_id %in% non_superbinder_promoters$gene_id]

# setting up superbinders metaplot_Df
superbinder_metaplot_df <- data.frame(x = integer(), dens = numeric(), dbp = character())
lncrna_mrna_promoters = super_lncrna_mrna_gr
# for loop to populate super binder _metaplot
for(i in 1:length(filtered_consensus_list)) {
  print(names(filtered_consensus_list)[[i]])
  tmp_df <- profile_tss(filtered_consensus_list[[i]], lncrna_mrna_promoters)
  tmp_df$dbp <- names(filtered_consensus_list)[[i]]
  superbinder_metaplot_df <- bind_rows(superbinder_metaplot_df, tmp_df)
}
```

    ## [1] "ADNP"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrY
    ##   - in 'y': GL000214.1, GL000218.1, GL000220.1, GL000252.2, KI270438.1, KI270442.1, KI270467.1, KI270733.1, KI270767.1, KI270857.1, KI270878.1, KN196487.1, KQ031389.1, KV575244.1, ML143377.1, ML143380.1, GL000194.1, GL000216.2, GL000219.1, GL000224.1, GL000225.1, GL383577.2, KI270712.1, KI270728.1, KI270729.1, KI270861.1, KQ458384.1, KZ208915.1, ML143372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "AFF4"
    ## [1] "AHDC1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000218.1, GL000219.1, GL000251.2, GL000252.2, GL000254.2, GL339449.2, GL383527.1, GL383556.1, GL383563.3, GL383577.2, GL877875.1, KI270467.1, KI270713.1, KI270721.1, KI270744.1, KI270782.1, KI270803.1, KI270829.1, KI270857.1, KI270861.1, KI270862.1, KI270878.1, KI270879.1, KI270908.1, KN196484.1, KQ458384.1, KV575244.1, KV880768.1, KZ208915.1, ML143371.1, ML143377.1, ML143381.1, GL000224.1, KI270712.1, KI270728.1, KI270733.1, KI270830.1, KI270853.1, KI270869.1, KN538372.1, KV766196.1, KV880764.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "AHR"
    ## [1] "AKAP8"
    ## [1] "ARID3A"
    ## [1] "ARID4A"
    ## [1] "ARID4B"
    ## [1] "ARID5B"
    ## [1] "ARNTL"
    ## [1] "ASH2L"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000219.1, GL383577.2, GL383581.2, KI270713.1, KI270729.1, KI270879.1, KI270904.1, KQ458384.1, KV880768.1, KZ208922.1, GL000218.1, GL000251.2, KZ208915.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ATAD3A"
    ## [1] "ATF2"
    ## [1] "ATF3"
    ## [1] "ATF5"
    ## [1] "ATF6"
    ## [1] "ATF7"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000009.2, GL000194.1, GL000208.1, GL000214.1, GL000216.2, GL000218.1, GL000219.1, GL000220.1, GL000221.1, GL000224.1, GL000225.1, GL000252.2, GL000255.2, GL000256.2, GL339449.2, GL383522.1, GL383542.1, GL383563.3, GL383577.2, GL383581.2, JH159146.1, KI270442.1, KI270706.1, KI270711.1, KI270712.1, KI270713.1, KI270714.1, KI270718.1, KI270719.1, KI270721.1, KI270722.1, KI270723.1, KI270725.1, KI270726.1, KI270728.1, KI270731.1, KI270732.1, KI270733.1, KI270734.1, KI270736.1, KI270742.1, KI270743.1, KI270744.1, KI270745.1, KI270750.1, KI270751.1, KI270754.1, KI270762.1, KI270764.1, KI270772.1, KI270805.1, KI270819.1, KI270830.1, KI270832.1, KI270846.1, KI270847.1, KI270850.1, KI270853.1, KI270861.1, KI270862.1, KI270865.1, KI270868.1, KI270869.1, KI270879.1, KI270894.1, KI270896.1, KI270903.1, KI270905.1, KI270907.1, KI270908.1, KI270913.1, KI270936.1, KN196484.1, KN196487.1, KN538372.1, KQ031389.1, KQ090026.1, KQ090028.1, KQ759762.1, KV575244.1, KV766193.1, KV880764.1, KV880768.1, KZ208908.1, KZ208915.1, KZ208921.1, KZ208922.1, KZ559100.1, ML143353.1, ML143364.1, ML143365.1, ML143371.1, ML143372.1, ML143375.1, ML143377.1, ML143380.1, GL000008.2, GL000195.1, GL000251.2, GL000254.2, GL383556.1, KI270438.1, KI270866.1, KI270878.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "BCL3"
    ## [1] "BCL6"
    ## [1] "BHLHE40"
    ## [1] "BRCA1"
    ## [1] "CAMTA2"
    ## [1] "CBFB"
    ## [1] "CBX5"
    ## [1] "CEBPA"
    ## [1] "CEBPB"
    ## [1] "CEBPD"
    ## [1] "CEBPG"
    ## [1] "CERS6"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrY
    ##   - in 'y': GL000224.1, GL877875.1, KI270438.1, KI270712.1, KI270713.1, KI270733.1, KI270744.1, KN196487.1, KV880768.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "CHD2"
    ## [1] "CREB3"
    ## [1] "CREM"
    ## [1] "CTCF"
    ## [1] "DBP"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrY
    ##   - in 'y': KI270438.1, KI270466.1, KI270467.1, KI270713.1, KI270830.1, KI270869.1, KI270879.1, KN196487.1, KV880768.1, ML143380.1, KI270712.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "DDIT3"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrY
    ##   - in 'y': GL000214.1, GL000216.2, GL000220.1, GL000224.1, GL000225.1, KI270330.1, KI270333.1, KI270336.1, KI270337.1, KI270438.1, KI270442.1, KI270466.1, KI270467.1, KI270713.1, KI270729.1, KI270733.1, KI270879.1, KN196487.1, KV880768.1, ML143380.1, GL383577.2
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "DLX6"
    ## [1] "DMAP1"
    ## [1] "DMTF1"
    ## [1] "DNMT1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000009.2, GL000194.1, GL000205.2, GL000213.1, GL000214.1, GL000216.2, GL000218.1, GL000219.1, GL000220.1, GL000221.1, GL000224.1, GL000225.1, GL000252.2, GL000253.2, GL000254.2, GL000255.2, GL000256.2, GL339449.2, GL383563.3, GL383566.1, GL383574.1, GL383577.2, JH159146.1, JH159147.1, KI270442.1, KI270512.1, KI270519.1, KI270538.1, KI270706.1, KI270709.1, KI270711.1, KI270712.1, KI270713.1, KI270714.1, KI270718.1, KI270719.1, KI270721.1, KI270722.1, KI270723.1, KI270724.1, KI270725.1, KI270728.1, KI270731.1, KI270732.1, KI270733.1, KI270734.1, KI270735.1, KI270736.1, KI270742.1, KI270743.1, KI270744.1, KI270745.1, KI270749.1, KI270750.1, KI270751.1, KI270753.1, KI270754.1, KI270762.1, KI270767.1, KI270772.1, KI270780.1, KI270782.1, KI270791.1, KI270792.1, KI270795.1, KI270797.1, KI270805.1, KI270815.1, KI270816.1, KI270819.1, KI270821.1, KI270830.1, KI270831.1, KI270832.1, KI270846.1, KI270847.1, KI270850.1, KI270853.1, KI270854.1, KI270856.1, KI270857.1, KI270860.1, KI270861.1, KI270862.1, KI270865.1, KI270868.1, KI270869.1, KI270870.1, KI270871.1, KI270872.1, KI270875.1, KI270877.1, KI270878.1, KI270879.1, KI270894.1, KI270896.1, KI270901.1, KI270903.1, KI270905.1, KI270908.1, KI270913.1, KI270937.1, KN196484.1, KN538361.1, KN538364.1, KN538368.1, KN538372.1, KN538373.1, KQ031389.1, KQ090026.1, KQ090028.1, KQ458383.1, KQ458386.1, KQ759759.1, KQ759762.1, KV575244.1, KV766196.1, KV880764.1, KV880768.1, KZ208904.1, KZ208906.1, KZ208915.1, KZ208918.1, KZ208919.1, KZ208921.1, KZ208922.1, KZ559100.1, KZ559109.1, ML143341.1, ML143353.1, ML143355.1, ML143358.1, ML143364.1, ML143365.1, ML143366.1, ML143371.1, ML143372.1, ML143375.1, ML143377.1, ML143378.1, ML143379.1, ML143380.1, GL000258.2, KI270589.1, KI270876.1, KI270907.1, KI270936.1, KQ090022.1, KZ208908.1, KZ208916.1, ML143370.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "DPF2"
    ## [1] "DRAP1"
    ## [1] "DZIP1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000214.1, GL000218.1, GL000219.1, GL000225.1, GL383577.2, KI270713.1, KI270721.1, KI270733.1, KI270869.1, KI270879.1, KI270908.1, KZ208915.1, GL000008.2, GL000009.2, GL000194.1, GL000208.1, GL000216.2, GL000220.1, GL000221.1, GL000224.1, GL000251.2, GL000252.2, GL000253.2, GL000254.2, GL000256.2, GL339449.2, GL383522.1, GL383533.1, GL383563.3, GL383566.1, GL383574.1, GL877875.1, KI270330.1, KI270442.1, KI270706.1, KI270709.1, KI270711.1, KI270712.1, KI270714.1, KI270718.1, KI270719.1, KI270720.1, KI270722.1, KI270723.1, KI270724.1, KI270725.1, KI270728.1, KI270731.1, KI270732.1, KI270734.1, KI270742.1, KI270743.1, KI270744.1, KI270745.1, KI270746.1, KI270750.1, KI270751.1, KI270753.1, KI270754.1, KI270762.1, KI270764.1, KI270772.1, KI270787.1, KI270791.1, KI270805.1, KI270819.1, KI270821.1, KI270824.1, KI270830.1, KI270831.1, KI270832.1, KI270845.1, KI270847.1, KI270853.1, KI270854.1, KI270857.1, KI270861.1, KI270868.1, KI270871.1, KI270872.1, KI270877.1, KI270878.1, KI270894.1, KI270902.1, KI270903.1, KI270905.1, KI270907.1, KN196484.1, KN196487.1, KN538372.1, KQ031389.1, KQ090020.1, KQ090026.1, KQ458386.1, KV575244.1, KV766196.1, KV880768.1, KZ208904.1, KZ208906.1, KZ208918.1, KZ208921.1, KZ208922.1, KZ559100.1, ML143353.1, ML143364.1, ML143365.1, ML143366.1, ML143367.1, ML143371.1, ML143372.1, ML143373.1, ML143374.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "E2F1"
    ## [1] "E2F2"
    ## [1] "E2F4"
    ## [1] "E2F5"
    ## [1] "EEA1"
    ## [1] "EED"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000219.1, GL000224.1, KI270712.1, KI270713.1, KI270721.1, KI270728.1, KI270733.1, KI270744.1, KI270869.1, KI270879.1, KI270894.1, KI270908.1, KZ208915.1, GL000009.2, GL000214.1, GL000216.2, GL000218.1, GL000225.1, GL000252.2, GL000255.2, GL339449.2, GL383533.1, GL383557.1, GL877875.1, JH159146.1, KI270442.1, KI270709.1, KI270711.1, KI270714.1, KI270718.1, KI270719.1, KI270722.1, KI270725.1, KI270734.1, KI270745.1, KI270750.1, KI270754.1, KI270772.1, KI270830.1, KI270847.1, KI270871.1, KI270877.1, KI270878.1, KI270903.1, KI270907.1, KN196484.1, KN538372.1, KQ090026.1, KV880768.1, KZ208908.1, KZ208921.1, KZ208922.1, ML143345.1, ML143364.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "EGR1"
    ## [1] "ELF1"
    ## [1] "ELF3"
    ## [1] "ELF4"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrY, chrM
    ##   - in 'y': GL000194.1, GL000219.1, GL339449.2, GL383577.2, KI270713.1, KI270879.1, GL000254.2
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ELK1"
    ## [1] "ERF"
    ## [1] "ESRRA"
    ## [1] "ETS1"
    ## [1] "ETV5"
    ## [1] "ETV6"
    ## [1] "EZH2"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrY, chrM
    ##   - in 'y': KI270712.1, KI270729.1, KI270824.1, KI270830.1, KI270894.1, KQ090026.1, ML143364.1, KI270733.1, KI270845.1, KI270871.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "FOSL2"
    ## [1] "FOXA1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000218.1, GL000219.1, GL000251.2, GL000254.2, GL000256.2, GL000257.2, GL383556.1, GL383577.2, GL383580.2, KI270713.1, KI270721.1, KI270744.1, KI270772.1, KI270827.1, KI270830.1, KI270853.1, KI270857.1, KI270861.1, KI270862.1, KI270878.1, KI270879.1, KI270904.1, KI270905.1, KQ031389.1, KQ090028.1, KQ458383.1, KV880768.1, KZ208915.1, KZ208918.1, KZ208919.1, ML143343.1, ML143381.1, GL339449.2, GL383527.1, GL383555.2, GL949742.1, KI270734.1, KI270801.1, KI270908.1, KN196474.1, KN196483.1, KN196484.1, KQ458384.1, KZ208922.1, KZ559105.1, ML143360.1, ML143366.1, ML143371.1, GL000253.2, KI270803.1, KN196487.1, GL383552.1, KI270733.1, KI270751.1, KI270850.1, KV575244.1, GL000214.1, GL000216.2, GL000224.1, GL000225.1, KI270330.1, KI270438.1, KI270712.1, KI270730.1, KI270754.1, ML143372.1, ML143377.1, ML143380.1, GL000194.1, KI270442.1, KI270709.1, KI270723.1, KI270735.1, KI270745.1, KI270782.1, KI270821.1, KI270896.1, KQ031384.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "FOXA2"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000205.2, GL000219.1, GL000251.2, GL000254.2, GL000257.2, GL339449.2, GL383526.1, GL383552.1, GL383556.1, GL383563.3, GL383571.1, GL383577.2, GL383580.2, GL383581.2, KI270713.1, KI270721.1, KI270723.1, KI270734.1, KI270744.1, KI270782.1, KI270801.1, KI270822.1, KI270827.1, KI270830.1, KI270853.1, KI270857.1, KI270861.1, KI270862.1, KI270878.1, KI270879.1, KI270896.1, KI270904.1, KI270905.1, KN196484.1, KQ031389.1, KQ090022.1, KQ090026.1, KQ090028.1, KQ458383.1, KV575244.1, KV880768.1, KZ208915.1, KZ208918.1, KZ208919.1, KZ559105.1, ML143343.1, ML143360.1, ML143371.1, ML143377.1, ML143381.1, KI270333.1, KI270337.1, KI270466.1, KI270467.1, KN196487.1, ML143380.1, GL000214.1, GL000216.2, GL000224.1, GL000225.1, GL000252.2, KI270330.1, KI270438.1, KI270709.1, KI270729.1, KI270730.1, KI270735.1, KI270751.1, KI270754.1, ML143366.1, GL000194.1, GL000195.1, GL000218.1, GL000253.2, GL000255.2, KI270442.1, KI270519.1, KI270720.1, KI270733.1, KI270803.1, KI270908.1, KQ458384.1, KQ759762.1, KZ208922.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "FOXA3"
    ## [1] "FOXJ3"
    ## [1] "FOXK1"
    ## [1] "FOXM1"
    ## [1] "FOXO1"
    ## [1] "FOXP1"
    ## [1] "FOXP4"
    ## [1] "FUBP1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000219.1, GL000252.2, GL000255.2, GL339449.2, GL383563.3, GL383577.2, GL383580.2, KI270713.1, KI270721.1, KI270734.1, KI270736.1, KI270744.1, KI270773.1, KI270830.1, KI270850.1, KI270857.1, KI270861.1, KI270862.1, KI270869.1, KI270878.1, KI270879.1, KI270905.1, KQ090026.1, KQ458383.1, KQ458384.1, KV880768.1, ML143365.1, ML143377.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "FUBP3"
    ## [1] "GABPA"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000219.1, GL000251.2, GL000255.2, GL339449.2, KI270816.1, KI270830.1, KI270879.1, KV880768.1, ML143371.1, ML143377.1, ML143380.1, GL000194.1, GL000218.1, GL000254.2, KZ208915.1, ML143364.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "GABPB1"
    ## [1] "GATA2"
    ## [1] "GATAD1"
    ## [1] "GATAD2A"
    ## [1] "GLI4"
    ## [1] "GMEB1"
    ## [1] "GMEB2"
    ## [1] "GPN1"
    ## [1] "GTF2F1"
    ## [1] "GTF3A"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': KI270713.1, KI270733.1, KI270879.1, KZ208906.1, GL000194.1, GL000251.2, GL339449.2, GL383577.2, KI270712.1, KI270721.1, KI270728.1, KI270905.1, KN196484.1, KV880768.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "GZF1"
    ## [1] "H2AFZ"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': KI270442.1, KI270729.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "H3K27ac"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000218.1, GL000251.2, GL000254.2, GL000255.2, GL000257.2, GL339449.2, GL383563.3, GL383574.1, KI270713.1, KI270721.1, KI270729.1, KI270744.1, KI270782.1, KI270857.1, KI270878.1, KI270879.1, KI270903.1, KI270904.1, KN196484.1, KQ458384.1, KV575244.1, KV880768.1, KZ208915.1, ML143366.1, ML143372.1, ML143377.1, ML143380.1, GL000194.1, GL000219.1, GL383577.2, GL383581.2, KI270861.1, KQ759762.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "H3K36me3"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000225.1, GL383557.1, KI270442.1, KI270711.1, KI270729.1, KI270736.1, KI270801.1, KI270853.1, KI270866.1, KI270879.1, KQ090028.1, KV575244.1, KZ559109.1, GL000218.1, GL000219.1, KI270467.1, KI270728.1, GL383522.1, KN196474.1, KV766194.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "H3K4me1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000218.1, GL000219.1, GL000257.2, GL383574.1, GL383577.2, GL383581.2, KI270442.1, KI270713.1, KI270729.1, KI270744.1, KI270829.1, KI270830.1, KI270847.1, KI270856.1, KI270861.1, KI270877.1, KI270878.1, KI270879.1, KI270903.1, KN196484.1, KQ090026.1, KQ458383.1, KQ458384.1, KV880764.1, KZ208919.1, KZ208922.1, ML143360.1, ML143366.1, ML143377.1, ML143380.1, KZ208915.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "H3K4me2"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000205.2, GL000214.1, GL000218.1, GL000219.1, GL000220.1, GL000251.2, GL000252.2, GL000253.2, GL000254.2, GL000255.2, GL000256.2, GL339449.2, GL383552.1, GL383563.3, GL383574.1, GL383577.2, GL383580.2, KI270706.1, KI270712.1, KI270713.1, KI270714.1, KI270729.1, KI270744.1, KI270782.1, KI270803.1, KI270809.1, KI270821.1, KI270830.1, KI270831.1, KI270832.1, KI270853.1, KI270856.1, KI270857.1, KI270860.1, KI270861.1, KI270878.1, KI270879.1, KI270903.1, KI270904.1, KI270905.1, KI270908.1, KN196481.1, KN196484.1, KN538364.1, KQ090026.1, KQ458383.1, KQ458384.1, KQ759762.1, KV575244.1, KV880764.1, KV880768.1, KZ208904.1, KZ208906.1, KZ208915.1, KZ208921.1, KZ208922.1, KZ559105.1, ML143360.1, ML143366.1, ML143367.1, ML143372.1, ML143377.1, ML143380.1, KV766196.1, ML143371.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "H3K4me3"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000218.1, GL000219.1, GL000220.1, GL000251.2, GL000252.2, GL000253.2, GL000254.2, GL000255.2, GL000256.2, GL339449.2, GL383563.3, KI270712.1, KI270713.1, KI270714.1, KI270717.1, KI270721.1, KI270744.1, KI270762.1, KI270816.1, KI270824.1, KI270830.1, KI270832.1, KI270847.1, KI270849.1, KI270850.1, KI270853.1, KI270857.1, KI270861.1, KI270872.1, KI270879.1, KI270894.1, KI270897.1, KI270904.1, KI270905.1, KI270908.1, KN538364.1, KQ031389.1, KQ090026.1, KQ759762.1, KV575244.1, KV880764.1, KV880768.1, KZ208906.1, KZ208915.1, KZ208922.1, KZ559105.1, KZ559109.1, ML143366.1, ML143371.1, ML143372.1, ML143377.1, ML143380.1, KI270765.1, KI270848.1, KI270856.1, KZ208917.1, KI270333.1, KI270337.1, KI270466.1, KI270467.1, KI270878.1, KI270336.1, KI270438.1, KI270442.1, KI270515.1, KI270729.1, KI270735.1, KI270736.1, KI270750.1, KN196487.1, GL000216.2, GL000224.1, GL000225.1, KI270435.1, KI270509.1, KI270519.1, KI270588.1, KI270591.1, KI270709.1, KI270746.1, KI270751.1, KI270757.1, GL383556.1, GL383574.1, KI270745.1, KI270815.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "H3K79me2"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL383557.1, KI270729.1, KI270744.1, KI270830.1, KI270832.1, KI270866.1, KI270879.1, KI270904.1, KN538364.1, KQ458384.1, KV575244.1, KV880768.1, KZ208915.1, ML143366.1, ML143377.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "H3K9ac"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000218.1, GL000219.1, GL000220.1, GL000251.2, GL000252.2, GL339449.2, GL383563.3, GL383574.1, KI270713.1, KI270721.1, KI270729.1, KI270736.1, KI270744.1, KI270745.1, KI270830.1, KI270857.1, KI270879.1, KI270904.1, KI270908.1, KN196484.1, KQ458384.1, KV575244.1, KV880768.1, KZ208915.1, ML143377.1, ML143380.1, GL000254.2, GL000255.2, KI270866.1, KI270872.1, KI270903.1, KQ090026.1, KV880764.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "H3K9me3"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000008.2, GL000009.2, GL000208.1, GL000216.2, GL000224.1, GL000225.1, KI270435.1, KI270438.1, KI270442.1, KI270468.1, KI270511.1, KI270512.1, KI270515.1, KI270516.1, KI270517.1, KI270519.1, KI270522.1, KI270538.1, KI270584.1, KI270587.1, KI270589.1, KI270590.1, KI270591.1, KI270709.1, KI270716.1, KI270718.1, KI270719.1, KI270720.1, KI270722.1, KI270723.1, KI270725.1, KI270729.1, KI270730.1, KI270732.1, KI270735.1, KI270736.1, KI270746.1, KI270747.1, KI270749.1, KI270750.1, KI270751.1, KI270753.1, KI270756.1, KI270757.1, KN196487.1, KN538372.1, KQ031387.1, KV766197.1, ML143354.1, ML143364.1, ML143378.1, GL000194.1, GL000218.1, GL000219.1, GL000221.1, KI270320.1, KI270322.1, KI270333.1, KI270337.1, KI270429.1, KI270448.1, KI270465.1, KI270466.1, KI270467.1, KI270507.1, KI270508.1, KI270509.1, KI270510.1, KI270518.1, KI270521.1, KI270580.1, KI270583.1, KI270588.1, KI270593.1, KI270710.1, KI270712.1, KI270724.1, KI270728.1, KI270733.1, KI270738.1, KI270742.1, KI270743.1, KI270775.1, KI270839.1, KI270881.1, KI270907.1, KN538360.1, KN538367.1, KN538368.1, KQ031384.1, KV766199.1, KZ208924.1, ML143359.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "H4K20me1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': KI270333.1, KI270729.1, KI270819.1, KI270866.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "HCFC1"
    ## [1] "HDAC1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000195.1, GL000218.1, GL000219.1, GL000251.2, GL000252.2, GL000253.2, GL000254.2, GL000255.2, GL000256.2, GL339449.2, GL383556.1, GL383563.3, GL383577.2, KI270466.1, KI270467.1, KI270712.1, KI270713.1, KI270714.1, KI270721.1, KI270733.1, KI270734.1, KI270744.1, KI270750.1, KI270751.1, KI270782.1, KI270816.1, KI270830.1, KI270832.1, KI270850.1, KI270853.1, KI270857.1, KI270860.1, KI270861.1, KI270862.1, KI270869.1, KI270872.1, KI270878.1, KI270879.1, KI270897.1, KI270904.1, KI270905.1, KI270908.1, KN196484.1, KN196487.1, KQ090026.1, KQ458384.1, KQ759762.1, KV575244.1, KV766196.1, KV880764.1, KV880768.1, KZ208915.1, KZ208921.1, KZ208922.1, KZ559109.1, ML143371.1, ML143377.1, ML143380.1, KI270745.1, KI270842.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "HDAC2"
    ## [1] "HINFP"
    ## [1] "HIVEP1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000009.2, GL000194.1, GL000214.1, GL000218.1, GL000219.1, GL000220.1, GL000221.1, GL000224.1, GL000225.1, GL000251.2, GL000252.2, GL000253.2, GL000254.2, GL000255.2, GL000257.2, GL339449.2, GL383522.1, GL383542.1, GL383563.3, GL383574.1, GL383577.2, GL383581.2, GL877875.1, JH159146.1, KI270442.1, KI270706.1, KI270709.1, KI270711.1, KI270712.1, KI270713.1, KI270719.1, KI270721.1, KI270722.1, KI270723.1, KI270728.1, KI270731.1, KI270733.1, KI270734.1, KI270743.1, KI270744.1, KI270745.1, KI270751.1, KI270754.1, KI270772.1, KI270791.1, KI270810.1, KI270816.1, KI270819.1, KI270827.1, KI270830.1, KI270832.1, KI270842.1, KI270846.1, KI270850.1, KI270853.1, KI270854.1, KI270857.1, KI270861.1, KI270862.1, KI270869.1, KI270878.1, KI270879.1, KI270894.1, KI270901.1, KI270903.1, KI270905.1, KI270907.1, KI270908.1, KI270913.1, KN196484.1, KN538364.1, KN538372.1, KQ090026.1, KQ759762.1, KV575244.1, KV766196.1, KV880764.1, KV880768.1, KZ208906.1, KZ208908.1, KZ208915.1, KZ208918.1, KZ208921.1, KZ208922.1, KZ559105.1, KZ559109.1, ML143343.1, ML143364.1, ML143371.1, ML143372.1, ML143377.1, ML143380.1, GL000205.2, GL000208.1, GL000216.2, GL000256.2, GL000258.2, GL383533.1, GL383550.2, GL383566.1, GL949742.1, KI270707.1, KI270714.1, KI270724.1, KI270725.1, KI270729.1, KI270732.1, KI270742.1, KI270749.1, KI270750.1, KI270753.1, KI270762.1, KI270765.1, KI270781.1, KI270792.1, KI270805.1, KI270824.1, KI270831.1, KI270856.1, KI270860.1, KI270868.1, KI270871.1, KI270872.1, KI270896.1, KI270902.1, KN196480.1, KN196487.1, KQ090022.1, KQ458383.1, KQ458386.1, KZ208914.1, KZ559103.1, ML143358.1, ML143367.1, ML143374.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "HMG20A"
    ## [1] "HMG20B"
    ## [1] "HMGXB3"
    ## [1] "HMGXB4"
    ## [1] "HNF1A"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000218.1, GL000219.1, GL000251.2, GL000252.2, GL000254.2, GL000257.2, GL339449.2, GL383563.3, GL383577.2, KI270442.1, KI270713.1, KI270721.1, KI270744.1, KI270782.1, KI270819.1, KI270830.1, KI270832.1, KI270853.1, KI270857.1, KI270861.1, KI270862.1, KI270872.1, KI270878.1, KI270879.1, KI270905.1, KI270907.1, KN196484.1, KQ031389.1, KQ983256.1, KV575244.1, KV880768.1, KZ208915.1, KZ208921.1, KZ208922.1, ML143343.1, ML143371.1, ML143372.1, ML143377.1, ML143380.1, GL000008.2, GL000009.2, GL000195.1, GL000205.2, GL000214.1, GL000216.2, GL000220.1, GL000224.1, GL000225.1, GL000258.2, GL383533.1, KI270438.1, KI270589.1, KI270591.1, KI270709.1, KI270712.1, KI270718.1, KI270719.1, KI270720.1, KI270722.1, KI270723.1, KI270725.1, KI270728.1, KI270730.1, KI270732.1, KI270733.1, KI270734.1, KI270735.1, KI270736.1, KI270745.1, KI270746.1, KI270750.1, KI270751.1, KI270753.1, KI270754.1, KI270756.1, KI270757.1, KI270762.1, KI270816.1, KI270821.1, KI270846.1, KI270869.1, KI270894.1, KI270913.1, KN538372.1, KV880764.1, ML143364.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "HNF1B"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000195.1, GL000218.1, GL000219.1, GL000224.1, GL000251.2, GL000252.2, GL000253.2, GL000254.2, GL000255.2, GL339449.2, GL383555.2, GL383556.1, GL383557.1, GL383563.3, GL383574.1, GL383577.2, GL383580.2, GL877875.1, KI270442.1, KI270711.1, KI270713.1, KI270714.1, KI270721.1, KI270728.1, KI270731.1, KI270734.1, KI270744.1, KI270751.1, KI270781.1, KI270782.1, KI270803.1, KI270809.1, KI270830.1, KI270832.1, KI270845.1, KI270846.1, KI270850.1, KI270853.1, KI270856.1, KI270857.1, KI270861.1, KI270869.1, KI270872.1, KI270878.1, KI270879.1, KI270894.1, KI270903.1, KI270904.1, KI270905.1, KI270907.1, KI270908.1, KN196483.1, KN196484.1, KQ090026.1, KQ458383.1, KQ458384.1, KQ759762.1, KQ983256.1, KV575244.1, KV880763.1, KV880764.1, KV880768.1, KZ208906.1, KZ208915.1, KZ208918.1, KZ208921.1, KZ208922.1, KZ559105.1, ML143360.1, ML143367.1, ML143372.1, ML143377.1, GL000205.2, GL000225.1, GL000256.2, GL383578.2, KI270762.1, KI270862.1, KI270875.1, KI270896.1, KV766197.1, ML143371.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "HNF4A"
    ## [1] "HNF4G"
    ## [1] "HOMEZ"
    ## [1] "HOXA3"
    ## [1] "HOXA5"
    ## [1] "HOXD1"
    ## [1] "HSF2"
    ## [1] "IKZF4"
    ## [1] "IKZF5"
    ## [1] "IRF1"
    ## [1] "IRF2"
    ## [1] "IRF5"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000218.1, GL000219.1, GL000220.1, GL339449.2, KI270712.1, KI270713.1, KI270721.1, KI270728.1, KI270733.1, KI270744.1, KI270821.1, KI270830.1, KI270861.1, KI270869.1, KI270879.1, KI270903.1, KI270908.1, KV880768.1, KZ208915.1, ML143377.1, GL383577.2, KI270442.1, KI270857.1, KI270878.1, KI270894.1, ML143372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ISL2"
    ## [1] "ISX"
    ## [1] "JUN"
    ## [1] "JUND"
    ## [1] "KAT7"
    ## [1] "KAT8"
    ## [1] "KDM1A"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000219.1, GL000252.2, GL383563.3, GL383577.2, KI270713.1, KI270853.1, KI270857.1, KI270861.1, KI270878.1, KV880768.1, ML143377.1, GL000216.2, GL000218.1, GL000224.1, GL000225.1, KI270330.1, KI270467.1, KI270733.1, KI270879.1, KI270442.1, KI270729.1, KN196480.1, KN196484.1, KN538361.1, KQ031389.1, KV880764.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "KDM2A"
    ## [1] "KDM3A"
    ## [1] "KDM4B"
    ## [1] "KDM5B"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000214.1, GL000218.1, GL000219.1, GL000255.2, GL000256.2, GL383563.3, KI270706.1, KI270712.1, KI270713.1, KI270714.1, KI270721.1, KI270728.1, KI270731.1, KI270733.1, KI270734.1, KI270744.1, KI270830.1, KI270832.1, KI270857.1, KI270862.1, KI270869.1, KI270879.1, KI270894.1, KI270903.1, KI270904.1, KI270905.1, KI270908.1, KN196484.1, KN538364.1, KQ090026.1, KV880768.1, KZ208915.1, KZ208921.1, KZ208922.1, ML143353.1, ML143372.1, GL000008.2, GL000195.1, GL000205.2, GL000220.1, GL000224.1, GL000225.1, GL000251.2, GL000252.2, GL000253.2, GL000254.2, GL000258.2, GL339449.2, GL383522.1, GL383533.1, GL383555.2, GL383556.1, GL383557.1, GL383566.1, GL383574.1, GL383577.2, GL383580.2, GL383581.2, GL877875.1, JH159146.1, JH159148.1, KI270442.1, KI270707.1, KI270709.1, KI270711.1, KI270719.1, KI270723.1, KI270726.1, KI270745.1, KI270746.1, KI270750.1, KI270751.1, KI270754.1, KI270762.1, KI270772.1, KI270780.1, KI270782.1, KI270791.1, KI270804.1, KI270809.1, KI270816.1, KI270819.1, KI270820.1, KI270821.1, KI270824.1, KI270831.1, KI270835.1, KI270842.1, KI270846.1, KI270847.1, KI270850.1, KI270853.1, KI270854.1, KI270856.1, KI270860.1, KI270861.1, KI270868.1, KI270870.1, KI270872.1, KI270877.1, KI270878.1, KI270896.1, KI270897.1, KI270907.1, KI270927.1, KI270936.1, KI270937.1, KN196480.1, KN196487.1, KN538372.1, KQ031389.1, KQ090021.1, KQ458383.1, KQ458386.1, KQ759762.1, KV575244.1, KV766193.1, KV766196.1, KV880764.1, KZ208904.1, KZ208906.1, KZ208908.1, KZ208919.1, KZ559105.1, KZ559109.1, ML143341.1, ML143358.1, ML143364.1, ML143365.1, ML143366.1, ML143367.1, ML143371.1, ML143377.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "KDM6A"
    ## [1] "KLF11"
    ## [1] "KLF12"
    ## [1] "KLF13"
    ## [1] "KLF16"
    ## [1] "KLF6"
    ## [1] "KLF9"
    ## [1] "KMT2A"
    ## [1] "KMT2B"
    ## [1] "LBX2"
    ## [1] "LCOR"
    ## [1] "LCORL"
    ## [1] "LIN54"
    ## [1] "LRRFIP1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrY
    ##   - in 'y': GL000216.2, GL000224.1, GL000225.1, KI270438.1, KI270442.1, KI270466.1, KI270467.1, KI270713.1, KI270729.1, KI270733.1, KI270744.1, KI270869.1, KN196487.1, KV880768.1, ML143380.1, GL000219.1, GL383563.3, GL383577.2, KI270330.1, KI270731.1, KI270879.1, KN196484.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "MAF1"
    ## [1] "MAFF"
    ## [1] "MAFG"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000224.1, GL000254.2, KI270712.1, KI270728.1, KI270733.1, KI270744.1, KI270762.1, KI270853.1, KI270908.1, KN196483.1, KN538372.1, KZ559100.1, ML143377.1, GL000220.1, KI270330.1, KI270754.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "MAFK"
    ## [1] "MATR3"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrY
    ##   - in 'y': GL000008.2, GL000194.1, GL000214.1, GL000216.2, GL000218.1, GL000219.1, GL000220.1, GL000224.1, GL000225.1, GL877875.1, KI270330.1, KI270438.1, KI270442.1, KI270589.1, KI270709.1, KI270712.1, KI270713.1, KI270718.1, KI270719.1, KI270723.1, KI270728.1, KI270729.1, KI270730.1, KI270733.1, KI270734.1, KI270735.1, KI270736.1, KI270744.1, KI270754.1, KI270756.1, KI270757.1, KI270762.1, KI270772.1, KI270830.1, KI270878.1, KI270879.1, KI270894.1, KN196487.1, KN538372.1, KV880768.1, ML143377.1, ML143380.1, KI270731.1, KI270819.1, KQ031384.1, KQ031389.1, KQ090026.1, KZ208915.1, KZ208922.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "MAX"
    ## [1] "MAZ"
    ## [1] "MED1"
    ## [1] "MED13"
    ## [1] "MEF2A"
    ## [1] "MEF2D"
    ## [1] "MIER2"
    ## [1] "MIER3"
    ## [1] "MIXL1"
    ## [1] "MLX"
    ## [1] "MNX1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000214.1, GL000218.1, GL000219.1, GL000224.1, GL000251.2, GL000252.2, GL000253.2, GL000254.2, GL000255.2, GL000256.2, GL339449.2, GL383556.1, GL383563.3, GL383574.1, GL383577.2, GL383581.2, KI270442.1, KI270712.1, KI270713.1, KI270714.1, KI270719.1, KI270721.1, KI270728.1, KI270733.1, KI270734.1, KI270744.1, KI270745.1, KI270754.1, KI270782.1, KI270816.1, KI270821.1, KI270830.1, KI270832.1, KI270850.1, KI270853.1, KI270857.1, KI270861.1, KI270862.1, KI270869.1, KI270872.1, KI270877.1, KI270878.1, KI270879.1, KI270894.1, KI270905.1, KI270907.1, KI270908.1, KN196484.1, KN538364.1, KN538372.1, KQ090026.1, KQ458384.1, KQ759762.1, KV575244.1, KV766196.1, KV880768.1, KZ208915.1, KZ208922.1, KZ559109.1, ML143371.1, ML143372.1, ML143377.1, ML143380.1, GL000009.2, GL000205.2, GL000216.2, GL000225.1, GL383527.1, GL383533.1, GL383550.2, GL383557.1, GL877875.1, JH159146.1, JH159147.1, KI270706.1, KI270709.1, KI270711.1, KI270717.1, KI270722.1, KI270723.1, KI270726.1, KI270729.1, KI270731.1, KI270748.1, KI270762.1, KI270765.1, KI270772.1, KI270780.1, KI270791.1, KI270802.1, KI270804.1, KI270813.1, KI270819.1, KI270827.1, KI270835.1, KI270845.1, KI270846.1, KI270847.1, KI270854.1, KI270856.1, KI270868.1, KI270871.1, KI270897.1, KI270903.1, KI270904.1, KI270913.1, KQ090021.1, KQ090022.1, KQ458383.1, KQ458386.1, KZ208904.1, KZ208906.1, KZ208909.1, KZ208914.1, KZ208921.1, KZ559105.1, ML143355.1, ML143358.1, ML143364.1, ML143366.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "MTA1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000216.2, GL000218.1, GL000219.1, GL000252.2, GL383563.3, GL383577.2, KI270438.1, KI270442.1, KI270713.1, KI270714.1, KI270718.1, KI270721.1, KI270754.1, KI270782.1, KI270857.1, KI270861.1, KI270869.1, KI270878.1, KI270905.1, KN196484.1, KV880768.1, KZ208915.1, ML143377.1, GL000009.2, GL000194.1, GL000208.1, GL000214.1, GL000220.1, GL000221.1, GL000224.1, GL000225.1, GL000251.2, GL000253.2, GL000254.2, GL000255.2, GL000256.2, GL000258.2, GL339449.2, GL383522.1, GL383533.1, GL383550.2, GL383555.2, GL383557.1, GL383574.1, GL383581.2, GL877875.1, GL949742.1, JH159146.1, JH159147.1, KI270508.1, KI270519.1, KI270706.1, KI270707.1, KI270709.1, KI270712.1, KI270719.1, KI270720.1, KI270722.1, KI270723.1, KI270724.1, KI270725.1, KI270726.1, KI270728.1, KI270729.1, KI270731.1, KI270732.1, KI270733.1, KI270734.1, KI270735.1, KI270742.1, KI270743.1, KI270744.1, KI270745.1, KI270746.1, KI270749.1, KI270750.1, KI270751.1, KI270753.1, KI270762.1, KI270764.1, KI270765.1, KI270767.1, KI270772.1, KI270791.1, KI270792.1, KI270795.1, KI270797.1, KI270804.1, KI270805.1, KI270810.1, KI270813.1, KI270816.1, KI270819.1, KI270821.1, KI270824.1, KI270827.1, KI270830.1, KI270831.1, KI270832.1, KI270835.1, KI270842.1, KI270845.1, KI270846.1, KI270847.1, KI270850.1, KI270853.1, KI270854.1, KI270856.1, KI270860.1, KI270862.1, KI270868.1, KI270870.1, KI270871.1, KI270872.1, KI270876.1, KI270877.1, KI270879.1, KI270894.1, KI270897.1, KI270901.1, KI270902.1, KI270903.1, KI270907.1, KI270908.1, KI270910.1, KI270913.1, KI270927.1, KI270936.1, KI270937.1, KN196481.1, KN196487.1, KN538364.1, KN538372.1, KQ031389.1, KQ090017.1, KQ090022.1, KQ090026.1, KQ458384.1, KQ458386.1, KQ759759.1, KQ759762.1, KV575244.1, KV766193.1, KV766196.1, KV880764.1, KZ208904.1, KZ208906.1, KZ208908.1, KZ208909.1, KZ208917.1, KZ208918.1, KZ208919.1, KZ208921.1, KZ208922.1, KZ559100.1, KZ559105.1, KZ559109.1, KZ559113.1, ML143353.1, ML143355.1, ML143358.1, ML143361.1, ML143364.1, ML143365.1, ML143366.1, ML143367.1, ML143370.1, ML143371.1, ML143372.1, ML143373.1, ML143374.1, ML143378.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "MTF2"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000195.1, GL000251.2, GL000252.2, GL000255.2, GL000256.2, GL339449.2, GL383563.3, KI270712.1, KI270713.1, KI270717.1, KI270731.1, KI270765.1, KI270824.1, KI270829.1, KI270830.1, KI270850.1, KI270857.1, KI270862.1, KI270869.1, KI270879.1, KI270894.1, KI270905.1, KI270908.1, KQ458383.1, KV880768.1, KZ208906.1, KZ208915.1, KZ559113.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "MXD1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000214.1, GL000216.2, GL000218.1, GL000219.1, GL000220.1, GL000224.1, GL000225.1, GL000252.2, GL000254.2, GL000255.2, GL339449.2, GL383533.1, GL383563.3, GL383574.1, GL383577.2, JH159147.1, KI270435.1, KI270438.1, KI270442.1, KI270589.1, KI270591.1, KI270706.1, KI270709.1, KI270712.1, KI270713.1, KI270714.1, KI270719.1, KI270721.1, KI270722.1, KI270723.1, KI270726.1, KI270728.1, KI270733.1, KI270734.1, KI270744.1, KI270745.1, KI270750.1, KI270754.1, KI270762.1, KI270781.1, KI270782.1, KI270813.1, KI270816.1, KI270821.1, KI270830.1, KI270832.1, KI270846.1, KI270850.1, KI270853.1, KI270857.1, KI270861.1, KI270862.1, KI270868.1, KI270869.1, KI270879.1, KI270894.1, KI270896.1, KI270897.1, KI270904.1, KI270905.1, KI270908.1, KI270927.1, KN196484.1, KN538364.1, KN538372.1, KQ458383.1, KQ458386.1, KQ759762.1, KV575244.1, KV766193.1, KV766196.1, KV880768.1, KZ208906.1, KZ208915.1, KZ208921.1, KZ559113.1, ML143364.1, ML143371.1, ML143377.1, ML143380.1, GL000256.2, KI270872.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "MXD3"
    ## [1] "MXD4"
    ## [1] "MXI1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000225.1, KI270438.1, KI270467.1, KI270713.1, KI270729.1, KI270733.1, KI270879.1, KN196487.1, GL000194.1, GL000219.1, GL000252.2, GL339449.2, GL383563.3, GL383577.2, KI270830.1, KI270857.1, KI270861.1, KI270862.1, KQ759762.1, KV575244.1, KV880768.1, KZ208915.1, ML143367.1, ML143377.1, ML143380.1, KI270466.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "MYNN"
    ## [1] "MYRF"
    ## [1] "NAIF1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000219.1, KI270713.1, KI270908.1, KN196487.1, KI270466.1, KI270467.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "NFAT5"
    ## [1] "NFATC3"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000214.1, GL000219.1, GL000255.2, GL383577.2, KI270713.1, KI270721.1, KI270733.1, KI270744.1, KI270745.1, KI270830.1, KI270862.1, KI270879.1, KN196484.1, KQ090026.1, KV880768.1, ML143377.1, KI270711.1, KQ458383.1, KZ208915.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "NFE2L1"
    ## [1] "NFE2L2"
    ## [1] "NFIA"
    ## [1] "NFIC"
    ## [1] "NFIL3"
    ## [1] "NFKB2"
    ## [1] "NFKBIZ"
    ## [1] "NFYA"
    ## [1] "NFYB"
    ## [1] "NFYC"
    ## [1] "NKX3-1"
    ## [1] "NR0B2"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000208.1, GL000214.1, GL000216.2, GL000218.1, GL000219.1, GL000220.1, GL000224.1, GL000225.1, GL339449.2, GL383566.1, KI270435.1, KI270442.1, KI270589.1, KI270712.1, KI270713.1, KI270721.1, KI270723.1, KI270725.1, KI270728.1, KI270731.1, KI270733.1, KI270744.1, KI270754.1, KI270816.1, KI270830.1, KI270832.1, KI270857.1, KI270869.1, KI270879.1, KI270894.1, KI270908.1, KN538360.1, KN538372.1, KV766196.1, KV880768.1, KZ208915.1, KZ208922.1, KZ559105.1, ML143343.1, ML143367.1, ML143377.1, ML143380.1, GL000008.2, GL000009.2, GL000205.2, GL000251.2, GL000252.2, GL000255.2, GL000256.2, GL383522.1, GL383527.1, GL383533.1, GL383563.3, GL383574.1, GL383577.2, GL383581.2, JH159147.1, KI270714.1, KI270719.1, KI270722.1, KI270724.1, KI270750.1, KI270753.1, KI270772.1, KI270791.1, KI270802.1, KI270821.1, KI270831.1, KI270842.1, KI270847.1, KI270850.1, KI270853.1, KI270856.1, KI270861.1, KI270868.1, KI270878.1, KI270903.1, KI270904.1, KI270907.1, KI270936.1, KN196484.1, KN538364.1, KQ090028.1, KQ458386.1, KV575244.1, KV766193.1, KZ208906.1, KZ208916.1, KZ208921.1, ML143353.1, ML143355.1, ML143360.1, ML143364.1, ML143372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "NR2C2"
    ## [1] "NR2F1"
    ## [1] "NR2F6"
    ## [1] "NR5A1"
    ## [1] "NRF1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrY
    ##   - in 'y': KI270438.1, KI270466.1, KI270467.1, KI270857.1, KN196487.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "NRL"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000009.2, GL000194.1, GL000195.1, GL000208.1, GL000214.1, GL000216.2, GL000218.1, GL000219.1, GL000220.1, GL000224.1, GL000225.1, GL000252.2, GL000253.2, GL000254.2, GL000255.2, GL000256.2, GL339449.2, GL383563.3, GL383574.1, GL383577.2, GL877875.1, JH159146.1, KI270442.1, KI270706.1, KI270709.1, KI270712.1, KI270713.1, KI270714.1, KI270717.1, KI270719.1, KI270721.1, KI270722.1, KI270723.1, KI270724.1, KI270728.1, KI270731.1, KI270732.1, KI270733.1, KI270734.1, KI270742.1, KI270743.1, KI270744.1, KI270745.1, KI270746.1, KI270748.1, KI270750.1, KI270751.1, KI270754.1, KI270772.1, KI270782.1, KI270791.1, KI270810.1, KI270816.1, KI270819.1, KI270830.1, KI270831.1, KI270832.1, KI270845.1, KI270846.1, KI270847.1, KI270850.1, KI270853.1, KI270856.1, KI270857.1, KI270861.1, KI270862.1, KI270868.1, KI270869.1, KI270872.1, KI270878.1, KI270879.1, KI270894.1, KI270903.1, KI270905.1, KI270907.1, KI270908.1, KN196484.1, KN196487.1, KN538364.1, KN538372.1, KQ090020.1, KQ090026.1, KQ458383.1, KQ458386.1, KQ759762.1, KV575244.1, KV766196.1, KV880768.1, KZ208915.1, KZ208918.1, KZ208921.1, KZ208922.1, KZ559105.1, KZ559109.1, ML143372.1, ML143373.1, ML143374.1, ML143377.1, ML143380.1, GL000221.1, GL000251.2, GL383527.1, GL383531.1, GL383542.1, GL383550.2, GL383566.1, JH159147.1, KI270438.1, KI270519.1, KI270711.1, KI270718.1, KI270720.1, KI270725.1, KI270726.1, KI270729.1, KI270735.1, KI270753.1, KI270762.1, KI270795.1, KI270804.1, KI270805.1, KI270824.1, KI270876.1, KI270877.1, KI270901.1, KI270904.1, KQ031389.1, KQ759759.1, KZ208908.1, KZ559100.1, ML143341.1, ML143353.1, ML143364.1, ML143365.1, ML143366.1, ML143367.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ONECUT1"
    ## [1] "ONECUT2"
    ## [1] "PAF1"
    ## [1] "PATZ1"
    ## [1] "PAWR"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': KI270713.1, KI270862.1, KV880768.1, KI270728.1, KI270733.1, KI270846.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "PAXIP1"
    ## [1] "PHF20"
    ## [1] "PHF21A"
    ## [1] "PHF5A"
    ## [1] "PHF8"
    ## [1] "PITX1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000195.1, GL000205.2, GL000214.1, GL000218.1, GL000219.1, GL000220.1, GL000251.2, GL000252.2, GL000253.2, GL000254.2, GL000255.2, GL339449.2, GL383556.1, GL383563.3, GL383574.1, GL383577.2, GL383580.2, KI270442.1, KI270711.1, KI270712.1, KI270713.1, KI270714.1, KI270721.1, KI270726.1, KI270728.1, KI270731.1, KI270733.1, KI270744.1, KI270745.1, KI270751.1, KI270754.1, KI270782.1, KI270803.1, KI270829.1, KI270830.1, KI270832.1, KI270842.1, KI270846.1, KI270847.1, KI270850.1, KI270853.1, KI270857.1, KI270861.1, KI270862.1, KI270869.1, KI270877.1, KI270878.1, KI270879.1, KI270894.1, KI270903.1, KI270904.1, KI270905.1, KI270907.1, KI270908.1, KN196484.1, KN538364.1, KQ090026.1, KQ090028.1, KQ458383.1, KQ458384.1, KQ759762.1, KV575244.1, KV880768.1, KZ208906.1, KZ208908.1, KZ208915.1, KZ208921.1, KZ208922.1, KZ559105.1, ML143367.1, ML143371.1, ML143377.1, ML143380.1, ML143381.1, GL000224.1, GL000256.2, KI270722.1, KI270734.1, KI270816.1, KI270821.1, ML143353.1, ML143366.1, ML143375.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "POGZ"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000218.1, GL000219.1, GL000251.2, GL000252.2, GL000254.2, GL000255.2, GL000257.2, GL339449.2, GL383563.3, GL383574.1, GL383577.2, KI270589.1, KI270712.1, KI270713.1, KI270721.1, KI270728.1, KI270733.1, KI270744.1, KI270745.1, KI270782.1, KI270821.1, KI270829.1, KI270830.1, KI270850.1, KI270853.1, KI270857.1, KI270861.1, KI270862.1, KI270879.1, KI270904.1, KI270905.1, KI270908.1, KN196484.1, KN538364.1, KQ458383.1, KQ458384.1, KV575244.1, KV880768.1, KZ208915.1, KZ208921.1, KZ208922.1, ML143371.1, ML143377.1, ML143380.1, ML143381.1, KI270805.1, KI270832.1, KI270842.1, KZ208906.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "POLR2A"
    ## [1] "POLR2AphosphoS2"
    ## [1] "POLR2AphosphoS5"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000218.1, GL000219.1, GL000220.1, GL000251.2, GL000252.2, GL000253.2, GL000254.2, GL000255.2, GL339449.2, GL383563.3, GL383574.1, GL383577.2, KI270713.1, KI270721.1, KI270744.1, KI270782.1, KI270819.1, KI270831.1, KI270847.1, KI270853.1, KI270857.1, KI270866.1, KI270878.1, KI270879.1, KI270903.1, KI270904.1, KI270905.1, KI270908.1, KQ458383.1, KQ759762.1, KV575244.1, KV880768.1, KZ208915.1, KZ208921.1, ML143377.1, ML143380.1, KI270438.1, KI270816.1, KI270830.1, KI270832.1, KI270850.1, KI270861.1, KI270907.1, KN538364.1, KQ090026.1, KQ090028.1, KV766196.1, KZ208904.1, ML143366.1, ML143367.1, ML143372.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "PPARG"
    ## [1] "PRDM10"
    ## [1] "PRDM15"
    ## [1] "PREB"
    ## [1] "RAD21"
    ## [1] "RARA"
    ## [1] "RARG"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrY
    ##   - in 'y': KI270438.1, KI270713.1, KN196487.1, ML143380.1, GL000194.1, GL000218.1, KI270721.1, KI270730.1, KI270733.1, KI270744.1, KI270869.1, KI270879.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "RBAK"
    ## [1] "RBPJ"
    ## [1] "RCOR1"
    ## [1] "RCOR2"
    ## [1] "RELA"
    ## [1] "RERE"
    ## [1] "REST"
    ## [1] "RFX3"
    ## [1] "RFX5"
    ## [1] "RFXANK"
    ## [1] "RFXAP"
    ## [1] "RREB1"
    ## [1] "RXRA"
    ## [1] "RXRB"
    ## [1] "SAFB2"
    ## [1] "SAP130"
    ## [1] "SATB2"
    ## [1] "SFPQ"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000219.1, GL000251.2, GL000252.2, GL000255.2, GL339449.2, GL383577.2, KI270713.1, KI270721.1, KI270731.1, KI270744.1, KI270782.1, KI270861.1, KI270869.1, KI270879.1, KI270904.1, KI270908.1, KN538364.1, KQ458383.1, KV880768.1, KZ208921.1, ML143377.1, GL000195.1, GL000254.2, GL383563.3, KI270711.1, KI270714.1, KI270809.1, KI270816.1, KI270830.1, KI270831.1, KI270832.1, KI270847.1, KI270850.1, KI270857.1, KI270862.1, KI270866.1, KI270878.1, KI270894.1, KI270903.1, KI270907.1, KN196484.1, KQ458384.1, KQ759762.1, KV575244.1, KZ208906.1, KZ208915.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "SIN3A"
    ## [1] "SIN3B"
    ## [1] "SIX1"
    ## [1] "SIX4"
    ## [1] "SMAD1"
    ## [1] "SMAD3"
    ## [1] "SMAD4"
    ## [1] "SMAD9"
    ## [1] "SMC3"
    ## [1] "SNAI1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000252.2, GL383563.3, KI270713.1, KI270714.1, KI270721.1, KI270728.1, KI270733.1, KI270745.1, KI270846.1, KI270857.1, KI270862.1, KI270879.1, KI270894.1, KI270904.1, ML143380.1, GL000219.1, KI270466.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "SNAPC4"
    ## [1] "SOX13"
    ## [1] "SOX18"
    ## [1] "SOX6"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000218.1, GL000219.1, GL000251.2, GL000252.2, GL000254.2, GL000255.2, GL000257.2, GL339449.2, GL383563.3, GL383574.1, GL383577.2, GL383580.2, GL383581.2, KI270442.1, KI270713.1, KI270714.1, KI270721.1, KI270731.1, KI270734.1, KI270744.1, KI270745.1, KI270782.1, KI270816.1, KI270829.1, KI270830.1, KI270842.1, KI270850.1, KI270853.1, KI270856.1, KI270857.1, KI270861.1, KI270862.1, KI270872.1, KI270875.1, KI270878.1, KI270879.1, KI270894.1, KI270903.1, KI270904.1, KI270905.1, KI270908.1, KN196483.1, KN196484.1, KQ031389.1, KQ458383.1, KQ458384.1, KV575244.1, KV766197.1, KV880768.1, KZ208906.1, KZ208915.1, KZ208921.1, KZ208922.1, KZ559105.1, KZ559109.1, ML143360.1, ML143371.1, ML143372.1, ML143377.1, ML143380.1, ML143381.1, GL000253.2, KI270712.1, KI270827.1, KI270832.1, KN538364.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "SP140L"
    ## [1] "SP5"
    ## [1] "SPEN"
    ## [1] "SRF"
    ## [1] "STAG1"
    ## [1] "STAT5B"
    ## [1] "STAT6"
    ## [1] "SUZ12"
    ## [1] "TAF1"
    ## [1] "TBL1XR1"
    ## [1] "TBP"
    ## [1] "TBX2"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000218.1, GL000219.1, GL000251.2, GL000252.2, GL000254.2, GL339449.2, GL383577.2, KI270442.1, KI270712.1, KI270713.1, KI270721.1, KI270728.1, KI270733.1, KI270744.1, KI270830.1, KI270862.1, KI270869.1, KI270878.1, KI270879.1, KI270894.1, KI270908.1, KV575244.1, KV880768.1, KZ208915.1, KZ208921.1, KZ208922.1, ML143377.1, ML143380.1, GL000214.1, GL000253.2, GL000255.2, GL383556.1, GL383563.3, GL383574.1, KI270714.1, KI270731.1, KI270782.1, KI270816.1, KI270827.1, KI270832.1, KI270850.1, KI270853.1, KI270857.1, KI270861.1, KI270877.1, KI270903.1, KI270905.1, KI270907.1, KN196484.1, KN538364.1, KQ090026.1, KQ458384.1, KQ759762.1, KZ208906.1, KZ559105.1, ML143367.1, ML143371.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "TCF3"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000220.1, GL000225.1, GL383577.2, KI270713.1, KI270723.1, KI270733.1, KI270861.1, KI270879.1, KN196487.1, GL000194.1, GL000216.2, GL000219.1, GL339449.2, KI270442.1, KI270709.1, KI270721.1, KI270744.1, KI270754.1, KI270772.1, KI270857.1, KZ559109.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "TCF7L2"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000225.1, KI270333.1, KI270587.1, KI270729.1, KI270829.1, KI270337.1, KI270466.1, KI270467.1, KI270746.1, KI270750.1, ML143380.1, GL000194.1, GL000218.1, GL000219.1, GL000251.2, GL000252.2, GL000254.2, GL383577.2, KI270712.1, KI270713.1, KI270714.1, KI270721.1, KI270728.1, KI270733.1, KI270744.1, KI270861.1, KI270862.1, KI270869.1, KI270878.1, KI270879.1, KI270894.1, KI270903.1, KI270905.1, KN196484.1, KQ458384.1, KV575244.1, KV880768.1, KZ208915.1, KZ559109.1, ML143371.1, ML143377.1, GL000255.2, GL339449.2, KI270734.1, KI270782.1, KI270830.1, KI270842.1, KI270857.1, KZ208921.1, KZ208922.1, ML143381.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "TEAD1"
    ## [1] "TEAD3"
    ## [1] "TEF"
    ## [1] "TFAP4"
    ## [1] "TFDP1"
    ## [1] "TFDP2"
    ## [1] "TFE3"
    ## [1] "THAP11"
    ## [1] "THAP8"
    ## [1] "THAP9"
    ## [1] "THRA"
    ## [1] "THRB"
    ## [1] "TIGD3"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': KI270713.1, GL000194.1, GL000219.1, GL000252.2, GL339449.2, GL383577.2, GL383580.2, KI270712.1, KI270721.1, KI270733.1, KI270744.1, KI270830.1, KI270869.1, KI270879.1, KI270908.1, KV880768.1, ML143377.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "TIGD6"
    ## [1] "TOE1"
    ## [1] "TP53"
    ## [1] "TRAFD1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000218.1, GL000219.1, GL000251.2, GL000254.2, GL000255.2, GL383563.3, GL383574.1, GL383577.2, GL383580.2, KI270711.1, KI270712.1, KI270713.1, KI270721.1, KI270728.1, KI270731.1, KI270734.1, KI270744.1, KI270830.1, KI270850.1, KI270853.1, KI270857.1, KI270869.1, KI270878.1, KI270879.1, KI270894.1, KI270903.1, KI270905.1, KN196484.1, KQ090026.1, KV880768.1, KZ208915.1, KZ208921.1, ML143377.1, ML143380.1, GL000214.1, GL000220.1, KI270733.1, KI270736.1, KI270753.1, KI270908.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "TSC22D2"
    ## [1] "UBTF"
    ## [1] "USF1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000205.2, GL000218.1, GL000251.2, GL000252.2, GL000255.2, GL383556.1, GL383557.1, KI270733.1, KI270841.1, KI270842.1, KI270878.1, KI270879.1, KI270908.1, KQ759762.1, KV880768.1, ML143377.1, ML143381.1, GL000194.1, GL000219.1, GL000220.1, GL339449.2, GL383555.2, GL383577.2, KI270713.1, KI270734.1, KI270743.1, KI270744.1, KI270830.1, KI270850.1, KI270853.1, KI270857.1, KI270861.1, KI270869.1, KI270872.1, KI270877.1, KQ090026.1, KZ208908.1, KZ208915.1, KZ208921.1, ML143364.1, ML143366.1, GL000214.1, GL000216.2, GL000224.1, GL000225.1, GL000254.2, GL000256.2, KI270438.1, KI270709.1, KI270735.1, KI270745.1, KI270748.1, KI270754.1, KI270782.1, KI270821.1, KI270832.1, KN538372.1, KV575244.1, KZ208913.1, ML143371.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "USF2"
    ## [1] "WIZ"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000218.1, GL000219.1, GL000220.1, GL000224.1, GL000225.1, GL000252.2, KI270442.1, KI270712.1, KI270713.1, KI270721.1, KI270728.1, KI270733.1, KI270744.1, KI270750.1, KI270753.1, KI270782.1, KI270830.1, KI270857.1, KI270869.1, KI270879.1, KV880768.1, KZ208921.1, ML143377.1, GL000255.2, GL383556.1, GL383574.1, GL383577.2, KI270734.1, KI270894.1, KI270908.1, KN196484.1, KZ208915.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "XBP1"
    ## [1] "YEATS2"
    ## [1] "YEATS4"
    ## [1] "YY1"
    ## [1] "ZBED4"
    ## [1] "ZBED5"
    ## [1] "ZBTB1"
    ## [1] "ZBTB10"
    ## [1] "ZBTB14"
    ## [1] "ZBTB21"
    ## [1] "ZBTB24"
    ## [1] "ZBTB25"
    ## [1] "ZBTB26"
    ## [1] "ZBTB38"
    ## [1] "ZBTB39"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000252.2, KI270713.1, KI270861.1, KI270862.1, GL000195.1, GL000251.2, GL000255.2, GL000256.2, GL339449.2, KI270721.1, KI270734.1, KI270857.1, KI270869.1, KI270879.1, KI270908.1, KV880768.1, KZ208906.1, KZ208915.1, KZ208921.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZBTB4"
    ## [1] "ZBTB42"
    ## [1] "ZBTB43"
    ## [1] "ZBTB44"
    ## [1] "ZBTB46"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, KI270712.1, KI270713.1, KI270721.1, KI270728.1, KI270733.1, KI270905.1, KV880768.1, GL000220.1, GL000225.1, KI270744.1, KI270754.1, KN538372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZBTB49"
    ## [1] "ZBTB7A"
    ## [1] "ZBTB7B"
    ## [1] "ZC3H4"
    ## [1] "ZC3H8"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrY, chrM
    ##   - in 'y': KI270713.1, KI270879.1, KI270908.1, GL000194.1, KI270733.1, KI270830.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZCCHC11"
    ## [1] "ZFP1"
    ## [1] "ZFP14"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrY
    ##   - in 'y': GL000194.1, GL000216.2, GL000220.1, GL000224.1, GL000225.1, GL000252.2, KI270438.1, KI270713.1, KI270728.1, KI270730.1, KI270733.1, KI270744.1, KI270869.1, KI270712.1, KI270754.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZFP36L2"
    ## [1] "ZFP37"
    ## [1] "ZFP41"
    ## [1] "ZFP64"
    ## [1] "ZFP82"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000218.1, GL000219.1, KI270466.1, KI270467.1, KI270713.1, KI270721.1, KI270733.1, KI270816.1, KI270830.1, KI270832.1, KI270869.1, KI270879.1, KI270908.1, KN538364.1, KV880768.1, KZ208915.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZFY"
    ## [1] "ZGPAT"
    ## [1] "ZHX3"
    ## [1] "ZIK1"
    ## [1] "ZKSCAN1"
    ## [1] "ZKSCAN5"
    ## [1] "ZKSCAN8"
    ## [1] "ZMYM3"
    ## [1] "ZNF12"
    ## [1] "ZNF124"
    ## [1] "ZNF138"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000219.1, GL000255.2, GL339449.2, KI270712.1, KI270713.1, KI270721.1, KI270728.1, KI270733.1, KI270821.1, KI270830.1, KI270869.1, KI270879.1, KN196484.1, KN538364.1, KV880768.1, KZ208922.1, GL000194.1, GL000218.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF142"
    ## [1] "ZNF143"
    ## [1] "ZNF205"
    ## [1] "ZNF217"
    ## [1] "ZNF224"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000214.1, GL000216.2, GL000219.1, GL000220.1, GL000224.1, GL000225.1, GL000252.2, GL000253.2, KI270442.1, KI270709.1, KI270711.1, KI270712.1, KI270713.1, KI270714.1, KI270726.1, KI270728.1, KI270731.1, KI270733.1, KI270735.1, KI270744.1, KI270750.1, KI270751.1, KI270754.1, KI270869.1, KI270879.1, KI270894.1, KI270903.1, KN196487.1, KQ031389.1, KV880768.1, ML143372.1, ML143377.1, GL383577.2, KI270721.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF225"
    ## [1] "ZNF230"
    ## [1] "ZNF232"
    ## [1] "ZNF234"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrY, chrM
    ##   - in 'y': GL000219.1, KI270733.1, KI270869.1, GL000216.2, GL000220.1, KI270442.1, KI270466.1, KI270713.1, KI270728.1, KI270744.1, KI270850.1, KI270894.1, KI270908.1, KN538372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF25"
    ## [1] "ZNF256"
    ## [1] "ZNF263"
    ## [1] "ZNF264"
    ## [1] "ZNF274"
    ## [1] "ZNF276"
    ## [1] "ZNF280B"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000008.2, GL000009.2, GL000194.1, GL000195.1, GL000205.2, GL000208.1, GL000214.1, GL000216.2, GL000218.1, GL000219.1, GL000220.1, GL000221.1, GL000224.1, GL000225.1, GL000251.2, GL000252.2, GL000253.2, GL000254.2, GL000255.2, GL000256.2, GL000258.2, GL339449.2, GL383530.1, GL383542.1, GL383555.2, GL383557.1, GL383563.3, GL383574.1, GL383577.2, GL383581.2, GL877875.1, GL949742.1, JH159146.1, KI270442.1, KI270706.1, KI270709.1, KI270711.1, KI270712.1, KI270713.1, KI270714.1, KI270717.1, KI270718.1, KI270719.1, KI270721.1, KI270722.1, KI270723.1, KI270724.1, KI270725.1, KI270726.1, KI270728.1, KI270731.1, KI270733.1, KI270734.1, KI270742.1, KI270743.1, KI270744.1, KI270745.1, KI270746.1, KI270748.1, KI270750.1, KI270751.1, KI270753.1, KI270754.1, KI270762.1, KI270765.1, KI270767.1, KI270772.1, KI270779.1, KI270780.1, KI270782.1, KI270792.1, KI270802.1, KI270805.1, KI270809.1, KI270813.1, KI270815.1, KI270816.1, KI270819.1, KI270821.1, KI270824.1, KI270827.1, KI270830.1, KI270831.1, KI270832.1, KI270845.1, KI270846.1, KI270847.1, KI270850.1, KI270853.1, KI270854.1, KI270856.1, KI270857.1, KI270860.1, KI270861.1, KI270862.1, KI270868.1, KI270869.1, KI270870.1, KI270871.1, KI270872.1, KI270878.1, KI270879.1, KI270894.1, KI270897.1, KI270900.1, KI270903.1, KI270904.1, KI270905.1, KI270907.1, KI270908.1, KI270913.1, KI270927.1, KI270934.1, KI270937.1, KN196484.1, KN196487.1, KN538364.1, KN538372.1, KQ031387.1, KQ031389.1, KQ090015.1, KQ090021.1, KQ090026.1, KQ090027.1, KQ458383.1, KQ458384.1, KQ458386.1, KQ759759.1, KQ759762.1, KV575244.1, KV766193.1, KV766196.1, KV880768.1, KZ208904.1, KZ208906.1, KZ208908.1, KZ208914.1, KZ208915.1, KZ208921.1, KZ208922.1, KZ559100.1, KZ559105.1, KZ559109.1, ML143353.1, ML143355.1, ML143358.1, ML143364.1, ML143365.1, ML143366.1, ML143371.1, ML143372.1, ML143373.1, ML143375.1, ML143377.1, ML143380.1, GL000213.1, GL383522.1, GL383526.1, GL383533.1, GL383550.2, JH159147.1, JH159148.1, KI270438.1, KI270519.1, KI270538.1, KI270707.1, KI270720.1, KI270729.1, KI270732.1, KI270735.1, KI270749.1, KI270764.1, KI270787.1, KI270791.1, KI270795.1, KI270797.1, KI270804.1, KI270810.1, KI270818.1, KI270823.1, KI270835.1, KI270842.1, KI270876.1, KI270896.1, KI270902.1, KI270909.1, KI270925.1, KI270936.1, KI270938.1, KN196472.1, KN538360.1, KN538361.1, KQ090020.1, KQ983256.1, KV880764.1, KZ208909.1, KZ208918.1, KZ208919.1, ML143343.1, ML143367.1, ML143374.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF280D"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000218.1, GL000219.1, GL000251.2, GL000254.2, GL000255.2, GL383577.2, GL383578.2, KI270712.1, KI270713.1, KI270721.1, KI270733.1, KI270744.1, KI270853.1, KI270862.1, KI270869.1, KI270878.1, KI270879.1, KI270905.1, KI270908.1, KQ458384.1, KV880768.1, KZ208915.1, KZ208922.1, KZ559105.1, ML143377.1, KI270726.1, KN196484.1, KV575244.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF281"
    ## [1] "ZNF296"
    ## [1] "ZNF30"
    ## [1] "ZNF317"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000218.1, GL000251.2, GL000255.2, GL383563.3, GL383577.2, GL949747.2, KI270713.1, KI270721.1, KI270746.1, KI270829.1, KI270853.1, KI270857.1, KI270860.1, KI270862.1, KI270879.1, KI270905.1, KN196484.1, KZ208906.1, KZ208915.1, GL000219.1, GL000252.2, GL000256.2, GL383574.1, KI270728.1, KI270830.1, KI270861.1, KI270869.1, KI270894.1, KI270903.1, KI270908.1, KN538364.1, KV880768.1, KZ208922.1, ML143341.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF318"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000008.2, GL000194.1, GL000195.1, GL000205.2, GL000214.1, GL000216.2, GL000218.1, GL000219.1, GL000220.1, GL000224.1, GL000225.1, GL000255.2, GL339449.2, GL383563.3, GL383577.2, KI270438.1, KI270442.1, KI270589.1, KI270706.1, KI270712.1, KI270713.1, KI270718.1, KI270719.1, KI270721.1, KI270722.1, KI270726.1, KI270728.1, KI270730.1, KI270731.1, KI270733.1, KI270734.1, KI270736.1, KI270742.1, KI270744.1, KI270745.1, KI270750.1, KI270751.1, KI270754.1, KI270762.1, KI270795.1, KI270805.1, KI270830.1, KI270846.1, KI270847.1, KI270856.1, KI270857.1, KI270861.1, KI270862.1, KI270869.1, KI270878.1, KI270879.1, KI270894.1, KI270903.1, KI270905.1, KI270908.1, KN196484.1, KN196487.1, KN538372.1, KQ031389.1, KQ090026.1, KV766193.1, KV766196.1, KV880768.1, KZ208915.1, KZ208922.1, ML143341.1, ML143353.1, ML143362.1, ML143364.1, ML143377.1, ML143379.1, ML143380.1, GL383574.1, KI270711.1, KI270714.1, KI270723.1, KI270772.1, KI270782.1, KI270816.1, KI270819.1, KI270829.1, KI270853.1, KI270866.1, KI270896.1, KN538364.1, KV880764.1, ML143355.1, ML143358.1, ML143365.1, ML143372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF326"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrY
    ##   - in 'y': GL000216.2, GL000219.1, GL000224.1, GL000225.1, GL383577.2, KI270330.1, KI270438.1, KI270442.1, KI270466.1, KI270467.1, KI270713.1, KI270721.1, KI270729.1, KI270733.1, KI270819.1, KI270847.1, KI270879.1, KN196487.1, KV880768.1, ML143377.1, ML143380.1, GL000214.1, KI270744.1, KI270830.1, KI270831.1, KI270866.1, KN196484.1, KQ031384.1, KQ090026.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF331"
    ## [1] "ZNF333"
    ## [1] "ZNF335"
    ## [1] "ZNF337"
    ## [1] "ZNF33B"
    ## [1] "ZNF34"
    ## [1] "ZNF343"
    ## [1] "ZNF350"
    ## [1] "ZNF362"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000195.1, GL000214.1, GL000218.1, GL000219.1, GL000250.2, GL000251.2, GL000252.2, GL000254.2, GL000255.2, GL000256.2, GL000257.2, GL339449.2, GL383527.1, GL383574.1, GL383577.2, GL383581.2, GL383582.2, KI270713.1, KI270714.1, KI270721.1, KI270733.1, KI270734.1, KI270744.1, KI270782.1, KI270785.1, KI270798.1, KI270806.1, KI270816.1, KI270827.1, KI270829.1, KI270830.1, KI270831.1, KI270832.1, KI270849.1, KI270853.1, KI270856.1, KI270857.1, KI270860.1, KI270861.1, KI270862.1, KI270868.1, KI270869.1, KI270875.1, KI270878.1, KI270879.1, KI270894.1, KI270903.1, KI270904.1, KI270905.1, KI270907.1, KI270908.1, KN196484.1, KN538369.1, KQ458384.1, KQ983256.1, KV575244.1, KV766198.1, KV880768.1, KZ208906.1, KZ208915.1, KZ208921.1, KZ208922.1, ML143353.1, ML143366.1, ML143372.1, ML143375.1, ML143377.1, ML143380.1, GL000220.1, KI270438.1, KI270442.1, KI270712.1, KI270718.1, KI270719.1, KI270728.1, KI270751.1, KI270754.1, KN538372.1, KQ031389.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF383"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': KI270713.1, KI270721.1, KI270879.1, KV880768.1, KI270733.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF384"
    ## [1] "ZNF407"
    ## [1] "ZNF414"
    ## [1] "ZNF44"
    ## [1] "ZNF446"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000219.1, GL000251.2, GL877875.1, KI270442.1, KI270713.1, KI270733.1, KI270830.1, KI270857.1, KI270869.1, KI270879.1, KN538372.1, KV880768.1, KZ208915.1, KZ208921.1, ML143367.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF451"
    ## [1] "ZNF460"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000218.1, GL000224.1, GL000254.2, GL339449.2, KI270442.1, KI270709.1, KI270712.1, KI270728.1, KI270745.1, KI270751.1, KI270832.1, KI270853.1, KI270856.1, KI270879.1, KI270905.1, KN196484.1, KQ458384.1, KV575244.1, KV880768.1, KZ208915.1, KZ208921.1, KZ559109.1, ML143377.1, GL000194.1, GL000219.1, GL000251.2, GL000252.2, GL000253.2, GL000255.2, GL000257.2, GL383556.1, GL383557.1, GL383577.2, GL383581.2, KI270333.1, KI270711.1, KI270713.1, KI270714.1, KI270721.1, KI270744.1, KI270782.1, KI270827.1, KI270849.1, KI270857.1, KI270861.1, KI270876.1, KI270878.1, KI270894.1, KI270908.1, KN196474.1, KQ031389.1, KV880764.1, KZ208904.1, KZ208908.1, KZ208919.1, ML143360.1, ML143374.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF483"
    ## [1] "ZNF485"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000214.1, GL000218.1, GL000219.1, GL000224.1, KI270442.1, KI270712.1, KI270719.1, KI270721.1, KI270728.1, KI270733.1, KI270734.1, KI270744.1, KI270751.1, KI270754.1, KI270762.1, KI270869.1, KI270879.1, KN196487.1, KN538372.1, KQ031389.1, KV880768.1, ML143377.1, GL000009.2, GL000195.1, GL000216.2, GL000220.1, GL000221.1, GL000225.1, GL000251.2, GL000252.2, GL000253.2, GL000254.2, GL000255.2, GL000256.2, GL339449.2, GL383550.2, GL383555.2, GL383574.1, GL383577.2, GL877875.1, KI270438.1, KI270538.1, KI270709.1, KI270713.1, KI270717.1, KI270718.1, KI270722.1, KI270725.1, KI270726.1, KI270742.1, KI270743.1, KI270745.1, KI270746.1, KI270750.1, KI270753.1, KI270772.1, KI270791.1, KI270815.1, KI270816.1, KI270821.1, KI270830.1, KI270850.1, KI270853.1, KI270854.1, KI270857.1, KI270862.1, KI270871.1, KI270872.1, KI270876.1, KI270877.1, KI270878.1, KI270894.1, KI270897.1, KI270903.1, KI270905.1, KI270908.1, KN538364.1, KV575244.1, KV766196.1, KZ208915.1, KZ208921.1, KZ208922.1, KZ559100.1, KZ559109.1, ML143353.1, ML143364.1, ML143366.1, ML143367.1, ML143371.1, ML143372.1, ML143375.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF490"
    ## [1] "ZNF501"
    ## [1] "ZNF503"
    ## [1] "ZNF510"
    ## [1] "ZNF511"
    ## [1] "ZNF512"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000219.1, GL000254.2, GL339449.2, GL383556.1, GL949746.1, KI270442.1, KI270712.1, KI270713.1, KI270725.1, KI270728.1, KI270733.1, KI270734.1, KI270762.1, KI270764.1, KI270827.1, KI270853.1, KI270879.1, KI270908.1, KV880764.1, KZ208915.1, ML143365.1, ML143372.1, ML143377.1, GL000214.1, GL000218.1, GL000220.1, GL000224.1, KI270330.1, KI270717.1, KI270718.1, KI270719.1, KI270722.1, KI270726.1, KI270744.1, KI270754.1, KI270821.1, KI270850.1, KI270857.1, KI270894.1, KI270905.1, KQ090027.1, KV880768.1, ML143367.1, ML143380.1, GL000194.1, GL000251.2, GL000252.2, GL000255.2, GL383577.2, GL383578.2, GL383580.2, KI270711.1, KI270714.1, KI270721.1, KI270782.1, KI270830.1, KI270842.1, KI270862.1, KI270869.1, KI270878.1, KN538364.1, KZ208906.1, KZ208921.1, KZ208922.1, KZ559105.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF512B"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000214.1, GL000219.1, GL000220.1, GL000251.2, GL000252.2, GL000254.2, GL000255.2, GL339449.2, GL383577.2, GL383578.2, GL383580.2, KI270711.1, KI270712.1, KI270713.1, KI270714.1, KI270721.1, KI270726.1, KI270733.1, KI270734.1, KI270744.1, KI270782.1, KI270821.1, KI270830.1, KI270842.1, KI270850.1, KI270857.1, KI270862.1, KI270869.1, KI270878.1, KI270879.1, KI270894.1, KI270905.1, KI270908.1, KN538364.1, KV880764.1, KV880768.1, KZ208906.1, KZ208915.1, KZ208921.1, KZ208922.1, KZ559105.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF513"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrY
    ##   - in 'y': GL000219.1, KI270438.1, KI270713.1, KI270733.1, KI270879.1, KN196487.1, KV880768.1, GL000220.1, GL000251.2, KI270712.1, KI270728.1, KI270744.1, KI270830.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF526"
    ## [1] "ZNF543"
    ## [1] "ZNF547"
    ## [1] "ZNF548"
    ## [1] "ZNF550"
    ## [1] "ZNF556"
    ## [1] "ZNF557"
    ## [1] "ZNF572"
    ## [1] "ZNF574"
    ## [1] "ZNF576"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000219.1, GL383563.3, KI270712.1, KI270713.1, KI270714.1, KI270721.1, KI270728.1, KI270879.1, KV880768.1, KZ208909.1, KZ208915.1, GL000194.1, GL000214.1, GL000254.2, GL000255.2, GL383533.1, GL383555.2, GL383574.1, KI270442.1, KI270734.1, KI270744.1, KI270750.1, KI270772.1, KI270802.1, KI270857.1, KI270894.1, KI270903.1, KI270908.1, KV766193.1, KZ208908.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF580"
    ## [1] "ZNF598"
    ## [1] "ZNF607"
    ## [1] "ZNF608"
    ## [1] "ZNF609"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000219.1, GL000254.2, GL383577.2, KI270712.1, KI270713.1, KI270721.1, KI270728.1, KI270733.1, KI270861.1, KI270862.1, KI270878.1, KI270879.1, KV880768.1, KZ208915.1, GL000194.1, GL000218.1, GL000220.1, GL000251.2, GL000252.2, GL000255.2, GL339449.2, GL383563.3, GL383574.1, KI270442.1, KI270714.1, KI270734.1, KI270744.1, KI270745.1, KI270754.1, KI270782.1, KI270830.1, KI270832.1, KI270847.1, KI270853.1, KI270857.1, KI270869.1, KI270872.1, KI270894.1, KI270903.1, KI270904.1, KI270905.1, KI270908.1, KN196484.1, KN538364.1, KN538372.1, KQ031389.1, KQ090026.1, KQ458384.1, KV575244.1, KZ208918.1, KZ208921.1, KZ208922.1, ML143371.1, ML143377.1, ML143380.1, ML143381.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF614"
    ## [1] "ZNF616"
    ## [1] "ZNF619"
    ## [1] "ZNF629"
    ## [1] "ZNF639"
    ## [1] "ZNF644"
    ## [1] "ZNF652"
    ## [1] "ZNF660"
    ## [1] "ZNF674"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000220.1, KI270713.1, KI270728.1, KI270733.1, GL000194.1, GL000219.1, KI270712.1, KI270869.1, KN196487.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF678"
    ## [1] "ZNF687"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000195.1, GL000214.1, GL000218.1, GL000219.1, GL000220.1, GL000251.2, GL000252.2, GL000253.2, GL000254.2, GL000255.2, GL000256.2, GL000257.2, GL000258.2, GL339449.2, GL383526.1, GL383555.2, GL383556.1, GL383557.1, GL383563.3, GL383566.1, GL383574.1, GL383577.2, GL383580.2, GL949752.1, JH159146.1, KI270706.1, KI270712.1, KI270713.1, KI270714.1, KI270717.1, KI270719.1, KI270721.1, KI270723.1, KI270726.1, KI270728.1, KI270731.1, KI270733.1, KI270734.1, KI270744.1, KI270745.1, KI270748.1, KI270751.1, KI270753.1, KI270754.1, KI270762.1, KI270772.1, KI270780.1, KI270782.1, KI270804.1, KI270809.1, KI270810.1, KI270813.1, KI270816.1, KI270819.1, KI270821.1, KI270827.1, KI270830.1, KI270831.1, KI270832.1, KI270842.1, KI270846.1, KI270847.1, KI270848.1, KI270850.1, KI270853.1, KI270854.1, KI270856.1, KI270857.1, KI270860.1, KI270861.1, KI270862.1, KI270868.1, KI270869.1, KI270872.1, KI270877.1, KI270878.1, KI270879.1, KI270894.1, KI270897.1, KI270903.1, KI270904.1, KI270905.1, KI270907.1, KI270908.1, KI270910.1, KI270913.1, KI270927.1, KI270938.1, KN196479.1, KN196484.1, KN538364.1, KN538372.1, KQ090026.1, KQ090028.1, KQ458383.1, KQ458384.1, KQ458386.1, KQ759762.1, KV575244.1, KV766196.1, KV880768.1, KZ208906.1, KZ208908.1, KZ208914.1, KZ208915.1, KZ208918.1, KZ208921.1, KZ208922.1, KZ559105.1, KZ559109.1, ML143353.1, ML143358.1, ML143364.1, ML143366.1, ML143367.1, ML143371.1, ML143372.1, ML143377.1, ML143380.1, GL383533.1, GL383550.2, GL383581.2, KI270709.1, KI270765.1, KI270788.1, KI270791.1, KI270792.1, KI270805.1, KI270824.1, KI270845.1, KI270849.1, KI270871.1, KN196472.1, KQ031389.1, KZ208904.1, KZ559113.1, ML143365.1, ML143378.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF691"
    ## [1] "ZNF697"
    ## [1] "ZNF709"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000216.2, GL000219.1, GL000251.2, GL339449.2, GL383563.3, GL383577.2, KI270442.1, KI270712.1, KI270713.1, KI270714.1, KI270725.1, KI270726.1, KI270728.1, KI270730.1, KI270733.1, KI270742.1, KI270744.1, KI270830.1, KI270832.1, KI270857.1, KI270869.1, KI270879.1, KN196484.1, KN538364.1, KQ759762.1, KV575244.1, KV880768.1, KZ208915.1, KZ208921.1, KI270721.1, KI270861.1, KI270908.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF710"
    ## [1] "ZNF713"
    ## [1] "ZNF720"
    ## [1] "ZNF737"
    ## [1] "ZNF740"
    ## [1] "ZNF761"
    ## [1] "ZNF766"
    ## [1] "ZNF768"
    ## [1] "ZNF772"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000195.1, GL000251.2, GL000252.2, GL000254.2, GL339449.2, GL383577.2, KI270713.1, KI270721.1, KI270733.1, KI270861.1, KI270862.1, KI270879.1, KN538372.1, KV880768.1, KZ208921.1, ML143377.1, GL000194.1, GL000218.1, GL000219.1, GL000224.1, GL383574.1, KI270712.1, KI270714.1, KI270728.1, KI270744.1, KI270780.1, KI270816.1, KI270821.1, KI270830.1, KI270832.1, KI270857.1, KI270869.1, KI270894.1, KI270905.1, KI270908.1, KN196484.1, KN538364.1, KQ031389.1, KQ458384.1, KV766196.1, KZ208915.1, KZ208922.1, ML143372.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF773"
    ## [1] "ZNF775"
    ## [1] "ZNF776"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000009.2, GL000194.1, GL000214.1, GL000218.1, GL000219.1, GL000220.1, GL000224.1, GL000225.1, GL000254.2, GL000255.2, GL339449.2, GL383563.3, GL383574.1, KI270442.1, KI270706.1, KI270712.1, KI270713.1, KI270714.1, KI270719.1, KI270721.1, KI270722.1, KI270723.1, KI270728.1, KI270731.1, KI270733.1, KI270734.1, KI270736.1, KI270744.1, KI270745.1, KI270750.1, KI270754.1, KI270772.1, KI270816.1, KI270830.1, KI270832.1, KI270850.1, KI270854.1, KI270860.1, KI270861.1, KI270869.1, KI270879.1, KI270894.1, KI270901.1, KI270905.1, KI270908.1, KN196484.1, KN538364.1, KN538372.1, KV575244.1, KV766196.1, KV880768.1, KZ208915.1, KZ208921.1, KZ208922.1, ML143360.1, ML143364.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF777"
    ## [1] "ZNF781"
    ## [1] "ZNF782"
    ## [1] "ZNF784"
    ## [1] "ZNF786"
    ## [1] "ZNF788"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000195.1, GL000205.2, GL000214.1, GL000216.2, GL000218.1, GL000219.1, GL000224.1, GL000251.2, GL000252.2, GL000253.2, GL000254.2, GL000255.2, GL339449.2, GL383563.3, GL383574.1, GL383577.2, GL383580.2, KI270442.1, KI270711.1, KI270712.1, KI270713.1, KI270714.1, KI270719.1, KI270721.1, KI270728.1, KI270731.1, KI270733.1, KI270734.1, KI270744.1, KI270772.1, KI270780.1, KI270782.1, KI270816.1, KI270821.1, KI270827.1, KI270829.1, KI270830.1, KI270832.1, KI270842.1, KI270850.1, KI270853.1, KI270857.1, KI270861.1, KI270862.1, KI270869.1, KI270875.1, KI270877.1, KI270878.1, KI270879.1, KI270894.1, KI270903.1, KI270905.1, KI270907.1, KI270908.1, KN196484.1, KN538364.1, KN538372.1, KQ090026.1, KQ458384.1, KQ458386.1, KQ759762.1, KV575244.1, KV880768.1, KZ208906.1, KZ208915.1, KZ208921.1, KZ208922.1, KZ559109.1, ML143358.1, ML143366.1, ML143367.1, ML143371.1, ML143372.1, ML143377.1, ML143380.1, GL000225.1, KI270718.1, KI270754.1, KQ031389.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF792"
    ## [1] "ZNF800"
    ## [1] "ZNF83"
    ## [1] "ZNF839"
    ## [1] "ZNF883"
    ## [1] "ZNF891"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000008.2, GL000009.2, GL000194.1, GL000195.1, GL000208.1, GL000214.1, GL000216.2, GL000218.1, GL000219.1, GL000220.1, GL000221.1, GL000224.1, GL000225.1, GL000251.2, GL000252.2, GL000253.2, GL000254.2, GL000255.2, GL000256.2, GL339449.2, GL383522.1, GL383533.1, GL383556.1, GL383563.3, GL383566.1, GL383574.1, GL383577.2, GL383581.2, GL877875.1, GL949742.1, JH159146.1, JH159147.1, JH159148.1, KI270438.1, KI270442.1, KI270538.1, KI270706.1, KI270707.1, KI270709.1, KI270712.1, KI270713.1, KI270714.1, KI270718.1, KI270719.1, KI270720.1, KI270721.1, KI270722.1, KI270723.1, KI270724.1, KI270726.1, KI270728.1, KI270731.1, KI270732.1, KI270733.1, KI270734.1, KI270742.1, KI270743.1, KI270744.1, KI270745.1, KI270746.1, KI270750.1, KI270751.1, KI270753.1, KI270754.1, KI270762.1, KI270772.1, KI270782.1, KI270804.1, KI270805.1, KI270810.1, KI270816.1, KI270819.1, KI270821.1, KI270823.1, KI270824.1, KI270830.1, KI270831.1, KI270832.1, KI270845.1, KI270846.1, KI270847.1, KI270850.1, KI270853.1, KI270854.1, KI270857.1, KI270860.1, KI270861.1, KI270869.1, KI270871.1, KI270872.1, KI270876.1, KI270877.1, KI270878.1, KI270879.1, KI270894.1, KI270896.1, KI270897.1, KI270903.1, KI270904.1, KI270905.1, KI270907.1, KI270908.1, KI270910.1, KI270913.1, KN196484.1, KN196487.1, KN538364.1, KN538372.1, KQ031389.1, KQ090026.1, KQ458386.1, KQ759762.1, KV575244.1, KV766196.1, KV880764.1, KV880768.1, KZ208904.1, KZ208906.1, KZ208908.1, KZ208915.1, KZ208918.1, KZ208921.1, KZ208922.1, KZ559109.1, ML143353.1, ML143364.1, ML143365.1, ML143366.1, ML143367.1, ML143371.1, ML143372.1, ML143377.1, ML143378.1, ML143380.1, GL000205.2, KI270711.1, KI270725.1, KI270764.1, KI270862.1, KQ458383.1, KV766197.1, ML143341.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZSCAN20"
    ## [1] "ZSCAN22"
    ## [1] "ZSCAN29"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000218.1, GL000251.2, GL000252.2, GL000254.2, GL000255.2, GL339449.2, GL383574.1, GL383577.2, KI270713.1, KI270721.1, KI270733.1, KI270829.1, KI270857.1, KI270879.1, KI270894.1, KI270905.1, KN196484.1, KV575244.1, KV880768.1, KZ208915.1, GL000194.1, GL000253.2, KI270908.1, KZ208922.1, ML143350.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZSCAN31"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000214.1, GL000218.1, GL000219.1, GL000220.1, GL000224.1, GL000225.1, GL000251.2, GL000252.2, GL000253.2, GL000254.2, GL000255.2, GL339449.2, GL383574.1, GL383577.2, GL383580.2, KI270713.1, KI270714.1, KI270721.1, KI270728.1, KI270733.1, KI270744.1, KI270782.1, KI270816.1, KI270830.1, KI270861.1, KI270869.1, KI270879.1, KI270894.1, KI270905.1, KI270908.1, KN196483.1, KN196484.1, KV880768.1, KZ208915.1, KZ208922.1, ML143377.1, ML143380.1, KI270904.1, KZ208921.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZSCAN9"
    ## [1] "ZXDC"
    ## [1] "ZZZ3"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000218.1, GL000219.1, GL000220.1, GL000224.1, GL000225.1, GL339449.2, GL383563.3, KI270589.1, KI270712.1, KI270713.1, KI270726.1, KI270728.1, KI270733.1, KI270742.1, KI270744.1, KI270762.1, KI270830.1, KI270857.1, KI270869.1, KI270894.1, KN538372.1, KQ031389.1, KZ208915.1, ML143353.1, ML143362.1, ML143372.1, GL000195.1, GL000216.2, GL000251.2, GL000255.2, GL383577.2, KI270438.1, KI270719.1, KI270721.1, KI270754.1, KI270756.1, KI270853.1, KI270862.1, KI270879.1, KI270904.1, KI270905.1, KV880768.1, ML143377.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

``` r
# non super binder meta_df
non_superbinder_metaplot_df <- data.frame(x = integer(), dens = numeric(), dbp = character())
lncrna_mrna_promoters=non_super_lncrna_mrna_gr
# for loop to populate mRNA_metaplot
for(i in 1:length(filtered_consensus_list)) {
 print(names(filtered_consensus_list)[[i]])
  tmp_df <- profile_tss(filtered_consensus_list[[i]], lncrna_mrna_promoters)
  tmp_df$dbp <- names(filtered_consensus_list)[[i]]
  non_superbinder_metaplot_df <- bind_rows(non_superbinder_metaplot_df, tmp_df)
}
```

    ## [1] "ADNP"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrY
    ##   - in 'y': GL000214.1, GL000218.1, GL000220.1, GL000252.2, KI270438.1, KI270442.1, KI270467.1, KI270733.1, KI270767.1, KI270857.1, KI270878.1, KN196487.1, KQ031389.1, KV575244.1, ML143377.1, ML143380.1, GL000194.1, GL000216.2, GL000219.1, GL000224.1, GL000225.1, GL383577.2, KI270712.1, KI270728.1, KI270729.1, KI270861.1, KQ458384.1, KZ208915.1, ML143372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "AFF4"
    ## [1] "AHDC1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000218.1, GL000219.1, GL000251.2, GL000252.2, GL000254.2, GL339449.2, GL383527.1, GL383556.1, GL383563.3, GL383577.2, GL877875.1, KI270467.1, KI270713.1, KI270721.1, KI270744.1, KI270782.1, KI270803.1, KI270829.1, KI270857.1, KI270861.1, KI270862.1, KI270878.1, KI270879.1, KI270908.1, KN196484.1, KQ458384.1, KV575244.1, KV880768.1, KZ208915.1, ML143371.1, ML143377.1, ML143381.1, GL000224.1, KI270712.1, KI270728.1, KI270733.1, KI270830.1, KI270853.1, KI270869.1, KN538372.1, KV766196.1, KV880764.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "AHR"
    ## [1] "AKAP8"
    ## [1] "ARID3A"
    ## [1] "ARID4A"
    ## [1] "ARID4B"
    ## [1] "ARID5B"
    ## [1] "ARNTL"
    ## [1] "ASH2L"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000219.1, GL383577.2, GL383581.2, KI270713.1, KI270729.1, KI270879.1, KI270904.1, KQ458384.1, KV880768.1, KZ208922.1, GL000218.1, GL000251.2, KZ208915.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ATAD3A"
    ## [1] "ATF2"
    ## [1] "ATF3"
    ## [1] "ATF5"
    ## [1] "ATF6"
    ## [1] "ATF7"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000009.2, GL000194.1, GL000208.1, GL000214.1, GL000216.2, GL000218.1, GL000219.1, GL000220.1, GL000221.1, GL000224.1, GL000225.1, GL000252.2, GL000255.2, GL000256.2, GL339449.2, GL383522.1, GL383542.1, GL383563.3, GL383577.2, GL383581.2, JH159146.1, KI270442.1, KI270706.1, KI270711.1, KI270712.1, KI270713.1, KI270714.1, KI270718.1, KI270719.1, KI270721.1, KI270722.1, KI270723.1, KI270725.1, KI270726.1, KI270728.1, KI270731.1, KI270732.1, KI270733.1, KI270734.1, KI270736.1, KI270742.1, KI270743.1, KI270744.1, KI270745.1, KI270750.1, KI270751.1, KI270754.1, KI270762.1, KI270764.1, KI270772.1, KI270805.1, KI270819.1, KI270830.1, KI270832.1, KI270846.1, KI270847.1, KI270850.1, KI270853.1, KI270861.1, KI270862.1, KI270865.1, KI270868.1, KI270869.1, KI270879.1, KI270894.1, KI270896.1, KI270903.1, KI270905.1, KI270907.1, KI270908.1, KI270913.1, KI270936.1, KN196484.1, KN196487.1, KN538372.1, KQ031389.1, KQ090026.1, KQ090028.1, KQ759762.1, KV575244.1, KV766193.1, KV880764.1, KV880768.1, KZ208908.1, KZ208915.1, KZ208921.1, KZ208922.1, KZ559100.1, ML143353.1, ML143364.1, ML143365.1, ML143371.1, ML143372.1, ML143375.1, ML143377.1, ML143380.1, GL000008.2, GL000195.1, GL000251.2, GL000254.2, GL383556.1, KI270438.1, KI270866.1, KI270878.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "BCL3"
    ## [1] "BCL6"
    ## [1] "BHLHE40"
    ## [1] "BRCA1"
    ## [1] "CAMTA2"
    ## [1] "CBFB"
    ## [1] "CBX5"
    ## [1] "CEBPA"
    ## [1] "CEBPB"
    ## [1] "CEBPD"
    ## [1] "CEBPG"
    ## [1] "CERS6"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrY
    ##   - in 'y': GL000224.1, GL877875.1, KI270438.1, KI270712.1, KI270713.1, KI270733.1, KI270744.1, KN196487.1, KV880768.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "CHD2"
    ## [1] "CREB3"
    ## [1] "CREM"
    ## [1] "CTCF"
    ## [1] "DBP"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrY
    ##   - in 'y': KI270438.1, KI270466.1, KI270467.1, KI270713.1, KI270830.1, KI270869.1, KI270879.1, KN196487.1, KV880768.1, ML143380.1, KI270712.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "DDIT3"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrY
    ##   - in 'y': GL000214.1, GL000216.2, GL000220.1, GL000224.1, GL000225.1, KI270330.1, KI270333.1, KI270336.1, KI270337.1, KI270438.1, KI270442.1, KI270466.1, KI270467.1, KI270713.1, KI270729.1, KI270733.1, KI270879.1, KN196487.1, KV880768.1, ML143380.1, GL383577.2
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "DLX6"
    ## [1] "DMAP1"
    ## [1] "DMTF1"
    ## [1] "DNMT1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000009.2, GL000194.1, GL000205.2, GL000213.1, GL000214.1, GL000216.2, GL000218.1, GL000219.1, GL000220.1, GL000221.1, GL000224.1, GL000225.1, GL000252.2, GL000253.2, GL000254.2, GL000255.2, GL000256.2, GL339449.2, GL383563.3, GL383566.1, GL383574.1, GL383577.2, JH159146.1, JH159147.1, KI270442.1, KI270512.1, KI270519.1, KI270538.1, KI270706.1, KI270709.1, KI270711.1, KI270712.1, KI270713.1, KI270714.1, KI270718.1, KI270719.1, KI270721.1, KI270722.1, KI270723.1, KI270724.1, KI270725.1, KI270728.1, KI270731.1, KI270732.1, KI270733.1, KI270734.1, KI270735.1, KI270736.1, KI270742.1, KI270743.1, KI270744.1, KI270745.1, KI270749.1, KI270750.1, KI270751.1, KI270753.1, KI270754.1, KI270762.1, KI270767.1, KI270772.1, KI270780.1, KI270782.1, KI270791.1, KI270792.1, KI270795.1, KI270797.1, KI270805.1, KI270815.1, KI270816.1, KI270819.1, KI270821.1, KI270830.1, KI270831.1, KI270832.1, KI270846.1, KI270847.1, KI270850.1, KI270853.1, KI270854.1, KI270856.1, KI270857.1, KI270860.1, KI270861.1, KI270862.1, KI270865.1, KI270868.1, KI270869.1, KI270870.1, KI270871.1, KI270872.1, KI270875.1, KI270877.1, KI270878.1, KI270879.1, KI270894.1, KI270896.1, KI270901.1, KI270903.1, KI270905.1, KI270908.1, KI270913.1, KI270937.1, KN196484.1, KN538361.1, KN538364.1, KN538368.1, KN538372.1, KN538373.1, KQ031389.1, KQ090026.1, KQ090028.1, KQ458383.1, KQ458386.1, KQ759759.1, KQ759762.1, KV575244.1, KV766196.1, KV880764.1, KV880768.1, KZ208904.1, KZ208906.1, KZ208915.1, KZ208918.1, KZ208919.1, KZ208921.1, KZ208922.1, KZ559100.1, KZ559109.1, ML143341.1, ML143353.1, ML143355.1, ML143358.1, ML143364.1, ML143365.1, ML143366.1, ML143371.1, ML143372.1, ML143375.1, ML143377.1, ML143378.1, ML143379.1, ML143380.1, GL000258.2, KI270589.1, KI270876.1, KI270907.1, KI270936.1, KQ090022.1, KZ208908.1, KZ208916.1, ML143370.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "DPF2"
    ## [1] "DRAP1"
    ## [1] "DZIP1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000214.1, GL000218.1, GL000219.1, GL000225.1, GL383577.2, KI270713.1, KI270721.1, KI270733.1, KI270869.1, KI270879.1, KI270908.1, KZ208915.1, GL000008.2, GL000009.2, GL000194.1, GL000208.1, GL000216.2, GL000220.1, GL000221.1, GL000224.1, GL000251.2, GL000252.2, GL000253.2, GL000254.2, GL000256.2, GL339449.2, GL383522.1, GL383533.1, GL383563.3, GL383566.1, GL383574.1, GL877875.1, KI270330.1, KI270442.1, KI270706.1, KI270709.1, KI270711.1, KI270712.1, KI270714.1, KI270718.1, KI270719.1, KI270720.1, KI270722.1, KI270723.1, KI270724.1, KI270725.1, KI270728.1, KI270731.1, KI270732.1, KI270734.1, KI270742.1, KI270743.1, KI270744.1, KI270745.1, KI270746.1, KI270750.1, KI270751.1, KI270753.1, KI270754.1, KI270762.1, KI270764.1, KI270772.1, KI270787.1, KI270791.1, KI270805.1, KI270819.1, KI270821.1, KI270824.1, KI270830.1, KI270831.1, KI270832.1, KI270845.1, KI270847.1, KI270853.1, KI270854.1, KI270857.1, KI270861.1, KI270868.1, KI270871.1, KI270872.1, KI270877.1, KI270878.1, KI270894.1, KI270902.1, KI270903.1, KI270905.1, KI270907.1, KN196484.1, KN196487.1, KN538372.1, KQ031389.1, KQ090020.1, KQ090026.1, KQ458386.1, KV575244.1, KV766196.1, KV880768.1, KZ208904.1, KZ208906.1, KZ208918.1, KZ208921.1, KZ208922.1, KZ559100.1, ML143353.1, ML143364.1, ML143365.1, ML143366.1, ML143367.1, ML143371.1, ML143372.1, ML143373.1, ML143374.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "E2F1"
    ## [1] "E2F2"
    ## [1] "E2F4"
    ## [1] "E2F5"
    ## [1] "EEA1"
    ## [1] "EED"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000219.1, GL000224.1, KI270712.1, KI270713.1, KI270721.1, KI270728.1, KI270733.1, KI270744.1, KI270869.1, KI270879.1, KI270894.1, KI270908.1, KZ208915.1, GL000009.2, GL000214.1, GL000216.2, GL000218.1, GL000225.1, GL000252.2, GL000255.2, GL339449.2, GL383533.1, GL383557.1, GL877875.1, JH159146.1, KI270442.1, KI270709.1, KI270711.1, KI270714.1, KI270718.1, KI270719.1, KI270722.1, KI270725.1, KI270734.1, KI270745.1, KI270750.1, KI270754.1, KI270772.1, KI270830.1, KI270847.1, KI270871.1, KI270877.1, KI270878.1, KI270903.1, KI270907.1, KN196484.1, KN538372.1, KQ090026.1, KV880768.1, KZ208908.1, KZ208921.1, KZ208922.1, ML143345.1, ML143364.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "EGR1"
    ## [1] "ELF1"
    ## [1] "ELF3"
    ## [1] "ELF4"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrY, chrM
    ##   - in 'y': GL000194.1, GL000219.1, GL339449.2, GL383577.2, KI270713.1, KI270879.1, GL000254.2
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ELK1"
    ## [1] "ERF"
    ## [1] "ESRRA"
    ## [1] "ETS1"
    ## [1] "ETV5"
    ## [1] "ETV6"
    ## [1] "EZH2"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrY, chrM
    ##   - in 'y': KI270712.1, KI270729.1, KI270824.1, KI270830.1, KI270894.1, KQ090026.1, ML143364.1, KI270733.1, KI270845.1, KI270871.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "FOSL2"
    ## [1] "FOXA1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000218.1, GL000219.1, GL000251.2, GL000254.2, GL000256.2, GL000257.2, GL383556.1, GL383577.2, GL383580.2, KI270713.1, KI270721.1, KI270744.1, KI270772.1, KI270827.1, KI270830.1, KI270853.1, KI270857.1, KI270861.1, KI270862.1, KI270878.1, KI270879.1, KI270904.1, KI270905.1, KQ031389.1, KQ090028.1, KQ458383.1, KV880768.1, KZ208915.1, KZ208918.1, KZ208919.1, ML143343.1, ML143381.1, GL339449.2, GL383527.1, GL383555.2, GL949742.1, KI270734.1, KI270801.1, KI270908.1, KN196474.1, KN196483.1, KN196484.1, KQ458384.1, KZ208922.1, KZ559105.1, ML143360.1, ML143366.1, ML143371.1, GL000253.2, KI270803.1, KN196487.1, GL383552.1, KI270733.1, KI270751.1, KI270850.1, KV575244.1, GL000214.1, GL000216.2, GL000224.1, GL000225.1, KI270330.1, KI270438.1, KI270712.1, KI270730.1, KI270754.1, ML143372.1, ML143377.1, ML143380.1, GL000194.1, KI270442.1, KI270709.1, KI270723.1, KI270735.1, KI270745.1, KI270782.1, KI270821.1, KI270896.1, KQ031384.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "FOXA2"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000205.2, GL000219.1, GL000251.2, GL000254.2, GL000257.2, GL339449.2, GL383526.1, GL383552.1, GL383556.1, GL383563.3, GL383571.1, GL383577.2, GL383580.2, GL383581.2, KI270713.1, KI270721.1, KI270723.1, KI270734.1, KI270744.1, KI270782.1, KI270801.1, KI270822.1, KI270827.1, KI270830.1, KI270853.1, KI270857.1, KI270861.1, KI270862.1, KI270878.1, KI270879.1, KI270896.1, KI270904.1, KI270905.1, KN196484.1, KQ031389.1, KQ090022.1, KQ090026.1, KQ090028.1, KQ458383.1, KV575244.1, KV880768.1, KZ208915.1, KZ208918.1, KZ208919.1, KZ559105.1, ML143343.1, ML143360.1, ML143371.1, ML143377.1, ML143381.1, KI270333.1, KI270337.1, KI270466.1, KI270467.1, KN196487.1, ML143380.1, GL000214.1, GL000216.2, GL000224.1, GL000225.1, GL000252.2, KI270330.1, KI270438.1, KI270709.1, KI270729.1, KI270730.1, KI270735.1, KI270751.1, KI270754.1, ML143366.1, GL000194.1, GL000195.1, GL000218.1, GL000253.2, GL000255.2, KI270442.1, KI270519.1, KI270720.1, KI270733.1, KI270803.1, KI270908.1, KQ458384.1, KQ759762.1, KZ208922.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "FOXA3"
    ## [1] "FOXJ3"
    ## [1] "FOXK1"
    ## [1] "FOXM1"
    ## [1] "FOXO1"
    ## [1] "FOXP1"
    ## [1] "FOXP4"
    ## [1] "FUBP1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000219.1, GL000252.2, GL000255.2, GL339449.2, GL383563.3, GL383577.2, GL383580.2, KI270713.1, KI270721.1, KI270734.1, KI270736.1, KI270744.1, KI270773.1, KI270830.1, KI270850.1, KI270857.1, KI270861.1, KI270862.1, KI270869.1, KI270878.1, KI270879.1, KI270905.1, KQ090026.1, KQ458383.1, KQ458384.1, KV880768.1, ML143365.1, ML143377.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "FUBP3"
    ## [1] "GABPA"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000219.1, GL000251.2, GL000255.2, GL339449.2, KI270816.1, KI270830.1, KI270879.1, KV880768.1, ML143371.1, ML143377.1, ML143380.1, GL000194.1, GL000218.1, GL000254.2, KZ208915.1, ML143364.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "GABPB1"
    ## [1] "GATA2"
    ## [1] "GATAD1"
    ## [1] "GATAD2A"
    ## [1] "GLI4"
    ## [1] "GMEB1"
    ## [1] "GMEB2"
    ## [1] "GPN1"
    ## [1] "GTF2F1"
    ## [1] "GTF3A"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': KI270713.1, KI270733.1, KI270879.1, KZ208906.1, GL000194.1, GL000251.2, GL339449.2, GL383577.2, KI270712.1, KI270721.1, KI270728.1, KI270905.1, KN196484.1, KV880768.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "GZF1"
    ## [1] "H2AFZ"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': KI270442.1, KI270729.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "H3K27ac"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000218.1, GL000251.2, GL000254.2, GL000255.2, GL000257.2, GL339449.2, GL383563.3, GL383574.1, KI270713.1, KI270721.1, KI270729.1, KI270744.1, KI270782.1, KI270857.1, KI270878.1, KI270879.1, KI270903.1, KI270904.1, KN196484.1, KQ458384.1, KV575244.1, KV880768.1, KZ208915.1, ML143366.1, ML143372.1, ML143377.1, ML143380.1, GL000194.1, GL000219.1, GL383577.2, GL383581.2, KI270861.1, KQ759762.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "H3K36me3"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000225.1, GL383557.1, KI270442.1, KI270711.1, KI270729.1, KI270736.1, KI270801.1, KI270853.1, KI270866.1, KI270879.1, KQ090028.1, KV575244.1, KZ559109.1, GL000218.1, GL000219.1, KI270467.1, KI270728.1, GL383522.1, KN196474.1, KV766194.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "H3K4me1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000218.1, GL000219.1, GL000257.2, GL383574.1, GL383577.2, GL383581.2, KI270442.1, KI270713.1, KI270729.1, KI270744.1, KI270829.1, KI270830.1, KI270847.1, KI270856.1, KI270861.1, KI270877.1, KI270878.1, KI270879.1, KI270903.1, KN196484.1, KQ090026.1, KQ458383.1, KQ458384.1, KV880764.1, KZ208919.1, KZ208922.1, ML143360.1, ML143366.1, ML143377.1, ML143380.1, KZ208915.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "H3K4me2"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000205.2, GL000214.1, GL000218.1, GL000219.1, GL000220.1, GL000251.2, GL000252.2, GL000253.2, GL000254.2, GL000255.2, GL000256.2, GL339449.2, GL383552.1, GL383563.3, GL383574.1, GL383577.2, GL383580.2, KI270706.1, KI270712.1, KI270713.1, KI270714.1, KI270729.1, KI270744.1, KI270782.1, KI270803.1, KI270809.1, KI270821.1, KI270830.1, KI270831.1, KI270832.1, KI270853.1, KI270856.1, KI270857.1, KI270860.1, KI270861.1, KI270878.1, KI270879.1, KI270903.1, KI270904.1, KI270905.1, KI270908.1, KN196481.1, KN196484.1, KN538364.1, KQ090026.1, KQ458383.1, KQ458384.1, KQ759762.1, KV575244.1, KV880764.1, KV880768.1, KZ208904.1, KZ208906.1, KZ208915.1, KZ208921.1, KZ208922.1, KZ559105.1, ML143360.1, ML143366.1, ML143367.1, ML143372.1, ML143377.1, ML143380.1, KV766196.1, ML143371.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "H3K4me3"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000218.1, GL000219.1, GL000220.1, GL000251.2, GL000252.2, GL000253.2, GL000254.2, GL000255.2, GL000256.2, GL339449.2, GL383563.3, KI270712.1, KI270713.1, KI270714.1, KI270717.1, KI270721.1, KI270744.1, KI270762.1, KI270816.1, KI270824.1, KI270830.1, KI270832.1, KI270847.1, KI270849.1, KI270850.1, KI270853.1, KI270857.1, KI270861.1, KI270872.1, KI270879.1, KI270894.1, KI270897.1, KI270904.1, KI270905.1, KI270908.1, KN538364.1, KQ031389.1, KQ090026.1, KQ759762.1, KV575244.1, KV880764.1, KV880768.1, KZ208906.1, KZ208915.1, KZ208922.1, KZ559105.1, KZ559109.1, ML143366.1, ML143371.1, ML143372.1, ML143377.1, ML143380.1, KI270765.1, KI270848.1, KI270856.1, KZ208917.1, KI270333.1, KI270337.1, KI270466.1, KI270467.1, KI270878.1, KI270336.1, KI270438.1, KI270442.1, KI270515.1, KI270729.1, KI270735.1, KI270736.1, KI270750.1, KN196487.1, GL000216.2, GL000224.1, GL000225.1, KI270435.1, KI270509.1, KI270519.1, KI270588.1, KI270591.1, KI270709.1, KI270746.1, KI270751.1, KI270757.1, GL383556.1, GL383574.1, KI270745.1, KI270815.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "H3K79me2"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL383557.1, KI270729.1, KI270744.1, KI270830.1, KI270832.1, KI270866.1, KI270879.1, KI270904.1, KN538364.1, KQ458384.1, KV575244.1, KV880768.1, KZ208915.1, ML143366.1, ML143377.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "H3K9ac"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000218.1, GL000219.1, GL000220.1, GL000251.2, GL000252.2, GL339449.2, GL383563.3, GL383574.1, KI270713.1, KI270721.1, KI270729.1, KI270736.1, KI270744.1, KI270745.1, KI270830.1, KI270857.1, KI270879.1, KI270904.1, KI270908.1, KN196484.1, KQ458384.1, KV575244.1, KV880768.1, KZ208915.1, ML143377.1, ML143380.1, GL000254.2, GL000255.2, KI270866.1, KI270872.1, KI270903.1, KQ090026.1, KV880764.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "H3K9me3"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000008.2, GL000009.2, GL000208.1, GL000216.2, GL000224.1, GL000225.1, KI270435.1, KI270438.1, KI270442.1, KI270468.1, KI270511.1, KI270512.1, KI270515.1, KI270516.1, KI270517.1, KI270519.1, KI270522.1, KI270538.1, KI270584.1, KI270587.1, KI270589.1, KI270590.1, KI270591.1, KI270709.1, KI270716.1, KI270718.1, KI270719.1, KI270720.1, KI270722.1, KI270723.1, KI270725.1, KI270729.1, KI270730.1, KI270732.1, KI270735.1, KI270736.1, KI270746.1, KI270747.1, KI270749.1, KI270750.1, KI270751.1, KI270753.1, KI270756.1, KI270757.1, KN196487.1, KN538372.1, KQ031387.1, KV766197.1, ML143354.1, ML143364.1, ML143378.1, GL000194.1, GL000218.1, GL000219.1, GL000221.1, KI270320.1, KI270322.1, KI270333.1, KI270337.1, KI270429.1, KI270448.1, KI270465.1, KI270466.1, KI270467.1, KI270507.1, KI270508.1, KI270509.1, KI270510.1, KI270518.1, KI270521.1, KI270580.1, KI270583.1, KI270588.1, KI270593.1, KI270710.1, KI270712.1, KI270724.1, KI270728.1, KI270733.1, KI270738.1, KI270742.1, KI270743.1, KI270775.1, KI270839.1, KI270881.1, KI270907.1, KN538360.1, KN538367.1, KN538368.1, KQ031384.1, KV766199.1, KZ208924.1, ML143359.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "H4K20me1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': KI270333.1, KI270729.1, KI270819.1, KI270866.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "HCFC1"
    ## [1] "HDAC1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000195.1, GL000218.1, GL000219.1, GL000251.2, GL000252.2, GL000253.2, GL000254.2, GL000255.2, GL000256.2, GL339449.2, GL383556.1, GL383563.3, GL383577.2, KI270466.1, KI270467.1, KI270712.1, KI270713.1, KI270714.1, KI270721.1, KI270733.1, KI270734.1, KI270744.1, KI270750.1, KI270751.1, KI270782.1, KI270816.1, KI270830.1, KI270832.1, KI270850.1, KI270853.1, KI270857.1, KI270860.1, KI270861.1, KI270862.1, KI270869.1, KI270872.1, KI270878.1, KI270879.1, KI270897.1, KI270904.1, KI270905.1, KI270908.1, KN196484.1, KN196487.1, KQ090026.1, KQ458384.1, KQ759762.1, KV575244.1, KV766196.1, KV880764.1, KV880768.1, KZ208915.1, KZ208921.1, KZ208922.1, KZ559109.1, ML143371.1, ML143377.1, ML143380.1, KI270745.1, KI270842.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "HDAC2"
    ## [1] "HINFP"
    ## [1] "HIVEP1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000009.2, GL000194.1, GL000214.1, GL000218.1, GL000219.1, GL000220.1, GL000221.1, GL000224.1, GL000225.1, GL000251.2, GL000252.2, GL000253.2, GL000254.2, GL000255.2, GL000257.2, GL339449.2, GL383522.1, GL383542.1, GL383563.3, GL383574.1, GL383577.2, GL383581.2, GL877875.1, JH159146.1, KI270442.1, KI270706.1, KI270709.1, KI270711.1, KI270712.1, KI270713.1, KI270719.1, KI270721.1, KI270722.1, KI270723.1, KI270728.1, KI270731.1, KI270733.1, KI270734.1, KI270743.1, KI270744.1, KI270745.1, KI270751.1, KI270754.1, KI270772.1, KI270791.1, KI270810.1, KI270816.1, KI270819.1, KI270827.1, KI270830.1, KI270832.1, KI270842.1, KI270846.1, KI270850.1, KI270853.1, KI270854.1, KI270857.1, KI270861.1, KI270862.1, KI270869.1, KI270878.1, KI270879.1, KI270894.1, KI270901.1, KI270903.1, KI270905.1, KI270907.1, KI270908.1, KI270913.1, KN196484.1, KN538364.1, KN538372.1, KQ090026.1, KQ759762.1, KV575244.1, KV766196.1, KV880764.1, KV880768.1, KZ208906.1, KZ208908.1, KZ208915.1, KZ208918.1, KZ208921.1, KZ208922.1, KZ559105.1, KZ559109.1, ML143343.1, ML143364.1, ML143371.1, ML143372.1, ML143377.1, ML143380.1, GL000205.2, GL000208.1, GL000216.2, GL000256.2, GL000258.2, GL383533.1, GL383550.2, GL383566.1, GL949742.1, KI270707.1, KI270714.1, KI270724.1, KI270725.1, KI270729.1, KI270732.1, KI270742.1, KI270749.1, KI270750.1, KI270753.1, KI270762.1, KI270765.1, KI270781.1, KI270792.1, KI270805.1, KI270824.1, KI270831.1, KI270856.1, KI270860.1, KI270868.1, KI270871.1, KI270872.1, KI270896.1, KI270902.1, KN196480.1, KN196487.1, KQ090022.1, KQ458383.1, KQ458386.1, KZ208914.1, KZ559103.1, ML143358.1, ML143367.1, ML143374.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "HMG20A"
    ## [1] "HMG20B"
    ## [1] "HMGXB3"
    ## [1] "HMGXB4"
    ## [1] "HNF1A"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000218.1, GL000219.1, GL000251.2, GL000252.2, GL000254.2, GL000257.2, GL339449.2, GL383563.3, GL383577.2, KI270442.1, KI270713.1, KI270721.1, KI270744.1, KI270782.1, KI270819.1, KI270830.1, KI270832.1, KI270853.1, KI270857.1, KI270861.1, KI270862.1, KI270872.1, KI270878.1, KI270879.1, KI270905.1, KI270907.1, KN196484.1, KQ031389.1, KQ983256.1, KV575244.1, KV880768.1, KZ208915.1, KZ208921.1, KZ208922.1, ML143343.1, ML143371.1, ML143372.1, ML143377.1, ML143380.1, GL000008.2, GL000009.2, GL000195.1, GL000205.2, GL000214.1, GL000216.2, GL000220.1, GL000224.1, GL000225.1, GL000258.2, GL383533.1, KI270438.1, KI270589.1, KI270591.1, KI270709.1, KI270712.1, KI270718.1, KI270719.1, KI270720.1, KI270722.1, KI270723.1, KI270725.1, KI270728.1, KI270730.1, KI270732.1, KI270733.1, KI270734.1, KI270735.1, KI270736.1, KI270745.1, KI270746.1, KI270750.1, KI270751.1, KI270753.1, KI270754.1, KI270756.1, KI270757.1, KI270762.1, KI270816.1, KI270821.1, KI270846.1, KI270869.1, KI270894.1, KI270913.1, KN538372.1, KV880764.1, ML143364.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "HNF1B"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000195.1, GL000218.1, GL000219.1, GL000224.1, GL000251.2, GL000252.2, GL000253.2, GL000254.2, GL000255.2, GL339449.2, GL383555.2, GL383556.1, GL383557.1, GL383563.3, GL383574.1, GL383577.2, GL383580.2, GL877875.1, KI270442.1, KI270711.1, KI270713.1, KI270714.1, KI270721.1, KI270728.1, KI270731.1, KI270734.1, KI270744.1, KI270751.1, KI270781.1, KI270782.1, KI270803.1, KI270809.1, KI270830.1, KI270832.1, KI270845.1, KI270846.1, KI270850.1, KI270853.1, KI270856.1, KI270857.1, KI270861.1, KI270869.1, KI270872.1, KI270878.1, KI270879.1, KI270894.1, KI270903.1, KI270904.1, KI270905.1, KI270907.1, KI270908.1, KN196483.1, KN196484.1, KQ090026.1, KQ458383.1, KQ458384.1, KQ759762.1, KQ983256.1, KV575244.1, KV880763.1, KV880764.1, KV880768.1, KZ208906.1, KZ208915.1, KZ208918.1, KZ208921.1, KZ208922.1, KZ559105.1, ML143360.1, ML143367.1, ML143372.1, ML143377.1, GL000205.2, GL000225.1, GL000256.2, GL383578.2, KI270762.1, KI270862.1, KI270875.1, KI270896.1, KV766197.1, ML143371.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "HNF4A"
    ## [1] "HNF4G"
    ## [1] "HOMEZ"
    ## [1] "HOXA3"
    ## [1] "HOXA5"
    ## [1] "HOXD1"
    ## [1] "HSF2"
    ## [1] "IKZF4"
    ## [1] "IKZF5"
    ## [1] "IRF1"
    ## [1] "IRF2"
    ## [1] "IRF5"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000218.1, GL000219.1, GL000220.1, GL339449.2, KI270712.1, KI270713.1, KI270721.1, KI270728.1, KI270733.1, KI270744.1, KI270821.1, KI270830.1, KI270861.1, KI270869.1, KI270879.1, KI270903.1, KI270908.1, KV880768.1, KZ208915.1, ML143377.1, GL383577.2, KI270442.1, KI270857.1, KI270878.1, KI270894.1, ML143372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ISL2"
    ## [1] "ISX"
    ## [1] "JUN"
    ## [1] "JUND"
    ## [1] "KAT7"
    ## [1] "KAT8"
    ## [1] "KDM1A"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000219.1, GL000252.2, GL383563.3, GL383577.2, KI270713.1, KI270853.1, KI270857.1, KI270861.1, KI270878.1, KV880768.1, ML143377.1, GL000216.2, GL000218.1, GL000224.1, GL000225.1, KI270330.1, KI270467.1, KI270733.1, KI270879.1, KI270442.1, KI270729.1, KN196480.1, KN196484.1, KN538361.1, KQ031389.1, KV880764.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "KDM2A"
    ## [1] "KDM3A"
    ## [1] "KDM4B"
    ## [1] "KDM5B"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000214.1, GL000218.1, GL000219.1, GL000255.2, GL000256.2, GL383563.3, KI270706.1, KI270712.1, KI270713.1, KI270714.1, KI270721.1, KI270728.1, KI270731.1, KI270733.1, KI270734.1, KI270744.1, KI270830.1, KI270832.1, KI270857.1, KI270862.1, KI270869.1, KI270879.1, KI270894.1, KI270903.1, KI270904.1, KI270905.1, KI270908.1, KN196484.1, KN538364.1, KQ090026.1, KV880768.1, KZ208915.1, KZ208921.1, KZ208922.1, ML143353.1, ML143372.1, GL000008.2, GL000195.1, GL000205.2, GL000220.1, GL000224.1, GL000225.1, GL000251.2, GL000252.2, GL000253.2, GL000254.2, GL000258.2, GL339449.2, GL383522.1, GL383533.1, GL383555.2, GL383556.1, GL383557.1, GL383566.1, GL383574.1, GL383577.2, GL383580.2, GL383581.2, GL877875.1, JH159146.1, JH159148.1, KI270442.1, KI270707.1, KI270709.1, KI270711.1, KI270719.1, KI270723.1, KI270726.1, KI270745.1, KI270746.1, KI270750.1, KI270751.1, KI270754.1, KI270762.1, KI270772.1, KI270780.1, KI270782.1, KI270791.1, KI270804.1, KI270809.1, KI270816.1, KI270819.1, KI270820.1, KI270821.1, KI270824.1, KI270831.1, KI270835.1, KI270842.1, KI270846.1, KI270847.1, KI270850.1, KI270853.1, KI270854.1, KI270856.1, KI270860.1, KI270861.1, KI270868.1, KI270870.1, KI270872.1, KI270877.1, KI270878.1, KI270896.1, KI270897.1, KI270907.1, KI270927.1, KI270936.1, KI270937.1, KN196480.1, KN196487.1, KN538372.1, KQ031389.1, KQ090021.1, KQ458383.1, KQ458386.1, KQ759762.1, KV575244.1, KV766193.1, KV766196.1, KV880764.1, KZ208904.1, KZ208906.1, KZ208908.1, KZ208919.1, KZ559105.1, KZ559109.1, ML143341.1, ML143358.1, ML143364.1, ML143365.1, ML143366.1, ML143367.1, ML143371.1, ML143377.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "KDM6A"
    ## [1] "KLF11"
    ## [1] "KLF12"
    ## [1] "KLF13"
    ## [1] "KLF16"
    ## [1] "KLF6"
    ## [1] "KLF9"
    ## [1] "KMT2A"
    ## [1] "KMT2B"
    ## [1] "LBX2"
    ## [1] "LCOR"
    ## [1] "LCORL"
    ## [1] "LIN54"
    ## [1] "LRRFIP1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrY
    ##   - in 'y': GL000216.2, GL000224.1, GL000225.1, KI270438.1, KI270442.1, KI270466.1, KI270467.1, KI270713.1, KI270729.1, KI270733.1, KI270744.1, KI270869.1, KN196487.1, KV880768.1, ML143380.1, GL000219.1, GL383563.3, GL383577.2, KI270330.1, KI270731.1, KI270879.1, KN196484.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "MAF1"
    ## [1] "MAFF"
    ## [1] "MAFG"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000224.1, GL000254.2, KI270712.1, KI270728.1, KI270733.1, KI270744.1, KI270762.1, KI270853.1, KI270908.1, KN196483.1, KN538372.1, KZ559100.1, ML143377.1, GL000220.1, KI270330.1, KI270754.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "MAFK"
    ## [1] "MATR3"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrY
    ##   - in 'y': GL000008.2, GL000194.1, GL000214.1, GL000216.2, GL000218.1, GL000219.1, GL000220.1, GL000224.1, GL000225.1, GL877875.1, KI270330.1, KI270438.1, KI270442.1, KI270589.1, KI270709.1, KI270712.1, KI270713.1, KI270718.1, KI270719.1, KI270723.1, KI270728.1, KI270729.1, KI270730.1, KI270733.1, KI270734.1, KI270735.1, KI270736.1, KI270744.1, KI270754.1, KI270756.1, KI270757.1, KI270762.1, KI270772.1, KI270830.1, KI270878.1, KI270879.1, KI270894.1, KN196487.1, KN538372.1, KV880768.1, ML143377.1, ML143380.1, KI270731.1, KI270819.1, KQ031384.1, KQ031389.1, KQ090026.1, KZ208915.1, KZ208922.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "MAX"
    ## [1] "MAZ"
    ## [1] "MED1"
    ## [1] "MED13"
    ## [1] "MEF2A"
    ## [1] "MEF2D"
    ## [1] "MIER2"
    ## [1] "MIER3"
    ## [1] "MIXL1"
    ## [1] "MLX"
    ## [1] "MNX1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000214.1, GL000218.1, GL000219.1, GL000224.1, GL000251.2, GL000252.2, GL000253.2, GL000254.2, GL000255.2, GL000256.2, GL339449.2, GL383556.1, GL383563.3, GL383574.1, GL383577.2, GL383581.2, KI270442.1, KI270712.1, KI270713.1, KI270714.1, KI270719.1, KI270721.1, KI270728.1, KI270733.1, KI270734.1, KI270744.1, KI270745.1, KI270754.1, KI270782.1, KI270816.1, KI270821.1, KI270830.1, KI270832.1, KI270850.1, KI270853.1, KI270857.1, KI270861.1, KI270862.1, KI270869.1, KI270872.1, KI270877.1, KI270878.1, KI270879.1, KI270894.1, KI270905.1, KI270907.1, KI270908.1, KN196484.1, KN538364.1, KN538372.1, KQ090026.1, KQ458384.1, KQ759762.1, KV575244.1, KV766196.1, KV880768.1, KZ208915.1, KZ208922.1, KZ559109.1, ML143371.1, ML143372.1, ML143377.1, ML143380.1, GL000009.2, GL000205.2, GL000216.2, GL000225.1, GL383527.1, GL383533.1, GL383550.2, GL383557.1, GL877875.1, JH159146.1, JH159147.1, KI270706.1, KI270709.1, KI270711.1, KI270717.1, KI270722.1, KI270723.1, KI270726.1, KI270729.1, KI270731.1, KI270748.1, KI270762.1, KI270765.1, KI270772.1, KI270780.1, KI270791.1, KI270802.1, KI270804.1, KI270813.1, KI270819.1, KI270827.1, KI270835.1, KI270845.1, KI270846.1, KI270847.1, KI270854.1, KI270856.1, KI270868.1, KI270871.1, KI270897.1, KI270903.1, KI270904.1, KI270913.1, KQ090021.1, KQ090022.1, KQ458383.1, KQ458386.1, KZ208904.1, KZ208906.1, KZ208909.1, KZ208914.1, KZ208921.1, KZ559105.1, ML143355.1, ML143358.1, ML143364.1, ML143366.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "MTA1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000216.2, GL000218.1, GL000219.1, GL000252.2, GL383563.3, GL383577.2, KI270438.1, KI270442.1, KI270713.1, KI270714.1, KI270718.1, KI270721.1, KI270754.1, KI270782.1, KI270857.1, KI270861.1, KI270869.1, KI270878.1, KI270905.1, KN196484.1, KV880768.1, KZ208915.1, ML143377.1, GL000009.2, GL000194.1, GL000208.1, GL000214.1, GL000220.1, GL000221.1, GL000224.1, GL000225.1, GL000251.2, GL000253.2, GL000254.2, GL000255.2, GL000256.2, GL000258.2, GL339449.2, GL383522.1, GL383533.1, GL383550.2, GL383555.2, GL383557.1, GL383574.1, GL383581.2, GL877875.1, GL949742.1, JH159146.1, JH159147.1, KI270508.1, KI270519.1, KI270706.1, KI270707.1, KI270709.1, KI270712.1, KI270719.1, KI270720.1, KI270722.1, KI270723.1, KI270724.1, KI270725.1, KI270726.1, KI270728.1, KI270729.1, KI270731.1, KI270732.1, KI270733.1, KI270734.1, KI270735.1, KI270742.1, KI270743.1, KI270744.1, KI270745.1, KI270746.1, KI270749.1, KI270750.1, KI270751.1, KI270753.1, KI270762.1, KI270764.1, KI270765.1, KI270767.1, KI270772.1, KI270791.1, KI270792.1, KI270795.1, KI270797.1, KI270804.1, KI270805.1, KI270810.1, KI270813.1, KI270816.1, KI270819.1, KI270821.1, KI270824.1, KI270827.1, KI270830.1, KI270831.1, KI270832.1, KI270835.1, KI270842.1, KI270845.1, KI270846.1, KI270847.1, KI270850.1, KI270853.1, KI270854.1, KI270856.1, KI270860.1, KI270862.1, KI270868.1, KI270870.1, KI270871.1, KI270872.1, KI270876.1, KI270877.1, KI270879.1, KI270894.1, KI270897.1, KI270901.1, KI270902.1, KI270903.1, KI270907.1, KI270908.1, KI270910.1, KI270913.1, KI270927.1, KI270936.1, KI270937.1, KN196481.1, KN196487.1, KN538364.1, KN538372.1, KQ031389.1, KQ090017.1, KQ090022.1, KQ090026.1, KQ458384.1, KQ458386.1, KQ759759.1, KQ759762.1, KV575244.1, KV766193.1, KV766196.1, KV880764.1, KZ208904.1, KZ208906.1, KZ208908.1, KZ208909.1, KZ208917.1, KZ208918.1, KZ208919.1, KZ208921.1, KZ208922.1, KZ559100.1, KZ559105.1, KZ559109.1, KZ559113.1, ML143353.1, ML143355.1, ML143358.1, ML143361.1, ML143364.1, ML143365.1, ML143366.1, ML143367.1, ML143370.1, ML143371.1, ML143372.1, ML143373.1, ML143374.1, ML143378.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "MTF2"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000195.1, GL000251.2, GL000252.2, GL000255.2, GL000256.2, GL339449.2, GL383563.3, KI270712.1, KI270713.1, KI270717.1, KI270731.1, KI270765.1, KI270824.1, KI270829.1, KI270830.1, KI270850.1, KI270857.1, KI270862.1, KI270869.1, KI270879.1, KI270894.1, KI270905.1, KI270908.1, KQ458383.1, KV880768.1, KZ208906.1, KZ208915.1, KZ559113.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "MXD1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000214.1, GL000216.2, GL000218.1, GL000219.1, GL000220.1, GL000224.1, GL000225.1, GL000252.2, GL000254.2, GL000255.2, GL339449.2, GL383533.1, GL383563.3, GL383574.1, GL383577.2, JH159147.1, KI270435.1, KI270438.1, KI270442.1, KI270589.1, KI270591.1, KI270706.1, KI270709.1, KI270712.1, KI270713.1, KI270714.1, KI270719.1, KI270721.1, KI270722.1, KI270723.1, KI270726.1, KI270728.1, KI270733.1, KI270734.1, KI270744.1, KI270745.1, KI270750.1, KI270754.1, KI270762.1, KI270781.1, KI270782.1, KI270813.1, KI270816.1, KI270821.1, KI270830.1, KI270832.1, KI270846.1, KI270850.1, KI270853.1, KI270857.1, KI270861.1, KI270862.1, KI270868.1, KI270869.1, KI270879.1, KI270894.1, KI270896.1, KI270897.1, KI270904.1, KI270905.1, KI270908.1, KI270927.1, KN196484.1, KN538364.1, KN538372.1, KQ458383.1, KQ458386.1, KQ759762.1, KV575244.1, KV766193.1, KV766196.1, KV880768.1, KZ208906.1, KZ208915.1, KZ208921.1, KZ559113.1, ML143364.1, ML143371.1, ML143377.1, ML143380.1, GL000256.2, KI270872.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "MXD3"
    ## [1] "MXD4"
    ## [1] "MXI1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000225.1, KI270438.1, KI270467.1, KI270713.1, KI270729.1, KI270733.1, KI270879.1, KN196487.1, GL000194.1, GL000219.1, GL000252.2, GL339449.2, GL383563.3, GL383577.2, KI270830.1, KI270857.1, KI270861.1, KI270862.1, KQ759762.1, KV575244.1, KV880768.1, KZ208915.1, ML143367.1, ML143377.1, ML143380.1, KI270466.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "MYNN"
    ## [1] "MYRF"
    ## [1] "NAIF1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000219.1, KI270713.1, KI270908.1, KN196487.1, KI270466.1, KI270467.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "NFAT5"
    ## [1] "NFATC3"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000214.1, GL000219.1, GL000255.2, GL383577.2, KI270713.1, KI270721.1, KI270733.1, KI270744.1, KI270745.1, KI270830.1, KI270862.1, KI270879.1, KN196484.1, KQ090026.1, KV880768.1, ML143377.1, KI270711.1, KQ458383.1, KZ208915.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "NFE2L1"
    ## [1] "NFE2L2"
    ## [1] "NFIA"
    ## [1] "NFIC"
    ## [1] "NFIL3"
    ## [1] "NFKB2"
    ## [1] "NFKBIZ"
    ## [1] "NFYA"
    ## [1] "NFYB"
    ## [1] "NFYC"
    ## [1] "NKX3-1"
    ## [1] "NR0B2"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000208.1, GL000214.1, GL000216.2, GL000218.1, GL000219.1, GL000220.1, GL000224.1, GL000225.1, GL339449.2, GL383566.1, KI270435.1, KI270442.1, KI270589.1, KI270712.1, KI270713.1, KI270721.1, KI270723.1, KI270725.1, KI270728.1, KI270731.1, KI270733.1, KI270744.1, KI270754.1, KI270816.1, KI270830.1, KI270832.1, KI270857.1, KI270869.1, KI270879.1, KI270894.1, KI270908.1, KN538360.1, KN538372.1, KV766196.1, KV880768.1, KZ208915.1, KZ208922.1, KZ559105.1, ML143343.1, ML143367.1, ML143377.1, ML143380.1, GL000008.2, GL000009.2, GL000205.2, GL000251.2, GL000252.2, GL000255.2, GL000256.2, GL383522.1, GL383527.1, GL383533.1, GL383563.3, GL383574.1, GL383577.2, GL383581.2, JH159147.1, KI270714.1, KI270719.1, KI270722.1, KI270724.1, KI270750.1, KI270753.1, KI270772.1, KI270791.1, KI270802.1, KI270821.1, KI270831.1, KI270842.1, KI270847.1, KI270850.1, KI270853.1, KI270856.1, KI270861.1, KI270868.1, KI270878.1, KI270903.1, KI270904.1, KI270907.1, KI270936.1, KN196484.1, KN538364.1, KQ090028.1, KQ458386.1, KV575244.1, KV766193.1, KZ208906.1, KZ208916.1, KZ208921.1, ML143353.1, ML143355.1, ML143360.1, ML143364.1, ML143372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "NR2C2"
    ## [1] "NR2F1"
    ## [1] "NR2F6"
    ## [1] "NR5A1"
    ## [1] "NRF1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrY
    ##   - in 'y': KI270438.1, KI270466.1, KI270467.1, KI270857.1, KN196487.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "NRL"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000009.2, GL000194.1, GL000195.1, GL000208.1, GL000214.1, GL000216.2, GL000218.1, GL000219.1, GL000220.1, GL000224.1, GL000225.1, GL000252.2, GL000253.2, GL000254.2, GL000255.2, GL000256.2, GL339449.2, GL383563.3, GL383574.1, GL383577.2, GL877875.1, JH159146.1, KI270442.1, KI270706.1, KI270709.1, KI270712.1, KI270713.1, KI270714.1, KI270717.1, KI270719.1, KI270721.1, KI270722.1, KI270723.1, KI270724.1, KI270728.1, KI270731.1, KI270732.1, KI270733.1, KI270734.1, KI270742.1, KI270743.1, KI270744.1, KI270745.1, KI270746.1, KI270748.1, KI270750.1, KI270751.1, KI270754.1, KI270772.1, KI270782.1, KI270791.1, KI270810.1, KI270816.1, KI270819.1, KI270830.1, KI270831.1, KI270832.1, KI270845.1, KI270846.1, KI270847.1, KI270850.1, KI270853.1, KI270856.1, KI270857.1, KI270861.1, KI270862.1, KI270868.1, KI270869.1, KI270872.1, KI270878.1, KI270879.1, KI270894.1, KI270903.1, KI270905.1, KI270907.1, KI270908.1, KN196484.1, KN196487.1, KN538364.1, KN538372.1, KQ090020.1, KQ090026.1, KQ458383.1, KQ458386.1, KQ759762.1, KV575244.1, KV766196.1, KV880768.1, KZ208915.1, KZ208918.1, KZ208921.1, KZ208922.1, KZ559105.1, KZ559109.1, ML143372.1, ML143373.1, ML143374.1, ML143377.1, ML143380.1, GL000221.1, GL000251.2, GL383527.1, GL383531.1, GL383542.1, GL383550.2, GL383566.1, JH159147.1, KI270438.1, KI270519.1, KI270711.1, KI270718.1, KI270720.1, KI270725.1, KI270726.1, KI270729.1, KI270735.1, KI270753.1, KI270762.1, KI270795.1, KI270804.1, KI270805.1, KI270824.1, KI270876.1, KI270877.1, KI270901.1, KI270904.1, KQ031389.1, KQ759759.1, KZ208908.1, KZ559100.1, ML143341.1, ML143353.1, ML143364.1, ML143365.1, ML143366.1, ML143367.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ONECUT1"
    ## [1] "ONECUT2"
    ## [1] "PAF1"
    ## [1] "PATZ1"
    ## [1] "PAWR"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': KI270713.1, KI270862.1, KV880768.1, KI270728.1, KI270733.1, KI270846.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "PAXIP1"
    ## [1] "PHF20"
    ## [1] "PHF21A"
    ## [1] "PHF5A"
    ## [1] "PHF8"
    ## [1] "PITX1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000195.1, GL000205.2, GL000214.1, GL000218.1, GL000219.1, GL000220.1, GL000251.2, GL000252.2, GL000253.2, GL000254.2, GL000255.2, GL339449.2, GL383556.1, GL383563.3, GL383574.1, GL383577.2, GL383580.2, KI270442.1, KI270711.1, KI270712.1, KI270713.1, KI270714.1, KI270721.1, KI270726.1, KI270728.1, KI270731.1, KI270733.1, KI270744.1, KI270745.1, KI270751.1, KI270754.1, KI270782.1, KI270803.1, KI270829.1, KI270830.1, KI270832.1, KI270842.1, KI270846.1, KI270847.1, KI270850.1, KI270853.1, KI270857.1, KI270861.1, KI270862.1, KI270869.1, KI270877.1, KI270878.1, KI270879.1, KI270894.1, KI270903.1, KI270904.1, KI270905.1, KI270907.1, KI270908.1, KN196484.1, KN538364.1, KQ090026.1, KQ090028.1, KQ458383.1, KQ458384.1, KQ759762.1, KV575244.1, KV880768.1, KZ208906.1, KZ208908.1, KZ208915.1, KZ208921.1, KZ208922.1, KZ559105.1, ML143367.1, ML143371.1, ML143377.1, ML143380.1, ML143381.1, GL000224.1, GL000256.2, KI270722.1, KI270734.1, KI270816.1, KI270821.1, ML143353.1, ML143366.1, ML143375.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "POGZ"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000218.1, GL000219.1, GL000251.2, GL000252.2, GL000254.2, GL000255.2, GL000257.2, GL339449.2, GL383563.3, GL383574.1, GL383577.2, KI270589.1, KI270712.1, KI270713.1, KI270721.1, KI270728.1, KI270733.1, KI270744.1, KI270745.1, KI270782.1, KI270821.1, KI270829.1, KI270830.1, KI270850.1, KI270853.1, KI270857.1, KI270861.1, KI270862.1, KI270879.1, KI270904.1, KI270905.1, KI270908.1, KN196484.1, KN538364.1, KQ458383.1, KQ458384.1, KV575244.1, KV880768.1, KZ208915.1, KZ208921.1, KZ208922.1, ML143371.1, ML143377.1, ML143380.1, ML143381.1, KI270805.1, KI270832.1, KI270842.1, KZ208906.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "POLR2A"
    ## [1] "POLR2AphosphoS2"
    ## [1] "POLR2AphosphoS5"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000218.1, GL000219.1, GL000220.1, GL000251.2, GL000252.2, GL000253.2, GL000254.2, GL000255.2, GL339449.2, GL383563.3, GL383574.1, GL383577.2, KI270713.1, KI270721.1, KI270744.1, KI270782.1, KI270819.1, KI270831.1, KI270847.1, KI270853.1, KI270857.1, KI270866.1, KI270878.1, KI270879.1, KI270903.1, KI270904.1, KI270905.1, KI270908.1, KQ458383.1, KQ759762.1, KV575244.1, KV880768.1, KZ208915.1, KZ208921.1, ML143377.1, ML143380.1, KI270438.1, KI270816.1, KI270830.1, KI270832.1, KI270850.1, KI270861.1, KI270907.1, KN538364.1, KQ090026.1, KQ090028.1, KV766196.1, KZ208904.1, ML143366.1, ML143367.1, ML143372.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "PPARG"
    ## [1] "PRDM10"
    ## [1] "PRDM15"
    ## [1] "PREB"
    ## [1] "RAD21"
    ## [1] "RARA"
    ## [1] "RARG"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrY
    ##   - in 'y': KI270438.1, KI270713.1, KN196487.1, ML143380.1, GL000194.1, GL000218.1, KI270721.1, KI270730.1, KI270733.1, KI270744.1, KI270869.1, KI270879.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "RBAK"
    ## [1] "RBPJ"
    ## [1] "RCOR1"
    ## [1] "RCOR2"
    ## [1] "RELA"
    ## [1] "RERE"
    ## [1] "REST"
    ## [1] "RFX3"
    ## [1] "RFX5"
    ## [1] "RFXANK"
    ## [1] "RFXAP"
    ## [1] "RREB1"
    ## [1] "RXRA"
    ## [1] "RXRB"
    ## [1] "SAFB2"
    ## [1] "SAP130"
    ## [1] "SATB2"
    ## [1] "SFPQ"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000219.1, GL000251.2, GL000252.2, GL000255.2, GL339449.2, GL383577.2, KI270713.1, KI270721.1, KI270731.1, KI270744.1, KI270782.1, KI270861.1, KI270869.1, KI270879.1, KI270904.1, KI270908.1, KN538364.1, KQ458383.1, KV880768.1, KZ208921.1, ML143377.1, GL000195.1, GL000254.2, GL383563.3, KI270711.1, KI270714.1, KI270809.1, KI270816.1, KI270830.1, KI270831.1, KI270832.1, KI270847.1, KI270850.1, KI270857.1, KI270862.1, KI270866.1, KI270878.1, KI270894.1, KI270903.1, KI270907.1, KN196484.1, KQ458384.1, KQ759762.1, KV575244.1, KZ208906.1, KZ208915.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "SIN3A"
    ## [1] "SIN3B"
    ## [1] "SIX1"
    ## [1] "SIX4"
    ## [1] "SMAD1"
    ## [1] "SMAD3"
    ## [1] "SMAD4"
    ## [1] "SMAD9"
    ## [1] "SMC3"
    ## [1] "SNAI1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000252.2, GL383563.3, KI270713.1, KI270714.1, KI270721.1, KI270728.1, KI270733.1, KI270745.1, KI270846.1, KI270857.1, KI270862.1, KI270879.1, KI270894.1, KI270904.1, ML143380.1, GL000219.1, KI270466.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "SNAPC4"
    ## [1] "SOX13"
    ## [1] "SOX18"
    ## [1] "SOX6"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000218.1, GL000219.1, GL000251.2, GL000252.2, GL000254.2, GL000255.2, GL000257.2, GL339449.2, GL383563.3, GL383574.1, GL383577.2, GL383580.2, GL383581.2, KI270442.1, KI270713.1, KI270714.1, KI270721.1, KI270731.1, KI270734.1, KI270744.1, KI270745.1, KI270782.1, KI270816.1, KI270829.1, KI270830.1, KI270842.1, KI270850.1, KI270853.1, KI270856.1, KI270857.1, KI270861.1, KI270862.1, KI270872.1, KI270875.1, KI270878.1, KI270879.1, KI270894.1, KI270903.1, KI270904.1, KI270905.1, KI270908.1, KN196483.1, KN196484.1, KQ031389.1, KQ458383.1, KQ458384.1, KV575244.1, KV766197.1, KV880768.1, KZ208906.1, KZ208915.1, KZ208921.1, KZ208922.1, KZ559105.1, KZ559109.1, ML143360.1, ML143371.1, ML143372.1, ML143377.1, ML143380.1, ML143381.1, GL000253.2, KI270712.1, KI270827.1, KI270832.1, KN538364.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "SP140L"
    ## [1] "SP5"
    ## [1] "SPEN"
    ## [1] "SRF"
    ## [1] "STAG1"
    ## [1] "STAT5B"
    ## [1] "STAT6"
    ## [1] "SUZ12"
    ## [1] "TAF1"
    ## [1] "TBL1XR1"
    ## [1] "TBP"
    ## [1] "TBX2"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000218.1, GL000219.1, GL000251.2, GL000252.2, GL000254.2, GL339449.2, GL383577.2, KI270442.1, KI270712.1, KI270713.1, KI270721.1, KI270728.1, KI270733.1, KI270744.1, KI270830.1, KI270862.1, KI270869.1, KI270878.1, KI270879.1, KI270894.1, KI270908.1, KV575244.1, KV880768.1, KZ208915.1, KZ208921.1, KZ208922.1, ML143377.1, ML143380.1, GL000214.1, GL000253.2, GL000255.2, GL383556.1, GL383563.3, GL383574.1, KI270714.1, KI270731.1, KI270782.1, KI270816.1, KI270827.1, KI270832.1, KI270850.1, KI270853.1, KI270857.1, KI270861.1, KI270877.1, KI270903.1, KI270905.1, KI270907.1, KN196484.1, KN538364.1, KQ090026.1, KQ458384.1, KQ759762.1, KZ208906.1, KZ559105.1, ML143367.1, ML143371.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "TCF3"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000220.1, GL000225.1, GL383577.2, KI270713.1, KI270723.1, KI270733.1, KI270861.1, KI270879.1, KN196487.1, GL000194.1, GL000216.2, GL000219.1, GL339449.2, KI270442.1, KI270709.1, KI270721.1, KI270744.1, KI270754.1, KI270772.1, KI270857.1, KZ559109.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "TCF7L2"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000225.1, KI270333.1, KI270587.1, KI270729.1, KI270829.1, KI270337.1, KI270466.1, KI270467.1, KI270746.1, KI270750.1, ML143380.1, GL000194.1, GL000218.1, GL000219.1, GL000251.2, GL000252.2, GL000254.2, GL383577.2, KI270712.1, KI270713.1, KI270714.1, KI270721.1, KI270728.1, KI270733.1, KI270744.1, KI270861.1, KI270862.1, KI270869.1, KI270878.1, KI270879.1, KI270894.1, KI270903.1, KI270905.1, KN196484.1, KQ458384.1, KV575244.1, KV880768.1, KZ208915.1, KZ559109.1, ML143371.1, ML143377.1, GL000255.2, GL339449.2, KI270734.1, KI270782.1, KI270830.1, KI270842.1, KI270857.1, KZ208921.1, KZ208922.1, ML143381.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "TEAD1"
    ## [1] "TEAD3"
    ## [1] "TEF"
    ## [1] "TFAP4"
    ## [1] "TFDP1"
    ## [1] "TFDP2"
    ## [1] "TFE3"
    ## [1] "THAP11"
    ## [1] "THAP8"
    ## [1] "THAP9"
    ## [1] "THRA"
    ## [1] "THRB"
    ## [1] "TIGD3"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': KI270713.1, GL000194.1, GL000219.1, GL000252.2, GL339449.2, GL383577.2, GL383580.2, KI270712.1, KI270721.1, KI270733.1, KI270744.1, KI270830.1, KI270869.1, KI270879.1, KI270908.1, KV880768.1, ML143377.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "TIGD6"
    ## [1] "TOE1"
    ## [1] "TP53"
    ## [1] "TRAFD1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000218.1, GL000219.1, GL000251.2, GL000254.2, GL000255.2, GL383563.3, GL383574.1, GL383577.2, GL383580.2, KI270711.1, KI270712.1, KI270713.1, KI270721.1, KI270728.1, KI270731.1, KI270734.1, KI270744.1, KI270830.1, KI270850.1, KI270853.1, KI270857.1, KI270869.1, KI270878.1, KI270879.1, KI270894.1, KI270903.1, KI270905.1, KN196484.1, KQ090026.1, KV880768.1, KZ208915.1, KZ208921.1, ML143377.1, ML143380.1, GL000214.1, GL000220.1, KI270733.1, KI270736.1, KI270753.1, KI270908.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "TSC22D2"
    ## [1] "UBTF"
    ## [1] "USF1"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000205.2, GL000218.1, GL000251.2, GL000252.2, GL000255.2, GL383556.1, GL383557.1, KI270733.1, KI270841.1, KI270842.1, KI270878.1, KI270879.1, KI270908.1, KQ759762.1, KV880768.1, ML143377.1, ML143381.1, GL000194.1, GL000219.1, GL000220.1, GL339449.2, GL383555.2, GL383577.2, KI270713.1, KI270734.1, KI270743.1, KI270744.1, KI270830.1, KI270850.1, KI270853.1, KI270857.1, KI270861.1, KI270869.1, KI270872.1, KI270877.1, KQ090026.1, KZ208908.1, KZ208915.1, KZ208921.1, ML143364.1, ML143366.1, GL000214.1, GL000216.2, GL000224.1, GL000225.1, GL000254.2, GL000256.2, KI270438.1, KI270709.1, KI270735.1, KI270745.1, KI270748.1, KI270754.1, KI270782.1, KI270821.1, KI270832.1, KN538372.1, KV575244.1, KZ208913.1, ML143371.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "USF2"
    ## [1] "WIZ"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000218.1, GL000219.1, GL000220.1, GL000224.1, GL000225.1, GL000252.2, KI270442.1, KI270712.1, KI270713.1, KI270721.1, KI270728.1, KI270733.1, KI270744.1, KI270750.1, KI270753.1, KI270782.1, KI270830.1, KI270857.1, KI270869.1, KI270879.1, KV880768.1, KZ208921.1, ML143377.1, GL000255.2, GL383556.1, GL383574.1, GL383577.2, KI270734.1, KI270894.1, KI270908.1, KN196484.1, KZ208915.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "XBP1"
    ## [1] "YEATS2"
    ## [1] "YEATS4"
    ## [1] "YY1"
    ## [1] "ZBED4"
    ## [1] "ZBED5"
    ## [1] "ZBTB1"
    ## [1] "ZBTB10"
    ## [1] "ZBTB14"
    ## [1] "ZBTB21"
    ## [1] "ZBTB24"
    ## [1] "ZBTB25"
    ## [1] "ZBTB26"
    ## [1] "ZBTB38"
    ## [1] "ZBTB39"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000252.2, KI270713.1, KI270861.1, KI270862.1, GL000195.1, GL000251.2, GL000255.2, GL000256.2, GL339449.2, KI270721.1, KI270734.1, KI270857.1, KI270869.1, KI270879.1, KI270908.1, KV880768.1, KZ208906.1, KZ208915.1, KZ208921.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZBTB4"
    ## [1] "ZBTB42"
    ## [1] "ZBTB43"
    ## [1] "ZBTB44"
    ## [1] "ZBTB46"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, KI270712.1, KI270713.1, KI270721.1, KI270728.1, KI270733.1, KI270905.1, KV880768.1, GL000220.1, GL000225.1, KI270744.1, KI270754.1, KN538372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZBTB49"
    ## [1] "ZBTB7A"
    ## [1] "ZBTB7B"
    ## [1] "ZC3H4"
    ## [1] "ZC3H8"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrY, chrM
    ##   - in 'y': KI270713.1, KI270879.1, KI270908.1, GL000194.1, KI270733.1, KI270830.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZCCHC11"
    ## [1] "ZFP1"
    ## [1] "ZFP14"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrY
    ##   - in 'y': GL000194.1, GL000216.2, GL000220.1, GL000224.1, GL000225.1, GL000252.2, KI270438.1, KI270713.1, KI270728.1, KI270730.1, KI270733.1, KI270744.1, KI270869.1, KI270712.1, KI270754.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZFP36L2"
    ## [1] "ZFP37"
    ## [1] "ZFP41"
    ## [1] "ZFP64"
    ## [1] "ZFP82"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000218.1, GL000219.1, KI270466.1, KI270467.1, KI270713.1, KI270721.1, KI270733.1, KI270816.1, KI270830.1, KI270832.1, KI270869.1, KI270879.1, KI270908.1, KN538364.1, KV880768.1, KZ208915.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZFY"
    ## [1] "ZGPAT"
    ## [1] "ZHX3"
    ## [1] "ZIK1"
    ## [1] "ZKSCAN1"
    ## [1] "ZKSCAN5"
    ## [1] "ZKSCAN8"
    ## [1] "ZMYM3"
    ## [1] "ZNF12"
    ## [1] "ZNF124"
    ## [1] "ZNF138"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000219.1, GL000255.2, GL339449.2, KI270712.1, KI270713.1, KI270721.1, KI270728.1, KI270733.1, KI270821.1, KI270830.1, KI270869.1, KI270879.1, KN196484.1, KN538364.1, KV880768.1, KZ208922.1, GL000194.1, GL000218.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF142"
    ## [1] "ZNF143"
    ## [1] "ZNF205"
    ## [1] "ZNF217"
    ## [1] "ZNF224"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000214.1, GL000216.2, GL000219.1, GL000220.1, GL000224.1, GL000225.1, GL000252.2, GL000253.2, KI270442.1, KI270709.1, KI270711.1, KI270712.1, KI270713.1, KI270714.1, KI270726.1, KI270728.1, KI270731.1, KI270733.1, KI270735.1, KI270744.1, KI270750.1, KI270751.1, KI270754.1, KI270869.1, KI270879.1, KI270894.1, KI270903.1, KN196487.1, KQ031389.1, KV880768.1, ML143372.1, ML143377.1, GL383577.2, KI270721.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF225"
    ## [1] "ZNF230"
    ## [1] "ZNF232"
    ## [1] "ZNF234"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrY, chrM
    ##   - in 'y': GL000219.1, KI270733.1, KI270869.1, GL000216.2, GL000220.1, KI270442.1, KI270466.1, KI270713.1, KI270728.1, KI270744.1, KI270850.1, KI270894.1, KI270908.1, KN538372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF25"
    ## [1] "ZNF256"
    ## [1] "ZNF263"
    ## [1] "ZNF264"
    ## [1] "ZNF274"
    ## [1] "ZNF276"
    ## [1] "ZNF280B"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000008.2, GL000009.2, GL000194.1, GL000195.1, GL000205.2, GL000208.1, GL000214.1, GL000216.2, GL000218.1, GL000219.1, GL000220.1, GL000221.1, GL000224.1, GL000225.1, GL000251.2, GL000252.2, GL000253.2, GL000254.2, GL000255.2, GL000256.2, GL000258.2, GL339449.2, GL383530.1, GL383542.1, GL383555.2, GL383557.1, GL383563.3, GL383574.1, GL383577.2, GL383581.2, GL877875.1, GL949742.1, JH159146.1, KI270442.1, KI270706.1, KI270709.1, KI270711.1, KI270712.1, KI270713.1, KI270714.1, KI270717.1, KI270718.1, KI270719.1, KI270721.1, KI270722.1, KI270723.1, KI270724.1, KI270725.1, KI270726.1, KI270728.1, KI270731.1, KI270733.1, KI270734.1, KI270742.1, KI270743.1, KI270744.1, KI270745.1, KI270746.1, KI270748.1, KI270750.1, KI270751.1, KI270753.1, KI270754.1, KI270762.1, KI270765.1, KI270767.1, KI270772.1, KI270779.1, KI270780.1, KI270782.1, KI270792.1, KI270802.1, KI270805.1, KI270809.1, KI270813.1, KI270815.1, KI270816.1, KI270819.1, KI270821.1, KI270824.1, KI270827.1, KI270830.1, KI270831.1, KI270832.1, KI270845.1, KI270846.1, KI270847.1, KI270850.1, KI270853.1, KI270854.1, KI270856.1, KI270857.1, KI270860.1, KI270861.1, KI270862.1, KI270868.1, KI270869.1, KI270870.1, KI270871.1, KI270872.1, KI270878.1, KI270879.1, KI270894.1, KI270897.1, KI270900.1, KI270903.1, KI270904.1, KI270905.1, KI270907.1, KI270908.1, KI270913.1, KI270927.1, KI270934.1, KI270937.1, KN196484.1, KN196487.1, KN538364.1, KN538372.1, KQ031387.1, KQ031389.1, KQ090015.1, KQ090021.1, KQ090026.1, KQ090027.1, KQ458383.1, KQ458384.1, KQ458386.1, KQ759759.1, KQ759762.1, KV575244.1, KV766193.1, KV766196.1, KV880768.1, KZ208904.1, KZ208906.1, KZ208908.1, KZ208914.1, KZ208915.1, KZ208921.1, KZ208922.1, KZ559100.1, KZ559105.1, KZ559109.1, ML143353.1, ML143355.1, ML143358.1, ML143364.1, ML143365.1, ML143366.1, ML143371.1, ML143372.1, ML143373.1, ML143375.1, ML143377.1, ML143380.1, GL000213.1, GL383522.1, GL383526.1, GL383533.1, GL383550.2, JH159147.1, JH159148.1, KI270438.1, KI270519.1, KI270538.1, KI270707.1, KI270720.1, KI270729.1, KI270732.1, KI270735.1, KI270749.1, KI270764.1, KI270787.1, KI270791.1, KI270795.1, KI270797.1, KI270804.1, KI270810.1, KI270818.1, KI270823.1, KI270835.1, KI270842.1, KI270876.1, KI270896.1, KI270902.1, KI270909.1, KI270925.1, KI270936.1, KI270938.1, KN196472.1, KN538360.1, KN538361.1, KQ090020.1, KQ983256.1, KV880764.1, KZ208909.1, KZ208918.1, KZ208919.1, ML143343.1, ML143367.1, ML143374.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF280D"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000218.1, GL000219.1, GL000251.2, GL000254.2, GL000255.2, GL383577.2, GL383578.2, KI270712.1, KI270713.1, KI270721.1, KI270733.1, KI270744.1, KI270853.1, KI270862.1, KI270869.1, KI270878.1, KI270879.1, KI270905.1, KI270908.1, KQ458384.1, KV880768.1, KZ208915.1, KZ208922.1, KZ559105.1, ML143377.1, KI270726.1, KN196484.1, KV575244.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF281"
    ## [1] "ZNF296"
    ## [1] "ZNF30"
    ## [1] "ZNF317"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000218.1, GL000251.2, GL000255.2, GL383563.3, GL383577.2, GL949747.2, KI270713.1, KI270721.1, KI270746.1, KI270829.1, KI270853.1, KI270857.1, KI270860.1, KI270862.1, KI270879.1, KI270905.1, KN196484.1, KZ208906.1, KZ208915.1, GL000219.1, GL000252.2, GL000256.2, GL383574.1, KI270728.1, KI270830.1, KI270861.1, KI270869.1, KI270894.1, KI270903.1, KI270908.1, KN538364.1, KV880768.1, KZ208922.1, ML143341.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF318"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000008.2, GL000194.1, GL000195.1, GL000205.2, GL000214.1, GL000216.2, GL000218.1, GL000219.1, GL000220.1, GL000224.1, GL000225.1, GL000255.2, GL339449.2, GL383563.3, GL383577.2, KI270438.1, KI270442.1, KI270589.1, KI270706.1, KI270712.1, KI270713.1, KI270718.1, KI270719.1, KI270721.1, KI270722.1, KI270726.1, KI270728.1, KI270730.1, KI270731.1, KI270733.1, KI270734.1, KI270736.1, KI270742.1, KI270744.1, KI270745.1, KI270750.1, KI270751.1, KI270754.1, KI270762.1, KI270795.1, KI270805.1, KI270830.1, KI270846.1, KI270847.1, KI270856.1, KI270857.1, KI270861.1, KI270862.1, KI270869.1, KI270878.1, KI270879.1, KI270894.1, KI270903.1, KI270905.1, KI270908.1, KN196484.1, KN196487.1, KN538372.1, KQ031389.1, KQ090026.1, KV766193.1, KV766196.1, KV880768.1, KZ208915.1, KZ208922.1, ML143341.1, ML143353.1, ML143362.1, ML143364.1, ML143377.1, ML143379.1, ML143380.1, GL383574.1, KI270711.1, KI270714.1, KI270723.1, KI270772.1, KI270782.1, KI270816.1, KI270819.1, KI270829.1, KI270853.1, KI270866.1, KI270896.1, KN538364.1, KV880764.1, ML143355.1, ML143358.1, ML143365.1, ML143372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF326"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrY
    ##   - in 'y': GL000216.2, GL000219.1, GL000224.1, GL000225.1, GL383577.2, KI270330.1, KI270438.1, KI270442.1, KI270466.1, KI270467.1, KI270713.1, KI270721.1, KI270729.1, KI270733.1, KI270819.1, KI270847.1, KI270879.1, KN196487.1, KV880768.1, ML143377.1, ML143380.1, GL000214.1, KI270744.1, KI270830.1, KI270831.1, KI270866.1, KN196484.1, KQ031384.1, KQ090026.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF331"
    ## [1] "ZNF333"
    ## [1] "ZNF335"
    ## [1] "ZNF337"
    ## [1] "ZNF33B"
    ## [1] "ZNF34"
    ## [1] "ZNF343"
    ## [1] "ZNF350"
    ## [1] "ZNF362"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000195.1, GL000214.1, GL000218.1, GL000219.1, GL000250.2, GL000251.2, GL000252.2, GL000254.2, GL000255.2, GL000256.2, GL000257.2, GL339449.2, GL383527.1, GL383574.1, GL383577.2, GL383581.2, GL383582.2, KI270713.1, KI270714.1, KI270721.1, KI270733.1, KI270734.1, KI270744.1, KI270782.1, KI270785.1, KI270798.1, KI270806.1, KI270816.1, KI270827.1, KI270829.1, KI270830.1, KI270831.1, KI270832.1, KI270849.1, KI270853.1, KI270856.1, KI270857.1, KI270860.1, KI270861.1, KI270862.1, KI270868.1, KI270869.1, KI270875.1, KI270878.1, KI270879.1, KI270894.1, KI270903.1, KI270904.1, KI270905.1, KI270907.1, KI270908.1, KN196484.1, KN538369.1, KQ458384.1, KQ983256.1, KV575244.1, KV766198.1, KV880768.1, KZ208906.1, KZ208915.1, KZ208921.1, KZ208922.1, ML143353.1, ML143366.1, ML143372.1, ML143375.1, ML143377.1, ML143380.1, GL000220.1, KI270438.1, KI270442.1, KI270712.1, KI270718.1, KI270719.1, KI270728.1, KI270751.1, KI270754.1, KN538372.1, KQ031389.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF383"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': KI270713.1, KI270721.1, KI270879.1, KV880768.1, KI270733.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF384"
    ## [1] "ZNF407"
    ## [1] "ZNF414"
    ## [1] "ZNF44"
    ## [1] "ZNF446"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000219.1, GL000251.2, GL877875.1, KI270442.1, KI270713.1, KI270733.1, KI270830.1, KI270857.1, KI270869.1, KI270879.1, KN538372.1, KV880768.1, KZ208915.1, KZ208921.1, ML143367.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF451"
    ## [1] "ZNF460"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000218.1, GL000224.1, GL000254.2, GL339449.2, KI270442.1, KI270709.1, KI270712.1, KI270728.1, KI270745.1, KI270751.1, KI270832.1, KI270853.1, KI270856.1, KI270879.1, KI270905.1, KN196484.1, KQ458384.1, KV575244.1, KV880768.1, KZ208915.1, KZ208921.1, KZ559109.1, ML143377.1, GL000194.1, GL000219.1, GL000251.2, GL000252.2, GL000253.2, GL000255.2, GL000257.2, GL383556.1, GL383557.1, GL383577.2, GL383581.2, KI270333.1, KI270711.1, KI270713.1, KI270714.1, KI270721.1, KI270744.1, KI270782.1, KI270827.1, KI270849.1, KI270857.1, KI270861.1, KI270876.1, KI270878.1, KI270894.1, KI270908.1, KN196474.1, KQ031389.1, KV880764.1, KZ208904.1, KZ208908.1, KZ208919.1, ML143360.1, ML143374.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF483"
    ## [1] "ZNF485"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000214.1, GL000218.1, GL000219.1, GL000224.1, KI270442.1, KI270712.1, KI270719.1, KI270721.1, KI270728.1, KI270733.1, KI270734.1, KI270744.1, KI270751.1, KI270754.1, KI270762.1, KI270869.1, KI270879.1, KN196487.1, KN538372.1, KQ031389.1, KV880768.1, ML143377.1, GL000009.2, GL000195.1, GL000216.2, GL000220.1, GL000221.1, GL000225.1, GL000251.2, GL000252.2, GL000253.2, GL000254.2, GL000255.2, GL000256.2, GL339449.2, GL383550.2, GL383555.2, GL383574.1, GL383577.2, GL877875.1, KI270438.1, KI270538.1, KI270709.1, KI270713.1, KI270717.1, KI270718.1, KI270722.1, KI270725.1, KI270726.1, KI270742.1, KI270743.1, KI270745.1, KI270746.1, KI270750.1, KI270753.1, KI270772.1, KI270791.1, KI270815.1, KI270816.1, KI270821.1, KI270830.1, KI270850.1, KI270853.1, KI270854.1, KI270857.1, KI270862.1, KI270871.1, KI270872.1, KI270876.1, KI270877.1, KI270878.1, KI270894.1, KI270897.1, KI270903.1, KI270905.1, KI270908.1, KN538364.1, KV575244.1, KV766196.1, KZ208915.1, KZ208921.1, KZ208922.1, KZ559100.1, KZ559109.1, ML143353.1, ML143364.1, ML143366.1, ML143367.1, ML143371.1, ML143372.1, ML143375.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF490"
    ## [1] "ZNF501"
    ## [1] "ZNF503"
    ## [1] "ZNF510"
    ## [1] "ZNF511"
    ## [1] "ZNF512"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000219.1, GL000254.2, GL339449.2, GL383556.1, GL949746.1, KI270442.1, KI270712.1, KI270713.1, KI270725.1, KI270728.1, KI270733.1, KI270734.1, KI270762.1, KI270764.1, KI270827.1, KI270853.1, KI270879.1, KI270908.1, KV880764.1, KZ208915.1, ML143365.1, ML143372.1, ML143377.1, GL000214.1, GL000218.1, GL000220.1, GL000224.1, KI270330.1, KI270717.1, KI270718.1, KI270719.1, KI270722.1, KI270726.1, KI270744.1, KI270754.1, KI270821.1, KI270850.1, KI270857.1, KI270894.1, KI270905.1, KQ090027.1, KV880768.1, ML143367.1, ML143380.1, GL000194.1, GL000251.2, GL000252.2, GL000255.2, GL383577.2, GL383578.2, GL383580.2, KI270711.1, KI270714.1, KI270721.1, KI270782.1, KI270830.1, KI270842.1, KI270862.1, KI270869.1, KI270878.1, KN538364.1, KZ208906.1, KZ208921.1, KZ208922.1, KZ559105.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF512B"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000214.1, GL000219.1, GL000220.1, GL000251.2, GL000252.2, GL000254.2, GL000255.2, GL339449.2, GL383577.2, GL383578.2, GL383580.2, KI270711.1, KI270712.1, KI270713.1, KI270714.1, KI270721.1, KI270726.1, KI270733.1, KI270734.1, KI270744.1, KI270782.1, KI270821.1, KI270830.1, KI270842.1, KI270850.1, KI270857.1, KI270862.1, KI270869.1, KI270878.1, KI270879.1, KI270894.1, KI270905.1, KI270908.1, KN538364.1, KV880764.1, KV880768.1, KZ208906.1, KZ208915.1, KZ208921.1, KZ208922.1, KZ559105.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF513"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrY
    ##   - in 'y': GL000219.1, KI270438.1, KI270713.1, KI270733.1, KI270879.1, KN196487.1, KV880768.1, GL000220.1, GL000251.2, KI270712.1, KI270728.1, KI270744.1, KI270830.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF526"
    ## [1] "ZNF543"
    ## [1] "ZNF547"
    ## [1] "ZNF548"
    ## [1] "ZNF550"
    ## [1] "ZNF556"
    ## [1] "ZNF557"
    ## [1] "ZNF572"
    ## [1] "ZNF574"
    ## [1] "ZNF576"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000219.1, GL383563.3, KI270712.1, KI270713.1, KI270714.1, KI270721.1, KI270728.1, KI270879.1, KV880768.1, KZ208909.1, KZ208915.1, GL000194.1, GL000214.1, GL000254.2, GL000255.2, GL383533.1, GL383555.2, GL383574.1, KI270442.1, KI270734.1, KI270744.1, KI270750.1, KI270772.1, KI270802.1, KI270857.1, KI270894.1, KI270903.1, KI270908.1, KV766193.1, KZ208908.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF580"
    ## [1] "ZNF598"
    ## [1] "ZNF607"
    ## [1] "ZNF608"
    ## [1] "ZNF609"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000219.1, GL000254.2, GL383577.2, KI270712.1, KI270713.1, KI270721.1, KI270728.1, KI270733.1, KI270861.1, KI270862.1, KI270878.1, KI270879.1, KV880768.1, KZ208915.1, GL000194.1, GL000218.1, GL000220.1, GL000251.2, GL000252.2, GL000255.2, GL339449.2, GL383563.3, GL383574.1, KI270442.1, KI270714.1, KI270734.1, KI270744.1, KI270745.1, KI270754.1, KI270782.1, KI270830.1, KI270832.1, KI270847.1, KI270853.1, KI270857.1, KI270869.1, KI270872.1, KI270894.1, KI270903.1, KI270904.1, KI270905.1, KI270908.1, KN196484.1, KN538364.1, KN538372.1, KQ031389.1, KQ090026.1, KQ458384.1, KV575244.1, KZ208918.1, KZ208921.1, KZ208922.1, ML143371.1, ML143377.1, ML143380.1, ML143381.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF614"
    ## [1] "ZNF616"
    ## [1] "ZNF619"
    ## [1] "ZNF629"
    ## [1] "ZNF639"
    ## [1] "ZNF644"
    ## [1] "ZNF652"
    ## [1] "ZNF660"
    ## [1] "ZNF674"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000220.1, KI270713.1, KI270728.1, KI270733.1, GL000194.1, GL000219.1, KI270712.1, KI270869.1, KN196487.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF678"
    ## [1] "ZNF687"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000195.1, GL000214.1, GL000218.1, GL000219.1, GL000220.1, GL000251.2, GL000252.2, GL000253.2, GL000254.2, GL000255.2, GL000256.2, GL000257.2, GL000258.2, GL339449.2, GL383526.1, GL383555.2, GL383556.1, GL383557.1, GL383563.3, GL383566.1, GL383574.1, GL383577.2, GL383580.2, GL949752.1, JH159146.1, KI270706.1, KI270712.1, KI270713.1, KI270714.1, KI270717.1, KI270719.1, KI270721.1, KI270723.1, KI270726.1, KI270728.1, KI270731.1, KI270733.1, KI270734.1, KI270744.1, KI270745.1, KI270748.1, KI270751.1, KI270753.1, KI270754.1, KI270762.1, KI270772.1, KI270780.1, KI270782.1, KI270804.1, KI270809.1, KI270810.1, KI270813.1, KI270816.1, KI270819.1, KI270821.1, KI270827.1, KI270830.1, KI270831.1, KI270832.1, KI270842.1, KI270846.1, KI270847.1, KI270848.1, KI270850.1, KI270853.1, KI270854.1, KI270856.1, KI270857.1, KI270860.1, KI270861.1, KI270862.1, KI270868.1, KI270869.1, KI270872.1, KI270877.1, KI270878.1, KI270879.1, KI270894.1, KI270897.1, KI270903.1, KI270904.1, KI270905.1, KI270907.1, KI270908.1, KI270910.1, KI270913.1, KI270927.1, KI270938.1, KN196479.1, KN196484.1, KN538364.1, KN538372.1, KQ090026.1, KQ090028.1, KQ458383.1, KQ458384.1, KQ458386.1, KQ759762.1, KV575244.1, KV766196.1, KV880768.1, KZ208906.1, KZ208908.1, KZ208914.1, KZ208915.1, KZ208918.1, KZ208921.1, KZ208922.1, KZ559105.1, KZ559109.1, ML143353.1, ML143358.1, ML143364.1, ML143366.1, ML143367.1, ML143371.1, ML143372.1, ML143377.1, ML143380.1, GL383533.1, GL383550.2, GL383581.2, KI270709.1, KI270765.1, KI270788.1, KI270791.1, KI270792.1, KI270805.1, KI270824.1, KI270845.1, KI270849.1, KI270871.1, KN196472.1, KQ031389.1, KZ208904.1, KZ559113.1, ML143365.1, ML143378.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF691"
    ## [1] "ZNF697"
    ## [1] "ZNF709"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000216.2, GL000219.1, GL000251.2, GL339449.2, GL383563.3, GL383577.2, KI270442.1, KI270712.1, KI270713.1, KI270714.1, KI270725.1, KI270726.1, KI270728.1, KI270730.1, KI270733.1, KI270742.1, KI270744.1, KI270830.1, KI270832.1, KI270857.1, KI270869.1, KI270879.1, KN196484.1, KN538364.1, KQ759762.1, KV575244.1, KV880768.1, KZ208915.1, KZ208921.1, KI270721.1, KI270861.1, KI270908.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF710"
    ## [1] "ZNF713"
    ## [1] "ZNF720"
    ## [1] "ZNF737"
    ## [1] "ZNF740"
    ## [1] "ZNF761"
    ## [1] "ZNF766"
    ## [1] "ZNF768"
    ## [1] "ZNF772"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000195.1, GL000251.2, GL000252.2, GL000254.2, GL339449.2, GL383577.2, KI270713.1, KI270721.1, KI270733.1, KI270861.1, KI270862.1, KI270879.1, KN538372.1, KV880768.1, KZ208921.1, ML143377.1, GL000194.1, GL000218.1, GL000219.1, GL000224.1, GL383574.1, KI270712.1, KI270714.1, KI270728.1, KI270744.1, KI270780.1, KI270816.1, KI270821.1, KI270830.1, KI270832.1, KI270857.1, KI270869.1, KI270894.1, KI270905.1, KI270908.1, KN196484.1, KN538364.1, KQ031389.1, KQ458384.1, KV766196.1, KZ208915.1, KZ208922.1, ML143372.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF773"
    ## [1] "ZNF775"
    ## [1] "ZNF776"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000009.2, GL000194.1, GL000214.1, GL000218.1, GL000219.1, GL000220.1, GL000224.1, GL000225.1, GL000254.2, GL000255.2, GL339449.2, GL383563.3, GL383574.1, KI270442.1, KI270706.1, KI270712.1, KI270713.1, KI270714.1, KI270719.1, KI270721.1, KI270722.1, KI270723.1, KI270728.1, KI270731.1, KI270733.1, KI270734.1, KI270736.1, KI270744.1, KI270745.1, KI270750.1, KI270754.1, KI270772.1, KI270816.1, KI270830.1, KI270832.1, KI270850.1, KI270854.1, KI270860.1, KI270861.1, KI270869.1, KI270879.1, KI270894.1, KI270901.1, KI270905.1, KI270908.1, KN196484.1, KN538364.1, KN538372.1, KV575244.1, KV766196.1, KV880768.1, KZ208915.1, KZ208921.1, KZ208922.1, ML143360.1, ML143364.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF777"
    ## [1] "ZNF781"
    ## [1] "ZNF782"
    ## [1] "ZNF784"
    ## [1] "ZNF786"
    ## [1] "ZNF788"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000195.1, GL000205.2, GL000214.1, GL000216.2, GL000218.1, GL000219.1, GL000224.1, GL000251.2, GL000252.2, GL000253.2, GL000254.2, GL000255.2, GL339449.2, GL383563.3, GL383574.1, GL383577.2, GL383580.2, KI270442.1, KI270711.1, KI270712.1, KI270713.1, KI270714.1, KI270719.1, KI270721.1, KI270728.1, KI270731.1, KI270733.1, KI270734.1, KI270744.1, KI270772.1, KI270780.1, KI270782.1, KI270816.1, KI270821.1, KI270827.1, KI270829.1, KI270830.1, KI270832.1, KI270842.1, KI270850.1, KI270853.1, KI270857.1, KI270861.1, KI270862.1, KI270869.1, KI270875.1, KI270877.1, KI270878.1, KI270879.1, KI270894.1, KI270903.1, KI270905.1, KI270907.1, KI270908.1, KN196484.1, KN538364.1, KN538372.1, KQ090026.1, KQ458384.1, KQ458386.1, KQ759762.1, KV575244.1, KV880768.1, KZ208906.1, KZ208915.1, KZ208921.1, KZ208922.1, KZ559109.1, ML143358.1, ML143366.1, ML143367.1, ML143371.1, ML143372.1, ML143377.1, ML143380.1, GL000225.1, KI270718.1, KI270754.1, KQ031389.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZNF792"
    ## [1] "ZNF800"
    ## [1] "ZNF83"
    ## [1] "ZNF839"
    ## [1] "ZNF883"
    ## [1] "ZNF891"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000008.2, GL000009.2, GL000194.1, GL000195.1, GL000208.1, GL000214.1, GL000216.2, GL000218.1, GL000219.1, GL000220.1, GL000221.1, GL000224.1, GL000225.1, GL000251.2, GL000252.2, GL000253.2, GL000254.2, GL000255.2, GL000256.2, GL339449.2, GL383522.1, GL383533.1, GL383556.1, GL383563.3, GL383566.1, GL383574.1, GL383577.2, GL383581.2, GL877875.1, GL949742.1, JH159146.1, JH159147.1, JH159148.1, KI270438.1, KI270442.1, KI270538.1, KI270706.1, KI270707.1, KI270709.1, KI270712.1, KI270713.1, KI270714.1, KI270718.1, KI270719.1, KI270720.1, KI270721.1, KI270722.1, KI270723.1, KI270724.1, KI270726.1, KI270728.1, KI270731.1, KI270732.1, KI270733.1, KI270734.1, KI270742.1, KI270743.1, KI270744.1, KI270745.1, KI270746.1, KI270750.1, KI270751.1, KI270753.1, KI270754.1, KI270762.1, KI270772.1, KI270782.1, KI270804.1, KI270805.1, KI270810.1, KI270816.1, KI270819.1, KI270821.1, KI270823.1, KI270824.1, KI270830.1, KI270831.1, KI270832.1, KI270845.1, KI270846.1, KI270847.1, KI270850.1, KI270853.1, KI270854.1, KI270857.1, KI270860.1, KI270861.1, KI270869.1, KI270871.1, KI270872.1, KI270876.1, KI270877.1, KI270878.1, KI270879.1, KI270894.1, KI270896.1, KI270897.1, KI270903.1, KI270904.1, KI270905.1, KI270907.1, KI270908.1, KI270910.1, KI270913.1, KN196484.1, KN196487.1, KN538364.1, KN538372.1, KQ031389.1, KQ090026.1, KQ458386.1, KQ759762.1, KV575244.1, KV766196.1, KV880764.1, KV880768.1, KZ208904.1, KZ208906.1, KZ208908.1, KZ208915.1, KZ208918.1, KZ208921.1, KZ208922.1, KZ559109.1, ML143353.1, ML143364.1, ML143365.1, ML143366.1, ML143367.1, ML143371.1, ML143372.1, ML143377.1, ML143378.1, ML143380.1, GL000205.2, KI270711.1, KI270725.1, KI270764.1, KI270862.1, KQ458383.1, KV766197.1, ML143341.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZSCAN20"
    ## [1] "ZSCAN22"
    ## [1] "ZSCAN29"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000218.1, GL000251.2, GL000252.2, GL000254.2, GL000255.2, GL339449.2, GL383574.1, GL383577.2, KI270713.1, KI270721.1, KI270733.1, KI270829.1, KI270857.1, KI270879.1, KI270894.1, KI270905.1, KN196484.1, KV575244.1, KV880768.1, KZ208915.1, GL000194.1, GL000253.2, KI270908.1, KZ208922.1, ML143350.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZSCAN31"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000214.1, GL000218.1, GL000219.1, GL000220.1, GL000224.1, GL000225.1, GL000251.2, GL000252.2, GL000253.2, GL000254.2, GL000255.2, GL339449.2, GL383574.1, GL383577.2, GL383580.2, KI270713.1, KI270714.1, KI270721.1, KI270728.1, KI270733.1, KI270744.1, KI270782.1, KI270816.1, KI270830.1, KI270861.1, KI270869.1, KI270879.1, KI270894.1, KI270905.1, KI270908.1, KN196483.1, KN196484.1, KV880768.1, KZ208915.1, KZ208922.1, ML143377.1, ML143380.1, KI270904.1, KZ208921.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## [1] "ZSCAN9"
    ## [1] "ZXDC"
    ## [1] "ZZZ3"

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000194.1, GL000218.1, GL000219.1, GL000220.1, GL000224.1, GL000225.1, GL339449.2, GL383563.3, KI270589.1, KI270712.1, KI270713.1, KI270726.1, KI270728.1, KI270733.1, KI270742.1, KI270744.1, KI270762.1, KI270830.1, KI270857.1, KI270869.1, KI270894.1, KN538372.1, KQ031389.1, KZ208915.1, ML143353.1, ML143362.1, ML143372.1, GL000195.1, GL000216.2, GL000251.2, GL000255.2, GL383577.2, KI270438.1, KI270719.1, KI270721.1, KI270754.1, KI270756.1, KI270853.1, KI270862.1, KI270879.1, KI270904.1, KI270905.1, KV880768.1, ML143377.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

``` r
non_superbinder_metaplot_df$gene_type <- "non_super_binder"
superbinder_metaplot_df$gene_type <- "superbinder"
combined_super_binder_metaplot_profile <- bind_rows(non_superbinder_metaplot_df, superbinder_metaplot_df)

#ggplot(combined_super_binder_metaplot_profile, aes(x = x, y = dens, color = gene_type)) +
  #geom_vline(xintercept = 0, lty = 2) + 
  #geom_line(size = 1.5) +
  #facet_wrap(~ dbp, scales = "free_y") +
  #ggtitle("Promoter Metaplot") + 
  #scale_x_continuous(breaks = c(-1000, 0, 1000),
                     #labels = c("-1kb", "TSS", "+1kb"),
                     #name = "") + 
  #ylab("Peak frequency") +
  #xlab("Distance from TSS") +  # add label to x-axis
  #scale_color_manual(values = c("#424242", "#a8404c"))
#ggsave("/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/final_project/analysis/metaplots/superbinder_metaplot.pdf")


# We can now write these out for safekeeping / use in other analyses
#write_csv(promoter_features_df, "/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/final_project/analysis/04_binding_vs_expression/results/promoter_features_df.csv")
#write_csv(tpm, "/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/final_project/analysis/04_binding_vs_expression/results/mean_tpm_per_condition.csv")
#write_csv(samplesheet, "/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/final_project/analysis/04_binding_vs_expression/results/samplesheet.csv")
```
