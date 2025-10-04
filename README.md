# Assessment-4  RNA-Seq Analysis, Growth Data Analysis, and Biological Sequence Diversity
Author: Pavan Manda
Project Overview
This repository contains a comprehensive bioinformatics analysis project divided into three main sections: gene expression profiling from RNA-seq data, statistical analysis of tree growth patterns, and comparative genomics between Escherichia coli and Deinococcus radiodurans. The analysis demonstrates proficiency in R programming, data manipulation, statistical testing, and biological sequence analysis.
Table of Contents
Installation and Dependencies
Data Sources
Part 1: Gene Expression Analysis
Part 2: Growth Data Analysis
Part 3: Biological Sequence Diversity
Results Summary
Usage Instructions
Contact Information
Installation and Dependencies
Required R Packages 
install.packages("R.utils")    # Version 2.13.0 - File compression/decompression 
install.packages("seqinr")     # Version 4.2-x - Sequence analysis toolkit 
install.packages("ggplot2")    # Version 3.x - Data visualization
System Requirements
R version 4.1 or higher
RStudio (recommended for R Markdown rendering)
Minimum 2GB RAM for sequence analysis
Internet connection for data download
Data Sources
Part 1: Gene Expression Data
File: gene_expression.tsv
URL: https://raw.githubusercontent.com/ghazkha/Assessment4/main/gene_expression.tsv
Format: Tab-separated values with gene identifiers and expression counts across three samples
Size: 60,000+ genes
Part 2: Growth Data
File: growth_data.csv
URL: https://raw.githubusercontent.com/ghazkha/Assessment4/main/growth_data.csv
Format: CSV with tree measurements from 2005-2020
Size: 100 trees (50 per site)
Part 3: Genomic Sequences
E. coli K-12 MG1655 (GCA_000005845): Ensembl Bacteria Release 53
D. radiodurans R1 (GCA_000008565): Ensembl Bacteria Release 53
Format: Compressed FASTA (.fa.gz) containing coding DNA sequences
Part 1: Gene Expression Analysis
Purpose
This section analyzes RNA-seq count data to identify expression patterns, highly expressed genes, and overall distribution of gene activity across samples. The workflow involves data import, statistical calculations, and visualization of expression levels.
Workflow
Question 1: Data Import
The gene expression data is imported from a tab-separated file using read.table(). Gene identifiers from the first column are assigned as row names, creating a matrix where rows represent genes and columns represent samples. This structure facilitates row-wise operations for calculating expression statistics.

Question 2: Mean Expression Calculation
For each gene, mean expression is calculated across all three samples using rowMeans(). This provides a single representative value for each gene's activity level, enabling comparison and ranking. The mean values are added as a new column to the data frame.

Question 3: Top Gene Identification
Genes are sorted by mean expression in descending order using order(), and the top 10 are extracted. Results show that mitochondrial genes (MT-CO1, MT-ND4, MT-CO3) dominate with expression levels ranging from 232,187 to 529,317, indicating high metabolic activity in the sampled tissue.

Question 4: Low Expression Threshold
Using a logical comparison (gene$mean < 10), the analysis identifies 35,988 genes with mean expression below 10. This represents approximately 60% of all genes, indicating that most genes show minimal or no expression in these particular samples.

Question 5: Distribution Visualization
A histogram with 30 bins visualizes the frequency distribution of mean expression values. The distribution is heavily right-skewed, with the vast majority of genes clustered at low expression levels and a long tail extending to high values, characteristic of typical RNA-seq data.

Key Outputs
Table of top 10 expressed genes with mean values
Count of lowly expressed genes: 35,988
Histogram showing right-skewed expression distribution

Part 2: Growth Data Analysis
Purpose
This section examines tree circumference measurements over 15 years at two geographical sites (northeast and southwest) to assess growth patterns and test for statistically significant differences between locations.
Workflow
Question 6: Data Import and Structure
The CSV file is loaded using read.csv(), containing six columns: Site, TreeID, and circumference measurements from 2005, 2010, 2015, and 2020. The dataset tracks 100 individual trees (50 per site) over four time points.

Question 7: Descriptive Statistics
Using aggregate() with the formula syntax, mean and standard deviation are calculated for both sites at baseline (2005) and endpoint (2020). Northeast trees started slightly larger (mean 5.29 cm vs 4.86 cm) and ended substantially larger (mean 54.23 cm vs 45.60 cm). Standard deviations increased dramatically from ~1 cm to 17-25 cm, indicating heterogeneous growth rates within each site.

Question 8: Visual Comparison
Side-by-side boxplots created with base R graphics display the distribution of circumferences at start and end points. The plots reveal substantial growth at both sites, with northeast showing higher median values and greater variability (longer whiskers and more outliers) by 2020.

Question 9: Growth Rate Calculation
Ten-year growth (2010-2020) is computed by subtracting 2010 values from 2020 values, creating a new variable. Mean growth aggregated by site shows northeast with 42.94 cm and southwest with 35.49 cm, a difference of 7.45 cm favoring the northeast location.

Question 10: Statistical Hypothesis Testing
A Welch's two-sample t-test compares growth between sites, accounting for unequal variances. The test yields t = 1.89 with p = 0.062, which exceeds the conventional α = 0.05 threshold. While the difference is not statistically significant at the 5% level, the p-value suggests marginal significance that might warrant further investigation with larger sample sizes.

Key Outputs
Mean circumference at start: Northeast 5.29 cm, Southwest 4.86 cm
Mean circumference at end: Northeast 54.23 cm, Southwest 45.60 cm
Mean 10-year growth: Northeast 42.94 cm, Southwest 35.49 cm
Statistical test: t = 1.89, p = 0.062 (not significant at α = 0.05)

Part 3: Biological Sequence Diversity
Purpose
This section performs comparative genomic analysis between E. coli (a model mesophilic bacterium) and D. radiodurans(an extremophile known for radiation resistance) to explore differences in genome size, sequence composition, codon usage, and protein motif prevalence.
Workflow
Question 1: CDS Enumeration
Compressed FASTA files are downloaded from Ensembl using download.file(), decompressed with gunzip(), and read with seqinr::read.fasta(). Counting sequence objects reveals E. coli contains 4,239 coding sequences while D. radiodurans has 3,101, indicating E. coli has 27% more genes despite similar lifestyle constraints.

Question 2: Total Coding Content
Summing lengths of all sequences using sapply() and sum() shows E. coli has 3,978,528 bp of coding DNA compared to D. radiodurans with 2,879,559 bp. This 1.1 Mb difference (38% more) reflects both the greater gene count and slightly longer average gene length in E. coli.

Question 3: Sequence Length Distribution
Calculating mean and median lengths reveals similar values between organisms (E. coli: mean 939 bp, median 831 bp; D. radiodurans: mean 929 bp, median 789 bp). Boxplots created with ggplot2 show comparable distributions with similar interquartile ranges and outliers, suggesting conserved gene length distributions despite genomic differences.

Question 4: Nucleotide and Amino Acid Composition
Custom functions calculate base frequencies by unlisting all sequences and creating frequency tables. D. radioduransexhibits dramatically higher GC content (~67%) compared to E. coli (~51%), a hallmark adaptation to radiation damage where GC pairs provide greater stability. After translating sequences with seqinr::translate(), amino acid frequency analysis shows both organisms favor leucine, alanine, and valine, though D. radiodurans shows subtle enrichment in amino acids encoded by GC-rich codons.

Question 5: Codon Usage Bias
A custom function extracts all codons by parsing sequences in triplets, creating frequency tables for all 64 codons. The top 20 codon comparison plot clearly demonstrates that E. coli preferentially uses AT-rich codons (like AAA for lysine) while D. radiodurans favors GC-rich alternatives (like GAG for glutamate, CTG for leucine). This codon bias reflects both genomic GC content and optimization for available tRNA pools in each organism.

Question 6: Protein Motif Analysis
By extracting 4-amino acid sliding windows from translated proteins, k-mer frequency analysis identifies over- and under-represented motifs. D. radiodurans shows 2-3× enrichment of hydrophobic, stability-promoting motifs (LAAL, AALA, AALL) composed primarily of alanine and leucine. These sequences likely contribute to protein structural stability under radiation stress. Conversely, both organisms show extremely low frequencies of tyrosine-rich motifs (YYVY, YYWD, YYWW), possibly due to oxidative susceptibility of aromatic residues or structural constraints. Bar plots visualize these differences, with D. radiodurans showing clear over-representation of stability motifs compared to E. coli.

Key Outputs
CDS count: E. coli 4,239, D. radiodurans 3,101
Total coding DNA: E. coli 3.98 Mb, D. radiodurans 2.88 Mb
GC content: E. coli 51%, D. radiodurans 67%
Codon bias: D. radiodurans strongly favors GC-rich codons
Motif enrichment: Hydrophobic stability motifs 2-3× higher in D. radiodurans

Results Summary
Gene Expression
Analysis of 60,000+ genes revealed that mitochondrial genes dominate expression profiles, with 9 of the top 10 genes being mitochondrial-encoded. The majority of genes (60%) show minimal expression, consistent with tissue-specific and condition-specific gene regulation.
Tree Growth
Northeast site trees showed 17% greater growth over 10 years compared to southwest trees. While not reaching statistical significance (p = 0.062), this trend suggests potential environmental or soil differences warranting further investigation.
Sequence Diversity
D. radiodurans demonstrates clear genomic adaptations to extreme environments through elevated GC content (67% vs 51%), corresponding codon bias toward GC-rich codons, and enrichment of hydrophobic protein motifs that enhance structural stability under radiation stress.
 
