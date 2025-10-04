Assessment 4: RNA-Seq Analysis, Growth Data Analysis, and Biological Sequence Diversity
Author: Pavan Manda
Table of Contents
Project Overview
Installation and Dependencies
Data Sources
Part 1: Gene Expression Analysis
Part 2: Growth Data Analysis
Part 3: Biological Sequence Diversity
Results Summary

Project Overview
This repository contains a comprehensive bioinformatics analysis project demonstrating proficiency in R programming, statistical analysis, and computational biology. The work is divided into three distinct sections, each addressing different aspects of biological data analysis. The first section focuses on gene expression profiling using RNA-seq data to identify highly expressed genes and understand expression patterns across samples. The second section performs statistical analysis of tree growth measurements collected over 15 years at two different geographical sites. The third section conducts comparative genomics analysis between two bacterial species with contrasting ecological niches to explore genomic adaptations and sequence diversity.

The analysis workflow incorporates multiple R packages and demonstrates skills in data import and manipulation, statistical testing, data visualization, and biological sequence analysis. All code is fully documented with clear explanations of methodology and interpretation of results. The project uses publicly available datasets and downloads sequence data directly from Ensembl, ensuring reproducibility across different computing environments.

Installation and Dependencies

Required R Packages
The analysis requires three main R packages for data processing and visualization. The R.utils package (version 2.13.0 or higher) provides functionality for file compression and decompression, essential for handling the compressed FASTA files downloaded from genomic databases. The seqinr package (version 4.2 or higher) serves as the primary toolkit for biological sequence analysis, offering functions to read FASTA files, translate DNA sequences to proteins, and perform various sequence manipulations. The ggplot2 package (version 3.0 or higher) enables creation of publication-quality data visualizations with consistent styling across all plots.

r
install.packages("R.utils") install.packages("seqinr") install.packages("ggplot2")

System Requirements
The analysis has been tested and runs successfully on R version 4.1 or higher. While the code can be executed in any R environment, RStudio is recommended for working with R Markdown files and generating HTML reports through the knit function. The sequence analysis components, particularly k-mer extraction and frequency calculations, require a minimum of 2GB RAM to process the complete genomic datasets efficiently. An active internet connection is necessary during initial setup to download gene expression data, growth measurements, and coding sequence files from remote repositories.

Data Sources

Part 1: Gene Expression Data
The gene expression dataset consists of RNA-seq count data stored in a tab-separated values file named gene_expression.tsv. This file is publicly available from the GitHub repository at https://raw.githubusercontent.com/ghazkha/Assessment4/main/gene_expression.tsv. The dataset contains expression measurements for over 60,000 genes across three biological samples, with gene identifiers in the first column and raw count values in subsequent columns. This format allows for straightforward import into R and facilitates row-wise operations for calculating gene-level statistics.

Part 2: Growth Data
Tree circumference measurements are provided in a CSV file named growth_data.csv, accessible at https://raw.githubusercontent.com/ghazkha/Assessment4/main/growth_data.csv. The dataset tracks 100 individual trees over 15 years, with 50 trees measured at a northeast site and 50 at a southwest site. Each tree has circumference measurements recorded at four time points: 2005, 2010, 2015, and 2020. This longitudinal design enables analysis of growth trajectories and comparison of growth patterns between the two geographical locations.

Part 3: Genomic Sequences
Coding DNA sequences for comparative analysis are obtained from the Ensembl Bacteria database, Release 53. The dataset for Escherichia coli strain K-12 MG1655 (assembly accession GCA_000005845) represents a well-characterized model organism used extensively in molecular biology research. The Deinococcus radiodurans R1 genome (assembly accession GCA_000008565) provides comparison data from an extremophile bacterium known for exceptional resistance to radiation damage. Both datasets are provided as compressed FASTA files containing all protein-coding sequences from each organism's genome.

Part 1: Gene Expression Analysis
Purpose and Approach
This analysis examines RNA-seq count data to characterize gene expression patterns, identify the most highly expressed genes, and visualize the overall distribution of gene activity. The workflow progresses systematically through data import, statistical calculations, ranking operations, and visualization to provide comprehensive insights into the transcriptional landscape of the sampled tissue.

Data Import and Processing
The gene expression data is imported using the read.table function, which reads tab-separated files and automatically parses the header row to extract column names. After import, gene identifiers from the first column are assigned as row names, transforming the data frame into a gene-by-sample matrix. This structure is particularly advantageous for RNA-seq analysis because it allows efficient application of row-wise functions to calculate gene-level statistics. The first column containing gene names is then removed since this information is now stored in the row names, leaving only the numerical expression values.

Mean Expression Calculation
For each gene in the dataset, mean expression is calculated across all three samples using the rowMeans function. This function efficiently computes row-wise averages, producing a single representative value for each gene's expression level. The mean provides a robust summary statistic that captures the typical expression level while being less sensitive to individual sample variation than any single measurement. These mean values are added as a new column to the data frame, enabling subsequent ranking and filtering operations.

Identification of Highly Expressed Genes
Genes are ranked by their mean expression levels using the order function with the decreasing parameter set to TRUE. This produces indices sorted from highest to lowest expression, which are then used to reorder the entire data frame. The top 10 genes are extracted and presented in a clean summary table. Analysis reveals that mitochondrial genes dominate the highest expression category, with nine of the top 10 genes being mitochondrial-encoded. Specifically, MT-CO1, MT-ND4, and MT-CO3 show expression levels ranging from 232,187 to 529,317, indicating exceptionally high metabolic activity in the tissue from which RNA was extracted. This pattern is consistent with tissues having high energy demands, such as muscle or heart tissue.

Low Expression Gene Counting
To characterize the proportion of genes showing minimal activity, a logical comparison identifies all genes with mean expression below 10. The sum function counts TRUE values in this logical vector, revealing that 35,988 genes fall below this threshold. This represents approximately 60% of all genes in the dataset, which is consistent with the biological principle that most genes in any given cell type or condition are either silent or expressed at very low levels. Only a subset of genes specific to the tissue type, developmental stage, and environmental conditions will be actively transcribed at high levels.

Expression Distribution Visualization
A histogram with 30 bins visualizes the frequency distribution of mean expression values across all genes. The resulting plot shows a heavily right-skewed distribution, where the vast majority of genes cluster at low expression levels with a long tail extending toward high values. This pattern is characteristic of RNA-seq data across diverse organisms and tissues, reflecting the fundamental organization of cellular gene expression where a small number of highly expressed genes account for most of the transcriptional activity while the majority of genes contribute minimally to the total RNA pool.

Part 2: Growth Data Analysis

Purpose and Approach
This section examines tree circumference measurements collected over 15 years at two geographically distinct sites to assess growth patterns and test whether observed differences between sites are statistically significant. The analysis incorporates descriptive statistics, visual comparisons, growth rate calculations, and formal hypothesis testing to provide comprehensive insights into tree development at each location.

Data Import and Structure
The CSV file is loaded using read.csv, which automatically detects comma delimiters and processes the header row. The resulting data frame contains six columns: Site (identifying northeast or southwest location), TreeID (unique identifier for each tree), and four circumference measurements recorded in 2005, 2010, 2015, and 2020. This structure represents a longitudinal dataset with 100 trees tracked across four time points, enabling both cross-sectional comparisons at specific times and within-tree growth trajectory analysis.

Descriptive Statistics
Using the aggregate function with formula syntax, mean and standard deviation are calculated separately for each site at the study's beginning (2005) and end (2020). The formula notation (Circumf_2005_cm ~ Site) specifies that circumference should be summarized by site, with the FUN parameter determining whether mean or standard deviation is calculated. Results show that northeast trees started slightly larger with a mean circumference of 5.29 cm compared to 4.86 cm for southwest trees. By 2020, this difference had amplified substantially, with northeast trees reaching 54.23 cm versus 45.60 cm for southwest trees. Standard deviations increased dramatically from approximately 1 cm at baseline to 17-25 cm at endpoint, indicating that individual trees within each site exhibited highly heterogeneous growth rates over the 15-year period.

Visual Comparison
Side-by-side boxplots created using base R graphics provide visual comparison of circumference distributions at the start and end of the study period. The par function with mfrow parameter arranges two plots horizontally, allowing direct visual comparison. Boxplots reveal the full distribution including median (center line), interquartile range (box), whiskers extending to extreme values, and outliers plotted as individual points. Visual inspection confirms substantial growth at both sites over 15 years. The northeast site shows higher median values and considerably greater variability by 2020, evidenced by longer whiskers and more outliers, suggesting that environmental conditions or genetic factors at this location produce more variable growth outcomes.

Growth Rate Calculation
Ten-year growth is calculated by subtracting 2010 circumference values from 2020 values, creating a new variable that captures growth during the second decade of the study. This approach focuses on the period when trees were more established and potentially showing different growth dynamics than the initial planting period. Mean growth aggregated by site using the aggregate function shows northeast trees growing 42.94 cm compared to 35.49 cm for southwest trees, a difference of 7.45 cm or approximately 17% greater growth at the northeast location. This substantial difference suggests potential environmental factors such as soil quality, water availability, temperature patterns, or sunlight exposure that favor growth at the northeast site.

Statistical Hypothesis Testing
A Welch's two-sample t-test formally tests whether the observed growth difference between sites is statistically significant or could reasonably occur by chance. The Welch variant is chosen over the standard Student's t-test because it does not assume equal variances between groups, making it more appropriate given the substantial difference in standard deviations observed between sites. The test yields a t-statistic of 1.89 with 87.98 degrees of freedom and a p-value of 0.062. Since this p-value exceeds the conventional significance threshold of 0.05, we cannot reject the null hypothesis of equal mean growth at the 5% significance level. However, the p-value is marginally significant (below 0.10), suggesting the data provide some evidence for a difference that might reach statistical significance with a larger sample size or longer observation period.

Part 3: Biological Sequence Diversity

Purpose and Approach
This section performs comparative genomic analysis between Escherichia coli, a well-characterized mesophilic bacterium that grows optimally at moderate temperatures, and Deinococcus radiodurans, an extremophile bacterium renowned for surviving extraordinarily high levels of ionizing radiation. The analysis explores differences in genome organization, base composition, codon usage preferences, and protein sequence motifs to understand the genomic basis of adaptation to extreme environments.

Coding Sequence Enumeration
Compressed FASTA files containing all coding DNA sequences are downloaded from the Ensembl Bacteria database using the download.file function. The files are automatically decompressed using gunzip from the R.utils package, then read into R using the read.fasta function from seqinr. Counting the number of sequence objects in each list reveals that E. coli contains 4,239 coding sequences while D. radiodurans has 3,101, indicating E. coli possesses 27% more genes. This difference is somewhat surprising given that both organisms are free-living bacteria with similar generation times, but it reflects E. coli's more complex metabolic capabilities and versatile adaptation to diverse environments including mammalian intestinal tracts.

Total Coding DNA Content
Summing the lengths of all sequences using sapply to extract lengths followed by sum provides the total number of base pairs dedicated to protein-coding regions. E. coli has 3,978,528 bp of coding DNA compared to D. radiodurans with 2,879,559 bp, representing a difference of approximately 1.1 megabases or 38% more coding content in E. coli. This difference reflects both the greater number of genes and slightly longer average gene length in E. coli. The coding content represents a substantial fraction of each organism's total genome size, consistent with the compact genomes typical of bacteria that minimize non-coding DNA.

Sequence Length Distribution
Calculating mean and median coding sequence lengths reveals remarkably similar values between the two organisms. E. coli genes average 939 bp with a median of 831 bp, while D. radiodurans genes average 929 bp with a median of 789 bp. The similarity of mean and median values, along with boxplots showing comparable interquartile ranges and outlier patterns, suggests that despite differences in gene number and total coding content, both organisms maintain similar distributions of gene lengths. This conservation likely reflects fundamental constraints on protein size and function that transcend ecological niche, with most proteins requiring similar numbers of amino acids to form functional domains regardless of the organism.

Nucleotide and Amino Acid Composition
Custom functions calculate base frequencies by unlisting all sequences to create a single vector of nucleotides, then generating frequency tables normalized by total length. This analysis reveals one of the most striking differences between the organisms: D. radiodurans exhibits dramatically higher GC content at approximately 67% compared to E. coli's 51%. This 16 percentage point difference is functionally significant because GC base pairs form three hydrogen bonds compared to two for AT pairs, providing greater thermal stability and resistance to radiation-induced strand breaks. The elevated GC content in D. radiodurans represents a key genomic adaptation enabling survival in high-radiation environments where DNA damage is frequent. After translating coding sequences to amino acids using the seqinr translate function, frequency analysis shows both organisms favor leucine, alanine, and valine as the most common amino acids. However, D. radiodurans shows subtle enrichment in amino acids encoded by GC-rich codons, a direct consequence of its elevated genomic GC content.

Codon Usage Bias
A custom function extracts all codons by parsing sequences into non-overlapping triplets, creating comprehensive frequency tables for all 64 possible codons. The top 20 codon comparison plot provides clear visual evidence of divergent codon preferences between organisms. E. coli preferentially uses AT-rich codons such as AAA for lysine, while D. radiodurans strongly favors GC-rich alternatives such as GAG for glutamate and CTG for leucine. This codon bias directly reflects each organism's genomic GC content and represents optimization for the available tRNA pools. Organisms evolve codon preferences that match their tRNA gene complement, ensuring efficient translation of the most frequently used codons. The strong GC bias in D. radiodurans codons means its tRNA genes have coevolved to recognize these codons preferentially, creating a matched system optimized for translating its GC-rich genes.

Protein Motif Analysis
By extracting all possible 4-amino acid windows from translated proteins using a sliding window approach, k-mer frequency analysis identifies motifs that are over-represented or under-represented in each organism. D. radiodurans shows 2-3 fold enrichment of hydrophobic motifs composed primarily of alanine and leucine, including LAAL, AALA, and AALL. These sequences are small, hydrophobic amino acids that promote tight packing of protein cores and enhance structural stability. The enrichment of these motifs in D. radiodurans likely contributes to protein stability under radiation stress, where maintaining proper protein folding is essential for cellular survival. Conversely, both organisms show extremely low frequencies of tyrosine-rich motifs such as YYVY, YYWD, and YYWW. These sequences containing multiple aromatic tyrosine residues may be disfavored due to oxidative susceptibility, structural constraints, or aggregation propensity. The bar plots visualizing these comparisons clearly demonstrate D. radiodurans' preference for stability-promoting motifs, consistent with adaptation to extreme environments requiring robust protein function under stress conditions.

Results Summary

Gene Expression Insights
Analysis of over 60,000 genes revealed that mitochondrial-encoded genes dominate the high expression category, with nine of the top 10 genes originating from the mitochondrial genome. This pattern indicates the sampled tissue has high energy demands requiring substantial oxidative phosphorylation capacity. The majority of genes, approximately 60%, show minimal expression below 10 counts, consistent with tissue-specific and condition-specific regulation where only relevant genes are transcriptionally active.

Tree Growth Patterns
Northeast site trees demonstrated 17% greater growth over the 10-year analysis period compared to southwest site trees, with mean growth of 42.94 cm versus 35.49 cm. While this difference did not reach statistical significance at the conventional 5% level (p = 0.062), the marginal significance suggests potential environmental or soil chemistry differences between sites that warrant further investigation. The increased variability observed at the northeast site indicates heterogeneous growing conditions or genetic diversity affecting individual tree responses.

Genomic Adaptations
D. radiodurans demonstrates clear genomic adaptations to extreme radiation environments through multiple coordinated changes. The elevated GC content of 67% compared to 51% in E. coli provides enhanced DNA stability through stronger base pairing. This GC bias extends to codon usage, with strong preference for GC-rich codons requiring coordinated evolution of the tRNA gene complement. At the protein level, enrichment of hydrophobic alanine-leucine motifs by 2-3 fold likely enhances protein stability under radiation stress. These multilevel adaptations illustrate how organisms evolve integrated solutions to environmental challenges through coordinated changes across genomic, transcriptional, and protein levels.
 
