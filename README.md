# Identification-of-human-ape-specific-splice-sites-and-the-associated-variants

1. Quality control of RNA-seq data
To ensure that only high-quality data were used for downstream analyses, RNA-seq reads with low quality (-q 30) were filtered out using fastp

2. RNA-seq reads alignment
We used HISAT2 to align the cleaned reads from each species to the corresponding reference genome (hg38 for human, panTro6 for chimpanzee, rheMac10 for rhesus macaque, and mm39 for mouse).

3. Identify candidate splice site by extracting jucntion reads
We extracted junction reads from the aligned BAM files using samtools and custom shell scripts to identify splice sites utilized in each species.

4. Cross-species genomic coordinate comparisons using LiftOver
A. We converted splice site positions from one species to their orthologous positions in another by using LiftOver.
B. We converted positions back to the original genome by using LiftOver.
C. We retained only the splice sites with consistent positions in both forward and reverse mappings.

5. Identify candidate species-specific splice sites
Based on the splice sites usage and cross-species positional correspondence obtained in the previous steps, we extracted candidate human-specific splice sites labeled as “human_nochimp_norhesus_nomouse”, and candidate ape-specific splice sites labeled as “human_chimp_norhesus_nomouse”.
In addition, species-specific loss of splice sites were also identified in this step: candidate human-specific loss of splice sites labeled as “nohuman_chimp_rhesus_mouse”, and candidate ape-specific loss of splice sites labeled as “nohuman_nochimp_rhesus_mouse”.

6. Pipeline for filtering species-specific splice sites
To ensure the accuracy of the identified candidate human/ape-specific splice sites, we performed filtering based on the following steps.
A. First, we retained junction reads in which at least one end aligned with annotated exon boundaries.
B. Next, we applied a series of filtering steps: Read coverage (number of junction reads ≥3), read quality (overhang ≥5 bp), the generality of splicing events (observed in at least three samples since there are six replicates),
C. Finally, we retained the sites differing by a single nucleotide within the dinucleotide motif compared to their orthologous sites in other species. 
