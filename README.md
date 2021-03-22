# primerg
Primer design for multiple gRNA hits within a continuous genome

# Background and problem definition
Members from the same gene cluster can perform redundant functions. Functional knock-out of all genes within a cluster using CRISPR-Cas9 is an approach to evaluate the role of the genes. There is already bioinformatic software to find the minimum number of gRNAs to target an entire gene cluster (Multitargeter: https://multicrispr.net/). Following a CRISPR-Cas9 knock-out experiment is the confirmation of mutation at Cas9 cleavage sites using NGS. However, there is currently no bioinformatic solution to resolve primer design for gene clusters by minimum gRNA targeting. Such a solution would involve designing unique primers per genome target site, if not ensuring the amplicon sequences derived from different genomic sites contain suitable polymorphism for NGS reads to be classified to each target. The primer design problem is complicated by targeted loci that are repeats or copies of each other and the exponential increase of complexity as the number of targeted loci increases.

# Solution
primerg designs primers per Cas9-targeted locus using primer3; subsequently, primer3 evaluates each primer and its amplicons for uniqueness amongst all primers of every site. The user can then pick primers that binds specifically to one locus or primers that bind to multiple loci but generates unique amplicons.
