# PRIMERg
Designs deep-sequencing/ NGS primers to assess editing outcomes from multiple gRNA targets within a continuous genome. 


## Background and problem definition
Members from the same gene cluster can perform redundant functions. Functional knock-out of genes within a cluster using CRISPR-Cas9 is an approach to evaluate the role of such genes. The bioinformatic support for the design of gRNAs to achieve multi-gene knock-out is currently available: e.g. Multitargeter: https://multicrispr.net/. The subsequent step after the CRISPR-Cas9 knock-out experiment is the confirmation of mutation at Cas9 cleavage sites by deep-sequencing. 

Deep-sequencing to evaluate the editing status of each gene requires unique amplicons per genomic target site. Hence, ideally, specific primers per genomic target site should be designed. When specific primers are impossible to generate, the next acceptable primer design solution at least ensures that the amplicon sequences derived from different genomic sites contain suitable polymorphism for NGS reads to be assigned to each target. 

Manual design of unique primers is decreasingly efficient as the number of targeted loci increases, especially when many gene members are (nearly) identical copies of each other. There are currently no bioinformatics solutions to resolve primer design for gene clusters CRISPR-Cas knockout experiments.


## PRIMERg as a solution
PRIMERg automatically designs deep-sequencing primers per Cas9-targeted locus using primer3 in gene clusters with members of high sequence similarity. Two sets of primers are produced: the first PCR product is large (default: 1200-2000 bp) and flanks the region to be sequenced; the second PCR produces the amplicon to be sequenced (default: 250-280 bp). The purpose of the first PCR is to restrict the template for the second PCR to increase primer specificity. 
PRIMERg will design specific primers, if possible. The uniqueness of the second PCR primer (whether the amplicon to be sequenced has “unique” SNPs to help differentiate the desired amplicon from others produced by non-specific binding) is checked and marked in the output.


## PRIMERg algorithm
![Alt text](algorithm_map.png?raw=true)


## Installations
PRIMERg runs on Linux. For Windows users, you can download Ubuntu (tested on Ubutun 20.04 but new versions should work. Installation guide: https://ubuntu.com/tutorials/ubuntu-on-windows#1-overview)

Python 3.8.5  on Linux 

Before running the rest of the codes here on Linux, update the package list:
```sudo apt update```

BLAST 2.9.0+ on Linux 
 ```sudo apt-get install ncbi-blast+```

Primer3 2.4.0-2 on Linux 
```sudo apt-get install primer3```

Primer3 wrapper for python3 
```sudo apt install python3-pip``` then ```pip install primer3-py```

Biopython for python3 
```pip install biopython```

pandas for python3 
```pip install pandas```

Openpyxl for python3
```pip install openpyxl```


## Quick start
1.	Put primerg.py and fasta file of your [1] genomic template and [2] gRNA sequences (5’-3’ direction) in the same directory.
2.	Open primerg.py with text editor or python IDE. Edit section on “Required input” (directory and input file names) and “Optional input”.
3.	Navigate to the directory in Ubuntu. E.g. if your directory named “PRIMERg” is on desktop: 
```cd /mnt/c/Users/cherw/Desktop/PRIMERg```
5.	```python3 primerg.py```

## Interpreting output
Columns:
 cleavage_pos	| Expected Cas9 cleavage position on user-supplied genomic template
 gRNA_pam	| gRNA + PAM sequence
 F1	| Forward primer for first PCR
 R1 | Reverse primer for first PCR
 F2	| Forward primer for second PCR
 R2	| Reverse primer for second PCR
 unique_F2	| If second PCR forward primer is unique (sequence from 3' end of primer to the 150th position cannot be found anywhere else in the genome), value = 1. Otherwise 0.
 unique_R2 | If second PCR reverse primer is unique (sequence from 3' end of primer to the 150th position cannot be found anywhere else in the genome), value = 1. Otherwise 0.


## Caveats:
  [1] Amplicons can be derived from copies of targeted genes elsewhere in the genome (not included as part of user-supplied genomic template fasta).
      Problem: filtered reads belong to two genes instead of just the targeted gene.
      Solution: Check if targeted genes have copies via BLAST. 
      Concatenate the copies to the continuous genome fasta file supplied to PRIMERg.

  [2] "Specificity" of primers is evaluated based on the user-supplied genomic template only without regard for the remaining genomic sequences. 
      Problem: Primers deemed as "specific" in PRIMERg might bind to other areas of the genome.
      Solution: To increase specificity of the deep-sequencing primers, run a first PCR to generate a 800-2000 bp amplicon that flanks the gRNA binding site. Then, run a second                 PCR with the deep-sequencing primers and use the 800-2000 bp amplicon as template. The two-step PCR first restricts the genomic template and increases hence 
                increase specificity for the second PCR that generates the NGS amplicon. Primers for both PCRs can be designed with PRIMERg (for first PCR, set  
                primer_product_size_range = [800, 2000]; for second PCR, set primer_product_size_range = [250, 280] or whichever values that is suitable for your NGS system).


## Troubleshooting 
  1. "OSError: Unrecognized base in input sequence" | Convert any non-ATCG bases (e.g. N/R) in the genomic template or gRNA fasta into A/T/C/G. 
  2. Primer designed produces multiple bands in PCR | See caveat [2] above. Also, try diluting the DNA template 100-fold before use in PCR.
  3. I want to include more parameters in primer3 | Refer to the arguments for Primer3 command line (http://primer3.org/manual.html) and implement them in the python script.
