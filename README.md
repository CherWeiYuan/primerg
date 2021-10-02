# PRIMERg
Designs deep-sequencing/ NGS primers to assess editing outcomes from multiple gRNA targets within a continuous genome. 


## Background and problem definition
Members from the same gene cluster can perform redundant functions. Functional knock-out of multiple genes within a gene cluster using CRISPR-Cas9 is an approach to evaluate the role of such genes. The bioinformatic support for the design of gRNAs to achieve multi-gene knock-out is currently available: e.g. Multitargeter: https://multicrispr.net/. The subsequent step after the CRISPR-Cas9 knock-out experiment is the confirmation of mutation at Cas9 cleavage sites by deep-sequencing. 

Deep-sequencing to evaluate the editing status of each gene requires specific primers per genomic target site. This is a problem for gene clusters, where each locus (and possibly their intergenic region) is highly similar to each other: primers designed for one locus binds to other sites, for example in RPW8/ HR4 cluster in TueWa1-2 (10 members in gene cluster):

![Alt text](https://raw.githubusercontent.com/CherWeiYuan/primerg/main/image/multiple_binding_sites.png)

Ideally, specific primers per genomic target site should be designed. When specific primers are impossible to generate, the next acceptable primer design solution at least ensures that the amplicon sequences derived from different genomic sites contain suitable polymorphism for NGS reads to be assigned to each target. 

Manual design of unique primers is increasingly difficult as the number of targeted loci increases, especially when many gene members are (nearly) identical copies of each other. There are currently no bioinformatic solutions to resolve primer design for gene clusters CRISPR-Cas knockout experiments.


## PRIMERg as a solution
PRIMERg can automatically design specific and/ or unique deep-sequencing primers per Cas9-targeted locus in a continuous genome (e.g. sequence consisting of all members of a gene cluster with high sequence similarity). 

Two sets of primers are produced: the first PCR product is large (default: 1200-2000 bp) and flanks the region to be sequenced, thus providing the template DNA for the second PCR; the second PCR produces the amplicon to be sequenced (default: 250-280 bp). The purpose of the first PCR is to restrict the template for the second PCR to increase primer specificity. 

All primer design starts with primer3, so they are theoretically optimized. PRIMERg will filter the primers to keep the specific ones, if it is possible. The uniqueness of the second PCR primer (whether the amplicon to be sequenced has “unique” SNPs to help differentiate the desired amplicon from others produced by non-specific binding) is checked and marked in the output.


## PRIMERg algorithm
![Alt text](https://raw.githubusercontent.com/CherWeiYuan/primerg/main/image/PRIMERg_algorithm.png)


## Installations
PRIMERg runs on Linux. For Windows users, you can download Ubuntu (tested on Ubutun 20.04 but new versions should work. Installation guide: https://ubuntu.com/tutorials/ubuntu-on-windows#1-overview)

Python 3.8.5  on Linux (Confirm that Python3 is already install with ```python3 --version```. If you need to update your python version, do: ```sudo apt update```, ```sudo apt upgrade```, then ```sudo apt upgrade python3```)

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
3.	Navigate to the directory in Ubuntu, e.g. if your directory named “PRIMERg” is on desktop: 
```cd /mnt/c/Users/cherw/Desktop/PRIMERg```
5.	```python3 primerg.py```

Note: PRIMERg can output a Pandas DataFrame instead of an excel sheet. If you wish to do so, change output_excel variable in primerg.py to False.

## Interpreting output
Columns of output excel sheet:

cleavage_pos	| Expected Cas9 cleavage position on user-supplied genomic template 

gRNA_pam	| gRNA + PAM sequence

F1	| Forward primer for first PCR

R1 | Reverse primer for first PCR

F2	| Forward primer for second PCR

R2	| Reverse primer for second PCR

unique_F2	| If second PCR forward primer is unique (sequence from 3' end of primer to the 150th position cannot be found anywhere else in the genome), value = 1. Otherwise 0.

unique_R2 | If second PCR reverse primer is unique (sequence from 3' end of primer to the 150th position cannot be found anywhere else in the genome), value = 1. Otherwise 0.


## Troubleshooting 
  1. "NA" instead of primers in excel cells | This occurs when there are no suitable primers under the user-defined conditions. Try relaxing “Optional input” conditions in python script. I personally find reducing "primer_min_GC" from the recommended 40% to 30% or changing Tm helps to get primers. In certain genomic regions, GC content is low (<30%) so primers under 30 bp is impossible to design. 
  2. "OSError: Unrecognized base in input sequence" | Convert any non-ATCG bases (e.g. N/R) in the genomic template or gRNA fasta into A/T/C/G. 
  3. Primer designed produces multiple bands in PCR | See caveat [2] above. Also, try diluting the DNA template 100-fold before use in PCR.
  4. I want to include more parameters in primer3 | Refer to the arguments for Primer3 command line (http://primer3.org/manual.html) and implement them in the python script.
