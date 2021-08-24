# primerg
Designs deep-sequencing/ NGS primers to assess editing outcomes from multiple gRNA targets within a continuous genome. 

## Background and problem definition
Members from the same gene cluster can perform redundant functions. Functional knock-out of genes within a cluster using CRISPR-Cas9 is an approach to evaluate the role of such genes. The bioinformatic support for the design of gRNAs to achieve multi-gene knock-out is currently available: e.g. Multitargeter: https://multicrispr.net/. The subsequent step after the CRISPR-Cas9 knock-out experiment is the confirmation of mutation at Cas9 cleavage sites by deep-sequencing. 

Deep-sequencing to evaluate editing status of each gene requires unique amplicons per genomic target site. Hence, ideally, unique primers per genomic target site should be designed. When unique primers is impossible to generate, the next acceptable primer design solution at least ensures that the amplicon sequences derived from different genomic sites contain suitable polymorphism for NGS reads to be assigned to each target. 

Manual design of unique primers are decreasingly efficient as the number of targeted loci increases, especially when many gene members are (nearly) identical copies of each other. There are currently no bioinformatic solution to resolve primer design for gene clusters CRISPR-Cas knockout experiments.

## Solution
We present the primerg software. primerg designs primers per Cas9-targeted locus using primer3; subsequently, primerg evaluates each primer for its specificity (binds to only one site in the user-supplied genomic template), and its amplicons for uniqueness amongst all primers of every site (if classified as "unique", it means no two amplicons from different sites share an identical sequence). The user can then pick primers that binds specifically to one locus and/ or primers that bind to multiple loci but generates unique amplicons.

## Outline of primerg 
Finds NGS primers for multiple gRNA targets within a continuous genomic region

Goal:  Find specific primer/ unique amplicon sequences per targeted site using primer3

Input: (Default settings) [1] gRNA list as fasta, [2] continuous genome seq as fasta and [3] directory to fasta files (to be specified within python script)

Output: Excel sheet of all primer3-generated primers and their amplicons per targeted site; specific primers and primers generating unique amplicons are marked with "1" in the respective columns

![Alt text](algorithm_map.png?raw=true)

Programme outline in default mode:
 
  Process gRNAs in fasta file into a list > List all possibilities of gRNA (+ PAM/PAM-less) >
  Find position of gRNA (+ PAM/PAM-less)  > Define region for paired-end sequencing (default: 150 bp up/downstream of cleavage site) >
  Design primers with Primer3 within region > BLAST each primer against continuous genome template to define unspecific primers >
  Generate list of all possible amplicons made by the primer, which can be aligned in the future to find signature motif of desired amplicon if necessary >
  BLAST to check if desired amplicon is unique or not > Output as excel sheet

Caveats:

  [1] Amplicons can be derived from copies of targeted genes elsewhere in the genome (not included as part of user-supplied genomic template fasta).
      Problem: filtered reads belong to two genes instead of just the targeted gene.
      Solution: Check if targeted genes have copies via BLAST. 
      Concatenate the copies to the continuous genome fasta file supplied to Primerg.

  [2] "Specificity" of primers is evaluated based on the user-supplied genomic template only without regard for the remaining genomic sequences. 
      Problem: Primers deemed as "specific" in primerg might bind to other areas of the genome.
      Solution: To increase specificity of the deep-sequencing primers, run a first PCR to generate a 800-2000 bp amplicon that flanks the gRNA binding site. Then, run a second                 PCR with the deep-sequencing primers and use the 800-2000 bp amplicon as template. The two-step PCR first restricts the genomic template and increases hence 
                increase specificity for the second PCR that generates the NGS amplicon. Primers for both PCRs can be designed with primerg (for first PCR, set  
                primer_product_size_range = [800, 2000]; for second PCR, set primer_product_size_range = [250, 280] or whichever values that is suitable for your NGS system).

# Guide

## System requirements 
  Ubuntu 20.04, with the following packages installed:
  Python 3.8.5 on Linux 
  BLAST 2.9.0+ on Linux    (sudo apt-get install ncbi-blast+)
  primer3 2.4.0-2 on Linux (sudo apt-get install primer3)

## Python packages
  Install the following packages with pip3 on Ubuntu terminal
  Biopython 1.7.8
  primer3 0.6.1
  pandas 0.25.3
  
## How to run primerg?
  1. Move script, genomic template fasta and fasta of gRNA sequences into the same folder
  2. Edit parameters in script ("Required input" and "Optional input")
  3. On the Ubuntu terminal, enter "python3 primerg.py"
  4. Retrieve the output in the same folder
 
## Troubleshooting 
  1. "OSError: Unrecognized base in input sequence" | Convert any non-ATCG bases (e.g. N/R) in genomic template or gRNA fasta into A/T/C/G. 
  2. Primer designed produces multiple bands in PCR | See caveat [2] above. Also, try diluting the DNA template by 100-fold prior to use in PCR.
  3. I want to include more parameters in primer3 | Refer to the arguments for Primer3 command line (http://primer3.org/manual.html) and implement them in the python script. 
