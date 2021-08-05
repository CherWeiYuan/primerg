# primerg
Primer design for multiple gRNA targets within a continuous genome

## Background and problem definition
Members from the same gene cluster can perform redundant functions. Functional knock-out of all genes within a cluster using CRISPR-Cas9 is an approach to evaluate the role of the genes. Bioinformatic software to find the minimum number of gRNAs to target an entire gene cluster currently exist (Multitargeter: https://multicrispr.net/). Following a CRISPR-Cas9 knock-out experiment is the confirmation of mutation at Cas9 cleavage sites using NGS. However, there is currently no bioinformatic solution to resolve primer design for gene clusters by minimum gRNA targeting. Such a solution would involve designing unique primers per genome target site, if not ensuring the amplicon sequences derived from different genomic sites contain suitable polymorphism for NGS reads to be classified to each target. The primer design problem is complicated by targeted loci that are repeats or copies of each other and the exponential increase of complexity as the number of targeted loci increases.

## Solution
primerg designs primers per Cas9-targeted locus using primer3; subsequently, primerg evaluates each primer for its specificity and its amplicons for uniqueness amongst all primers of every site. The user can then pick primers that binds specifically to one locus or primers that bind to multiple loci but generates unique amplicons.

## Primerg (primer3 for gRNAs)
Finds NGS primers for multiple gRNA hits within a continuous genomic region

Goal:  Find unique NGS primer/ amplicon sequences per targeted site using primer3

Input: (Default settings) [1] gRNA list as fasta, [2] continuous genome seq as fasta and [3] directory to fasta files 

Output: Excel sheet of primers and amplicons per targeted site; unique ones are marked with "1" in adjacent column

Programme outline in default mode: 
  Process gRNAs in fasta file into a list > List all possibilities of gRNA (+ PAM/PAM-less) >
  Find position of gRNA (+ PAM/PAM-less)  > Define region for paired-end sequencing (default: 150 bp up/downstream of cleavage site) >
  Design primers with Primer3 within region > BLAST each primer against continuous genome template to define unspecific primers >
  Generate list of all possible amplicons made by the primer, which can be aligned in the future to find signature motif of desired amplicon if necessary >
  BLAST to check if desired amplicon is unique or not > Output as excel sheet

Caveats
  [1] Amplicons can be derived from copies of targeted genes elsewhere in the genome (not included as input).
      Problem: filtered reads belong to two genes instead of just the targeted gene.
      Solution: Check if targeted genes have copies via BLAST. 
      Concatenate the copies to the continuous genome fasta file supplied to Primerg.

  [2] Specificity of primers assume to be within genomic template because for NGS we generally perform two PCRs.
      First PCR (~1000-1200 bp amplicon) increases specificity of second PCR (~250-280 bp); primers of second PCR is designed in Primerg. 

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
  
## Quick start 
  1. Move script, genomic template fasta and fasta of gRNA sequences into the same folder
  2. Edit parameters in script ("Required input" and "Optional input")
  3. On the Ubuntu terminal, enter "python3 primerg.py"
  4. Retrieve the output in the same folder
 
## Troubleshooting 
  1. "OSError: Unrecognized base in input sequence" | Convert any non-ATCG bases (e.g. N/R) in genomic template or gRNA fasta into A/T/C/G. 

## Useful manual:
  Arguments for Primer3 command line: http://primer3.org/manual.html
