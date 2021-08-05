## Primerg (primer3 for gRNAs)
#  Finds NGS primers for multiple gRNA hits within a continuous genomic region
#  Goal:  Find unique NGS primer/ amplicon sequences per targeted site using primer3
#  Input: (Default settings) [1] gRNA list as fasta, [2] continuous genome seq as fasta and 
#                            [3] directory to fasta files 

#  Output: Excel sheet of primers and amplicons per targeted site; 
#          unique ones are marked with "1" in adjacent column

#  Programme outline in default mode: 
#          Process gRNAs in fasta file into a list > List all possibilities of gRNA (+ PAM/PAM-less) >
#          Find position of gRNA (+ PAM/PAM-less)  > Define region for paired-end sequencing (default: 150 bp up/downstream of cleavage site) >
#          Design primers with Primer3 within region > BLAST each primer against continuous genome template to define unspecific primers >
#	   Generate list of all possible amplicons made by the primer, which can be aligned in the future to find signature motif of desired amplicon if necessary >
#          BLAST to check if desired amplicon is unique or not > Output as excel sheet

#  Caveats: [1] Amplicons can be derived from copies of targeted genes elsewhere in the genome (not included as input).
#               Problem: filtered reads belong to two genes instead of just the targeted gene.
#               Solution: Check if targeted genes have copies via BLAST. 
#               Concatenate the copies to the continuous genome fasta file supplied to Primerg.

#           [2] Specificity of primers assume to be within genomic template because for NGS we generally perform two PCRs.
#               First PCR (~1000-1200 bp amplicon) increases specificity of second PCR (~250-280 bp); primers of second PCR is designed in Primerg. 

## Requirement(s):
#  Ubuntu
#  Python3 on Linux (primer3 is not available on Windows)
#  BLAST on Linux (sudo apt-get install ncbi-blast+)
#  primer3 on Linux (sudo apt-get install primer3)
#  Also see 'Import packages' below for required python packages, primarily: primer3, biopython, pandas

## Useful manuals:
#  Arguments for Primer3 command line: http://primer3.org/manual.html

## Troubleshooting
#  1. "OSError: Unrecognized base in input sequence" | Are there bases other than ATCG in the genomic template or gRNA fasta files?

## Variables for users

#  Required input
directory = '/mnt/c/Users/cherw/Desktop/primerg'
gRNA_fasta = "TueWa1-2_gRNA1022,1023,1027.fasta"       #list of gRNAs in 5'-3' direction
genomic_DNA_fasta = "TueWa1-2_RPW8_HR_genomic.fasta"

#  Optional input

#  Primer3
primer_opt_size = 20                  #recommend 20 bp
primer_min_size = 18                  #recommend 18 bp
primer_max_size = 30                  #recommend 30 bp
primer_opt_tm = 65                    #recommend 65
primer_min_tm = 50                    #recommend 50
primer_max_tm = 70                    #recommend 70
primer_min_GC = 40                    #recommend 40%        
primer_max_GC = 60                    #recommend 60%
primer_product_size_range = [250, 280] #recommend [250, 280] for NGS sequencing amplicons
primer_dna_conc = 50                  #in nm
primer_dNTP_conc = 0.5                #in mM
primer_pair_max_diff_TM = 5           #recommend 5
primer_salt_divalent = 1.5            #mM of divalent cations e.g. Mg2+
primer_num_return = 50                #maximum number of primers to return per primer3 run 

#  CRISPR-Cas parameters
pam = "NGG"                           #e.g. "NGG", let pam be '' if pamless 
gRNA_len = 20
cleavage_pos = 3                      #nth base from 3' end of gRNA where Cas9 cleaves genomic template,
                                      #default = 3 where NNNNNNNNNNNNNNNNN|NNN
#  Pseudo-Primer-BLAST
min_primer_len = 10                    # The minimum length of alignment between primer and template to consider a primer as unspecific; default 10 bp is minimum for annealing prior to extension
total_mismatch = 6                     # The minimum number of mismatches between primer and template to consider a primer as unspecific; default is 6
three_prime_match = 3                  # At least X matches within the last 5 bp at the 3' end; default is 3
valid_amplicon_size = 500              # Size of amplicon produced by unspecific primers to be considered valid amplicon, default is 500 ***************** 500
NGS_amplicon_size = 750                # Length of amplicon of EACH primer in paired-end sequencing; default is set at 150 for iSeq ********** 

#  Primerg
output_excel = "output.xlsx"          #if an excel name is specified (e.g. filename.xlsx), output will be an excel sheet. If False, output will be a pandas dataframe.

#######################################################################################################################
## Source code
#  Do not edit unnecessarily

## Import packages
import os
import primer3
from Bio import SeqIO
from Bio.Seq import Seq
import re
import pandas as pd
import sys

from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

## Reading parameters
os.chdir(directory)
template = str(SeqIO.read(genomic_DNA_fasta, "fasta").seq).upper()    #template read as string


## Read gRNA fasta file into a list
def gRNA_listing(gRNA_fasta):
    gRNA_list = []
    try:
        gRNA_fasta_list = SeqIO.parse(gRNA_fasta, "fasta")
        for record in gRNA_fasta_list:
            gRNA_list += [str(record.seq).upper()]      
    except ValueError:                #to catch fasta files with single gRNA (use SeqIO.read instead of SeqIO.parse)
        gRNA_fasta_list = SeqIO.read(gRNA_fasta, "fasta")
        gRNA_list += [str(gRNA_fasta_list.seq.upper())] 
    return gRNA_list


## List all possibilities of gRNA + pam
def gRNA_pam_listing(gRNA_list, pam, template):
    gRNA_pam_list_F = []
    gRNA_pam_list_R = []
    gRNA_pam_list = []    
    
    #for CRISPR with pam
    if pam != '':
        #find all possibilities of pam    
        pam_list = [pam.replace('N', 'A'), pam.replace('N', 'T'), pam.replace('N', 'C'), pam.replace('N', 'G')]
        
        #generate list of gRNA + pam (forward)
        temp_list_F = [a + s for s in pam_list for a in gRNA_list] 
        
        #generate list of gRNA + pam (reverse complement)
        temp_list_R = [str(Seq(s).reverse_complement()).upper() + str(Seq(a).reverse_complement()).upper() for s in pam_list for a in gRNA_list]
        
        #create new list of gRNA + pam that is found in template
        for gRNA_pam in temp_list_F:
            if gRNA_pam in template:
                gRNA_pam_list_F += [ gRNA_pam ]
        for gRNA_pam in temp_list_R:
            if gRNA_pam in template:
                #converts gRNA + pam into reverse complement before adding to list
                gRNA_pam_list_R += [str(Seq(gRNA_pam).reverse_complement())]
        #output list will consist of [ [forward gRNA + pam], [reverse complement gRNA + pam]]
        gRNA_pam_list = [gRNA_pam_list_F, gRNA_pam_list_R]
        return gRNA_pam_list
    
    #for pamless CRISPR
    elif pam == '': 
        #generate list of gRNA (forward)
        temp_list_F = [a for a in gRNA_list] 
        
        #generate list of gRNA (reverse complement)
        temp_list_R = [str(Seq(a).reverse_complement()).upper() for a in gRNA_list]
        
        #create new list of gRNA that is found in template
        for gRNA_pam in temp_list_F:
            if gRNA_pam in template:
                gRNA_pam_list_F += [ gRNA_pam ]
        for gRNA_pam in temp_list_R:
            if gRNA_pam in template:
                #converts gRNA + pam into reverse complement before adding to list
                gRNA_pam_list_R += [str(Seq(gRNA_pam).reverse_complement())]
        #output list will consist of [ [forward gRNA + pam], [reverse complement gRNA + pam]]
        gRNA_pam_list = [gRNA_pam_list_F, gRNA_pam_list_R]
        return gRNA_pam_list        

## Primer design

#  Create error class
class SequenceError(Exception):
    pass

#  Find gRNA in genomic DNA: returns list of positions of each gRNAs; only one position per gRNA will be returned
def gRNA_finder(gRNA_pam_list, template):
    temp = []
    lst_F = []
    lst_RC = []
    lst_grand = []
    
    #note: data structure of gRNA_pam_list as [ [forward gRNA + pam], [reverse complement gRNA + pam]]
    
    #find position of forward gRNA cleavage site
    for gRNA_pam in gRNA_pam_list[0]:
        #find overlapping matches for forward gRNA AND keep list of Cas9 cleavage site position
        temp  += [ (m.start() + len(gRNA_pam) - len(pam) - cleavage_pos) for m in re.finditer('(?=' + gRNA_pam + ')', template) ] 
        lst_F += [(gRNA_pam, temp)]
        temp = []
        
    #find position of reverse complement gRNA cleavage site
    for gRNA_pam in gRNA_pam_list[1]:
        #find overlapping matches with reverse comp gRNA AND keep list Cas9 cleavage site position
        temp += [ (m.start() + len(pam) + cleavage_pos) for m in re.finditer('(?=' + str(Seq(gRNA_pam).reverse_complement()) + ')', template) ] 
        lst_RC += [(gRNA_pam, temp)]
        temp = []
        
    #structure of lst_grand is e.g. [('GGCGTTGACAGATGAGGGGCAGG', [403, 501]), ('AATGCTGGATTTTCTGCCTGTGG', [643])]
    lst_grand += lst_F + lst_RC
    return lst_grand

#  Design primers      
                
def primer_design(gRNA_pos, template):
    #initiate dataframe
    df = pd.DataFrame(columns = ['gRNA_index', 'region_index', 'combined_index', 'gRNA_pam', 'region', 'primer_F', 'primer_R',
                                                      'specific_F_primer', 'specific_R_primer', 'F_amplicons', 'R_amplicons', 'unique_F_amplicon', 'unique_R_amplicon'])                             

    #find targeted regions where primers are designed
    half_len = round(primer_product_size_range[1]/2) #half the length of amplicon
    gRNA_index = 0
    region_index = 1
    
    for unit in gRNA_pos:
        gRNA_seq = unit[0]
        pos_lst = unit[1]
        gRNA_index += 1
        #iterate through the position of cleavage sites
        for pos in range(len(pos_lst)):
            
            #reset regions for this cleavage site 
            regions = []
            
            #initiate template seqs per cleavage site
            regions += [template[pos_lst[pos] - half_len: pos_lst[pos] + half_len]]

            #design primers for current position
            for i in regions:
                try:
                    design_primer = primer3.bindings.designPrimers(
                        {
                            'SEQUENCE_TEMPLATE': i
                        },
                        {
                            'PRIMER_OPT_SIZE': primer_opt_size,
                            'PRIMER_MIN_SIZE': primer_min_size,
                            'PRIMER_MAX_SIZE': primer_max_size,
                            'PRIMER_OPT_TM': primer_opt_tm,
                            'PRIMER_MIN_TM': primer_min_tm,
                            'PRIMER_MAX_TM': primer_max_tm,
                            'PRIMER_MIN_GC': primer_min_GC,
                            'PRIMER_MAX_GC': primer_max_GC,
            
                            'PRIMER_DNA_CONC': primer_dna_conc,
                            'PRIMER_PRODUCT_SIZE_RANGE': primer_product_size_range,
                            'PRIMER_DNTP_CONC': primer_dNTP_conc,
                            'PRIMER_PAIR_MAX_DIFF_TM': primer_pair_max_diff_TM,
                            'PRIMER_SALT_DIVALENT': primer_salt_divalent,
                            'PRIMER_NUM_RETURN': primer_num_return,
                            
                        })
                except(OSError):
                    print("Error: PCR product size parameter demand primers to be designed beyond genomic template. Include longer genomic template sequence or 'decrease primer_product_size_range'.")
                    sys.exit()
                
                #retrieve all primers from design 
                for j in range(len(design_primer)):
                    try:
                        left = None
                        right = None
                        left = design_primer['PRIMER_LEFT_' + str(j) + '_SEQUENCE']
                        right = design_primer['PRIMER_RIGHT_' + str(j) + '_SEQUENCE']
                        primer_R_pos = template.find(str(Seq(right).reverse_complement())) + len(right)
                        df2 = pd.DataFrame([ [gRNA_index, region_index, f"{gRNA_index}_{region_index}", gRNA_seq, i, left, right, '', '', '', '', '', ''],],
                                           columns = ['gRNA_index', 'region_index', 'combined_index', 'gRNA_pam', 'region', 'primer_F', 'primer_R',
                                                      'specific_F_primer', 'specific_R_primer', 'F_amplicons', 'R_amplicons', 'unique_F_amplicon', 'unique_R_amplicon'])
                        df = df.append(df2) 
                    except KeyError:
                        break
                        
                    
                region_index += 1  
    return df

## BLAST function to count the number of non-unique primers the selected primer produces
## Input: list of amplicons produced by primer (amplicon library: library_F or library_R)

def amplicon_library_blast(library, genomic_DNA_fasta):
    cline = NcbiblastnCommandline(query = genomic_DNA_fasta, 
                              subject = library, 
                              out = "amplicon-blast.xml", 
                              outfmt = 5)                               #format of output as XML means outfmt = 5
    cline()
    result_handle = open("amplicon-blast.xml")
    blast_records = NCBIXML.read(result_handle)                         #use .read, not .parse, because there is only one query seq
    
    # Calculate amplicon length
    amplicon_length = None
    libraryFasta = SeqIO.parse(library, "fasta")
    for i in libraryFasta:
        amplicon_length = len(i.seq)
        break

    # Start counting number of non-unique amplicons generated for this primer
    count = 0    
    
    for alignment in blast_records.alignments:
        for hsp in alignment.hsps:
            if hsp.align_length == amplicon_length:                     #check if entire amplicon sequence is aligned; 
                                                                        #if not, then it indicates differences between aligned sequence
                if round((hsp.match.count('|')/hsp.align_length) * 100) == 100: 
                    count += 1
    
    return count



def pseudo_primer_blast(df, genomic_DNA_fasta):
    
    # Iterate through every primer pair in df
    len_df = len(df)
    for index, row in df.iterrows():
        print(f"Checking primer specificity at combined_index {df.loc[index,'combined_index']}, row {index+1}/{len_df+1}")
        primer_F = row['primer_F']
        primer_R = row['primer_R']
        primer_FR = primer_F + 'N'*20 + primer_R 
        
        # Write primer_F + NNNNNNNNNNNNNNNNNNNN + primer_R into fasta file for BLAST input
        SeqIO.write(SeqRecord(Seq(primer_FR), id = "primer_pair", annotations={"molecule_type": "DNA"}, description = ""), 
            "primer_pair.fasta", "fasta")
        
        ## BLAST to find binding sites for both primers in genomic template
        ## Note: The 20 Ns between F and R will allow each primer to be treated as separate entities in BLAST
        ##       i.e. We will find binding sites for F and R individually
        
        # BLAST set-up
        cline = NcbiblastnCommandline(query = "primer_pair.fasta", 
                                      subject = genomic_DNA_fasta, 
                                      evalue= 30000,                            #high e-value to allow more chance hits so we can evaluate all possible binding sites
                                      word_size = 7,                            #small seed is required as primer seq is short; allows more possible hits
                                      task = "blastn-short",                    #BLASTN program optimized for sequences shorter than 50 bases
                                      dust = "'no'",                            #No masking of low complexity sequence
                                      soft_masking = "'false'",                 #No masking of low complexity sequence
                                      max_target_seqs = 50000,                  #Allow more hits to be shown
                                      out = "primer-blast.xml",                 #store BLAST result in XML file
                                      outfmt = 5 )                              #format of output as XML means outfmt = 5
        
        # Run BLAST
        cline()
        
        # Read BLAST results
        result_handle = open("primer-blast.xml")
        blast_records = NCBIXML.read(result_handle)
        
        ## Iterate through BLAST results and mark out unspecific primers
        ## Valid unspecific primers are survivors of the following filters (parameters adjustable as variables above):
        ## [1] Length of alignment must exceed 10
        ## [2] Total number of mismatches does NOT exceed 6
        ## [3] Last 5 bp of 3' end of aligned sequences has at least 3 matches
        
        ## Unspecific primers will have their 3' positions stored in data structure
        ## i.e. query_start (if reverse primer), query_end (if forward primer); query refers to genomic PCR template
        
        # Data structure:
        # Generate database for positions of primer binding; record forward and reverse primer binding separately
        # i.e. [ [FORWARD PRIMER], [REVERSE PRIMER] ]
        
        pos_binding = [ [], [] ]
        for alignment in blast_records.alignments:
            for hsp in alignment.hsps:
                if len(hsp.query) >= min_primer_len:
                    if hsp.match.count(' ') <= total_mismatch:
                        if hsp.strand[1] == 'Plus':
                            if hsp.match[-5:].count('|') >= three_prime_match:
                                pos_binding[0] += [hsp.sbjct_end]     # store position at 3' end of forward primer
                        elif hsp.strand[1] == 'Minus':
                            if hsp.match[0:4].count('|') >= three_prime_match:
                                pos_binding[1] += [hsp.sbjct_start]   # store position at 3' end of reverse primer
                    
        # Remove positions if they cannot form valid amplicons
        # i.e. Predicted amplicons are too long to be amplified    
        pos_valid = [ [], [] ]
        for posF in pos_binding[0]:
            for posR in pos_binding[1]:
                if posR - posF > 0 and posR - posF >= valid_amplicon_size: #if posR is at 5' end of posF on plus strand of genomic template, no amplicon can be formed;
                                                                           #thereby the filter posR - posF > 0 where posR must be bigger than posF
                    pos_valid[0] += [posF]
                    pos_valid[1] += [posR]
        
        pos_valid = [ set(pos_valid[0]), set(pos_valid[1]) ]
        
        # Mark out unspecific primers
        if len(pos_valid[0]) > 1:
            df.loc[index, 'specific_F_primer'] = 0
        elif len(pos_valid[0]) == 1:
            df.loc[index, 'specific_F_primer'] = 1
        
        if len(pos_valid[1]) > 1:
            df.loc[index, 'specific_R_primer'] = 0
        elif len(pos_valid[1]) == 1:
            df.loc[index, 'specific_R_primer'] = 1
            
        ## Generate amplicons of F and R primers separately
        ## Generate amplicon sequence from each position
        ## Note: amplicon here does not include primer sequence because we are checking 
        ##       if amplified sequences can be differentiated from each other
        ##       (amplified sequences can differ but their primers has fixed sequence)
        
        amplicon_F_list = []
        amplicon_R_list = []
        
        for posF in pos_valid[0]:
            amplicon_F_list += [template[posF: posF + NGS_amplicon_size - len(primer_F)]]     #template = genomic template in string
        
        for posR in pos_valid[1]:
            amplicon_R_list += [template[posR - NGS_amplicon_size + len(primer_F): posR - 1]] #template = genomic template in string
        
        # Update df with amplicons
        # df.at with the column dtype=object will allow the df cells to accept a list
        df['F_amplicons'] = df['F_amplicons'].astype('object')
        df['R_amplicons'] = df['R_amplicons'].astype('object')
        df.at[index, 'F_amplicons'] = [amplicon_F_list]        #nested list is accepted into cell easily
        df.at[index, 'R_amplicons'] = [amplicon_R_list]        #nested list is accepted into cell easily
        
        ## BLAST template against F & R amplicon library separately
        ## If there are >1 full length, 100% identical matches, then amplicon is not unique (cannot be differentiated from each other) 
        ## The primer (F and/or R) thus does not generate unique amplicon 
        ## If either primer F or R generated unique amplicons, the primer pair can still be used
        
        # Retrieve expected amplicon sequence
        region = row['region']
        desired_amplicon_F = region[region.find(primer_F) + len(primer_F): region.find(primer_F) + 150]
        library_F = SeqRecord(Seq(desired_amplicon_F), annotations={"molecule_type": "DNA"}, description = "", id = 'desired F amplicon')
        SeqIO.write(library_F, "library_F.fasta", "fasta")
        
        primer_R_rc = str(Seq(primer_R).reverse_complement())
        desired_amplicon_R = region[region.find(primer_R_rc)-150+len(primer_R):region.find(primer_R_rc)]        
        library_R = SeqRecord(Seq(desired_amplicon_R), annotations={"molecule_type": "DNA"}, description = "", id = 'desired R amplicon')
        SeqIO.write(library_R, "library_R.fasta", "fasta")

        # BLAST against library_F
        same_F_amplicon_count = amplicon_library_blast("library_F.fasta", genomic_DNA_fasta)
        if same_F_amplicon_count > 1:
            df.loc[index, 'unique_F_amplicon'] = 0
        elif same_F_amplicon_count == 1:
            df.loc[index, 'unique_F_amplicon'] = 1
        
        # BLAST against library_R
        same_R_amplicon_count = amplicon_library_blast("library_R.fasta", genomic_DNA_fasta)
        if same_R_amplicon_count > 1:
            df.loc[index, 'unique_R_amplicon'] = 0
        elif same_R_amplicon_count == 1:
            df.loc[index, 'unique_R_amplicon'] = 1        
        
        
'''        
        ## BLAST template against F & R amplicon library separately
        ## If there are >1 full length, 100% identical matches, then amplicon is not unique (cannot be differentiated from each other) 
        ## The primer (F and/or R) thus does not generate unique amplicon 
        ## If either primer F or R generated unique amplicons, the primer pair can still be used
     
        # Generate amplicon library
        library_F = [SeqRecord(Seq(seq), annotations={"molecule_type": "DNA"}, description = "", id = str(index)) for index,seq in enumerate(amplicon_F_list) ]
        SeqIO.write(library_F, "library_F.fasta", "fasta")
        
        library_R = [SeqRecord(Seq(seq), annotations={"molecule_type": "DNA"}, description = "", id = str(index)) for index,seq in enumerate(amplicon_R_list) ]
        SeqIO.write(library_R, "library_R.fasta", "fasta")
        
        # BLAST against library_F
        same_F_amplicon_count = amplicon_library_blast("library_F.fasta", genomic_DNA_fasta)
        if same_F_amplicon_count > 1:
            df.loc[index, 'unique_F_amplicon'] = 0
        elif same_F_amplicon_count == 1:
            df.loc[index, 'unique_F_amplicon'] = 1
        
        # BLAST against library_R
        same_R_amplicon_count = amplicon_library_blast("library_R.fasta", genomic_DNA_fasta)
        if same_R_amplicon_count > 1:
            df.loc[index, 'unique_R_amplicon'] = 0
        elif same_R_amplicon_count == 1:
            df.loc[index, 'unique_R_amplicon'] = 1
'''          

def export(df):
    #export dataframe as excel sheet
    if output_excel == False:
        return df
    else:
        df.to_excel(output_excel, index = False)
    pass

def execute_all():
    gRNA_list = gRNA_listing(gRNA_fasta)
    print("gRNA list obtained.")
    
    gRNA_pam_list = gRNA_pam_listing(gRNA_list, pam, template)
    if pam != '':
        print("gRNA + PAM list obtained.")
    elif pam == '':
        print("PAM-less gRNA mode active.")
        
    gRNA_pos = gRNA_finder(gRNA_pam_list, template)
    if pam != '':
        print("All gRNA + PAM positions are found in template.")
    elif pam == '':
        print("All PAM-less gRNA positions are found in template.")
    
    print("")
    print("Designing primers...")
    df = primer_design(gRNA_pos, template).reset_index(drop = True)
    
    print("")
    print("Checking primer specificity and amplicon uniqueness")
    pseudo_primer_blast(df, genomic_DNA_fasta)
    
    print("")
    export(df)
    print("Execution completed.")
    print(f"Output is placed at {directory} as {output_excel}.")

##  Run programme
execute_all()

#Feedback
#2. Filter to keep only rows with 1(s)
#3. Primer BLAST!
#4. template.count for amplicon should consider overlapping but not for primer
#5. find signature motif for read sorting: only for unique amplicon = 1 but unique primer = 0






