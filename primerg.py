## PRIMERg
# Deep-sequencing primer design (with primerg) for multiple targets on a continuous genome

## Source code
#  Do not edit unnecessarily

#  Required input
directory = None
gRNA_fasta = None                 #list of gRNAs in 5'-3' direction (PAM is not required here)
genomic_DNA_fasta = None

#  Optional input

#  CRISPR-Cas parameters
pam = None                          #e.g. "NGG", let pam be '' if pamless, follows IUPAC DNA (e.g. R = A/G)
gRNA_len = None
cleavage_pos = None                    #nth base from 3' end of gRNA where Cas9 cleaves genomic template,
                                      #default = 3 where NNNNNNNNNNNNNNNNN|NNN
#  Primer3
primer_opt_size = None                 #recommend 20 bp
primer_min_size = None                 #recommend 18 bp
primer_max_size = None                  #recommend 30 bp
primer_opt_tm = None                    #recommend 65
primer_min_tm = None                    #recommend 60
primer_max_tm = None                   #recommend 70
primer_min_GC = None                    #recommend 40%; 30% for leniency ###check###
primer_max_GC = None                   #recommend 60%
primer_dna_conc = None                 #in nm
primer_dNTP_conc = None                #in mM
primer_pair_max_diff_TM = None           #recommend 5
primer_salt_divalent = None            #mM of divalent cations e.g. Mg2+
primer_num_return = None                #maximum number of primers to return per primer3 run 

#  Pseudo-Primer-BLAST
min_primer_len = None                   # The minimum length of alignment between primer and template to consider a primer as unspecific; default 10 bp is minimum for annealing prior to extension
total_mismatch = None                    # The minimum number of mismatches between primer and template to consider a primer as unspecific; default is 6
three_prime_match = None                  # At least X matches within the last 5 bp at the 3' end; default is 3
valid_amplicon_size = None             # Size of amplicon produced by unspecific primers to be considered valid amplicon, default is 500 ***************** 500
NGS_amplicon_size = None               # Length of amplicon of EACH primer in paired-end sequencing; default is set at 150 for iSeq ********** 

#  Primerg
output_excel = None          #if an excel name is specified (e.g. filename.xlsx), output will be an excel sheet. If False, output will be a pandas dataframe.

#  New inputs
primary_PCR_amplicon_size = None # First PCR product size (size of amplicon as template for second PCR)
secondary_PCR_amplicon_size = None# Second PCR product size (size of amplicon for sequencing)
num_primer_pair_per_cleavage_pos = None     # Generate x number of primer pairs per cleavage site
increment = None                         # Increment of allowed primary amplicon size when no primer can be generated
max_primary_amplicon_size = None      # Maximum primary amplicon size

## Import packages
from primerg_parameters import *
from os import chdir
from sys import exit
from re import findall, finditer
from pandas import DataFrame
import primer3
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

## Reading parameters
chdir(directory)
genomic_template = str(SeqIO.read(genomic_DNA_fasta, "fasta").seq).upper()    #template read as string                      


## Read gRNA fasta file into a list
def gRNA_listing(gRNA_fasta):
    # Time: O(n) where n is length of gRNA list
    # Space: O(n) where n is length of gRNA list
    
    gRNA_list = []
    
    try:
        gRNA_fasta_list = SeqIO.parse(gRNA_fasta, "fasta")
        for record in gRNA_fasta_list:
            gRNA_list += [str(record.seq).upper()]      
            
    except ValueError:    #to catch fasta files with single gRNA (use SeqIO.read instead of SeqIO.parse)
        gRNA_fasta_list = SeqIO.read(gRNA_fasta, "fasta")
        gRNA_list += [str(gRNA_fasta_list.seq.upper())] 
        
    return gRNA_list
              

## Find all gRNA + pam sequences
def gRNA_pam_listing(gRNA_list, pam, genomic_template):
    # Time and space depends on Regex's re.findall
    
    IUPAC_dict = {'A': 'A', 'T': 'T', 'C': 'C', 'G': 'G',
                  'N': '[ATCG]', 'R': '[AG]', 'Y': '[CT]',
                  'S': '[GC]', 'W':'[AT]', 'K':'[GT]',
                  'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
                  'H': '[ACT]', 'V': '[ACG]'}
    
    reg_pam = ''
    F = []
    R = []
    for letter in pam:
        reg_pam += IUPAC_dict[letter]
    
    gRNA_pam_list = [x + reg_pam for x in gRNA_list]
    
    for gRNA_pam in gRNA_pam_list:
        F += findall(gRNA_pam, genomic_template)
        R += findall(gRNA_pam, str(Seq(genomic_template).reverse_complement()))
        
    return [list(set(F)), list(set(R))]

## Primer design

#  Create error class
class SequenceError(Exception):
    pass

#  Find gRNA in genomic DNA: returns list of positions of each gRNA's cleavage site; only one position per gRNA will be returned
def gRNA_finder(gRNA_pam_list, genomic_template):
    # Time: Depends on Regex's re.finditer
    # Space: O(n) where n is len(gRNA_pam_list)
    
    temp = []
    lst_grand = []
    unpacked_lst = []
    
    #note: data structure of gRNA_pam_list is [ [forward gRNA + pam], [reverse complement gRNA + pam]]
    
    #find position of forward gRNA cleavage site
    for gRNA_pam in gRNA_pam_list[0]:
        #find overlapping matches for forward gRNA AND keep list of Cas9 cleavage site position
        temp  += [ (m.start() + len(gRNA_pam) - len(pam) - cleavage_pos) for m in finditer('(?=' + gRNA_pam + ')', genomic_template) ] 
        lst_grand += [(gRNA_pam, temp)]
        temp = []
        
    #find position of reverse complement gRNA cleavage site
    for gRNA_pam in gRNA_pam_list[1]:
        #find overlapping matches with reverse comp gRNA AND keep list Cas9 cleavage site position
        temp += [ (m.start() + len(pam) + cleavage_pos) for m in finditer('(?=' + str(Seq(gRNA_pam).reverse_complement()) + ')', genomic_template) ] 
        lst_grand += [(gRNA_pam, temp)]
        temp = []
    
    #structure of lst_grand is e.g. [('GGCGTTGACAGATGAGGGGCAGG', [403, 501]), ('AATGCTGGATTTTCTGCCTGTGG', [643])]
    for i in lst_grand:
        for j in i[1]:
            unpacked_lst += [tuple((i[0], j))]
    
    #reset list to save space
    lst_grand = []
    
    return unpacked_lst

# Find tuple of primers F and R for ONE cleavage site
def primer_design(pos, template, amplicon_size):

    #find targeted regions where primers will be designed
    half_len = round(amplicon_size[1]/2) #half the length of amplicon
    region = template[pos - half_len: pos + half_len]        

    #design primers for current position
    try:
        design_primer = primer3.bindings.designPrimers(
            {
                'SEQUENCE_TEMPLATE': region
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
                'PRIMER_PRODUCT_SIZE_RANGE': amplicon_size,
                'PRIMER_DNTP_CONC': primer_dNTP_conc,
                'PRIMER_PAIR_MAX_DIFF_TM': primer_pair_max_diff_TM,
                'PRIMER_SALT_DIVALENT': primer_salt_divalent,
                'PRIMER_NUM_RETURN': primer_num_return,
                
            })
        
    except(OSError):
        print("Error: PCR product size parameter demand primers to be designed beyond genomic template. Include longer genomic template sequence or 'decrease primer_product_size_range'.")
        exit()
    
    #retrieve all primers from design 
    F_lst = []
    R_lst = []
    
    for j in range(len(design_primer)):
        try:
            F_lst += [design_primer['PRIMER_LEFT_' + str(j) + '_SEQUENCE']]
            R_lst += [design_primer['PRIMER_RIGHT_' + str(j) + '_SEQUENCE']]
        
        except KeyError:
            break
    
    if len(F_lst) < primer_num_return:
        print(f"User-defined number of forward primers to return = {primer_num_return} is not met. Explanation from Primer3:")
        print(design_primer['PRIMER_LEFT_EXPLAIN'])

    if len(F_lst) < primer_num_return:
        print(f"User-defined number of reverse primers to return = {primer_num_return} is not met. Explanation from Primer3:")
        print(design_primer['PRIMER_RIGHT_EXPLAIN'])
    
    F_lst = tuple(set(F_lst))
    R_lst = tuple(set(R_lst))
    
    # Return 
    return (pos, F_lst, R_lst)

# Filter primer_list (primer_list = ((list of F primers), (list of R primers)) by pseudo_primer_blast 
# Returns primer_list with only specific primers left
def pseudo_primer_blast(primer_list, template):
    survivors = [[], [], []] # primer_list[0] is pos from primer_design(), in turn from gRNA_finder()
    for k in range(1,3):
        for i in primer_list[k]:     
            bind_count = 0
            # Write primer in fasta format for BLAST
            SeqIO.write(SeqRecord(Seq(i), 
                                  id = "primer", 
                                  annotations={"molecule_type": "DNA"}, 
                                  description = ""), 
                        "primerg_specificity_check_primer.fasta", "fasta")

            SeqIO.write(SeqRecord(Seq(template), 
                                  id = "template", 
                                  annotations={"molecule_type": "DNA"}, 
                                  description = ""), 
                        "primerg_specificity_check_template.fasta", "fasta")
        
            # BLAST set-up
            cline = NcbiblastnCommandline(query = "primerg_specificity_check_primer.fasta", 
                                          subject = "primerg_specificity_check_template.fasta", 
                                          evalue= 30000,                            #high e-value to allow more chance hits so we can evaluate all possible binding sites
                                          word_size = 7,                            #small seed is required as primer seq is short; allows more possible hits
                                          task = "blastn-short",                    #BLASTN program optimized for sequences shorter than 50 bases
                                          dust = "'no'",                            #No masking of low complexity sequence
                                          soft_masking = "'false'",                 #No masking of low complexity sequence
                                          strand = "plus",                          #Only plus strand is used as query
                                          max_target_seqs = 50000,                  #Allow more hits to be shown
                                          out = "primer-blast.xml",                 #store BLAST result in XML file
                                          outfmt = 5 )                              #format of output as XML means outfmt = 5
            
            # Run BLAST
            cline()
            
            # Read BLAST results
            result_handle = open("primer-blast.xml")
            blast_records = NCBIXML.read(result_handle)
            
            ## Iterate through BLAST results and add specific primers to returned list
            ## Pseudo-primer-blast checks if primer binds to site if they meet the requirements below (parameters adjustable as variables above):
            ## [1] Length of alignment must exceed 10
            ## [2] Total number of mismatches does NOT exceed 6
            ## [3] Last 5 bp of 3' end of aligned sequences has at least 3 matches
            
            # bind_count adds one when primer binds to a site
            # if bind_count exceeds one, primer is unspecific
            
            for alignment in blast_records.alignments:
                for hsp in alignment.hsps:
                    if len(hsp.query) >= min_primer_len:
                        if hsp.match.count(' ') <= total_mismatch:
                            if hsp.strand[1] == 'Plus':
                                if hsp.match[-5:].count('|') >= three_prime_match:
                                    bind_count += 1

                            elif hsp.strand[1] == 'Minus':
                                if hsp.match[0:4].count('|') >= three_prime_match:
                                    bind_count += 1
                                    
            if bind_count > 1 or bind_count == 0:
                pass
            
            elif bind_count == 1:
                if k == 1:
                    survivors[1] += [i]
                elif k == 2:
                    survivors[2] += [i]

    survivors = [primer_list[0], list(set(survivors[1])), list(set(survivors[2])) ]
    return survivors


# Generate sets of primers given primer_list ((list of F primers), (list of R primers))
# N is the number of primer sets desired
def generate_primer_pairs(primer_list, N, template, amplicon_size):
    
    original_N = N # keeps track of starting N value
    n = N          # keeps track of alternative list length
    
    F_lst = []     # Generate F and R primer lists where positions in each list is corresponding e.g. lst[0][0] corresponds to list[1][0]  
    R_lst = []
    alt_F_lst = [] # alternate lists keep specific + unspecific primer pairs
    alt_R_lst = []
    
    local_template = template[int(primer_list[0]-(0.5*amplicon_size[1])): int(primer_list[0]+(0.5*amplicon_size[1])) ]
    
    for i in range(1,3):
        if i == 1: # Means we are looking at forward (left) primer
            if N <= 0:
                break
            
            for primer in primer_list[i]:
                if N <= 0:
                    break
                
                try:

                    design_primer = primer3.bindings.designPrimers(
                    {
                        'SEQUENCE_TEMPLATE': local_template
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
                        'PRIMER_PRODUCT_SIZE_RANGE': amplicon_size,
                        'PRIMER_DNTP_CONC': primer_dNTP_conc,
                        'PRIMER_PAIR_MAX_DIFF_TM': primer_pair_max_diff_TM,
                        'PRIMER_SALT_DIVALENT': primer_salt_divalent,
                        'PRIMER_NUM_RETURN': primer_num_return,
                        
                        'SEQUENCE_PRIMER': primer
                        
                    })
                
                except OSError:
                    print("Error: PCR product size parameter demand primers to be designed beyond genomic template. Please provide longer genomic template sequence or decrease 'primer_product_size_range'.")
                    exit()
                
                except IndexError:
                    print("Error: local_template cannot be generated because cleavage site position +/- 0.5*desired_amplicon_size is beyond user-supplied genomic template. Please provide longer genomic template sequence.")
                    exit()                    
                    
                except: # List is empty
                    pass
        
                
                for j in range(len(design_primer)):
                    if N <= 0:
                        break
                    try:
                        if design_primer['PRIMER_RIGHT_' + str(j) + '_SEQUENCE'] in primer_list[2]:
                            F_lst += [primer]
                            R_lst += [design_primer['PRIMER_RIGHT_' + str(j) + '_SEQUENCE']]
                            N -= 1
                            break # ensures we only use the same primer only one time
                        else:
                            if n > 0:
                                alt_F_lst += [primer]
                                alt_R_lst += [design_primer['PRIMER_RIGHT_' + str(j) + '_SEQUENCE']] 
                                n -= 1
                                break
           
                    except KeyError:
                        break

        elif i == 2: # Means we are looking at reverse (right) primer
            if N <= 0:
                break
            
            for primer in primer_list[i]:
                if N <= 0:
                    break
                
                try:
                    design_primer = primer3.bindings.designPrimers(
                    {
                        'SEQUENCE_TEMPLATE': local_template
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
                        'PRIMER_PRODUCT_SIZE_RANGE': amplicon_size,
                        'PRIMER_DNTP_CONC': primer_dNTP_conc,
                        'PRIMER_PAIR_MAX_DIFF_TM': primer_pair_max_diff_TM,
                        'PRIMER_SALT_DIVALENT': primer_salt_divalent,
                        'PRIMER_NUM_RETURN': primer_num_return,
                        
                        'SEQUENCE_PRIMER_REVCOMP': primer
                        
                    })
                
                except(OSError):
                    print("Error: PCR product size parameter demand primers to be designed beyond genomic template. Include longer genomic template sequence or 'decrease primer_product_size_range'.")
                    exit()
                    
                except: # List is empty
                    pass
        
                
                for j in range(len(design_primer)):
                    if N <= 0:
                        break
                    try:
                        if design_primer['PRIMER_LEFT_' + str(j) + '_SEQUENCE'] in primer_list[1]:
                            F_lst += [design_primer['PRIMER_LEFT_' + str(j) + '_SEQUENCE']]
                            R_lst += [primer]
                            N -= 1
                            break # ensures we only use the same primer only one time
                        else:
                            if n > 0:
                                alt_F_lst += [design_primer['PRIMER_LEFT_' + str(j) + '_SEQUENCE']]
                                alt_R_lst += [primer]
                                n -= 1
                                break # ensures we only use the same primer only one time
           
                    except KeyError:
                        break                
                
                
    if N > 0: # when specific F cannot correspond to specific R
        try:
            F_lst += alt_F_lst[0: original_N - (original_N-N)]
            R_lst += alt_R_lst[0: original_N - (original_N-N)]
        except: # in case alt_F or alt_R is empty, which is unlikely
            pass
     
    # Returns [pos, [forward primers], [reverse primers]], where position of pairing primers are corresponding
    # e.g. list[0][0] corresponds to list[1][0]
    return [primer_list[0], F_lst, R_lst]

# Generates PCR amplicon made by F and R primer (without primer sequence within the returned amplicon)
def generate_amplicon(pos, F, R, template, amplicon_size):
    local_template = template[int(pos-(0.5*amplicon_size[1])): int(pos+(0.5*amplicon_size[1])) ]
    pos_start = str(local_template).find(str(F)) + len(F)
    pos_end = str(local_template).find(str(Seq(R).reverse_complement())) 
    return local_template[pos_start: pos_end]

# Test
#Input:
#x = "aatttcagaaatctgtttttttttctctctctatcttcttccgcgtaatgattgaaaacctttttttcaatctatgtgtaattttattaaaataagaattaataataatctttttatttcgtttttttgagatctatgtttcatgaattttttgcataa"
#F = "tcttcttccgcgtaa"
#R = "caaaaaattcatga"
#generate_amplicon(F, R, x)

#Output: 'tgattgaaaacctttttttcaatctatgtgtaattttattaaaataagaattaataataatctttttatttcgtttttttgagatctatgtt'


### Execution

## Set up dataframe
df = DataFrame(columns = ['cleavage_pos', 'gRNA_pam', 
                             'F1', 'R1', 
                             'F2', 'R2', 
                             'unique_F2', 'unique_R2'])

## Ensure excel file is not opened
try: 
    df.to_excel(output_excel, index = False)\

except PermissionError:
    print(f"Please close {output_excel} and try again")
    exit()

# Initialize empty variables
amplicon = None
filtered_primer_list = None

# Obtain list of gRNAs and cleavage site positions
gRNA_list = gRNA_listing(gRNA_fasta)
gRNA_pam_list = gRNA_pam_listing(gRNA_list, pam, genomic_template)
gRNA_pos = gRNA_finder(gRNA_pam_list, genomic_template)

index = 0
start_index = 0
original_value = primary_PCR_amplicon_size[1]

for pos in gRNA_pos:    
    start_index = index
    # Update meta-data
    df2 = DataFrame([ [pos[1], pos[0], '', '', '', '', '', ''] ,], 
                       columns = ['cleavage_pos', 'gRNA_pam', 'F1', 'R1', 'F2', 'R2', 'unique_F2', 'unique_R2'])
    
    for i in range(num_primer_pair_per_cleavage_pos):
        df = df.append(df2).reset_index(drop = True)
    
    print(f"Processing gRNA + PAM {df.loc[index, 'gRNA_pam']} at position {df.loc[index, 'cleavage_pos']}")
    
    # Update primary PCR primer (first PCR)    
    print("Generating primary primer list")
    primer_list = primer_design(pos[1], genomic_template, primary_PCR_amplicon_size)
    
    print("Filtering primary primer list")
    filtered_primer_list = pseudo_primer_blast(primer_list, genomic_template)
    
    while filtered_primer_list[1] == [] and filtered_primer_list[2] == []:
        
        # Break if max amplicon size is reached
        if primary_PCR_amplicon_size[1] > max_primary_amplicon_size: 
           break 
       
        primary_PCR_amplicon_size[1] += increment
        # Generate new primary primer list with increased allowed amplicon size
        print("Generating primary primer list")
        primer_list = primer_design(pos[1], genomic_template, primary_PCR_amplicon_size)
    
        print("Filtering primary primer list")
        filtered_primer_list = pseudo_primer_blast(primer_list, genomic_template)
        
        print(f"No primary primer found. Increased allowed primary amplicon size to {primary_PCR_amplicon_size}")
    
    primary_PCR_amplicon_size[1] = original_value
    
    print("Generating primary primer pairs")
    primer_pairs = generate_primer_pairs(filtered_primer_list, num_primer_pair_per_cleavage_pos, genomic_template, primary_PCR_amplicon_size)
    ###check###
    print(primer_pairs)
    for i in range(num_primer_pair_per_cleavage_pos):
        #df.loc[index, 'cleavage_pos'] = pos[1] 
        #df.loc[index, 'gRNA_pam'] = df.loc[index, 'gRNA_pam']
        
        # Add primer to df except when there are no specific primary PCR primers
        try:
            df.loc[index, 'F1'] = primer_pairs[1][i]
        except:
            df.loc[index, 'F1'] = "NA"
        try:
            df.loc[index, 'R1'] = primer_pairs[2][i]
        except:
            df.loc[index, 'R1'] = "NA"
        index += 1
    
    # Update secondary PCR primer (second PCR)
    print("Generating secondary primer list")
    max_amplicon2_len = int(secondary_PCR_amplicon_size[1]/2)
    local_genomic_template = genomic_template[int(pos[1]-max_amplicon2_len) : int(pos[1] + max_amplicon2_len)]
    
    primer_list = primer_design(int(0.5*len(local_genomic_template)), local_genomic_template, secondary_PCR_amplicon_size)    
    
    print("Screening secondary primer list against genomic template to get globally specific primers")
    # Filter by genomic template to get globally specific primers
    filtered_primer_list = pseudo_primer_blast(primer_list, genomic_template)
    
    print("Generating globally specific secondary primer pairs")
    primer_pairs = generate_primer_pairs(filtered_primer_list, num_primer_pair_per_cleavage_pos, local_genomic_template, secondary_PCR_amplicon_size)
    
    if primer_pairs[1] == [] or primer_pairs[2] == []:        
        print("Globally specific primers are non-existent")
        print("Screening secondary primer list against first PCR amplicon to get locally specific primers")
        # Filter by genomic template to get globally specific primers
        filtered_primer_list = pseudo_primer_blast(primer_list, local_genomic_template)
        
        print("Generating locally specific secondary primer pairs")
        primer_pairs = generate_primer_pairs(filtered_primer_list, num_primer_pair_per_cleavage_pos, local_genomic_template, secondary_PCR_amplicon_size)

    for i in range(num_primer_pair_per_cleavage_pos):
        try:
            df.loc[start_index, 'F2'] = primer_pairs[1][i]
            df.loc[start_index, 'R2'] = primer_pairs[2][i]
            amplicon = generate_amplicon(primer_pairs[0], df.loc[start_index, 'F2'], df.loc[start_index, 'R2'], local_genomic_template, secondary_PCR_amplicon_size)
            ampliconF = amplicon[len(df.loc[start_index, 'F2']):secondary_PCR_amplicon_size[1]-len(df.loc[start_index, 'F2'])]
            ampliconR = amplicon[-(secondary_PCR_amplicon_size[1]) : -len(df.loc[start_index, 'R2']):]
            
            countF = 0
            countR = 0
            
            for i in finditer('(?=' + ampliconF + ')', genomic_template):
                countF += 1
            for i in finditer('(?=' + ampliconR + ')', genomic_template):
                countR += 1
                
            if countF == 1:
                df.loc[start_index, 'unique_F2'] = 1
            elif countF > 1:
                df.loc[start_index, 'unique_F2'] = 0
            else:
                df.loc[start_index, 'unique_F2'] = "Error: Amplicon F count in genome is negative"
            
            if countR == 1:
                df.loc[start_index, 'unique_R2'] = 1
            elif countR > 1:
                df.loc[start_index, 'unique_R2'] = 0
            else:
                df.loc[start_index, 'unique_R2'] = "Error: Amplicon R count in genome is negative"
                
            start_index += 1
            
        except IndexError: # when no secondary primer is specific to even second PCR amplicon
            df.loc[start_index, 'F2'] = "NA"
            df.loc[start_index, 'R2'] = "NA"
            start_index += 1            

    print("Dataframe updated:")    
    print(df)
    
# Export to excel
df = df.sort_values(by=['cleavage_pos']) # sort by ascending order
try: 
    if output_excel:
        df.to_excel(output_excel, index = False)
        print(f"Complete. Excel sheet exported as {directory}/{output_excel}")
    elif output_excel == False:
        print("Complete. Pandas dataframe is returned. For output as excel sheet, please give a name to variable 'output_excel' in python script.")
        df
    
except PermissionError:
    print(f"Please close {output_excel} and try again")












