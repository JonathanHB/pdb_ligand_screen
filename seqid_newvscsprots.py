import numpy as np
import os
import sys
import csv
import glob

import Bio
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist
from operator import itemgetter

#Jonathan Borowsky
#1/24/22
#This progarm compares the sequence identity of every protein used in
#the GVP cryptic pocket detection project to every protein in cryptosite.

################################################################################
#get FASTA amino acid sequences for all proteins
################################################################################

#lists of proteins
#-------------------------------------------------------------------------------

#cryptosite structures
cryptosite_multiseq = [["1HAG","E","1GHY","H","121"],["1JBU","H","1WUN","H","P5B"],["1KZ7","D","1GRN","A","AF3"],["1XCG","B","1OW3","B","GDP"],\
                       ["1XMG","B","1XVC","A","5BR"],["1Z92","A","1PY2","A","FRH"],["2AIR","H","1ZA1","D","CTP"],["2AKA","A","1YV3","A","BIT"],\
                       ["2BF3","A","3DHH","E","BML"],["2GFC","A","2JDS","A","L20"],["3CHE","A","2IUZ","B","D1H"],["3FDL","A","2YXJ","A","N3C"],\
                       ["3H5R","A","3H9J","D","APC"],["3MN9","A","3EKS","A","CY9"],["3PUW","E","1FQC","A","GLO"],["4HB2","C","4HAT","C","LMB"]]
cryptosite_multiseq_apo = [e[0].lower() for e in cryptosite_multiseq]

cryptosite_1seq = [["1RTC","A","1BR6","A","PT1"],["1CLL","A","1CTR","A","TFP"],["1B6B","A","1KUV","A","CA5"],["1JWP","A","1PZO","A","CBT"],["1BP5","A","1RYO","A","OXL"],["2F6V","A","1T49","A","892"],["1TQO","A","1TR5","A","THP"],["3CJ0","A","2BRL","A","POO"],["1SWX","A","2EUM","A","LAT"],["2CM2","A","2H4K","A","509"],["1H09","A","2IXU","A","MU2"],["2IYT","A","2IYQ","A","ADP/SKM"],["2ZB1","A","2NPQ","A","BOG"],["1KS9","A","2OFP","A","PAF"],["2AX9","A","2PIQ","A","RB1"],["2QFO","B","2WI7","A","2KL"],["2YQC","A","2YQS","A","UD1"],["1R1W","A","3F82","A","353"],["3NNU","A","3HL7","A","I46"],["1FXX","A","3HL8","A","BBP"],["3DXN","A","3HZT","A","J60"],["1HKA","A","3IP0","A","HHS"],["3KQA","B","3LTH","A","UD1"],["1G4E","B","1G67","B","POP/TZP"],["1EXM","A","1HA3","B","MAU"],["3HQD","A","1Q0B","B","NAT"],["3CJ0","A","3FQK","B","79Z"],["2BRK","A","2GIR","B","NN3"],["1RHB","A","2W5K","B","NDP"],["1SU4","A","3FGO","B","ACP"],["1NI6","D","3HOK","B","Q80"],["1G24","D","1GZF","C","NIR"],["1NEP","A","2HKA","C","C3S"],["3F74","C","3BQM","C","BQM"],["1W50","A","3IXJ","C","586"],["1FVR","A","2OO8","X","RAJ"],["1BSQ","A","1GX8","A","RTL"],["1ALB","A","1LIC","A","HDS"],["1PZT","A","1PZY","D","UDP"],["1ADE","A","1CIB","A","IMP"],["1HOO","B","1CIB","A","HDA"],["1EX6","A","1GKY","A","5GP"],["1RDW","X","1J6Z","A","RHO"],["1BNC","B","2V5A","A","LZL"],["2BLS","B","3GQZ","A","GF7"],["4AKE","B","1ANK","B","ANP"],["3GXD","B","2WCG","A","MT5"],["1PKL","B","3HQP","P","ATP/FDP/OXL"],["1NUW","A","1EYJ","B","AMP"],["1MY1","C","1FTL","A","DNQ"],["1E2X","A","1H9G","A","MYR"],["1ALV","A","1NX3","A","ISA"],["1RRG","A","1S9D","A","AFB"],["2BU8","A","2BU2","A","TF1"],["1UK2","A","2GZ7","A","D3F"],["2OHG","A","2OHV","A","NHL"],["2Q8F","A","2Q8H","A","TF4"],["2WGB","A","2V57","A","PRL"],["3BL9","B","3BL7","A","DD1"],["2WGQ","B","1D6Y","B","HY1"],["1IMF","A","1IMB","B","LIP"],["1K5H","C","2EGH","B","FOM"],["1A8I","A","2IEG","B","FRY"],["1QLW","B","2WKW","B","W22"],["1MY0","B","1N0T","D","AT1"],["2QLR","C","3DC1","A","AKG"],["2GPO","A","1S9Q","B","CHD"],["1OK8","A","1OKE","B","BOG"],["3L7U","C","2HVD","C","ADP"],["1DUB","D","1EY3","F","DAK"],["1K3F","B","1U1D","F","181"],["3PEO","G","2BYS","J","LOB"],["1ECJ","D","1ECC","B","PCP"],["1FA9","A","1L5S","B","URC"],["2H4E","B","3CFN","B","2AN"],["1ZAH","B","2OT1","D","N3P"],["2CGA","B","1AFQ","C","0FG"],["3B7D","E","2AL4","F","CX6"],]
cryptosite_1seq_apo = [e[0].lower() for e in cryptosite_1seq]

proteins_cs_all = cryptosite_1seq + cryptosite_multiseq
proteins_cs = cryptosite_1seq_apo + cryptosite_multiseq_apo

#new proteins
test_val_set = ['1ezm','1j8f','1kmo','1kx9','1s2o','1tvq','1urp','1y1a','2cey','2fjy','2hq8','2lao','2oy4','2w9t','2zku','3fvj','3kje','3nx1','3p53','3ppn','3qxw','3rwv','3ugk','4i92','4ic4','4p0i','4r72','4v38','4w51','5g1m','5h9a','5nia','5nzm','5uxa','5za4','6e5d','6hb0','6rvm','6ypk']
bricks = ['1igd','2fd7','1ofv','1qys','5bvl','4hjk','4tql','1amm','2alp']
moad_negatives = [e.lower() for e in ['1YDA','2HIZ','6CD4','1OVE','1OIT','1V2N','2ZJW','1E5O','3B0W','1UIB']]

#proteins used in five fold cross validation
filenames_ffcv = glob.glob("/project/bowmanlab/ameller/gvp/msft-share/*.pdb")
training = [fname.split("/")[-1].split(".")[0].lower() for fname in filenames_ffcv]
training_noncs = [e for e in training if (e not in proteins_cs)]
#this set of proteins overlaps with the proteins above
#Some of these proteins are SARS proteins with structures from swiss model
#homology modelling rather than PDB

#lists of proteins for the main loop:
proteins_gvp = test_val_set + bricks + moad_negatives + training_noncs
#proteins_cs = cryptosite_1seq_apo + cryptosite_multiseq_apo

proteins_all = proteins_cs + proteins_gvp

#get FASTA sequences for proteins without pdb ids (in practice sars-cov-2 swissmodel structures)
#-------------------------------------------------------------------------------

fasta_dir = "/project/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/fasta_files"

existing_ids = [i.split(".")[0] for i in os.listdir(fasta_dir)]

#get FASTA sequences from pdb files
for fname, pname in zip(filenames_ffcv, training_noncs):
    if pname not in existing_ids:
        with open(fname, 'r') as pdb_file:
            #running SeqIO.parse repeatedly causes it to crash with an empty file error
            SeqIO.write(SeqIO.parse(pdb_file, 'pdb-atom'), f"{fasta_dir}/{pname}.fasta", "fasta")
            existing_ids.append(i)

#download cryptosite/brick/new set FASTA files
proteins_rcsb = test_val_set + bricks + moad_negatives + proteins_cs
for i in proteins_rcsb:
#new_bricks_cryptosite+cryptosite_added_prots+moad_negatives_simulated:
    if i not in existing_ids:
        os.system(f"wget https://www.rcsb.org/fasta/entry/{i}/display -O {fasta_dir}/{i}.fasta")
        existing_ids.append(i)

#-------------------------------------------------------------------------------
#print utility

def getclass(prot1):
    if prot1 in test_val_set:
        p1_class = " (validation/test)"
    elif prot1 in bricks:
        p1_class = " (brick)"
    elif prot1 in moad_negatives:
        p1_class = " (moad negative)"
    elif prot1 in cryptosite_1seq_apo:
        p1_class = " (cryptosite 1seq)"
    elif prot1 in cryptosite_multiseq_apo:
        p1_class = " (cryptosite multiseq)"
    elif prot1 in training_noncs:
        p1_class = " (5fcv/training)"
    else:
        print("error: undefined set")
        p1_class = " error: undefined set "

    return p1_class

################################################################################
#check sequence identity of each pair of proteins
################################################################################

id_cutoff = 0.4 #maximum sequence identity

#eliminate low-identity hits which commonly occur by chance when aligning a
#short sequence to a much longer one
#set to 1 by default, set to 0 to disable this adjustment
#len_scale = 0 #1

zmax = len(proteins_cs)*len(proteins_gvp)
z = 0

#collect a list of matches above 40% identity
seqid_out = []

print(f"aligning {len(proteins_cs)} cryptosite proteins with {len(proteins_gvp)} gvp proteins; {int(zmax)} protein pairs")

for x, prot1 in enumerate(proteins_cs):

    print(f"{x}: {prot1}, {'{:.4f}'.format(z/zmax)} complete")

    try:
        fasta_seq1 = SeqIO.read(open(f"{fasta_dir}/{prot1}.fasta"), 'fasta').seq
    except: #parse multi-sequence FASTA files

        fasta_strs = [i for i in open(f"{fasta_dir}/{prot1}.fasta")]

        for strx, e in enumerate(fasta_strs):
            if strx%2 == 0:
                if e.split("|")[1][-1] == "]":
                    chain_auth = e.split("|")[1].split("[")[1].split(" ")[1].split("]")[0]
                    if chain_auth == proteins_cs_all[x][1]:
                        fasta_seq1 = fasta_strs[strx+1].strip()
                else:
                    chain_auth = e.split("|")[1].split(" ")[1]
                    if chain_auth == proteins_cs_all[x][1]:
                        fasta_seq1 = fasta_strs[strx+1].strip()

    for y, prot2 in enumerate(proteins_gvp):
        if prot1 == prot2:
            print(f"duplicate: {x}:{prot1}, {y}:{prot2}")

        fasta_seq2 = SeqIO.read(open(f"{fasta_dir}/{prot2}.fasta"), 'fasta').seq

        #matrix = matlist.blosum62
        alignment = pairwise2.align.globalds(fasta_seq1, fasta_seq2, matlist.blosum62, -11, -1) #parameters matching blastp
        #setting d uses the matrix, s uses the gap creation and extension penalties
        #note that format_alignment() can be used to display alignments nicely as long as they aren't wider than the terminal

        #print(alignment)

        minlen = min(len(fasta_seq1), len(fasta_seq2))
        #use the shorter of the two sequences as the denominator to determine the percent identity
        #as a fraction of the maximum achieveable percent identity given the sequence lengths
        #this approach is more conservative

        maxlen = max(len(fasta_seq1), len(fasta_seq2))

        identical_resis = 0

        for k in range(len(alignment[0].seqA)):
            if alignment[0].seqA[k] == alignment[0].seqB[k]:
                identical_resis += 1

        idres2 = sum([1*(alignment[0].seqA[k] == alignment[0].seqB[k]) for k in range(len(alignment[0].seqA))])
        if idres2 != identical_resis:
            sys.exit(0)

        seqid = identical_resis/len(alignment[0].seqA) #fraction of identical residues in the alignment

        #only works if the sequence space is only sparsely populated, which it appears to be for the set of filtered pairs
        #the second condition eliminates low-identity hits which commonly occur by chance when aligning a short sequence to a much longer one
        if seqid > id_cutoff:
            p1_class = getclass(prot1)
            p2_class = getclass(prot2)

            seqid_out.append([seqid, prot1, p1_class, prot2, p2_class])

            print(f"{prot1}{p1_class} and {prot2}{p2_class} have {seqid*100}% identity")

        z+=1 ##counter

print(sorted(seqid_out, key=itemgetter(0)))

