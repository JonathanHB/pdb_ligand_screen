import numpy as np
import itertools

import os
import csv

import mdtraj as md

directory = "/project/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/moad_negatives/iofiles"
serial = 2
align = False

blasthits = np.load(f"{directory}/blast_output_v{serial}.npy")

coord_threshold = 0.5 #try 3.5 or 4, although this risks including cryptic pocket residues as negatives
#.6 gets residue 243 on 1jwp; neither .6 nor .7 gets resi 249

testindex = 0

testprots = blasthits[testindex]

#-----------------------------------MUSCLE multiple sequence alignment code from Artur-----------------------------------#
import os
import mdtraj as md
import numpy as np
import skbio
from skbio import Protein, TabularMSA

def generate_msa(trajectories, protein_names):
    '''
        trajectories : mdtraj trajectories for which you wish to do the sequence alignment.
            Code assumes that you will supply only one protein chain

    '''
    with open('alignment.fasta', 'w') as f:
        for p, pdb in zip(protein_names, trajectories):
            prot = Protein(
                ''.join(r.code for r in pdb.top.residues if r.code),
                metadata={'id': p})
            prot.write(f, format='fasta')

    musclepath = "/project/bowmore/ameller/"
    #print(f'{musclepath}muscle -in alignment.fasta -out alignment-aligned.fasta')
    os.system(f'{musclepath}muscle -align alignment.fasta -output alignment-aligned.fasta')

    msa = TabularMSA.read('alignment-aligned.fasta', constructor=Protein)
    msa.reassign_index(minter='id')
    pdb2msa = {
        protein: np.where(~msa[msa.index == protein][0].gaps())[0]
        for protein in msa.index
    }
    msa2rind = {protein: {r: i for i, r in enumerate(*np.where(~msa[msa.index == protein][0].gaps()))}
                for protein in msa.index}

    return pdb2msa, msa2rind

#------------------------------------------------------------------------------------------------------------------------#

#get structure without downloading files repeatedly, will check resolution
def getstruct(struct, extant):

    #download the pdb structure if necessary
    if struct not in extant:
        try: #download and load .pdb file
            os.system(f"wget https://files.rcsb.org/download/{struct}.pdb -O {directory}/rcsb_pdb/{struct}.pdb")

        except IndexError:
            print("error downloading file; it may be unavailable in .pdb format")

            #large structures cannot be downloaded as pdb files and are too big to be desirable for our present analysis
            os.system(f"rm {directory}/rcsb_pdb/{struct}.pdb")
            #remove the empty file generated in the failed download attempt so that subsequent runs don't try to load it in mdtraj and subsequently crash
            return False

def getcoordresseqs(holo_xtal, ligand_resns, maxind):

    nonstandard_resis = ['B','J','O','X','Z']

    ligand_select_str = "resname '"+"' or resname '".join(ligand_resns)+"'"

    #Distance computations to identify ligand-coordinating holo atoms
    ligands = holo_xtal.top.select(f"{ligand_select_str} and not element H")
    protein = holo_xtal.top.select("protein and not element H")

    #print(holo_xtal.top.atom(protein[0]))

    lig_prot_pair_iis = np.array(list(itertools.product(ligands, protein)))
    prot_lig_dists = md.compute_distances(holo_xtal,lig_prot_pair_iis, periodic = False).flatten()

    lig_coord_iis = np.where(prot_lig_dists<coord_threshold)[0] #within 5 angstroms of ligand
    prot_iis = np.unique([lig_prot_pair_iis[i][-1] for i in lig_coord_iis])

    lining_resseqs = np.unique([holo_xtal.top.atom(k).residue.resSeq for k in prot_iis]) #get residue containing the atom
    lining_res_prot = [holo_xtal.top.atom(k).residue for k in prot_iis if holo_xtal.top.atom(k).residue.code not in nonstandard_resis]
    lining_resis = np.unique([k.index for k in lining_res_prot])

    all_resseqs = list([i.code for i in holo_xtal.top.residues if i.is_protein])

    return [lining_resseqs, lining_resis]


################################################################################
#----------------------------------Main loop-----------------------------------#
################################################################################

queryrseq = [] #93, 95, 97, 226

nonstandard_prots = ["6QEB", "6XWZ","5WG7"]#,"6DWU","6SWV"]
#6qeb has multiple states,
#6xwz has a disulfide bond error (but no such bond),
#5wg7 has ligand hydrogens with the wrong element annotation due to some kind of mdtraj bug
#7bs1 has indel mutations that complicate residue assignment <--remove from this list once MSA is working

existing_xids_moad = [i[0:4] for i in os.listdir(f"{directory}/moad_xml/")] #get a list of already-downloaded ligand lists to avoid downloading extra copies
existing_pdbids = [i[0:4] for i in os.listdir(f"{directory}/rcsb_pdb/")] #get a list of already-downloaded structures to avoid downloading extra copies

lining_resis_all = []
lining_resis_byprot = [] #rseqs by source protein for manual inspection

#------------------------------------------------test muscle msa-----------------------------------------
#b,j,o,x,z

if align == True:

    print("aligning sequences")

    xtal_all = []
    for x, prot in enumerate(testprots[0]):
        #print(f"{x}: {prot}")

        if prot in nonstandard_prots:
            continue

        getstruct(prot, existing_pdbids)
        xtal = md.load(f"{directory}/rcsb_pdb/{prot}.pdb")
        xtal_A = xtal.atom_slice(xtal.top.select("chainid 0"))
        xtal_all.append(xtal_A)

    msagen = generate_msa(xtal_all, testprots[0])

    np.save(f"{directory}/msagen-{testprots[0][0]}", msagen)

else:
    msagen = np.load(f"{directory}/msagen.npy")
pdb2msa = msagen[0]
msa2rind = msagen[1]

maxind = pdb2msa[testprots[0][0]][-1]

#print(msa2rind["1AXB"])
#print(pdb2msa["1ERQ"])

#------------------------------------------------test muscle msa end-------------------------------------

for x, prot in enumerate(testprots[0]):
    print(f"{x}: {prot}")
    prot = "6SAY"

    if prot in nonstandard_prots:
        continue

    getstruct(prot, existing_pdbids)

    xtal = md.load(f"{directory}/rcsb_pdb/{prot}.pdb")

    ligandsolvent = [[j for j in i.residues] for i in xtal.top.chains] #may not be robust for multiple chains

    #note that aa_resns includes d-amino acid residues such as d-phenylalanine (DPN)
    aa_resns = ["ASN", "ASP", "GLN", "GLU", "THR", "SER", "LYS", "ARG", "HIS", "PRO", "GLY", "CYS", "MET", "ALA", "VAL", "LEU", "ILE", "PHE", "TYR", "TRP", "DHI", "DPN"]
    solvent = ["HOH","DOD","NH3","NH4","CO3","NO3","PO4","SO4","LI","BE","K","NA","MG","CA","F","CL","BR","I","ZN","MN","CU","FE","O2","HDZ","CO2","CO2","XE","UNX"]

    ligands = [str(i) for i in ligandsolvent[-1]]
    ligand_resseqs = [i.resSeq for i in ligandsolvent[-1]]
    ligand_resns = [i[0:-len(str(ligand_resseqs[x]))] for x, i in enumerate(ligands)]
    ligand_resns_nonaq = [i for i in ligand_resns if i not in solvent and i not in aa_resns]

    if ligand_resns_nonaq != []:
        lining_rseqs_resis = getcoordresseqs(xtal, ligand_resns_nonaq, maxind)
        lining_rseqs = lining_rseqs_resis[0]
        lining_resis = lining_rseqs_resis[1]

        #find which proteins have the residues in queryrseq exposed
        for qsq in queryrseq:
            if qsq in lining_rseqs:
                print(f"{prot}---{qsq}------------------------------------------------")

        #print(lining_resis)
        print(pdb2msa[prot])
        #print(msa2rind[testprots[0][0]])
        lining_aligned = [msa2rind[testprots[0][0]][pdb2msa[prot][i]] for i in lining_resis if pdb2msa[prot][i] in msa2rind[testprots[0][0]]]

        lining_resis_all += lining_aligned
        lining_resis_byprot.append(lining_aligned)

    else:
        lining_resis_byprot.append([])

lining_resis_all = np.unique(lining_resis_all)

#-------------------------------outputs-----------------------------------------

#print(lining_rseqs_byprot)

print(testprots[0][0])

xtal = md.load(f"{directory}/rcsb_pdb/{testprots[0][0]}.pdb")
all_resseqs = list([i.resSeq for i in xtal.top.residues])

print("+".join([str(all_resseqs[i]) for i in lining_resis_all]))

################################################################################
# things to try:

#try tem and see if we get the horn site and don't get the residues known not to label
#try the sars-2 macrodomain
#use muscle multiple sequence alignment

################################################################################
#---------------------Greg says to try simulating the output-------------------#
################################################################################
