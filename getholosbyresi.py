import sys
import pickle
import numpy as np
import mdtraj as md

################################################################################
#return PDB IDs of the proteins in which a particular residue
#(numbered based on the centroid/reference structure) are within 5A or a ligand
################################################################################

#-------------------------------arguments---------------------------------------
#take the first argument as the centroid protein of interest
refprot = sys.argv[1]
#take the second argument as the query pdb residue number
queryrseq = sys.argv[2]

#-------------------------------static parameters-------------------------------
#directory containing protein lists
directory = "/project/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/moad_negatives/iofiles"

#--------------------------load processed blast outputs-------------------------
#list of proteins in the cluster of the argument centroid
aligned_prot = np.load(f"{refprot}-aligned-proteins.npy")

#ligand-lining residues by protein
resisbyprot = np.load(f"{refprot}-lining-resis-byprot.npy")

#dictionaries needed to convert between resseqs and residue indices (via msa indices)
with open(f"{refprot}-pdb2msa", 'rb') as pickle_file:
    pdb2msa = pickle.load(pickle_file)

with open(f"{refprot}-msa2rind", 'rb') as pickle_file:
    msa2rind = pickle.load(pickle_file)

#centroid crystal structure
xtal = md.load(f"{directory}/rcsb_pdb/{refprot}.pdb")

#-------------------------------------------------------------------------------
#convert the argument pdb residue to a centroid structure residue index
resseq2resi = xtal.top.atom(xtal.top.select(f"resSeq {queryrseq}")[0]).residue.index
query = msa2rind[refprot][pdb2msa[refprot][resseq2resi]]

#-------------------------------------------------------------------------------
#look for proteins with the query residue index and print their PDB IDs for inspection

#proteins with a ligand near the query residue
prots_with_resi_ligand = []

for x, i in enumerate(resisbyprot):
    #print(i)
    if query in i:
        prots_with_resi_ligand.append(aligned_prot[x])

print(prots_with_resi_ligand)
