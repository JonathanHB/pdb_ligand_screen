#get the proteins in which a particular residue (numbered based on the reference structure) are within 5A or a ligand
import sys
queryrseq = sys.argv[1]

refprot = "12CA"

import numpy as np
directory = "/project/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/moad_negatives/iofiles"
resisbyprot = np.load(f"{directory}/outputs-{refprot}.npy")

msagenprots = np.load(f"{directory}/msagen-{refprot}.npy")

prot_all = msagenprots[0]
pdb2msa = msagenprots[1][0]
msa2rind = msagenprots[1][1]

import mdtraj as md
xtal = md.load(f"{directory}/rcsb_pdb/{refprot}.pdb")
resseq2resi = xtal.top.atom(xtal.top.select(f"resSeq {queryrseq}")[0]).residue.index
print(resseq2resi)

query = msa2rind[refprot][pdb2msa[refprot][resseq2resi]]
print(query)

prots_with_resi_ligand = [] #proteins with a ligand near the query residue

for x, i in enumerate(resisbyprot[1]):
    #print(i)
    if query in i:
        prots_with_resi_ligand.append(resisbyprot[0][x])

print(prots_with_resi_ligand)
