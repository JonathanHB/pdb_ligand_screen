import numpy as np
import os

#find the centroid ids of the usearch clusters from which a set of proteins came

directory = "/project/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/moad_negatives/iofiles"
#serial number of input files
serial_in = 3

#get a list of files containing the usable proteins for each cluster/centroid
output_fns = os.listdir(f"{directory}/processed-blast-v{serial_in}")
processed_pdb = [i[0:4] for i in output_fns if i[4:] == f"-v{serial_in}-canonical-proteins.npy"]

#loop through each cluster
for pdbid in processed_pdb:

    protein_canonical = np.load(f"{directory}/processed-blast-v{serial_in}/{pdbid}-v{serial_in}-canonical-proteins.npy")
    prot_holo_monomer = np.load(f"{directory}/processed-blast-v{serial_in}/{pdbid}-v{serial_in}-holo-monomer.npy")
    prot_apo_monomer = np.load(f"{directory}/processed-blast-v{serial_in}/{pdbid}-v{serial_in}-apo-monomer.npy")

    print(f"{pdbid}: {(len(prot_holo_monomer)+len(prot_apo_monomer))/len(protein_canonical)}")

    #check if each query protein is in that cluster
    #for sprot in simulated_prots:
    #    if sprot in aligned_prots:
    #        print(f"{sprot} is in the cluster of centroid {afn[0:4]}")
