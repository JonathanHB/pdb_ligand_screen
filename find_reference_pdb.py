import numpy as np
import os

#find the centroid ids of the usearch clusters from which a set of proteins came

directory = "/project/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/moad_negatives/iofiles"
#input proteins associated with an unknown centroid
simulated_prots = ["7LFO","5UOJ","3GZ0","5E4P","4APE","1BTP","2ZHV","1HCL"]

#get a list of files containing the usable proteins for each cluster/centroid
output_fns = os.listdir(f"{directory}/processed-blast-v3")
aligned_fns = [i for i in output_fns if i[4:] == "-v3-aligned-proteins.npy"]

#loop through each cluster
for afn in aligned_fns:
    aligned_prots = np.load(f"{directory}/processed-blast-v3/{afn}")
    #check if each query protein is in that cluster
    for sprot in simulated_prots:
        if sprot in aligned_prots:
            print(f"{sprot} is in the cluster of centroid {afn[0:4]}")
