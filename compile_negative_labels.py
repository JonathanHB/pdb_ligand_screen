import numpy as np
#import mdtraj as md
import csv

#directory from which to load and save files
directory = "/project/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/moad_negatives/iofiles"
#moad file directory
blast_directory = "/project/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/new_pockets/iofiles"
#serial number of the blast output file to load
serial_in = 6
#the usearch cluster indices for which to compile results
start = 16
stop = 22

blasthits = np.load(f"{directory}/blast_output_v{serial_in}.npy")

with open(f'{directory}/usearch-{start}-{stop}-negatives.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter=',') #, quotechar='', quoting=csv.QUOTE_MINIMAL

    for testindex in range(start, stop):

        testprots = blasthits[testindex]
        centroid = testprots[0]
        print("===============================================================================")
        print(f"cluster {testindex}, centroid {centroid}")

        try:
            holo_prots = np.load(f"{directory}/{centroid}-holo-proteins.npy")
            resid_negative = np.load(f"{directory}/{centroid}-resid-negative.npy")
            rseqs_negative = np.load(f"{directory}/{centroid}-rseqs-negative.npy")
            rseqs_positive = np.load(f"{directory}/{centroid}-rseqs-positive.npy")

        except:
            print("output files not found; skipping to next usearch cluster")
            continue
        #usearch cluster index, centroid PDB ID, number of proteins with moad ligands, number of [resolved] centroid residues, number of centroid negative residues, centroid negative resseqs for pymol visualization
        writer.writerow([testindex, centroid, len(holo_prots), len(rseqs_negative)+len(rseqs_positive), len(rseqs_negative), "+".join(str(i) for i in rseqs_negative)]) #resid_negative, rseqs_negative, holo_prots
