import numpy as np

#directory from which to load and save files
directory = "/project/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/moad_negatives/iofiles"

proteins = ['1YDA','1UIB','1E5O','1V2N','2HIZ','1OIT','6CD4','1OVE','2ZJW','3B0W']

input = np.load(f"{directory}/blast_output_v1nobreak.npy")

output = [i for i in input if i[0] in proteins]

np.save(f"{directory}/blast_output_v1nobreak_simonly.npy", output)
