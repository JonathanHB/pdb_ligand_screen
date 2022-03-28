import os
import subprocess

dir = "/project/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/moad_negatives"

fnames = os.listdir(f"{dir}/usearch_msa")

lengths = []

for i in fnames:

    proc=subprocess.Popen(f"grep '>' usearch_msa/{i} | wc -l", shell=True, stdout=subprocess.PIPE)
    output=proc.communicate()[0]
    nprots = output.decode("utf-8").strip()

    lengths.append(i, nprots)

    #os.system(f"grep '>' usearch_msa/{i} | wc -l")
    #break
    #grep '>' moad_fasta_centroids | cut -d_ -f1 | cut -d'>' -f2 > centroid_pdb_ids

#for i in open("/project/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/moad_negatives/centroid_pdb_ids", 'r'):

#grep '>' moad_fasta_centroids | cut -d_ -f1 | cut -d'>' -f2
