import re
import numpy as np
import os
import sys
import Bio.Blast
from Bio.Blast import NCBIXML

#---------------------------------usearch commands------------------------------
#these commands were run locally over the mount to generate the input files for this script

#combine FASTA files from MOAD (https://www.bindingmoad.org/Home/download)
#find ../new_pockets/iofiles/rcsb_fasta_all -type f -exec cat {} \; > moad_fasta_all

#cluster MOAD proteins in sequence space to 95% identity and return the centroids
#make sure to add usearch to the PATH first (
#export PATH=$PATH:/home/jonathanb/software/usearch11.0.667_i86linux32
#)
#usearch11.0.667_i86linux32 -cluster_fast moad_fasta_all -id 0.95 -centroids moad_fasta_centroids
#to save an msa file for each cluster add -msaout [filename stem]
#this is still very fast

#make a list of pdb ids of the cluster centroids
#grep '>' moad_fasta_centroids | cut -d_ -f1 | cut -d'>' -f2 > centroid_pdb_ids

#-------------------------------------------------------------------------------

#for each centroid:
    #download the fasta sequence
    #blastp it against PDB structures
    #save the resulting structures with single hsps 100% sequence identity and coverage

#path to directory
directory = "/project/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/new_pockets/iofiles"

part = "all" #part a or b or both ("all") of the MOAD database, which is broken in two on account of its large size

idthresh = 0.90

#get the PDB IDs of all the MOAD structures
structures = np.hstack([np.unique([i[0:4].upper() for i in os.listdir(f"{directory}/every_part_a/BindingMOAD_2020/")]), np.unique([i[0:4].upper() for i in os.listdir(f"{directory}/every_part_b/BindingMOAD_2020/")])])
print(f"{len(structures)} structures")

#structures for which fasta files have already been downloaded
existing_ids = [i[0:4] for i in os.listdir(f"{directory}/rcsb_fasta_{part}/")]

hits_all = []
hits_alldata = []

debugcounter = 0

#loop through centroid structures
for i in open("/project/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/moad_negatives/centroid_pdb_ids_p90", 'r'):
#for i in structures:

    i = i[0:4] #remove \n from each line

    print("-----------------------------------------------------------------------------")
    print(debugcounter)
    print(f"query: {i}")

    #if debugcounter > 200:
    #    print("done")
    #    break
    debugcounter+=1

    if i in hits_all: #skip MOAD structures previously encountered as BLAST hits
        print("structure previously encountered as a BLAST hit; error")
        #continue

    hits = []

    #download fasta files if necessary
    if i not in existing_ids:
        os.system(f"wget https://www.rcsb.org/fasta/entry/{i}/display -O {directory}/rcsb_fasta_{part}/{i}.fasta")

    #blast the MOAD sequence against PDB
    #note that the ftp version of blastp must be used because the apt-get version is out of date and can't load the pdbaa database.
    #-outfmt 5 is required to produce a biopython-readable .xml file
    os.system(f"/project/bowmanlab/borowsky.jonathan/installations/ncbi-blast-2.12.0+/bin/blastp -db /project/bowmanlab/borowsky.jonathan/installations/pdbaa -query {directory}/rcsb_fasta_{part}/{i}.fasta -outfmt 5 -out {directory}/blast_xml_out/{i}-results.xml")

    #load the blastp results
    result_handle = open(f"{directory}/blast_xml_out/{i}-results.xml")

    #skip structures containing multiple different protein sequences
    #As multiple proteins are usually added to crystallization solutions for the purpose of studying complexes thereof,
    #there are likely relatively few multi protein PDB structures with monomeric bioassemblies
    #so I'm skipping them for the time being for ease of processing
    #The following can be used to read xml results for multimers, but further downstream adjustments would be required
    #  "blast_records = NCBIXML.parse(result_handle)
    #   blast_record = next(blast_records)"
    try:
        blast_record = NCBIXML.read(result_handle)
    except ValueError as err:
        print(f"skipping query structure containing multiple sequences: {err}")
        continue
    # see https://biopython.org/docs/1.75/api/Bio.Blast.Record.html for documentation

    #print(dir(blast_record)) #print a list of the record's contents

    #for each structure aligned to the moad structure
    for alignment in blast_record.alignments:

        #print(alignment.title)
        hsp = alignment.hsps[0] #there should be only one hsp for any alignment good enough to be useful; this is checked below
        #print(hsp)

        #check that there is one hsp with 100% identity and coverage
        if len(alignment.hsps) == 1 and hsp.identities/hsp.align_length > idthresh: # == hsp.align_length and hsp.align_length == alignment.length:
            print(alignment)
            print(hsp.identities) #equals the others here by transitivity

            entries = re.split("pdb\||>pdb\|", alignment.title) #separate out all pdb files with a given sequence
            pdb_ids = [title[0:4] for title in entries[1:]] #the first entry is "" from the left of the first pdb| header
            #print(len(pdb_ids))
            pdb_ids = [i for i in pdb_ids if i.upper() in structures]
            #print(len(pdb_ids))
            #print("#########################################################################")


            hits = hits+pdb_ids
        else:
            #print(len(alignment.hsps))
            #print(hsp.identities/alignment.length)
            print(alignment.title)
            print(hsp.identities)
            print(hsp.align_length)
            print(alignment.length)

            break # <-- heuristic that assumes that once a nonidentical alignment is encountered none of the remainder are any good either

            #print("---------------------")
            #print(hsp.identities)
            #print(hsp.align_length)
            #print(alignment.length)

    hits = list(np.unique(hits))
    hits_all += hits
    hits_alldata.append([i, hits, len(hits)])

    print(hits)

hits_alldata.sort(key = lambda x: x[2], reverse=True)

serial = 9 #<-- keep updated
np.save(f"/project/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/moad_negatives/iofiles/blast_output_v{serial}", hits_alldata)

print(hits_alldata)
