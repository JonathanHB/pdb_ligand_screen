import numpy as np
#import mdtraj as md
import csv

#directory from which to load and save files
directory = "/project/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/moad_negatives/iofiles"
#moad file directory
blast_directory = "/project/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/new_pockets/iofiles"
#serial number of the blast output file to load
serial_in = 10
#set output subfolder
serial_proc = 3 #seria output number from processing script
output_dir = f"{directory}/processed-blast-v{serial_proc}"
#the usearch cluster indices for which to compile results
start = 0
stop = 5

blasthits = np.load(f"{directory}/blast_output_v{serial_in}.npy")

with open(f'{directory}/usearch-{start}-{stop}-negatives.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter=',') #, quotechar='', quoting=csv.QUOTE_MINIMAL

    for testindex in range(start, stop):

        testprots = blasthits[testindex]
        centroid = testprots[0]
        print("===============================================================================")
        print(f"cluster {testindex}, centroid {centroid}")

        #variables to save and variable-specific filenames
        ##------------------------------reference-------------------------------
        savedict = {
        #lists of proteins in different categories;
        #apo and holo are mutually exclusive subsets of canonical, which is a subset of holo
         "aligned-proteins":prot_aligned,
         "canonical-proteins":protein_canonical,
         "holo-monomer":prot_holo_monomer,
         "apo-monomer":prot_apo_monomer,
        #the ligand-lining residues, projected onto the reference structure and msa respectively
        #"byprot" versions contain a dictionary of residues next to the ligand in each protein and
        #"all" versions contain a list of all residues next to a ligand in any protein in the cluster
         "lining-resis-byprot":lining_resis_byprot,
         "lining-resis-byprot":lining_resis_byprot_msa,
         "lining-resis-all":lining_resis_all,
         "lining-resis-all_msa":lining_resis_all_msa,
        #positive and negative resseqs and negative resids projected onto the reference protein
        #positive resids are saved as lining_resis_all (see above)
         "rseqs-positive_ref":ref_rseqs_p,
         "rseqs-negative_ref":ref_rseqs_n,
         "resid-negative_ref":ref_resid_n,
        #list of positive resids for all monomeric proteins
         "resid-positive-byprot":resids_byprot_p,
        #list of valid MOAD ligands and information about the fraction of unique ligands
         f"ligand-summary":ligand_summary,
        #information about the fraction of proteins which are multimers
         f"multimer-summary":multimer_summary}

        try:
            = np.load(f"{directory}/{centroid}-v{serial_proc}-.npy")


            #holo_prots = np.load(f"{directory}/{centroid}-v{serial_proc}-holo-proteins.npy")

            # = np.load(f"{directory}/{centroid}-v{serial_proc}-.npy")

            #resid_negative = np.load(f"{directory}/{centroid}-resid-negative.npy")
            #rseqs_negative = np.load(f"{directory}/{centroid}-rseqs-negative.npy")
            #rseqs_positive = np.load(f"{directory}/{centroid}-rseqs-positive.npy")

        except:
            print(f"output files for {centroid} not found; skipping to next usearch cluster")
            continue

        #things to return:
        # centroid PDB ID
        # usearch cluster index
        # number of monomeric holo proteins with moad ligands
        # number of [resolved] centroid residues
        # number of centroid negative residues
        # centroid negative resseqs for pymol visualization
        #



        outputlist = [centroid, testindex, len(aligned_proteins), len(apo_monomer)+len(holo_monomer), ]

        writer.writerow([testindex, centroid, len(holo_prots), len(rseqs_negative)+len(rseqs_positive), len(rseqs_negative), "+".join(str(i) for i in rseqs_negative)]) #resid_negative, rseqs_negative, holo_prots
