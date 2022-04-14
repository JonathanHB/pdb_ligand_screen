import numpy as np
import itertools

import os
import csv
import pickle
import re
import csv

import mdtraj as md
from skbio import Protein, TabularMSA

#directory from which to load and save files
directory = "/project/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/moad_negatives/iofiles"
#moad file directory
blast_directory = "/project/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/new_pockets/iofiles"
#serial number of the blast output file to load
serial_in = 10 #this is up to date for things containing apo structures
#whether to to MSA here or load a previously saved one
align = True
#the distance below which ligands are considered adjacent to residues
cutoff = 0.5
#used to distinguish different versions of the results
#(i.e. results computed with different cluster sequence identity cutoffs)
serial_out = 3
#folder to keep output organized
output_dir = f"{directory}/processed-blast-v{serial_out}"

#load BLAST hits and pull out the set of interest
blasthits = np.load(f"{directory}/blast_output_v{serial_in}.npy")

#canonical residues reference
aa_resns = ["ASN", "ASP", "GLN", "GLU", "THR",
            "SER", "LYS", "ARG", "HIS", "PRO",
            "GLY", "CYS", "MET", "ALA", "VAL",
            "LEU", "ILE", "PHE", "TYR", "TRP",
            "SEC", "MSE"]

#-----------------------------------MUSCLE multiple sequence alignment code from Artur-----------------------------------#
#import this if possible rather than copying it here; imports

#import os
#import mdtraj as md
#import numpy as np
#from skbio import Protein, TabularMSA

def generate_msa(trajectories, protein_names, centroid_id, version='new'):
    '''
        trajectories : mdtraj trajectories for which you wish to do the sequence alignment.
            Code assumes that you will supply only one protein chain

    '''

    #check for selenocysteine here

    input_sequences = {}
    with open(f'{output_dir}/{centroid_id}-alignment.fasta', 'w') as f:
        for p, pdb in zip(protein_names, trajectories):
            prot = Protein(
                ''.join(r.code for r in pdb.top.residues if r.code),
                metadata={'id': p})

            #print(p)
            #print(''.join(r.code for r in pdb.top.residues if r.code))
            prot.write(f, format='fasta')
            input_sequences[p] = str(prot)

    muscledir = "/project/bowmore/ameller/"

    if version == 'old': #deprecated; could be removed
        os.system(f'{muscledir}muscle -in {output_dir}/alignment.fasta -out {output_dir}/alignment-aligned.fasta')
        os.system(f'{muscledir}muscle -refine -in {output_dir}/alignment-aligned.fasta -out {output_dir}/alignment-refined.fasta')
        msa = TabularMSA.read(f'{output_dir}/alignment-refined.fasta',
                              constructor=Protein)
    else:
        os.system(f'{muscledir}muscle -align {output_dir}/{centroid_id}-alignment.fasta -output {output_dir}/{centroid_id}-alignment-aligned.fasta')
        msa = TabularMSA.read(f'{output_dir}/{centroid_id}-alignment-aligned.fasta', constructor=Protein)

    msa.reassign_index(minter='id')

    starting_offsets = {}
    for p in protein_names:
        gap_filter = msa[msa.index == p][0].gaps()
        aligned_sequence = str(msa[msa.index == p][0][~gap_filter])
        starting_offsets[p] = aligned_sequence.find(input_sequences[p])
        # we should always find the aligned sequence in the input sequence
        assert starting_offsets[p] >= 0

    pdb2msa = {
        protein: {
            i + starting_offsets[protein]: r
            for i, r in enumerate(*np.where(~msa[msa.index == protein][0].gaps()))
        }
        for protein in msa.index
    }
    msa2rind = {
        protein: {
            r: i + starting_offsets[protein]
            for i, r in enumerate(*np.where(~msa[msa.index == protein][0].gaps()))
        }
        for protein in msa.index
    }

    return pdb2msa, msa2rind, msa

#------------------------------------------------------------------------------------------------------------------------#

#download pdb structures from rcsb iff they are not already present
def getstruct(struct, extant):

    if struct not in extant:
        try:
            os.system(f"wget https://files.rcsb.org/download/{struct}.pdb -O {directory}/rcsb_pdb/{struct}.pdb")

        except IndexError:
            #large structures cannot be downloaded as pdb files and are too big to be desirable for our present analysis
            print("error downloading file; it may be unavailable in .pdb format")

            #remove the empty file generated in the failed download attempt so that subsequent runs don't try to load it in mdtraj and subsequently crash
            os.system(f"rm {directory}/rcsb_pdb/{struct}.pdb")
            return False

#get the protein residues coordinating the specified ligands
#the prot argument is for debugging
def getcoordresseqs_moad(xtal, ligand_resns, ligand_rseqs):

    #construct ligand selection query
    molecule_queries = []
    for x, e in enumerate(ligand_resns):

        #normal 1-residue ligands
        if len(e.split(" ")) == 1:
            molecule_queries.append(f"(resname '{e}' and resSeq {ligand_rseqs[x]})")

        #composite multi-residue ligands for which MOAD does not provide adequate residue numbers
        else:
            #print(e)
            for resi in e.split(" "): #split moad contents into residues
                if resi not in aa_resns:
                    #print(resi)
                    molecule_queries.append(f"resname '{resi}'")
                else:
                    print(f"amino acid residue {resi} detected in composite ligand; \
                          cannot include in query due to the lack of adequate residue indices from moad")

            #print("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||")
            #print(ligand_resns)
            #print(ligand_rseqs)

    if len(molecule_queries) == 0:
        print("no valid ligands found")
        return False

    ligand_select_str = " or ".join(molecule_queries)

    #get indices of ligand and protein atoms
    ligand_ai = xtal.top.select(f"({ligand_select_str}) and not element H")
    protein_ai = xtal.top.select("protein and not element H")

    #get indices of ligand-coordinating protein atoms; index 0 indicates the 0th trajectory frame
    lining_ai = md.compute_neighbors(xtal, cutoff, ligand_ai, haystack_indices=protein_ai, periodic = False)[0]

    #get residue indices and numbers
    lining_resids = np.unique([xtal.top.atom(i).residue.index for i in lining_ai])    #0-indices of the lining residues in the mdtraj structure
    lining_resseqs = np.unique([xtal.top.atom(i).residue.resSeq for i in lining_ai])  #rcsb pdb residue numbers of the lining residues

    return [lining_resseqs, lining_resids]

#get biological ligands from MOAD
def getligs(struct, chain=""):

    #WARNING: The use of --no-check-certificate compromises security
    #ligand information (contains validity information not found in pdb structure)
    if struct not in existing_xids_moad:
        try:
            os.system(f"wget https://www.bindingmoad.org/files/csv/{struct.lower()}.csv -O {blast_directory}/moad_xml/{struct.lower()}.csv --no-check-certificate")
            existing_xids_moad.append(struct)
            #keeps the list of existing xml files up to date for when the code encounters apo candidates which are in moad and were previously loaded as holo candidates
        except:
            print(f"{struct} not in moad")
            return []

    bio_ligands = [] #get all the valid ligands from the moad structure

    chain_elems = chain.split("_")

    with open(f'{blast_directory}/moad_xml/{struct.lower()}.csv') as csvfile:
        reader = csv.reader(csvfile)
        firstrow = True
        for row in reader:
            if firstrow: #skip first row, which lacks ligand information
                firstrow = False
                continue

            contents = row[3].split(":")
            if (contents[1] in chain_elems or chain == "") and row[4]=="valid":
            #the latter case (chain == "") is for use of this method in get_struct(),
            #where it is used to obtain a list of all chains in a pdb structure containing biologically relevant ligands
                bio_ligands.append(contents)

    return bio_ligands

#separate biological monomers from multimers
def checkoligomerization(pdbid):

    author_found = False #prevent the software from reading any software-determined-only biomolecules
    look_for_mer = False #read only the first line after each biomolecule line to avoid reading the

    for line in open(f"{directory}/rcsb_pdb/{pdbid}.pdb"): #get pdb apo indices

        #determine if the biological unit is a monomer, look for the chains only if it is
        if look_for_mer and (line[0:45] == "REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT:" or (line[0:52] == "REMARK 350 SOFTWARE DETERMINED QUATERNARY STRUCTURE:" and (not author_found))):
            if line[0:45] == "REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT:":
                author_found = True

            asm_type = re.split(": | ", line)[6]

            if asm_type != "MONOMERIC":
                print(f"skipping {asm_type} assembly")
                return False
            else:
                #print(asm_type)
                return True

        elif look_for_mer and not (line[0:52] == "REMARK 350 SOFTWARE DETERMINED QUATERNARY STRUCTURE:" and author_found):
            print("missing biological unit")
            return False

        if line[0:23] == f"REMARK 350 BIOMOLECULE:":
            look_for_mer = True

################################################################################
#----------------------------------Main loop-----------------------------------#
################################################################################

#residues that cause problems:
#[DMO?? (not the PDB name), MLY, SUI]

#get a list of already-downloaded structures to avoid downloading extra copies
existing_xids_moad = [i[0:4] for i in os.listdir(f"{blast_directory}/moad_xml/")] #get a list of already-downloaded ligand lists to avoid downloading extra copies
existing_pdbids = [i[0:4] for i in os.listdir(f"{directory}/rcsb_pdb/")]

#usearch cluster indices (sorted by descending size) to process
start = 24
stop = len(blasthits)

for testindex in range(start, stop):
    testprots = blasthits[testindex]
    print("===============================================================================")
    print(f"cluster {testindex}, centroid {testprots[0]}")

#------------------------------------------------muscle msa-----------------------------------------

    #compute the msa and save the results
    if align == True:

        print("aligning sequences")

        #get lists of crystal structures and pdb ids for MSA
        xtal_aligned = []
        prot_aligned = []

        #loop through all proteins collected by BLAST
        for x, prot in enumerate(testprots[1]):
            print(f"{x}: {prot}")

            #download the protein structure, skip further processing if this fails
            if getstruct(prot, existing_pdbids) != False:

                #excluded_na_complexes = ["1T0K", "6DKS"]
                #if prot in excluded_na_complexes:
                #    print(f"skipping {prot}, which has RNA as the first chain")
                #    continue

                try:
                    xtal = md.load(f"{directory}/rcsb_pdb/{prot}.pdb")
                    xtal_0 = xtal.atom_slice(xtal.top.select("chainid 0"))

                    #skip protein-na complexes stating with na chains
                    if [r.code for r in xtal_0.top.residues if r.code] == []:
                        print(f"skipping {prot}, which is probably a protein-nucleic acid complex with a nucleic acid as the first chain")
                        continue

                    xtal_aligned.append(xtal_0)
                    prot_aligned.append(prot)
                except:
                    print(f"error loading {prot}; investigate manually")

        #throw out clusters with no valid structures
        if len(prot_aligned) == 0:
            print("no valid structures found, skipping to next usearch cluster")
            continue

        #generate and save the msa
        msagen = generate_msa(xtal_aligned, prot_aligned, testprots[0])

        savedict_pickle_1 = {
         f"pdb2msa":msagen[0],
         f"msa2rind":msagen[1],
         f"msa":msagen[2]}

        for e in savedict_pickle_1.keys():
            with open(f"{output_dir}/{testprots[0]}-v{serial_out}-{e}", "wb") as output_file:
                pickle.dump(savedict_pickle_1[e], output_file)

        np.save(f"{output_dir}/{testprots[0]}-v{serial_out}-aligned-proteins", prot_aligned)

        #separate the returned objects for further processing
        pdb2msa = msagen[0]
        msa2rind = msagen[1]

    #load the results of a saved msa
    else:

        with open(f"{output_dir}/{testprots[0]}-v{serial_out}-pdb2msa", 'rb') as pickle_file:
            pdb2msa = pickle.load(pickle_file)

        with open(f"{output_dir}/{testprots[0]}-v{serial_out}-msa2rind", 'rb') as pickle_file:
            msa2rind = pickle.load(pickle_file)

        prot_aligned = np.load(f"{output_dir}/{testprots[0]}-v{serial_out}-aligned-proteins.npy")

    #------------------------------------------------muscle msa end-------------------------------------

    #set reference protein to the centroid,
    #or the first protein (in alphabetic order) if the centroid is unavailable
    if testprots[0] in prot_aligned:
        refprot = testprots[0]
    else:
        refprot = prot_aligned[0]

    #proteins with canonical residues or a single n-terminal acetyl group only
    protein_canonical = []

    #apo = no moad ligand, holo = moad ligand
    #monomer/multimer is assigned from PDB remark 350
    prot_holo_monomer = []
    prot_apo_monomer = []
    prot_holo_multimer = []
    prot_apo_multimer = []

    #note: lining residues = residues near ligands
    #reference structure indices
    lining_resis_all = [] #pooled list of all rind2msa reference structure lining residues
    lining_resis_byprot = {} #rcsb pdb lining residue indices in the reference structure by source protein
    #msa indices
    lining_resis_all_msa = [] #pooled list of all msa lining residue indices
    lining_resis_byprot_msa = {} #msa lining residue indices by source protein

    #track all ligands encountered to determine the number of unique ligands
    all_moad_ligands = []

    #track the fraction of proteins which are monomers
    physio_monomer = 0

    #loop through all proteins collected by BLAST and included in the MSA
    for x, prot in enumerate(prot_aligned):
        print(f"{x}: {prot}")

        xtal = md.load(f"{directory}/rcsb_pdb/{prot}.pdb")

        #skip proteins with nonstandard residues because these mess up the MSA

        aa_resns = ["ASN", "ASP", "GLN", "GLU", "THR",
                    "SER", "LYS", "ARG", "HIS", "PRO",
                    "GLY", "CYS", "MET", "ALA", "VAL",
                    "LEU", "ILE", "PHE", "TYR", "TRP",
                    "SEC", "MSE"]

        #get the residues of the first mdtraj chain, which should be protein
        #this is not robust against weird chain orderings
        proteinresidues = [[j for j in i.residues] for i in xtal.top.chains][0]

        #full residue names
        residues = [str(i) for i in proteinresidues]
        #residue pdb numbers
        prot_rseqs = [i.resSeq for i in proteinresidues]
        #get residue names
        prot_resns = [i[0:-len(str(prot_rseqs[x]))] for x, i in enumerate(residues)]
        #noncanonical residues
        prot_resns_noncanonical = [i for i in prot_resns if i not in aa_resns]

        #check if the protein has noncanonical residues other than one n-terminal acetyl group, which is handled below
        if len(prot_resns_noncanonical) > 0 and not (len(prot_resns_noncanonical) == 1 and str(xtal.top.residue(0))[0:3] == "ACE"):
            print("skipping protein due to noncanonical residues:")
            print(prot_resns_noncanonical)
            continue

        protein_canonical.append(prot)

        #extract moad ligand residue types and ids for neighbor calculations
        moad_ligands = getligs(prot)
        moad_resns = [i[0] for i in moad_ligands]
        moad_rseqs = [int(i[2]) for i in moad_ligands]

        all_moad_ligands += moad_resns

        peptideligandonly = False

        #skip structures with no ligands of interest
        if moad_ligands != []:

            #get the residues within cutoff distance of the ligand
            lining_rseqs_resis = getcoordresseqs_moad(xtal, moad_resns, moad_rseqs)

            #handle cases where the ligands turn out to be peptides,
            #in which case MOAD does not supply reliable numbering
            if lining_rseqs_resis != False:

                #separate the returned objects for further processing
                lining_rseqs = lining_rseqs_resis[0] #rcsb pdb numbers
                lining_resis = lining_rseqs_resis[1] #mdtraj residue indices

                #acetylated n-termini cause off-by-one errors when converting residue
                #numbering using the multiple sequence alignment. This if statement
                #detects them and adjusts the numbering by 1 to correct this
                if str(xtal.top.residue(0))[0:3] == "ACE":
                    print("shifting residue numbering by one to correct for N-terminal acetylation")
                    lining_resis = [i+1 for i in lining_resis]

                #convert indices using msa

                #get the indices of the lining residues in the msa
                lining_aligned_msa = [pdb2msa[prot][i] for i in lining_resis if i in pdb2msa[prot]]

                #get the indices of the lining residues in the reference structure
                #these are the mdtraj indices
                lining_aligned = [msa2rind[refprot][pdb2msa[prot][i]] for i in lining_resis
                                  if (i in pdb2msa[prot] and pdb2msa[prot][i] in msa2rind[refprot])]

                #append indices for output
                #reference structure indices
                lining_resis_all += lining_aligned
                lining_resis_byprot[prot] = lining_aligned
                #msa indices
                lining_resis_all_msa += lining_aligned_msa
                lining_resis_byprot_msa[prot] = lining_aligned_msa

                if checkoligomerization(prot):
                    physio_monomer+=1
                    prot_holo_monomer.append(prot)
                else:
                    prot_holo_multimer.append(prot)

            else:
                peptideligandonly = True

        #append empty arrays for structures with no ligands of interest
        if moad_ligands == [] or peptideligandonly == True:
            #print("apo structure found, indicating filtering error; please inspect-------------------------------------------------")

            #prot_apo.append(prot)

            lining_resis_byprot[prot] = []
            lining_resis_byprot_msa[prot] = []

            if checkoligomerization(prot):
                physio_monomer+=1
                prot_apo_monomer.append(prot)
            else:
                prot_apo_multimer.append(prot)

    #remove duplicate residue indices
    lining_resis_all = np.unique(lining_resis_all)
    lining_resis_all_msa = np.unique(lining_resis_all_msa)


    #-------------------------------outputs---------------------------------------------------------
    print("-------------------------results----------------------------------------------")

    #-------------------------------oligomerization info--------------------------------------------

    n_canon = len(protein_canonical)
    multimer_summary = [n_canon, physio_monomer, physio_monomer/n_canon]

    #-------------------------------ligand composition and uniqueness-------------------------------

    ligand_summary = [all_moad_ligands, len(all_moad_ligands), len(np.unique(all_moad_ligands)), len(np.unique(all_moad_ligands))/len(all_moad_ligands)]

    #-------------------------------positive examples projected onto monomers for training----------------------------------

    #get the mdtraj residue indices of positive residues for each protein from the MSA
    resids_byprot_p = {}

    prot_monomer = prot_apo_monomer+prot_holo_monomer
    print(prot_monomer)

    for x, prot in enumerate(prot_monomer):
        lining_resids_prot = [msa2rind[prot][e] for e in lining_resis_all_msa if e in msa2rind[prot]]
        resids_byprot_p[prot] = [prot, lining_resids_prot]

    #-------------------------------negative examples projected onto for validation and testing--------------------

    #load reference structure and extract the 0th chain,
    #which ought to be protein for any normal structure
    xtal = md.load(f"{directory}/rcsb_pdb/{refprot}.pdb")
    xtal_0 = xtal.atom_slice(xtal.top.select("chainid 0"))

    #calculate pdb resseqs of ligand-lining residues for pymol display
    ref_rseqs_p = [xtal_0.top.residue(i).resSeq for i in lining_resis_all]
    ref_rseqs_n = [
        r.resSeq
        for r in xtal_0.top.residues
        if r.index not in lining_resis_all and r.is_protein
    ]

    #calculate indices of negative resiudes in reference structure
    ref_resid_n = [
        r.index
        for r in xtal_0.top.residues
        if r.index not in lining_resis_all and r.is_protein
    ]

    #--------------------------------save results---------------------------------------------------

    savedict_pickle_2 = {
    #the ligand-lining residues, projected onto the reference structure and msa respectively
    #"byprot" versions contain a dictionary of residues next to the ligand in each protein
     "lining-resis-byprot":lining_resis_byprot,
     "lining-resis-byprot":lining_resis_byprot_msa,
    #list of positive resids for all monomeric proteins
     "resid-positive-byprot":resids_byprot_p}

    for e in savedict_pickle_2.keys():
        with open(f"{output_dir}/{testprots[0]}-v{serial_out}-{e}", "wb") as output_file:
            pickle.dump(savedict_pickle_2[e], output_file)

    #variables to save and variable-specific filenames
    savedict_npy = {
    #lists of proteins in different categories;
    #apo and holo are mutually exclusive subsets of canonical, which is a subset of holo
     "aligned-proteins":prot_aligned,
     "canonical-proteins":protein_canonical,
     "holo-monomer":prot_holo_monomer,
     "apo-monomer":prot_apo_monomer,
     "holo-multimer":prot_holo_multimer,
     "apo-multimer":prot_apo_multimer,
    #the ligand-lining residues, projected onto the reference structure and msa respectively
    #"all" versions contain a list of all residues next to a ligand in any protein in the cluster
     "lining-resis-all":lining_resis_all,
     "lining-resis-all-msa":lining_resis_all_msa,
    #positive and negative resseqs and negative resids projected onto the reference protein
    #positive resids are saved as lining_resis_all (see above)
     "rseqs-positive-ref":ref_rseqs_p,
     "rseqs-negative-ref":ref_rseqs_n,
     "resid-negative-ref":ref_resid_n,
    #list of valid MOAD ligands and information about the fraction of unique ligands
     f"ligand-summary":ligand_summary,
    #information about the fraction of proteins which are multimers
     f"multimer-summary":multimer_summary}

    for e in savedict_npy.keys():
        np.save(f"{output_dir}/{testprots[0]}-v{serial_out}-{e}", savedict_npy[e])

    #save csv output
    with open(f'{output_dir}/usearch-v{serial_out}-negatives.csv', 'a', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')

        #centroid, usearch index, reference protein, number of aligned proteins,
        writer.writerow([testprots[0], testindex, refprot, len(prot_aligned),
        # number of aligned proteins with only canonical residues, number of holo proteins,
        len(protein_canonical), len(prot_holo_monomer)+len(prot_holo_multimer),
        # fraction of physiological monomers among canonical proteins
        (len(prot_holo_monomer)+len(prot_apo_monomer))/len(protein_canonical),
        #number of moad ligands, number of unique moad ligands,
        len(all_moad_ligands), len(np.unique(all_moad_ligands)),
        #length of reference protein, number of negative residues,
        len(ref_rseqs_p)+len(ref_rseqs_n), len(ref_rseqs_n),
        #format negatives for visualization
        "+".join(str(i) for i in ref_rseqs_n)
        ])

    #print negative residues for manual inspection
    print(f"resseqs of negative residues in pdb id {testprots[0]}: {'+'.join(str(i) for i in ref_rseqs_n)}\n")

################################################################################
#                               things to try
#------------------------------------------------------------------------------#

#try tem and see if we get the horn site and don't get the residues known not to label - done
#try the sars-2 macrodomain
#use muscle multiple sequence alignment - done
#Greg says to try simulating the output - partially done


################################################################################
#                                 trimmings
#------------------------------------------------------------------------------#
