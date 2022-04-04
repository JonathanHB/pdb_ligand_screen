import numpy as np
import itertools

import os
import csv
import pickle
import re

import mdtraj as md

#directory from which to load and save files
directory = "/project/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/moad_negatives/iofiles"
#moad file directory
blast_directory = "/project/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/new_pockets/iofiles"
#serial number of the blast output file to load
serial_in = 6
#the index of the protein set in the blast file to use
#testindex = 0
#whether to to MSA here or load a previously saved one
align = False
#the distance below which ligands are considered adjacent to residues
cutoff = 0.5

#load BLAST hits and pull out the set of interest
blasthits = np.load(f"{directory}/blast_output_v{serial_in}.npy")

#-----------------------------------MUSCLE multiple sequence alignment code from Artur-----------------------------------#
#import this if possible rather than copying it here

import os
import mdtraj as md
import numpy as np
from skbio import Protein, TabularMSA

def generate_msa(trajectories, protein_names, version='new'):
    '''
        trajectories : mdtraj trajectories for which you wish to do the sequence alignment.
            Code assumes that you will supply only one protein chain

    '''
    input_sequences = {}
    with open('alignment.fasta', 'w') as f:
        for p, pdb in zip(protein_names, trajectories):
            prot = Protein(
                ''.join(r.code for r in pdb.top.residues if r.code),
                metadata={'id': p})
            prot.write(f, format='fasta')
            input_sequences[p] = str(prot)

    muscledir = "/project/bowmore/ameller/"
    if version == 'old':
        os.system(f'{muscledir}muscle -in alignment.fasta -out alignment-aligned.fasta')
        os.system(f'{muscledir}muscle -refine -in alignment-aligned.fasta -out alignment-refined.fasta')
        msa = TabularMSA.read('alignment-refined.fasta',
                              constructor=Protein)
    else:
        os.system(f'{muscledir}muscle -align alignment.fasta -output alignment-aligned.fasta')
        msa = TabularMSA.read('alignment-aligned.fasta',
                              constructor=Protein)

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
def getcoordresseqs(xtal, ligand_resns, ligand_rseqs, ligand_nonaq_indices):

    #construct ligand selection
    molecule_queries = []
    for x in ligand_nonaq_indices:
        molecule_queries.append(f"(resname '{ligand_resns[x]}' and resSeq {ligand_rseqs[x]})")

    ligand_select_str = " or ".join(molecule_queries)

    #get indices of ligand and protein atoms
    ligand_ai = xtal.top.select(f"{ligand_select_str} and not element H")
    protein_ai = xtal.top.select("protein and not element H")

    #get indices of ligand-coordinating protein atoms; index 0 indicates the 0th trajectory frame
    lining_ai = md.compute_neighbors(xtal, cutoff, ligand_ai, haystack_indices=protein_ai, periodic = False)[0]

    #get residue indices and numbers
    lining_resids = [xtal.top.atom(i).residue.index for i in lining_ai]    #0-indices of the lining residues in the mdtraj structure
    lining_resseqs = [xtal.top.atom(i).residue.resSeq for i in lining_ai]  #rcsb pdb residue numbers of the lining residues

    return [lining_resseqs, lining_resids]

#get the protein residues coordinating the specified ligands
def getcoordresseqs_moad(xtal, ligand_resns, ligand_rseqs):

    #construct ligand selection query
    molecule_queries = []
    for x, e in enumerate(ligand_resns):
        molecule_queries.append(f"(resname '{e}' and resSeq {ligand_rseqs[x]})")

    ligand_select_str = " or ".join(molecule_queries)

    #get indices of ligand and protein atoms
    ligand_ai = xtal.top.select(f"{ligand_select_str} and not element H")
    protein_ai = xtal.top.select("protein and not element H")

    #get indices of ligand-coordinating protein atoms; index 0 indicates the 0th trajectory frame
    lining_ai = md.compute_neighbors(xtal, cutoff, ligand_ai, haystack_indices=protein_ai, periodic = False)[0]

    #get residue indices and numbers
    lining_resids = np.unique([xtal.top.atom(i).residue.index for i in lining_ai])    #0-indices of the lining residues in the mdtraj structure
    lining_resseqs = np.unique([xtal.top.atom(i).residue.resSeq for i in lining_ai])  #rcsb pdb residue numbers of the lining residues

    return [lining_resseqs, lining_resids]

#get structural and biological ligands for getprot()
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

    return bio_ligands #[keep_ligands, bio_ligands, all_ligands, nst_ligands]

def checkoligomerization(pdbid):

    #--------------------------------------------check resolution and break structure into monomeric chains, if any---------------------------------------------------------------------------------------------

    author_found = False #prevent the software from reading any software-determined-only biomolecules
    look_for_mer = False #read only the first line after each biomolecule line to avoid reading the

    for line in open(f"{directory}/rcsb_pdb/{pdbid}.pdb"): #get pdb apo indices

        #print(line)

        #determine if the biological unit is a monomer, look for the chains only if it is
        if look_for_mer and (line[0:45] == "REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT:" or (line[0:52] == "REMARK 350 SOFTWARE DETERMINED QUATERNARY STRUCTURE:" and (not author_found))):
            if line[0:45] == "REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT:":
                author_found = True

            asm_type = re.split(": | ", line)[6]
            #return(asm_type)
            if asm_type != "MONOMERIC":
                print(f"skipping {asm_type} assembly")
                return False
            else:
                print(asm_type)
                return True
            #look_for_chains = True

        elif look_for_mer and not (line[0:52] == "REMARK 350 SOFTWARE DETERMINED QUATERNARY STRUCTURE:" and author_found):
            print("missing biological unit")
            return False

        if line[0:23] == f"REMARK 350 BIOMOLECULE:": # {biomolecule+1} #did not handle variable numbers of spaces well
            look_for_mer = True

################################################################################
#----------------------------------Main loop-----------------------------------#
################################################################################

#get a list of already-downloaded structures to avoid downloading extra copies
existing_xids_moad = [i[0:4] for i in os.listdir(f"{blast_directory}/moad_xml/")] #get a list of already-downloaded ligand lists to avoid downloading extra copies
existing_pdbids = [i[0:4] for i in os.listdir(f"{directory}/rcsb_pdb/")]

for testindex in range(0, 20):
    testprots = blasthits[testindex]

#------------------------------------------------test muscle msa-----------------------------------------

#compute the msa and save the results
    if align == True:

        print("aligning sequences")

        #get lists of crystal structures and pdb ids for MSA
        xtal_all = []
        prot_all = []

        #loop through all proteins collected by BLAST
        for x, prot in enumerate(testprots[1]):
            print(f"{x}: {prot}")

            #download the protein structure, skip further processing if this fails
            if getstruct(prot, existing_pdbids) != False:

                #load the protein structure and get the first chain; skip the protein if this fails
                try:
                    xtal = md.load(f"{directory}/rcsb_pdb/{prot}.pdb")
                    xtal_0 = xtal.atom_slice(xtal.top.select("chainid 0"))
                    xtal_all.append(xtal_0)
                    prot_all.append(prot)
                except:
                    print(f"error loading {prot}")

        #generate and save the msa
        msagen = generate_msa(xtal_all, prot_all)

        savedict = {
         f"{testprots[0]}-pdb2msa":msagen[0],
         f"{testprots[0]}-msa2rind":msagen[1],
         f"{testprots[0]}-msa":msagen[2]}

        for e in savedict.keys():
            with open(e, "wb") as output_file:
                pickle.dump(savedict[e], output_file)

        np.save(f"{testprots[0]}-aligned-proteins", prot_all)



        #np.save(f"{directory}/msagen-{testprots[0]}", [prot_all, msagen])

        #separate the returned objects for further processing
        pdb2msa = msagen[0]
        msa2rind = msagen[1]

    #load the results of a saved msa
    else:

        with open(f"{testprots[0]}-pdb2msa", 'rb') as pickle_file:
            pdb2msa = pickle.load(pickle_file)

        with open(f"{testprots[0]}-msa2rind", 'rb') as pickle_file:
            msa2rind = pickle.load(pickle_file)
        #pdb2msa = pickle.load(f"{testprots[0]}-pdb2msa")
        #msa2rind = pickle.load(f"{testprots[0]}-msa2rind")
        prot_all = np.load(f"{testprots[0]}-aligned-proteins.npy")

        #msagenprots = np.load(f"{directory}/msagen-{testprots[0]}.npy")

        #separate the returned objects for further processing
        #prot_all = msagenprots[0]
        #pdb2msa = msagenprots[1][0]
        #msa2rind = msagenprots[1][1]

    #------------------------------------------------test muscle msa end-------------------------------------

    #residues near ligands (= lining residues)
    prot_holo = [] #all proteins with relevant ligands
    lining_resis_all = [] #pooled list of all rind2msa reference structure lining residues
    lining_resis_byprot = [] #rcsb pdb lining residue numbers by source protein

    #loop through all proteins collected by BLAST and included in the MSA

    #prot_all = ["12CA","1CAN","6EQU"] #for debugging

    for x, prot in enumerate(prot_all):
        print(f"{x}: {prot}")

        #skip multimers
        if not checkoligomerization(prot):
            lining_resis_byprot.append([-999]) #to distinguish multimers from aligned structures with no MOAD ligands
            continue

        xtal = md.load(f"{directory}/rcsb_pdb/{prot}.pdb")

        moad_ligands = getligs(prot)
        moad_resns = [i[0] for i in moad_ligands]
        moad_rseqs = [int(i[2]) for i in moad_ligands]

        #print(moad_rseqs)
        #print(moad_ligands)

        #skip structures with more than two mdtraj chains
        # if xtal.top.n_chains != 2:
        #     #add an entry to the list so that protein numbers stay synchronized with their residue lists
        #     #-99 serves to distinguish the entry from empty lists for proteins with no ligands which were still processed successfully
        #     lining_resis_byprot.append([-99])
        #     print(f"skipping {prot} with {xtal.top.n_chains} chains")
        #     continue

        #get the residues of the second mdtraj chain, which should contain the ligand and solvent
        #ligandsolvent = [[j for j in i.residues] for i in xtal.top.chains][1]


        #aqueous solvents, gases, and ions which are unlikely to open cryptic pockets,
        #likely to cover most of the surface, or unrepresentative of biological ligands
        #solvent = ["HOH","DOD","NH3","NH4","CO3","NO3","PO4","SO4","LI","BE","K","NA","MG","CA","F","CL","BR","I","ZN","MN","CU","FE","O2","HDZ","CO2","CO2","XE","UNX"]

        #full residue names
        #ligands = [str(i) for i in ligandsolvent]
        #residue pdb numbers
        #ligand_rseqs = [i.resSeq for i in ligandsolvent]
        #get residue names
        #ligand_resns = [i[0:-len(str(ligand_rseqs[x]))] for x, i in enumerate(ligands)]
        #indices of residues which are of biological interest (not in the list above)
        #ligand_nonaq_indices = [x for x, i in enumerate(ligand_resns) if i in moad_ligands] # if i not in solvent]

        #skip structures with no ligands of interest
        if moad_ligands != []:

            prot_holo.append(prot)

            #get the residues within cutoff distance of the ligand
            lining_rseqs_resis = getcoordresseqs_moad(xtal, moad_resns, moad_rseqs)

            #separate the returned objects for further processing
            lining_rseqs = lining_rseqs_resis[0] #rcsb pdb numbers
            lining_resis = lining_rseqs_resis[1] #mdtraj residue indices

            #print(lining_resis)
            #print(lining_rseqs)

            #acetylated n-termini cause off-by-one errors when converting residue
            #numbering using the multiple sequence alignment. This if statement
            #detects them and adjusts the numbering by 1 to correct this
            if str(xtal.top.residue(0))[0:3] == "ACE":
                print("shifting residue numbering by one to correct for N-terminal acetylation")
                lining_resis = [i+1 for i in lining_resis]

            #get the indices of the lining residues in the reference structure
            #these are the mdtraj indices
            #print(pdb2msa.keys())
            #print(prot)
            #print(pdb2msa[prot])
            #print([pdb2msa[prot][i] for i in lining_resis])

            lining_aligned = [msa2rind[testprots[0]][pdb2msa[prot][i]] for i in lining_resis
                              if (i in pdb2msa[prot] and pdb2msa[prot][i] in msa2rind[testprots[0]])]

            #append indices for output
            lining_resis_all += lining_aligned
            lining_resis_byprot.append(lining_aligned)

        #append empty arrays for structures with no ligands of interest
        else:
            lining_resis_byprot.append([])

    #remove duplicate residue indices
    lining_resis_all = np.unique(lining_resis_all)

    #-------------------------------outputs-----------------------------------------

    #reference structure pdb id
    print("-----------------------------------------------------------------------")

    #load reference structure and extract the 0th chain,
    #which ought to be protein for any normal structure
    xtal = md.load(f"{directory}/rcsb_pdb/{testprots[0]}.pdb")
    xtal_0 = xtal.atom_slice(xtal.top.select("chainid 0"))

    #calculate pdb resseqs of ligand-lining residues for pymol display
    output_rseqs_p = [xtal_0.top.residue(i).resSeq for i in lining_resis_all]
    output_rseqs_n = [
        r.resSeq
        for r in xtal_0.top.residues
        if r.index not in lining_resis_all and r.is_protein
    ]

    #calculate indices of negative resiude in reference structure
    output_resid_n = [
        r.index
        for r in xtal_0.top.residues
        if r.index not in lining_resis_all and r.is_protein
    ]

    #save results
    savedict = {
     f"{testprots[0]}-aligned-proteins":prot_all,
     f"{testprots[0]}-holo-proteins":prot_holo,
     f"{testprots[0]}-lining-resis-byprot":lining_resis_byprot,
     f"{testprots[0]}-lining-resis-all":lining_resis_all,
     f"{testprots[0]}-rseqs-positive":output_rseqs_p,
     f"{testprots[0]}-rseqs-negative":output_rseqs_n,
     f"{testprots[0]}-resid-negative":output_resid_n}

    for e in savedict.keys():
        np.save(e, savedict[e])

    #print negative residues for manual inspection
    print(f"resseqs of lining residues in pdb id {testprots[0]}: {'+'.join(str(i) for i in output_rseqs_n)}")

################################################################################
#                               things to try
#------------------------------------------------------------------------------#

#try tem and see if we get the horn site and don't get the residues known not to label - done
#try the sars-2 macrodomain
#use muscle multiple sequence alignment - done
#Greg says to try simulating the output - partially done

################################################################################


#-------------------------------trimmings---------------------------------------


    #"{testprots[0]}-centroid_pdb_id":testprots[0], # <--unnecessary since this is in the filename


#output_resid_p = lining_resis_all
#"+".join(str(xtal.top.residue(i).resSeq) for i in lining_resis_all if xtal.top.residue(i).is_protein)]
#np.save(f"{directory}/outputs-{testprots[0]}", [testprots[0], prot_all, lining_resis_byprot, lining_resis_all, output_rseqs_p, output_rseqs_n, output_resid_n])

    #if type == "holo":
    #    holo_ligs = getligs(struct, "")
    #    lig_chains = [l[1] for l in holo_ligs]

        #
        # #check the resolution, exit if it's inadequate
        # if line[0:22] == "REMARK   2 RESOLUTION.":
        #     #print(line.split(" "))
        #     pdblist = line.split(" ")
        #     if pdblist[8] == "NOT" and (pdblist[9] == "AVAILABLE." or pdblist[9] == "APPLICABLE."):
        #         print("resolution unavailable")
        #         #return False
        #     else:
        #         try:
        #             if float(pdblist[8]) > 2.5: #if the resolution is inadequate
        #                 print("inadequate resolution of %s A" % float(pdblist[8]))
        #                 #return False
        #             else: #if the resolution is adequate
        #                 #print("resolution of %s A" % float(pdblist[8])) #not really useful for debugging
        #                 res_checked = True
        #         except ValueError:
        #             print(f"nonstandard REMARK 2: {line}")
        #             #return False
        #

#print aligned 1can resseqs only for debugging:
#print(f"{prot_all[1]} lining residue resseqs in 12ca: " + "+".join(str(xtal.top.residue(i).resSeq) for i in np.unique(lining_resis_byprot[1])))
#print(f"{prot_all[2]} lining residue resseqs in 12ca: " + "+".join(str(xtal.top.residue(i).resSeq) for i in np.unique(lining_resis_byprot[2])))


#aa = [str(xtal.top.residue(i).resSeq) for i in lining_resis_all]
#bb = [i for i in lining_resis_all]

#print([f"{aa[x]}, {bb[x]}" for x in range(len(aa))])

# #>1CAN
# ------SHHWGYGKHNGPEHWHKDFPIAKGERQSPVDIDTHTAKYDPSLKPLSVSYDQATSLRILNNGHAFNVEFDDSQD
# KAVLKGGPLDGTYRLIQFHFHWGSLDGQGSEHTVDKKKYAAELHLVHWNTKYGDFGKAVQQPDGLAVLGIFLKVGSAKPG
# LQKVVDVLDSIKTKGKSADFTNFDPRGLLPESLDYWTYPGSLTTPPLLECVTWIVLKEPISVSSEQVLKFRKLNFNGEGE
# PEELMVDNWRPAQPLKNRQIKASFK
# >12CA
# ---------WGYGKHNGPEHWHKDFPIAKGERQSPVDIDTHTAKYDPSLKPLSVSYDQATSLRILNNGHAFNVEFDDSQD
# KAVLKGGPLDGTYRLIQFHFHWGSLDGQGSEHTVDKKKYAAELHLAHWNTKYGDFGKAVQQPDGLAVLGIFLKVGSAKPG
# LQKVVDVLDSIKTKGKSADFTNFDPRGLLPESLDYWTYPGSLTTPPLLECVTWIVLKEPISVSSEQVLKFRKLNFNGEGE
# PEELMVDNWRPAQPLKNRQIKASF-

    # #inspect the first residue of each protein to determine which structures contain acetylated n-termini which cause off-by-1 errors in the msa
    # aa_resns = ["ASN", "ASP", "GLN", "GLU", "THR",
    #             "SER", "LYS", "ARG", "HIS", "PRO",
    #             "GLY", "CYS", "MET", "ALA", "VAL",
    #             "LEU", "ILE", "PHE", "TYR", "TRP"]
    #
    # res0 = str(xtal.top.residue(0))
    # if res0[0:3] not in aa_resns:
    #     print(f"{res0}-------------------------------------")
    #


        #mark where a particular residue is found next to the ligand
        #testresi = 115
        #if testresi in lining_aligned:
        #    print(f"--------------------{testresi} in aligned lining residues-----------------------------")


#maxind = np.max(msa2rind[prot_all[0]]) #the maximum reference index; used to generate the dictionary below
#pdbs_byresi = {[x, [] for x in range(0, maxind)]} #PDB ids with ligands near each rind2msa reference structure residue
