import numpy as np
import os
import csv
import pickle
import sys

#human carbonic anhydrases from MOAD with 100% sequence identity
#this could also be loaded from /project/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/moad_negatives/iofiles/blast_output_v100pct.npy
#hch_moad = ['1CA3', '1CNW', '1CNX', '1CNY', '1EOU', '1G6V', '1HCA', '1KWQ', '1KWR', '1T9N', '1TB0', '1TBT', '1TE3', '1TEQ', '1TEU', '1XEG', '1XEV', '2AX2', '2EU2', '2EU3', '2EZ7', '2FMG', '2FMZ', '2GEH', '2H15', '2HKK', '2O4Z', '2POU', '2POW', '2Q1B', '2Q1Q', '2Q38', '2VVA', '2VVB', '3B4F', '3BL0', '3BL1', '3C7P', '3CAJ', '3CYU', '3D8W', '3D92', '3D93', '3D9Z', '3DAZ', '3DD0', '3DD8', '3EFI', '3EFT', '3F4X', '3F8E', '3FFP', '3HFP', '3HKN', '3HKQ', '3HKT', '3HKU', '3HLJ', '3HS4', '3IEO', '3IGP', '3K2F', '3K34', '3K7K', '3KKX', '3KOI', '3KOK', '3KON', '3KS3', '3KWA', '3L14', '3M1J', '3M2N', '3M3X', '3M40', '3M5E', '3M67', '3M96', '3M98', '3MHC', '3MHI', '3MHL', '3MHM', '3MHO', '3ML2', '3MMF', '3MNA', '3MWO', '3MYQ', '3MZC', '3N0N', '3N2P', '3N3J', '3N4B', '3NB5', '3NI5', '3NJ9', '3OIM', '3OKU', '3OY0', '3OYQ', '3OYS', '3P3H', '3P3J', '3P44', '3P4V', '3P55', '3P58', '3P5A', '3P5L', '3QYK', '3S8X', '3S9T', '3SAP', '3SAX', '3SBH', '3SBI', '3T82', '3T83', '3T84', '3T85', '3V2J', '3V2M', '3V5G', '3ZP9', '4CA2', '4CQ0', '4DZ7', '4DZ9', '4E3D', '4E3F', '4E3G', '4E3H', '4E49', '4E4A', '4E5Q', '4FPT', '4FRC', '4FU5', '4FVN', '4FVO', '4HT0', '4ITO', '4KNI', '4KNJ', '4LHI', '4LP6', '4MLT', '4MLX', '4MO8', '4MTY', '4N0X', '4N16', '4PQ7', '4PYX', '4PYY', '4PZH', '4Q49', '4Q6D', '4Q6E', '4Q7P', '4Q7S', '4Q7V', '4Q7W', '4Q81', '4Q83', '4Q87', '4Q8X', '4Q8Y', '4Q8Z', '4Q90', '4Q99', '4Q9Y', '4QIY', '4QJM', '4QSA', '4QSB', '4QSI', '4QTL', '4R5B', '4RFC', '4RFD', '4RH2', '4RN4', '4RUX', '4RUY', '4RUZ', '4WL4', '4WW6', '4XE1', '4YX4', '4YXI', '4YXO', '4YXU', '4YYT', '5A6H', '5BYI', '5C8I', '5DOG', '5DOH', '5DRS', '5DSK', '5DSL', '5DSM', '5DSO', '5DSP', '5DSQ', '5DSR', '5E2R', '5EH5', '5EH7', '5EH8', '5EHE', '5EHV', '5EHW', '5FDC', '5FDI', '5FLO', '5FLP', '5FLQ', '5FLR', '5FLS', '5FLT', '5FNG', '5FNH', '5FNI', '5FNJ', '5FNK', '5FNL', '5FNM', '5GMN', '5J8Z', '5JQ0', '5JQT', '5L3O', '5L6K', '5L6T', '5L70', '5L9E', '5LL4', '5LLC', '5LLE', '5LLG', '5LLH', '5LMD', '5N0D', '5N0E', '5N1R', '5N1S', '5N24', '5N25', '5NEA', '5NEE', '5O07', '5TFX', '5THJ', '5THN', '5TXY', '5TY1', '5TY8', '5TY9', '5TYA', '5U0D', '5U0E', '5U0F', '5U0G', '5ULN', '5UMC', '5VGY', '5WEX', '5Y2R', '5Y2S', '5YUI', '5YUJ', '5YUK', '6C7W', '6C7X', '6D1L', '6D1M', '6EQU', '6G3Q', '6G6T', '6GOT', '6H29', '6H2Z', '6H33', '6H34', '6H3Q', '6HD2', '6HX5', '6IC2', '6KM3', '6KM4', '6KM5', '6KM6', '6LUU', '6LUV', '6LUW', '6LUX', '6LUY', '6LUZ', '6LV1', '6LV2', '6LV3', '6LV4', '6LV5', '6LV6', '6LV7', '6LV8', '6LV9', '6LVA', '6ODZ', '6OE0', '6OE1', '6OTO', '6PDV', '6PEA', '6PGX', '6QEB', '6R6F', '6R6J', '6RG3', '6RG4', '6RHJ', '6RHK', '6RVF', '6RVK', '6RVL', '6RW1', '6SDS', '6SX9', '6SYB', '6SYS', '6T4N', '6T4O', '6T4P', '6T5C', '6T81', '6UFB', '6UFC', '6UFD', '6UX1', '6WKA', '6WQ4', '6WQ5', '6WQ7', '6WQ8', '6WQ9', '6XWZ', '6XXT', '6YMA', '6YMB', '6YPW', '6YQT', '6YQU', '6YRI', '6YZQ', '6YZR', '6YZS', '6YZT', '6YZU', '6YZV', '6YZW', '6YZX', '6ZR8', '7BI5', '7JNR', '7JNV', '7JNW', '7JNX', '7JNZ', '7K6I', '7K6J', '7K6K', '7K6L']

serial = 4

pdbid = '1yda'
loadpath = f"/project/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/moad_negatives/iofiles/processed-blast-v{serial}" #blast_output_v100pct-{pdbid}.npy"

#resid_positive_byprot
with open(f"{loadpath}/{pdbid.upper()}-v{serial}-lining-resis-byprot", 'rb') as pickle_file:
    resbyprot = pickle.load(pickle_file)

moad_directory = "/project/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/new_pockets/iofiles/moad_xml"

existing_xids_moad = [i[0:4] for i in os.listdir(moad_directory)] #get a list of already-downloaded ligand lists to avoid downloading extra copies

#get biological ligands for a given protein from MOAD
def moad_valid_ligands(struct, chain=""):

    #WARNING: The use of --no-check-certificate compromises security
    #ligand information (contains validity information not found in pdb structure)
    if struct not in existing_xids_moad:
        try:
            os.system(f"wget https://www.bindingmoad.org/files/csv/{struct.lower()}.csv -O {moad_directory}/{struct.lower()}.csv --no-check-certificate")
            existing_xids_moad.append(struct)
            #keeps the list of existing xml files up to date for when the code encounters apo candidates which are in moad and were previously loaded as holo candidates
        except:
            print(f"{struct} not in moad")
            return []

    valid_ligands = []

    chain_elems = chain.split("_")

    with open(f'{moad_directory}/{struct.lower()}.csv') as csvfile:
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
                valid_ligands.append(contents)

    return valid_ligands
    #if len(valid_ligands) > 0:
    #    return True
    #else:
    #    return False

#run through all the proteins to compile a list of ligands----------------------

resid_all = []
proteins_out = []
n_moad = 0

for prot in list(resbyprot.keys()):
    vligs = moad_valid_ligands(prot)
    if len(vligs)>0:
        n_moad += 1
        if not set(resbyprot[prot]).issubset(set(resid_all)):
            resid_all += list(resbyprot[prot])
            resid_all = list(np.unique(resid_all))
            proteins_out.append([prot, vligs, resbyprot[prot]])

#eliminate further redundancies:

#maxlen_init = len(proteins_out)
#proteins_out_2 = []

for x in range(len(proteins_out)):
    if x >= len(proteins_out):
        break
    for y in range(x+1, len(proteins_out)):
        if y >= len(proteins_out):
            break
        if set(proteins_out[x][2]).issubset(set(proteins_out[y][2])):
            proteins_out.pop(x)
            x -= 1
            y -= 1
        elif set(proteins_out[y][2]).issubset(set(proteins_out[x][2])):
            proteins_out.pop(y)
            x -= 1
            y -= 1

serial_out = 2
savedir = f"/project/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/moad_negatives/iofiles/processed-blast-v{serial}" #100pct-processed-blast"
np.save(f"{savedir}/{pdbid}_allresi_ligand_subset_v{serial_out}.npy", proteins_out)

#outputs------------------------------------------------------------------------

#print(resbyprot.keys()) #looks fine
print("-----------------------")
print(proteins_out)
print()
print(len(resbyprot.keys()))
print(n_moad)
print(len(proteins_out))
print(len(proteins_out))
print(len(proteins_out)/n_moad)
