#manually inspect things in pymol

import numpy as np

mtdir = "/home/jonathanb/mount"
directory = f"{mtdir}/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/moad_negatives/iofiles"
serial = 3

blast_outputs = np.load(f"{directory}/blast_output_v{serial}.npy", allow_pickle=True)

@cmd.extend
def checkmoad(x):
    x = int(x)

    cmd.delete("all")

    print(blast_outputs[x][0])

    cmd.fetch(blast_outputs[x][0][0])
    for i in blast_outputs[x][0][1:]:
        cmd.fetch(i)
        cmd.align(i, blast_outputs[x][0][0])
    cmd.center(blast_outputs[x][0][0])
