import sys
import string

upperpath = "/home/jonathanb/mount"
directory = f"{upperpath}/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/moad_negatives/iofiles"

serial = 3
centroid='2HIZ'

alphabet = list(string.ascii_uppercase)

input = np.load(f"{directory}/processed-blast-v{serial}/{centroid}-v{serial}-holo-monomer.npy")

cmd.delete("all")
#cmd.fetch(centroid)

for j in alphabet:
    try:
        cmd.fetch(centroid+j)
        cmd.align(centroid+j, centroid)
    except:
        break

for i in input:
    print(i)
    if i != centroid: # and i not in ["1CNW", "1CNX","1CNY"]:

        for j in alphabet:
            try:
                cmd.fetch(i+j)
                cmd.align(i+j, centroid+"A")
            except:
                break
