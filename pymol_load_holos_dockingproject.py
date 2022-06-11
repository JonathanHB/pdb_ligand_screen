import sys
import string

upperpath = "/home/jonathanb/mount"
directory = f"{upperpath}/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/moad_negatives/iofiles/100pct-processed-blast"
protname = "human-brd4-bromodomain" #human-carbonic-anhydrase

serial = 4

protlist = np.load(f"{directory}/{protname}_kd.npy")

#centroid='1OIT'

#alphabet = list(string.ascii_uppercase)

#input = list(np.load(f"{directory}/processed-blast-v{serial}/{centroid}-v{serial}-holo-monomer.npy"))\
#+list(np.load(f"{directory}/processed-blast-v{serial}/{centroid}-v{serial}-holo-multimer.npy"))

cmd.delete("all")
#cmd.fetch(centroid)

#for j in alphabet:
#    try:
#        cmd.fetch(centroid+j)
#        cmd.align(centroid+j, centroid)
#    except:
#        break

prot0 = protlist[0][0]+protlist[0][2]

for x, prot in enumerate(protlist):

    print(prot)

    pname = prot[0]+prot[2]

    cmd.fetch(pname)

    if x != 0:
        cmd.align(pname, prot0)

    cmd.hide("sticks", f"{pname} and not resn {prot[1]}")
    # and i not in ["1CNW", "1CNX","1CNY"]:
    #    for j in alphabet:
    #        try:
    #            cmd.fetch(i+j)
    #            cmd.align(i+j, centroid+"A")
    #        except:
    #            break


#graphics settings and hiding highly soluble ligands
print(f"hide sticks, resn GOL+EDO+FMT+DMS+ACT+NO3 or elem H; hide spheres, resn NA+CL+K+PO4+SO4+NO3; \
hide nonbonded; hide cartoon, not {prot0}; util.cbag; center {prot0}")

#less soluble but still MOAD-invalid substances
print("hide sticks, resn PEG+PGE+BU3")

#more important ions
print("hide spheres, resn CA+NI+IOD")
