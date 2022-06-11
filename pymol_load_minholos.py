import sys
import string
import numpy as np

upperpath = "/Users/jonathanborowsky/mount"#"/home/jonathanb/mount"
directory = f"{upperpath}/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/moad_negatives/iofiles"

serial = 4
serial_2 = 2
centroid='1YDA'

alphabet = list(string.ascii_uppercase)

#loaddir = savedir = f"/project/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/moad_negatives/iofiles/processed-blast-v{serial}" #100pct-processed-blast"
input = np.load(f"{directory}/processed-blast-v{serial}/{centroid.lower()}_allresi_ligand_subset_v{serial_2}.npy", allow_pickle=True)

#list(np.load(f"{directory}/processed-blast-v{serial}/{centroid}-v{serial}-holo-monomer.npy"))\
#+list(np.load(f"{directory}/processed-blast-v{serial}/{centroid}-v{serial}-holo-multimer.npy"))

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
    if i[0] != centroid: # and i not in ["1CNW", "1CNX","1CNY"]:

        for j in alphabet:
            try:
                cmd.fetch(i[0]+j)
                cmd.align(i[0]+j, centroid+"A")

            except:
                break

        resquery = [f"(resn {k[0]} and chain {k[1]} and resi {k[2]})" for k in i[1]]

        if len(resquery) != 0:
            #if i[0] == "6S9Z":
            #    print(resquery)

            cmd.hide("sticks", f"{i[0]}*")
            cmd.show("sticks", f"{i[0]}* and {' or '.join(resquery)}")




#graphics settings and hiding highly soluble ligands
#print(f"hide sticks, resn GOL+EDO+FMT+DMS+ACT or elem H; hide spheres, resn NA+CL+K+PO4+SO4+NO3; \
#hide nonbonded; hide cartoon, not {centroid}; util.cbag; center {centroid}")

cmd.hide("cartoon", f"not {centroid}A")

cmd.center(centroid+"A")

#util.cbag(centroid+"A")
cmd.color("grey80")
cmd.color("red", f"{centroid}A and not resi 24+25+26+28+29+30+31+32+33+41+42+43+44+45+46+47+48+49+50+56+59+61+74+75+77+78+79+80+81+82+83+84+85+86+87+88+89+90+93+97+99+100+101+104+105+107+115+116+117+118+120+122+124+125+134+140+142+144+145+146+148+150+151+152+153+154+157+160+163+166+167+169+171+172+173+174+175+176+177+184+185+186+187+188+189+190+191+193+196+210+212+216+217+218+219+220+221+222+223+224+226+227+230+238+240+241+248+249+250+251+252+253+254+255+256+257+259")
util.cbac("not polymer.protein")

cmd.hide("sticks", "elem H")
cmd.hide("spheres", "not resn Zn")
cmd.hide("nonbonded")

cmd.set_view("\
     0.329743475,   -0.410843283,    0.849979579,\
     0.114306711,    0.911082625,    0.396031797,\
    -0.937119544,   -0.033424273,    0.347387850,\
     0.000151612,    0.000083841, -129.259933472,\
    -9.353231430,   -0.943585992,   14.482344627,\
    88.654739380,  169.870666504,  -20.000000000 ")

cmd.ray()

#cmd.set_view("\
#     0.526804566,    0.000159980,    0.849979579,\
#    -0.639185250,    0.659226179,    0.396031797,\
#    -0.560277700,   -0.751930475,    0.347387850,\
#    -0.000020377,    0.000261210, -127.018424988,\
#    -8.480367661,   -4.574350834,   16.481094360,\
#    86.411666870,  167.627578735,  -20.000000000")



#less soluble but still MOAD-invalid substances
#print("hide sticks, resn PEG+PGE+BU3")

#more (potentially) important ions
#print("hide spheres, resn CA+NI+IOD")
