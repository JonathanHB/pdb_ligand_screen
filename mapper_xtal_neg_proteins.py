# IMPORTS
import numpy as np
from glob import glob

# filepaths
upperpath = "/home/jonathanb/mount" #"/Users/ameller/projects/mnt/bowmore" #"X:/project"
structurepath = f"{upperpath}/bowmore/ameller/projects/pocket_prediction/val_structures"

# val_apo_path = f"{upperpath}/bowmore/ameller/gvp/data/validation_apo_ids.npy"
val_apo_path = f'{upperpath}/bowmore/ameller/projects/pocket_prediction/data/val_apo_ids.npy'

# val_apo_path = f"{upperpath}/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/new_pockets/labels/validation_apo_ids_all.npy"

predlabelpath = f"{upperpath}/bowmanlab/ameller/gvp/task2"

# parameters
binarize = False
binthresh = 0.5

# list of training schemes/network versions
label_trainschemes = [\
"train-with-4-residue-batches-no-balancing-intermediates-in-training/net_8-50_1-32_16-100_dr_0.1_nl_4_hd_100_lr_2e-05_b4resis_b1proteins_20epoch_feat_method_gp-to-nearest-resi-procedure_rank_7_stride_1_window_40_pos_30",
"train-with-4-residue-batches-constant-size-balanced-640-resi-draws-intermediates-in-training/net_8-50_1-32_16-100_dr_0.1_nl_4_hd_100_lr_2e-05_b4resis_b4proteins_30epoch_feat_method_nearby-pv-procedure_rank_7_stride_1_window_40_pos_116",
"train-with-4-residue-batches-constant-size-balanced-640-resi-draws-intermediates-in-training/net_8-50_1-32_16-100_dr_0.1_nl_4_hd_100_lr_2e-05_b4resis_b4proteins_30epoch_feat_method_nearby-pv-procedure_rank_7_stride_1_window_40_pos_87",
"cryptic-labels-train-with-1-protein-batches-no-balancing/net_8-50_1-32_16-100_dr_0.1_nl_4_hd_100_lr_2e-05_b1proteins_20epoch_feat_method_cryptic_labels_closed_structures",
"train-with-4-residue-batches-no-balancing-intermediates-in-training/net_8-50_1-32_16-100_dr_0.1_nl_4_hd_100_lr_2e-05_b4resis_b1proteins_20epoch_feat_method_gp-to-nearest-resi-procedure_rank_7_stride_1_window_40_pos_20",
"fpocket-drug-scores-unbinarized-train-with-4-residue-batches/net_8-50_1-32_16-100_dr_0.1_nl_4_hd_100_lr_2e-05_b4resis_b1proteins_20epoch_feat_method_fpocket_drug_scores_max_window_40_cutoff_0.3_stride_25",
"fpocket-drug-scores-unbinarized-train-with-4-residue-batches/net_8-50_1-32_16-100_dr_0.1_nl_4_hd_100_lr_2e-05_b4resis_b1proteins_20epoch_feat_method_fpocket_drug_scores_difference_window_40_cutoff_0.3_stride_25",
f"train-with-1-protein-batches-no-balancing-intermediates-in-training/net_8-50_1-32_16-100_dr_0.1_nl_4_hd_100_lr_2e-05_b1proteins_30epoch_feat_method_nearby-pv-procedure_rank_7_stride_1_window_40_pos_87",\
f"train-with-1-protein-batches-no-balancing-intermediates-in-training/net_8-50_1-32_16-100_dr_0.1_nl_4_hd_100_lr_2e-05_b1proteins_30epoch_feat_method_nearby-pv-procedure_rank_7_stride_1_window_40_pos_116",\
f"train-with-1-protein-batches-no-balancing-intermediates-in-training/net_8-50_1-32_16-100_dr_0.1_nl_4_hd_100_lr_2e-05_b1proteins_30epoch_feat_method_nearby-pv-procedure_rank_7_stride_1_window_40_pos_145",\
f"train-with-1-protein-batches-undersample-intermediates-in-training/net_8-50_1-32_16-100_dr_0.1_nl_4_hd_100_lr_2e-05_b1proteins_30epoch_feat_method_nearby-pv-procedure_rank_7_stride_1_window_40_pos_87",\
f"train-with-1-protein-batches-undersample-intermediates-in-training/net_8-50_1-32_16-100_dr_0.1_nl_4_hd_100_lr_2e-05_b1proteins_30epoch_feat_method_nearby-pv-procedure_rank_7_stride_1_window_40_pos_116",\
f"train-with-1-protein-batches-undersample-intermediates-in-training/net_8-50_1-32_16-100_dr_0.1_nl_4_hd_100_lr_2e-05_b1proteins_30epoch_feat_method_nearby-pv-procedure_rank_7_stride_1_window_40_pos_145"]

# global load commands

  # proteins in the validation set
val_apo_ids = np.load(val_apo_path)

# true labels
positive_examples = np.load(f"{upperpath}/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/all-GVP-project-ligand-resis.npy", allow_pickle=True)
truelabels_all = np.load(f'/{upperpath}/bowmore/ameller/projects/pocket_prediction/data/val_label_dictionary.npy',
                         allow_pickle=True).item()
truelabels_all_apo = [p for p in truelabels_all.keys()]
val_set_apo_ids_with_chainids = np.load(f'/{upperpath}/bowmore/ameller/projects/'
                                        'pocket_prediction/data/val_apo_ids_with_chainids.npy')

  # predicted labels
#identify the best epoch for each training parameter set; deprecated as Artur's code now takes care of this
#best_epochs = []

#for ts in label_trainschemes:

#    best_epoch = -1
#    best_pr_auc = 0

#    for epoch in range(0,30):
#        predprauc_mean = np.mean(np.load(f"{predlabelpath}/{ts}/val_protein_pr_aucs_{epoch}.npy"))
#        if predprauc_mean > best_pr_auc:
#            best_pr_auc = predprauc_mean
#            best_epoch = epoch

#    best_epochs.append(best_epoch)

#load the predicted labels and pr aucs for the epoch with the best average pr auc
#predlabels_all = [np.load(f"{predlabelpath}/{i}/val_y_pred_{best_epochs[x]}.npy") for x, i in enumerate(label_trainschemes)]
#predpraucs_all = [np.load(f"{predlabelpath}/{i}/val_protein_pr_aucs_{best_epochs[x]}.npy") for x, i in enumerate(label_trainschemes)]

predlabels_all = [np.load(f"{predlabelpath}/{i}/new_val_set_y_pred_best_epoch.npy") for x, i in enumerate(label_trainschemes)]
predpraucs_all = [[np.load(f"{predlabelpath}/{i}/val_new_auc_{epoch}.npy")]
                  for x, i in enumerate(label_trainschemes)
                  for epoch in range(len(glob(f"{predlabelpath}/{i}/val_new_auc_*.npy")))]

# graphics settings; presently unused

def set_graphics():
    cmd.do("bg grey50")
    cmd.do("set orthoscopic,1")
    cmd.do("set depth_cue,0")
    cmd.do("set auto_zoom, off")
    cmd.do("set sphere_quality, 5")
    cmd.do("set opaque_background, off")
    cmd.do("set ray_opaque_background, off")
    cmd.do("set antialias, 2")
    cmd.do("set ray_trace_mode, 1")
    cmd.do("set ray_trace_color, black")

#summary information to print when creating the macro

print(val_apo_ids)
print("Input:\nt2mapper ['t' to color labels on a [0, 1] scale, anything else to color them on a [0, max_label] scale], [index of the training scheme of interest; see above], [protein apo pdb id, or index in validation set protein list]")

#define method
@cmd.extend
def t2mapper(abscale, label_ts_ind, protin, cluster_labels):

    #identify and process inputs
    if abscale == "t":
        absolute_scale = True
    elif abscale == "f":
        absolute_scale = False
    else:
        absolute_scale = True
        print("invalid scale argument, enter 't' or 'f', defaulting to True")

    #identify and process inputs
    if cluster_labels == "t":
        cluster_labels = True
    elif cluster_labels == "f":
        cluster_labels = False
    else:
        cluster_labels = True
        print("invalid scale argument, enter 't' or 'f', defaulting to True")

    try:
        protx = int(protin)
        prot = val_apo_ids[protx]
        prot_with_chainid = val_set_apo_ids_with_chainids[protx]
    except ValueError:
        prot = protin
        protx = val_apo_ids.index(prot.lower())
        prot_with_chainid = val_set_apo_ids_with_chainids[protx]

    #get protein information
    truelabels = truelabels_all[prot]

    # determine if this is a cryptic pocket example
    if 1 in truelabels:
        load_holo = True
        positive_example_apoids = [i[0][0].lower() for i in positive_examples[0]]
        prot_index = positive_example_apoids.index(prot)
        prot_info = positive_examples[0][prot_index][0]
        print(prot_info)
    else:
        load_holo = False


    # prot_index = truelabels_all_apo.index(prot)
    # prot_info = truelabels_all[0][prot_index][0]

    # truelabels = truelabels_all[0][prot_index][4]

    #print summary information
    print(label_trainschemes[int(label_ts_ind)])
    # if protx < len(predpraucs_all):
    #     print(f"PR-AUC = {predpraucs_all[int(label_ts_ind)][protx]}")
    # print(prot_info)

    #delete all and load apo and holo structures
    cmd.delete("all")

    aponame = prot_with_chainid
    apo_path = f"{structurepath}/{aponame}_clean_h.pdb"
    if os.path.exists(apo_path):
        cmd.load(apo_path)
    else:
        apo_path = f"{structurepath}/{aponame.upper()}_clean_h.pdb"
        cmd.load(apo_path)

    if load_holo:
        holoname = prot_info[2]+prot_info[3]+"_clean_ligand_h"
        cmd.load(f"{structurepath}/{holoname}.pdb")
        cmd.align(aponame, holoname, cycles = 0)
        #misc graphics settings
        util.cbac(holoname)
        cmd.hide("cartoon", holoname)

    #show predicted labels by color

    #get initial pdb residue number
    init_resi = -999
    for line in open(apo_path):
        if init_resi == -999 and line[0:4] == "ATOM":
            init_resi = int(line[22:26])
            break

    #show true labels as sticks
    for i in np.where(truelabels == 1)[0]:
        cmd.show("sticks", f"{aponame} and resi {i + init_resi} and not elem H")

    #show true labels as sticks
    for i in np.where(truelabels == 0)[0]:
        cmd.show("lines", f"{aponame} and resi {i + init_resi} and not elem H")


    if not cluster_labels:
        predlabels = predlabels_all[int(label_ts_ind)][protx]

            #set label mapping to maximum color value
        if absolute_scale:
            scalemax = 1
        else:
            scalemax = max(predlabels)

        #color residues
        for resi, label in enumerate(predlabels):

            if binarize:
                if label > binthresh:
                    label = 1
                else:
                    label = 0

            grp = f"{aponame} and resi {int(resi)+init_resi}"
            s = f"b={label}"
            cmd.alter(grp, s)
    else:

        cluster_labels = np.load(f'../data/{aponame[0:4].lower() + aponame[4]}-resi-labels.npy')

        for resi, label in enumerate(cluster_labels):
            grp = f"{aponame} and resi {int(resi)+init_resi}"
            s = f"b={label}"
            cmd.alter(grp, s)

        scalemax = np.max(cluster_labels)

    # cmd.spectrum("b", "blue_white_red", aponame, 0, scalemax)
    cmd.spectrum("b", "gray_red", aponame, 0, scalemax)
