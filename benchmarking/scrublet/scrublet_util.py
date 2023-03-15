import scrublet as scr
import scipy.io
import numpy as np
import os
import sklearn.metrics as metrics
import gzip
import pandas as pd
from scipy.sparse import vstack
import random

def run_scrublet_bchmk_data(data_dir):
    # ValueError: n_components=30 must be between 1 and min(n_samples, n_features)=9 with svd_solver='arpack'
    cur_scrub, doublet_labels = load_bchmk_into_scrublet(data_dir)
    doublet_scores, predicted_doublets = cur_scrub.scrub_doublets()
    fpr, tpr, threshold = metrics.roc_curve(y_true=doublet_labels, y_score=doublet_scores, pos_label="doublet")
    roc_auc = metrics.auc(fpr, tpr)
    precision, recall, thresholds = metrics.precision_recall_curve(doublet_labels, doublet_scores, pos_label="doublet")
    auc_precision_recall = metrics.auc(recall, precision)
    out_df = pd.DataFrame({
        "roc": [roc_auc],
        "pr": [auc_precision_recall]
    })
    return out_df


def run_bchmk_on_bchmk_data(data_dir, data_sample_ID, num_to_test):
    df_ls = []
    for idx in range(num_to_test):
        cur_res = run_scrublet_bchmk_data(data_dir)
        df_ls.append(cur_res)

    all_stats = pd.concat(df_ls)
    all_stats['ID'] = data_sample_ID
    return all_stats


def load_bchmk_into_scrublet(data_dir):
    dataset_id = data_dir.split("/")[-1]
    count_file = data_dir + "/" + dataset_id + ".mtx"
    labels_file = data_dir + "/" + dataset_id + "_labels.csv"
    labels = list(pd.read_csv(labels_file, header=None)[0])
    count_mtx = scipy.io.mmread(count_file).T
    cur_scrub = scr.Scrublet(counts_matrix=count_mtx)
    return cur_scrub, labels


def load_fatemap_into_scrublet(data_dir, cell_labels_prefix, exp_dbl, doublets_pct=0.08):
    cell_id_file = data_dir + "/barcodes.tsv.gz"
    gene_id_file = data_dir + "/features.tsv.gz"
    mtx = data_dir + "/matrix.mtx.gz"
    singlets_file = cell_labels_prefix + "_singlets.txt"
    singlets_df = pd.read_csv(singlets_file, header=None)
    singlets_df[0] = singlets_df[0] + "-1"
    singlets_pair_for_doublet_simulation_file = data_dir + "/corrected_singlet_pairs.csv"
    singlets_pairs_df = pd.read_csv(singlets_pair_for_doublet_simulation_file, header=None, skipinitialspace=True)
    singlets_pairs_df[0] = singlets_pairs_df[0]

    with gzip.open(cell_id_file) as fp:
        cell_ids = np.array([x.decode('UTF-8').strip("\n") for x in fp.readlines()])

    singlets_in_seurat = [x for x in singlets_df[0] if x in cell_ids]
    idx_to_keep = [np.where(cell_ids == x)[0][0] for x in singlets_in_seurat]
    # cell id by gene id matrix
    counts_matrix = scipy.io.mmread(mtx).T.tocsr()[idx_to_keep, :]
    #     print("Total singlets: {}".format(str(counts_matrix.shape[0])))

    with gzip.open(gene_id_file) as fp:
        gene_IDs = np.array([x.decode('UTF-8').strip("\n").split("\t") for x in fp.readlines()])
    gene_IDs = gene_IDs[:, 1]

    # take average of expression to produce doublets
    singlets_pt1_idx = [singlets_in_seurat.index(x) for x in singlets_pairs_df[0]]
    singlets_pt2_idx = [singlets_in_seurat.index(x) for x in singlets_pairs_df[1]]
    avg_assay = (counts_matrix[singlets_pt1_idx, :] + counts_matrix[singlets_pt2_idx, :]) / 2
    doublet_id = pd.DataFrame(singlets_pairs_df[0] + "--" + singlets_pairs_df[1])
    doublet_id = doublet_id[0].str.strip()
    num_doublets = int((len(singlets_in_seurat) * doublets_pct) / (1 - doublets_pct))
    avg_assay = avg_assay[:num_doublets, :]
    doublet_labels = ["doublet" for i in range(num_doublets)]
    singlet_labels = ["singlet" for i in range(len(singlets_in_seurat))]
    doublet_labels.extend(singlet_labels)
    doublet_assay = vstack((avg_assay, counts_matrix))
    estimated_multiplet_rate = len(doublet_labels) / (len(doublet_labels) + len(singlet_labels))
    #     print("Total doublets: {}".format(str(avg_assay.shape[0])))
    #     print("Total dataset: {}".format(str(doublet_assay.shape[0])))
    # random initiation
    cur_random_state = random.randint(1, int(1e5))
    print("INFO: Random stat is {}".format(str(cur_random_state)))
    cur_scrub = scr.Scrublet(doublet_assay, expected_doublet_rate=exp_dbl, random_state = cur_random_state)

    # Unlike Seurat, we must return the labels separately instead of as a part of metadata
    return cur_scrub, doublet_labels


def run_scrublet(data_dir, cell_labels_prefix, exp_dbl, act_dbl):
    # ValueError: n_components=30 must be between 1 and min(n_samples, n_features)=9 with svd_solver='arpack'
    cur_scrub, doublet_labels = load_fatemap_into_scrublet(data_dir, cell_labels_prefix,
                                                           exp_dbl=exp_dbl, doublets_pct=act_dbl)
    doublet_scores, predicted_doublets = cur_scrub.scrub_doublets()
    fpr, tpr, threshold = metrics.roc_curve(y_true=doublet_labels, y_score=doublet_scores, pos_label="doublet")
    roc_auc = metrics.auc(fpr, tpr)
    precision, recall, thresholds = metrics.precision_recall_curve(doublet_labels, doublet_scores, pos_label="doublet")
    auc_precision_recall = metrics.auc(recall, precision)
    out_df = pd.DataFrame({
        "roc": [roc_auc],
        "pr": [auc_precision_recall]
    })
    return out_df


def run_bchmk(data_dir, cell_labels_prefix, data_sample_ID, num_to_test, exp_dbl, act_dbl):
    df_ls = []
    for idx in range(num_to_test):
        cur_res = run_scrublet(data_dir, cell_labels_prefix, exp_dbl, act_dbl)
        df_ls.append(cur_res)

    all_stats = pd.concat(df_ls)
    all_stats['ID'] = data_sample_ID
    return all_stats


# need to reload scrublet object everytime it is ran because it does not reset
def run_scrublet_on_dataset(data_dir, METHOD_ID, exp_dbl, act_dbl, num_to_test=5):
    print("Working on dataset: {}".format(data_dir))
    dir_10X = data_dir + "/10X"
    dataset_id = data_dir.split("/")[-1]
    sample_dirs = [dir_10X + "/" + x for x in os.listdir(dir_10X)]
    out_df = pd.DataFrame(columns=["roc", "pr"])
    labels_prefix = "{}/fatemapID/{}".format(data_dir, dataset_id)
    for cur_sample_dir in sample_dirs:
        #         cur_scrub_object, cur_labels = load_fatemap_into_scrublet(data_dir, cell_labels_prefix)
        cur_sample_id = cur_sample_dir.split("/")[-1]
        cur_data_sample_id = "{}_{}".format(dataset_id, cur_sample_id)
        cur_df = run_bchmk(cur_sample_dir, labels_prefix, cur_data_sample_id, num_to_test, exp_dbl, act_dbl)
        out_df = pd.concat([out_df, cur_df])
    return out_df
