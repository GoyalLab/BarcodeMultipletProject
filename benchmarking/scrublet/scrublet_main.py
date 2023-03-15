import matplotlib.pyplot as plt
from pathlib import Path
import pickle
from scrublet_util import *
import os
################################
#                              #
#     Artistic Parameters      #
#                              #
################################
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42

################################
#                              #
#       Input Parameters       #
#                              #
################################
root_dir = "/projects/p31666/zzhang/doublet-bchmk/data/fatemap_data"
out_dir = "/projects/p31666/zzhang/doublet-bchmk/output/run3"
METHOD_ID = "scrublet"
root_bchmk_dataset = "/projects/p31666/zzhang/doublet-bchmk/data/bchmk_paper_real_data"
num_test = 2
doublet_rates_to_test = [0.05, 0.08, 0.1, 0.15, 0.2, 0.25]



################################
#                              #
#             Main             #
#                              #
################################

def main(root_dir, METHOD_ID, out_dir):
    out_dfs = []
    out_stats = []
    all_bchmk_datasets_dirs = [os.path.join(root_bchmk_dataset, x) for x in next(os.walk(root_bchmk_dataset))[1]]
    print("INFO: Benchmarking datasets from benchmarking paper!")
    for cur_dataset_dir in all_bchmk_datasets_dirs:
        dataset_id = cur_dataset_dir.split("/")[-1]
        cur_out_df = run_bchmk_on_bchmk_data(cur_dataset_dir, dataset_id, num_to_test=5)
        cur_out_df["dbl_exp"] = "NA"
        cur_out_df["dbl_act"] = "NA"
        out_dfs.append(cur_out_df)
    print("INFO: Finished benchmarking datasets from benchmarking paper!")
    all_datasets_dirs = [os.path.join(root_dir, x) for x in next(os.walk(root_dir))[1]]
    for actual_doublet_rate in doublet_rates_to_test:
        for expected_doublet_rate in doublet_rates_to_test:
            print("INFO: Currently testing expected doublet rate:{}, actual doublet rate:{}".
                  format(str(expected_doublet_rate), str(actual_doublet_rate)))
            for cur_dataset_dir in all_datasets_dirs:
                cur_out_df = run_scrublet_on_dataset(cur_dataset_dir, METHOD_ID, num_to_test=num_test,
                                                     exp_dbl=expected_doublet_rate, act_dbl=actual_doublet_rate)
                out_dfs.append(cur_out_df)
                cur_out_df["dbl_exp"] = expected_doublet_rate
                cur_out_df["dbl_exp"] = actual_doublet_rate

    all_detection_rates_file = os.path.join(out_dir, METHOD_ID, "stats", "all_detection_rates.tsv")
    all_df = pd.concat(out_dfs)
    all_df.to_csv(all_detection_rates_file, sep="\t", index=False)

    return True


if __name__ == "__main__":
    main(root_dir, METHOD_ID, out_dir)
