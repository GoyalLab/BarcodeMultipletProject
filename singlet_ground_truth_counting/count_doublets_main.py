from count_doublets_utils import *
from scipy import io
import pandas as pd
import os
from tqdm import tqdm
import numpy as np
import glob
from pathlib import Path
import pathlib

################################
#                              #
#       Input Parameters       #
#                              #
################################

data_root = "/projects/p31666/zzhang/doublet-bchmk/data/fatemap_data"
keyword = "stepFourStarcodeShavedReads50"
pattern = '**/{}*'.format(keyword)
overwrite = True


def grab_input_files(root_dir, file_format):
    grabbed_files = []
    for filepath in pathlib.Path(root_dir).glob(pattern):
        grabbed_files.append(str(filepath.absolute()))
    return grabbed_files


if __name__ == "__main__":
    target_files = grab_input_files(data_root, pattern)
    for cur_file in target_files:
        # dataset name always at 2 directories above
        # i.e. ../../
        cur_dataset_ID = cur_file.split('/')[-3]
        cur_out_prefix = cur_file.split('/')[:-1]
        cur_out_prefix.append(cur_dataset_ID)
        cur_out_prefix = "/".join(cur_out_prefix)

        print("INFO: Work on on {}".format(cur_dataset_ID))
        count_doublets(cur_file, cur_out_prefix, overwrite=overwrite)

