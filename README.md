# 01. Fatemap Data Structure
```bash
├── FM01
│   ├── 10X
│   │   ├── sample1
│   │   │   │   ├── barcodes.tsv.gz
│   │   │   │   ├── features.tsv.gz
│   │   │   │   ├── matrix.mtx.gz
│   │   ├── sample2
│   │   ├──  ...
│   │   ├── sampleX
│   ├── fatemapID
│   │   ├── FM01_good_data.tsv
│   │   ├── FM01_multiplets.txt
│   │   ├── FM01_singlets.txt
│   │   ├── stepFourStarcodeShavedReads50.txt
├── FM02
├── FM03
├── ...
└── FM06
```
---------

# 02. Ground Truth Singlets Counting
## 02.1 Caveats with Preprocessing
The counting of doublets uses only the `cellID`, `BC50StarcodeD8` and `sampleNum` fields from the original files shared
from OneDrive. For dataset **FM02**, the file doesn't have header so it's added with the following command.
> `sed -i '1i cellID\tUMI\tOriginalBarcode\tBC50StarcodeD8\tBC40StarcodeD8\tBC30StarcodeD8\tBC30StarcodeD6\tSampleNum' stepThreeStarcodeShavedReads.txt`

In addition, the file appears to have a delimiter of `space`. This is converted to `tab` with the following command.
> `sed -i 's/ /\t/g' stepThreeStarcodeShavedReads.txt`

Lastly, the three columns used are separately stored with the following command. This applies to all dataset except for **FM01**.
>`cut -f 1,2,4,8 stepThreeStarcodeShavedReads.txt | awk '!seen[$0]++' > stepFourStarcodeShavedReads50.txt`

## 02.2 Singlet Identification
1. Input file description
   1. Four columns
      1. cellID
      2. UMI
      3. BC50StarcodeD8: i.e. fatemap barcode
      4. sampleNum
2. Select all the good data for each sample within each dataset
   1. Calculate a minimum UNI cutoff, which defaults to `3/4e5 * num_cellID`
   2. Perform the `groupby` command to count for each sample, how many `UMIs` are associated with a unique `cellID` and `fatemap barcode` combination
   3. Cells below the cutoff calculated in step 1 are dropped, and those kept are considered `good_data`
3. Pre-filtering data processing
   1. Construct a dictionary called `cellID2fatemap_dict` with `good_data` that documents how many fatemap barcode are associated with each `cellID`
   2. Construct a dictoinary called `fatemapID_dict` with multiplets that account for the number of times a combination of multiple fatemap barcodes are associated with `cellIDs`
   3. Identify `multilane_barcodes` that are seen within more than 1 sample
4. Filtering for true singlets
   1. If a `cellID` is seen in `multilane_barcodes`, this cellID is discarded. 
   2. If the `cellID` is associated with only 1 `fatemap barcode`, it is a singlet. 
   3. If multiple `cellIDs` are associated with the same combination of `fatemap barcode`, they are all singlets. 
   4. If a `cellID` is associated with multiple `fatemap barcodes`, but one of them is dominant over the others in terms of `nUMI`, then keep the dominant one as a singlet. 

## 02.3 Running Singlet Counting
1. Scripts Used
   1. `singlet_ground_truth_counting/count_doublets_utils.py`
   2. `singlet_ground_truth_counting/count_doublets_main.py`
   3. `singlet_ground_truth_counting/run_count_doublets.sh`
2. Expected output
   1. `data/fatemap_data/{Dataset_ID}/fatemapID/{Dataset_ID}_singlets.txt`
   2. `data/fatemap_data/{Dataset_ID}/fatemapID/{Dataset_ID}_multiplets.txt`
   3. `data/fatemap_data/{Dataset_ID}/fatemapID/{Dataset_ID}_{Sample_Num}_singlet_pairs.csv`
   4. `data/fatemap_data/{Dataset_ID}/fatemapID/{Dataset_ID}_singlets_stats.csv`

## 02.4 Sanity Test
1. Test if good data is generated correctly
   1. Good data is the nUMI count table for each `cellID + fatemap_barcode + sample_num combo`
   2. Each `cellID + fatemap_barcode + sample_num` combo is subset from the raw data
   3. The umi count is determined by evaluating how many rows this subset-ed dataframe has
   4. Check if this umi count is less then this sample's umi cutoff
   5. Cehck if this umi count is the same as the `nUMI` entry in `good_data`
2. Test if cellID to fatemap_barcode dictionary is generated correctly 
   1. Starting at the `good_data`
   2. Iterate through each of the samples, then each of the sample's cellIDs
   3. Check if the `good_data`, for that `cellID` and `sampleNum` combo, has the same number of fatemap_barcodes as the the number in the dictionary

# 03. Running Benchmark


# Dataset Usage Spreadsheet

| Dataset ID    | Number of Samples | Type                                      | Total Cell Count | Total Singlets | Singlet Proportion | Link                                                                           |
|---------------|-------------------|-------------------------------------------|------------------|----------------|--------------------|--------------------------------------------------------------------------------|
| FM01          | 4                 | Melanoma Patient Resistant Cell line      | 32189            | 20505          | 0.637              | /projects/p31666/zzhang/doublet-bchmk/data/fatemap_data/FM01                   |
| FM02          | 4                 | Melanoma Patient Resistant Cell line      | 39884            | 23886          | 0.598              | /projects/p31666/zzhang/doublet-bchmk/data/fatemap_data/FM02                   |
| FM03          | 2                 | Melanoma Patient Resistant Cell line      | 16961            | 8842           | 0.521              | /projects/p31666/zzhang/doublet-bchmk/data/fatemap_data/FM03                   |
| FM04          | 2                 | Breast Cancer Patient Resistant Cell line | 17069            | 10145          | 0.594              | /projects/p31666/zzhang/doublet-bchmk/data/fatemap_data/FM04                   |
| FM05          | 2                 | Melanoma Patient Resistant Cell line      | 26728            | 13476          | 0.504              | /projects/p31666/zzhang/doublet-bchmk/data/fatemap_data/FM05                   |
| FM06          | 2                 | Melanoma Patient Untreated Cell line      | 18106            | 13037          | 0.720              | /projects/p31666/zzhang/doublet-bchmk/data/fatemap_data/FM06                   |
| FM08          | 2                 | Primary Melanocytes                       | 6805             | 4438           | 0.652              | /projects/p31666/zzhang/doublet-bchmk/data/fatemap_data/FM08                   |
| non_cancer    | 2                 | HIPS differentiation                      | 10973            | 3457           | 0.315              | /projects/p31666/zzhang/doublet-bchmk/data/fatemap_data/non_cancer             |
| hm-12k        | NA                | Synthetic dataset                         | 12820            | 12090          | 0.943              | /projects/p31666/zzhang/doublet-bchmk/data/bchmk_paper_real_data/hm-12k        |
| J293t         | NA                | Potentially faulty annotation             | 500              | 458            | 0.916              | /projects/p31666/zzhang/doublet-bchmk/data/bchmk_paper_real_data/J293t         |
| pbmc-2ctrl-dm | NA                | Patient PBMCs                             | 13913            | 12315          | 0.885              | /projects/p31666/zzhang/doublet-bchmk/data/bchmk_paper_real_data/pbmc-2ctrl-dm |

---------
# Archive

## Demuxlet
Bam files are QC-ed according to GATK's instructions: https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery

Reference genome also requires a dict file, which can be created with:
>`gatk CreateSequenceDictionary -R genome.fa -O genome.dict`

The command ran is: 
>`gatk MarkDuplicatesSpark -I input.bam -O marked_duplicates.bam -M marked_dup_metrics.txt`

The UCSC fasta file has **chr** in front of the chromosomes 
>- `cat hg19.fa | sed 's/>chr/>/g' > hg19_new.fa`
>- `samtools faidx hg19_updated.fa`
>- `gatk CreateSequenceDictionary -R hg19_updated.fa`
