# Run classify on a selection of gappedPeaks
import os
from bedboss.bedclassifier.bedclassifier import get_bed_type

DATA_PATH = "/home/drc/test/gappedPeaks/"

all_files = []
for root, _, files in os.walk(DATA_PATH):
    for file in files:
        full_path = os.path.join(root, file)
        all_files.append(full_path)

count_neg = 0
count_pos = 0
for file in all_files:
    result = get_bed_type(file)
    if result[1] != 'gappedpeak':
        # print(f"This one is not classified as broadpeak: {file_path}")
        count_neg += 1
        print(f"{result}: {file}")
    else:

        # print("FOUND broadpeak")
        count_pos += 1
    #print(f"{result}: {file}")
print(f"gappedPeak: {count_pos} \nnot gappedPeak:{count_neg}")

# One-off testing

#result = get_bed_type("/home/drc/test/test_gappedPeaks_geofetched/data/GSE192575/GSM5751922_ATAC_resis_1_peaks.gappedPeak.gz")
#result = get_bed_type("/home/drc/test/test_gappedPeaks_geofetched/data/GSE192575/GSM5751923_ATAC_resis_2_peaks.gappedPeak.gz")

#print(result)