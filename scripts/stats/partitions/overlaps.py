import pybedtools


# Load your BED files
def main():
    # bed1 = pybedtools.BedTool('/home/bnt4me/virginia/repos/bedboss/scripts/stats/data/04dee6a07f7db768cba73481e7e8401a.bed.gz')
    bed1 = pybedtools.BedTool(
        "/home/bnt4me/virginia/repos/bedboss/scripts/stats/data/b72912f88165f2f252f273f4c3043fe3.bed.gz"
    )

    parts_hg38 = (
        "/home/bnt4me/virginia/repos/bedboss/scripts/stats/partitions/hg38_parts/"
    )
    exon = pybedtools.BedTool(f"{parts_hg38}/hg_38.exon.csv")
    intron = pybedtools.BedTool(f"{parts_hg38}hg_38.intron.csv")
    prm_cor = pybedtools.BedTool(f"{parts_hg38}hg_38.promoterCore.csv")
    prm_proc = pybedtools.BedTool(f"{parts_hg38}hg_38.promoterProx.csv")
    prm_5 = pybedtools.BedTool(f"{parts_hg38}hg_38.fiveUTR.csv")
    prm_3 = pybedtools.BedTool(f"{parts_hg38}hg_38.threeUTR.csv")

    f = 0.001
    r = True

    # Find the overlaps
    overlap_exon = bed1.intersect(exon, f=f, r=r, u=True)
    overlap_intron = bed1.intersect(intron, f=f, r=r, u=True)
    overlap_prm_cor = bed1.intersect(prm_cor, f=f, r=r, u=True)
    overlap_prm_proc = bed1.intersect(prm_proc, f=f, r=r, u=True)
    overlap_5 = bed1.intersect(prm_5, f=f, r=r, u=True)
    overlap_3 = bed1.intersect(prm_3, f=f, r=r, u=True)

    combined = overlap_exon.cat(
        overlap_intron,
        overlap_prm_cor,
        overlap_prm_proc,
        overlap_5,
        overlap_3,
        postmerge=False,
    )
    total_nujmber_of_regions = len(bed1)
    intragen = len(bed1) - len(combined)
    # print len for all overlaps:

    print("total: ", total_nujmber_of_regions)
    print("3: ", len(overlap_3))
    print("5: ", len(overlap_5))
    print("exon: ", len(overlap_exon))
    print("intron: ", len(overlap_intron))
    print("intergenic", intragen)
    print("prm_proc: ", len(overlap_prm_proc))
    print("prm_cor: ", len(overlap_prm_cor))

    # calculate the percentage of the overlaps
    print("3: ", len(overlap_3) / total_nujmber_of_regions)
    print("5: ", len(overlap_5) / total_nujmber_of_regions)
    print("exon: ", len(overlap_exon) / total_nujmber_of_regions)
    print("intron: ", len(overlap_intron) / total_nujmber_of_regions)
    print("intergenic", intragen / total_nujmber_of_regions)
    print("prm_proc: ", len(overlap_prm_proc) / total_nujmber_of_regions)
    print("prm_cor: ", len(overlap_prm_cor) / total_nujmber_of_regions)


if __name__ == "__main__":
    main()
