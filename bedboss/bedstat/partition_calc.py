import pybedtools


def calc_partitions(bed_id: str, bed_file: str, genome: str, out_folder: str) -> dict:

    bed_file = pybedtools.BedTool(bed_file)

    parts_hg38 = (
        "/home/bnt4me/virginia/repos/bedboss/scripts/stats/partitions/hg38_parts/"
    )
    exon = pybedtools.BedTool(f"{parts_hg38}/hg_38.exon.csv")
    intron = pybedtools.BedTool(f"{parts_hg38}hg_38.intron.csv")
    prm_cor = pybedtools.BedTool(f"{parts_hg38}hg_38.promoterCore.csv")
    prm_proc = pybedtools.BedTool(f"{parts_hg38}hg_38.promoterProx.csv")
    prm_5 = pybedtools.BedTool(f"{parts_hg38}hg_38.fiveUTR.csv")
    prm_3 = pybedtools.BedTool(f"{parts_hg38}hg_38.threeUTR.csv")

    f = 0.02
    r = True
    u = True

    # Find the overlaps
    overlap_exon = bed_file.intersect(exon, f=f, r=r, u=u)
    overlap_intron = bed_file.intersect(intron, f=f, r=r, u=u)
    overlap_prm_cor = bed_file.intersect(prm_cor, f=f, r=r, u=u)
    overlap_prm_proc = bed_file.intersect(prm_proc, f=f, r=r, u=u)
    overlap_5 = bed_file.intersect(prm_5, f=f, r=r, u=u)
    overlap_3 = bed_file.intersect(prm_3, f=f, r=r, u=u)

    overlap_exon_len = len(overlap_exon)
    overlap_intron_len = len(overlap_intron)
    overlap_prm_cor_len = len(overlap_prm_cor)
    overlap_prm_proc_len = len(overlap_prm_proc)
    overlap_5_len = len(overlap_5)
    overlap_3_len = len(overlap_3)

    total_number_of_regions = len(bed_file)

    # Overlap percentage
    exon_percentage = overlap_exon_len / total_number_of_regions
    intron_percentage = overlap_intron_len / total_number_of_regions
    prm_cor_percentage = overlap_prm_cor_len / total_number_of_regions
    prm_proc_percentage = overlap_prm_proc_len / total_number_of_regions
    prm_5_percentage = overlap_5_len / total_number_of_regions
    prm_3_percentage = overlap_3_len / total_number_of_regions

    # TODO: calculate log10(Obs/Exp) for each partition
    # TODO: how to do it?

    combined = overlap_exon.cat(
        overlap_intron,
        overlap_prm_cor,
        overlap_prm_proc,
        overlap_5,
        overlap_3,
        postmerge=False,
    )
    intragen = len(bed_file) - len(combined)
    intragen_percentage = intragen / total_number_of_regions

    import matplotlib.pyplot as plt

    # Sample data (replace these with your actual data)
    categories = [
        "exon",
        "fiveUTR",
        "intergenic",
        "intron",
        "promoterCore",
        "promoterProx",
        "threeUTR",
    ]
    frequencies = [
        exon_percentage,
        prm_5_percentage,
        intragen_percentage,
        intron_percentage,
        prm_cor_percentage,
        prm_proc_percentage,
        prm_3_percentage,
    ]

    # Create the bar plot
    plt.figure(figsize=(8, 6))  # Adjust figure size as needed
    plt.bar(categories, frequencies, color="gray")

    # Add labels and title
    plt.xlabel("Genomic partition", fontsize=12)
    plt.ylabel("Frequency (%)", fontsize=12)
    plt.title("Distribution across genomic partitions", fontsize=14)

    # Rotate x-axis labels for better readability
    plt.xticks(rotation=90)

    # Display the plot
    plt.tight_layout()
    plt.savefig(f"{out_folder}/{bed_id}_partitions.png")

    # return {
    #     "exon": exon_percentage,
    #     "fiveUTR": prm_5_percentage,
    #     "intergenic": intron_percentage,
    #     "intron": prm_cor_percentage,
    #     "promoterCore": prm_proc_percentage,
    #     "promoterProx": prm_3_percentage,
    # }
    exit(0)


# 2 -> Calculate gc content

# 3 -> Distance between neighbor regions

# 4 ->
