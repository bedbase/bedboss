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


def main2():
    import pybedtools
    import pandas as pd

    def calc_partitions(
        query, partition_list, remainder="intergenic", bp_proportion=False
    ):
        """
        Calculates the distribution of overlaps between query regions and an ordered list of genomic partitions.

        Parameters:
        - query (pybedtools.BedTool): The query regions to classify.
        - partition_list (dict): An ordered and named dictionary of genomic partitions (as BedTool objects).
                                 Keys are partition names (e.g., "promoter", "exon", etc.), values are BedTool objects.
        - remainder (str): Label for regions that do not overlap with any partition in partition_list. Defaults to "intergenic".
        - bp_proportion (bool): If True, calculates overlaps based on the number of base pairs. If False, assigns regions to partitions in priority order.

        Returns:
        - pd.DataFrame: Dataframe assigning each query region to a partition or the remainder category.
        """

        # Proportional overlap mode
        if bp_proportion:
            # Calculate base pair overlap with each partition
            total_overlap = {
                name: query.intersect(partition).total_coverage()
                for name, partition in partition_list.items()
            }

            # Calculate the base pairs that don't overlap any partition
            total_bp_query = sum([interval.length for interval in query])
            remainder_bp = total_bp_query - sum(total_overlap.values())

            # Prevent negative values in remainder
            if remainder_bp < 0:
                remainder_bp = 0

            # Create a DataFrame for proportional overlaps
            prop_partitions = pd.DataFrame(
                {
                    "partition": list(total_overlap.keys()) + [remainder],
                    "bp_overlap": list(total_overlap.values()) + [remainder_bp],
                }
            )
            prop_partitions["frequency"] = (
                prop_partitions["bp_overlap"] / prop_partitions["bp_overlap"].sum()
            )
            return prop_partitions

        # Priority overlap mode
        else:
            partition_assignments = []

            # Track assigned regions
            assigned = set()

            # Loop through partitions in priority order
            for partition_name, partition in partition_list.items():
                # Identify regions in query that have not yet been assigned
                unassigned = query.filter(lambda x: x.name not in assigned)

                # Find overlaps with the current partition
                overlaps = unassigned.intersect(partition, wa=True, u=True)

                # Assign regions to the current partition
                partition_assignments.extend(
                    [(interval.name, partition_name) for interval in overlaps]
                )
                assigned.update([interval.name for interval in overlaps])

            # Assign remaining regions to the remainder category
            remainder_regions = query.filter(lambda x: x.name not in assigned)
            partition_assignments.extend(
                [(interval.name, remainder) for interval in remainder_regions]
            )

            # Convert to a DataFrame and tally results
            df = pd.DataFrame(partition_assignments, columns=["region", "partition"])
            result = df["partition"].value_counts().reset_index()
            result.columns = ["partition", "Freq"]

            # Ensure all partitions, including the remainder, are represented
            all_partitions = list(partition_list.keys()) + [remainder]
            missing_partitions = [
                p for p in all_partitions if p not in result["partition"].values
            ]
            if missing_partitions:
                missing_df = pd.DataFrame(
                    {
                        "partition": missing_partitions,
                        "Freq": [0] * len(missing_partitions),
                    }
                )
                result = pd.concat([result, missing_df], ignore_index=True)

            return result

    bed1 = pybedtools.BedTool(
        "/home/bnt4me/virginia/repos/bedboss/scripts/stats/data/b72912f88165f2f252f273f4c3043fe3.bed.gz"
    )

    parts_hg38 = (
        "/home/bnt4me/virginia/repos/bedboss/scripts/stats/partitions/hg38_parts/"
    )
    parts = dict(
        exon=pybedtools.BedTool(f"{parts_hg38}/hg_38.exon.csv"),
        intron=pybedtools.BedTool(f"{parts_hg38}hg_38.intron.csv"),
        prm_cor=pybedtools.BedTool(f"{parts_hg38}hg_38.promoterCore.csv"),
        prm_proc=pybedtools.BedTool(f"{parts_hg38}hg_38.promoterProx.csv"),
        prm_5=pybedtools.BedTool(f"{parts_hg38}hg_38.fiveUTR.csv"),
        prm_3=pybedtools.BedTool(f"{parts_hg38}hg_38.threeUTR.csv"),
    )

    calc_partitions(bed1, parts, bp_proportion=True)


if __name__ == "__main__":
    main2()
