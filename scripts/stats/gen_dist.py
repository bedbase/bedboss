import pyranges as pr
import pandas as pd


# Function to calculate cumulative partitions with reference genome annotation
def calc_cumulative_partitions_ref(query, ref_assembly):
    # Validate input
    if not isinstance(query, pr.PyRanges):
        raise ValueError("query must be a PyRanges object")
    if not isinstance(ref_assembly, str):
        raise ValueError("refAssembly must be a string")

    # Get gene models (this part is pseudo-code, you need to define get_gene_models)
    gene_models = get_gene_models(ref_assembly)

    # Create partition list
    partition_list = genome_partition_list(
        gene_models["genesGR"],
        gene_models["exonsGR"],
        gene_models.get("threeUTRGR", None),
        gene_models.get("fiveUTRGR", None),
    )

    # Calculate cumulative partitions
    return calc_cumulative_partitions(query, partition_list)


# Function to create genome partition list
def genome_partition_list(
    genesGR,
    exonsGR,
    threeUTRGR=None,
    fiveUTRGR=None,
    get_core_promoter=True,
    get_prox_promoter=True,
    core_prom_size=100,
    prox_prom_size=2000,
):
    # Validate input
    if not isinstance(genesGR, pr.PyRanges) or not isinstance(exonsGR, pr.PyRanges):
        raise ValueError("genesGR and exonsGR must be PyRanges objects")
    if threeUTRGR and not isinstance(threeUTRGR, pr.PyRanges):
        raise ValueError("threeUTRGR must be a PyRanges object or None")
    if fiveUTRGR and not isinstance(fiveUTRGR, pr.PyRanges):
        raise ValueError("fiveUTRGR must be a PyRanges object or None")

    # # Helper function to generate promoter ranges
    # def get_promoters(granges, upstream, downstream=0):
    #     return granges.slack(upstream, -downstream)
    #
    # # Core and proximal promoters
    # prom_core = get_promoters(genesGR, upstream=core_prom_size) if get_core_promoter else None
    # prom_prox = get_promoters(genesGR, upstream=prox_prom_size) if get_prox_promoter else None
    #
    # # Remove core promoter regions from proximal promoter regions if both exist
    # if prom_core is not None and prom_prox is not None:
    #     prom_prox = prom_prox.subtract(prom_core)

    # Handle UTRs and exons
    if threeUTRGR is not None and fiveUTRGR is not None:
        # Remove UTR overlaps from exons
        fiveUTRGR = fiveUTRGR.subtract(threeUTRGR)
        exonsGR = exonsGR.subtract(threeUTRGR).subtract(fiveUTRGR)

        # Introns = genes - (5'UTR, 3'UTR, exons)
        non_three = genesGR.subtract(threeUTRGR)
        non_three_five = non_three.subtract(fiveUTRGR)
        intronGR = non_three_five.subtract(exonsGR)

    elif threeUTRGR is None and fiveUTRGR is not None:
        exonsGR = exonsGR.subtract(fiveUTRGR)
        intronGR = genesGR.subtract(fiveUTRGR).subtract(exonsGR)

    elif threeUTRGR is not None and fiveUTRGR is None:
        exonsGR = exonsGR.subtract(threeUTRGR)
        intronGR = genesGR.subtract(threeUTRGR).subtract(exonsGR)

    else:
        # No UTRs, introns = genes - exons
        intronGR = genesGR.subtract(exonsGR)

    # Create a partition list
    partition_list = {
        # 'promoterCore': prom_core,
        # 'promoterProx': prom_prox,
        "threeUTR": threeUTRGR,
        "fiveUTR": fiveUTRGR,
        "exon": exonsGR,
        "intron": intronGR,
    }

    # Remove None entries
    return {key: value for key, value in partition_list.items() if value is not None}


# Placeholder for get_gene_models
def get_gene_models():
    # Fetch gene models for the reference assembly, this part should be replaced with actual implementation
    # This could involve loading pre-processed annotations (e.g., from GTF or BED files)
    df = pd.read_csv(
        "/home/bnt4me/Downloads/files.hg38.ensembl_gtf", delimiter="\t", skiprows=5
    )

    df.columns = ["Chromosome", "b", "c", "Start", "End", "r", "a", "z", "y"]

    genes_df = df[df.iloc[:, 2] == "gene"].iloc[:, [0, 3, 4]]
    exons_df = df[df.iloc[:, 2] == "exon"].iloc[:, [0, 3, 4]]
    threeUTR_df = df[df.iloc[:, 2] == "three_prime_utr"].iloc[:, [0, 3, 4]]
    fiveUTR_df = df[df.iloc[:, 2] == "five_prime_utr"].iloc[:, [0, 3, 4]]

    full_dict = {
        "genesGR": pr.PyRanges(genes_df),  # Replace with actual gene ranges
        "exonsGR": pr.PyRanges(exons_df),  # Replace with actual exon ranges
        "threeUTRGR": pr.PyRanges(
            threeUTR_df
        ),  # Replace with 3' UTR ranges (if available)
        "fiveUTRGR": pr.PyRanges(
            fiveUTR_df
        ),  # Replace with 5' UTR ranges (if available)
    }
    # return full_dict

    genome_partition_list(
        full_dict["genesGR"],
        full_dict["exonsGR"],
        full_dict["threeUTRGR"],
        full_dict["fiveUTRGR"],
    )


# Placeholder for calc_cumulative_partitions
def calc_cumulative_partitions(query, partition_list):
    # Implement cumulative partition calculation using partition_list
    pass


#
# if __name__ == "__main__":
#     # # Example usage
#     # df = pd.read_csv("/home/bnt4me/Downloads/files.hg38.ensembl_gtf", delimiter="\t", skiprows=5)
#     # query = pr.PyRanges(df)  # Replace with actual query regions
#     # ref_assembly = "hg38"  # Replace with actual reference assembly
#     # calc_cumulative_partitions_ref(query, ref_assembly)
#     f = get_gene_models()


# ====

import os


def genom_distribution():
    import pybedtools
    import pandas as pd

    # df = pd.read_csv("/home/bnt4me/Downloads/rust_gc/hg38.ensGene.gtf", delimiter="\t")

    annotations = pybedtools.BedTool("/home/bnt4me/Downloads/rust_gc/hg38.ensGene.gtf")

    bed_file = pybedtools.BedTool(
        "/home/bnt4me/virginia/repos/bedboss/scripts/stats/8f37f6b696eb1862328efe7237ddccdf.bed.gz"
    )

    intersected = bed_file.intersect(annotations, wa=True, wb=True)

    if os.stat(intersected.fn).st_size == 0:
        print("No overlapping regions found between BED and annotations.")
    else:
        # Proceed with reading the output into a Pandas DataFrame
        df = pd.read_csv(intersected.fn, sep="\t", header=None)
        print(df.head())


if __name__ == "__main__":
    genom_distribution()
    # region_distribution()
    # import pandas as pd
    # z = pd.read_csv("/home/bnt4me/Downloads/ensembl_gtf__default/default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4.gtf", nrows=20, skiprows=5, sep='\t', header=None)
