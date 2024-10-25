"""
cd /home/bnt4me/virginia/repos/gdrs/bindings/

maturin build --release

uv pip install "gdrs @ /home/bnt4me/virginia/repos/gdrs/bindings/target/wheels/gdrs-0.1.0-cp310-cp310-manylinux_2_34_x86_64.whl"

"""

from gdrs import calc_gc_content, GenomeAssembly


def calc_gc(hg38, with_plot=False):

    gc_contents = calc_gc_content(
        "/home/bnt4me/virginia/repos/bedboss/scripts/stats/data/04dee6a07f7db768cba73481e7e8401a.bed.gz",
        hg38,
        ignore_unk_chroms=True,
    )

    import statistics

    gc_mean = statistics.mean(gc_contents)
    print(f"GC mean: {gc_mean}")

    if with_plot:
        import matplotlib.pyplot as plt
        import seaborn as sns
        from matplotlib.ticker import MaxNLocator

        plt.rcParams["font.size"] = 10
        sns.kdeplot(gc_contents, linewidth=0.8, color="black")
        plt.gca().xaxis.set_major_locator(MaxNLocator(nbins=5))

        plt.axvline(
            gc_mean,
            color="r",
            linestyle="--",
            linewidth=0.8,
            label=f"Mean: {gc_mean:.2f}",
        )
        sns.despine()
        plt.xlabel("GC Content")
        # Add a legend
        plt.legend()

        plt.title("GC Content Distribution")

        # Show the plot
        plt.show()


if __name__ == "__main__":

    import time

    start1 = time.time()
    hg38 = GenomeAssembly("/home/bnt4me/.refgenie/alias/hg38/fasta/default/hg38.fa")
    start = time.time()
    print(f"Time taken: {start1 - start}")
    calc_gc(hg38=hg38, with_plot=True)
    end = time.time()
    print(f"Time taken: {end-start}")
