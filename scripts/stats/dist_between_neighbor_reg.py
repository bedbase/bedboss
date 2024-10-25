"""
cd /home/bnt4me/virginia/repos/gdrs/bindings/

maturin build --release

uv pip install "gdrs @ /home/bnt4me/virginia/repos/gdrs/bindings/target/wheels/gdrs-0.1.0-cp310-cp310-manylinux_2_34_x86_64.whl"

"""

from gdrs import calc_neighbor_distances
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import seaborn as sns


def calc_dist():
    neighbor_dist = calc_neighbor_distances(
        "/home/bnt4me/virginia/repos/bedboss/scripts/stats/data/04dee6a07f7db768cba73481e7e8401a.bed.gz",
    )

    neighbor_dist_not_null = [x for x in neighbor_dist if x > 0]

    sns.kdeplot(
        data=neighbor_dist_not_null, log_scale=True, color="black", linewidth=0.8
    )

    plt.rcParams["font.size"] = 10
    plt.gca().yaxis.set_major_locator(MaxNLocator(nbins=5))

    sns.despine()
    plt.xlabel("bp distance")

    plt.title("Neighboring regions distance distribution")

    # Show the plot
    plt.show()


if __name__ == "__main__":
    calc_dist()
