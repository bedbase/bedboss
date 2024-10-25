"""
cd /home/bnt4me/virginia/repos/gdrs/bindings/

maturin build --release

uv pip install "gdrs @ /home/bnt4me/virginia/repos/gdrs/bindings/target/wheels/gdrs-0.1.0-cp310-cp310-manylinux_2_34_x86_64.whl"

"""

from gdrs import calc_widths
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import seaborn as sns
import statistics
import numpy as np
import pandas as pd


def width_hist():
    widths = calc_widths(
        "/home/bnt4me/virginia/repos/bedboss/scripts/stats/data/04dee6a07f7db768cba73481e7e8401a.bed.gz",
    )
    mean_width = statistics.mean(widths)
    print(widths[0:40])
    print(mean_width)
    data = widths

    # lower_quantile = 0.05  # Trim the lowest 5%
    # upper_quantile = 0.95
    #
    # lower_bound = np.percentile(data, lower_quantile * 100)
    # upper_bound = np.percentile(data, upper_quantile * 100)

    # Trim the data
    # trimmed_data = [x for x in data if lower_bound <= x <= upper_bound]

    # Create the histogram
    plt.figure(figsize=(12, 8))
    hist = sns.histplot(
        data,
        bins=list(range(20, 120, 3)),
        kde=False,
        stat="percent",
        color="darkgray",
        edgecolor="black",
    )

    x_ticks = np.arange(0, 80, 5)  # Adjust based on your data range
    plt.show()


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def calcDivisions(x, quantThresh=None, bins=None):
    # Calculate divisions using quantiles or specified number of bins
    if quantThresh is not None:
        divisions = np.quantile(x, quantThresh)
    else:
        divisions = np.linspace(
            min(x), max(x), bins + 1
        )  # Create uniform bins if no quantile thresholds provided
    return {
        "divisions": [-np.inf] + list(divisions) + [np.inf],
        "bins": len(divisions) - 1,
    }


def cutDists(x, divisions):
    # Create bins and frequencies
    cuts = pd.cut(x, bins=divisions, include_lowest=True)
    df = pd.DataFrame(cuts.value_counts().reset_index())
    df.columns = ["cuts", "Freq"]
    return df


def plotQTHist(x, quantThresh=None, bins=None, indep=False, numbers=False):
    if isinstance(x, list):
        if indep:
            plots = []
            for i, data in enumerate(x):
                g = plotQTHist(data, quantThresh, bins, indep=False, numbers=numbers)
                g.set_title(f"Plot {i + 1}")
                plots.append(g)
            return plots

    # Calculate divisions and handle recalculations
    output = calcDivisions(x, quantThresh=quantThresh, bins=bins)
    divisionCheck = output["divisions"]
    if len(divisionCheck) > len(set(divisionCheck)):
        if len(set(divisionCheck)) == 3:
            output["divisions"] = [
                -np.inf,
                divisionCheck[1],
                divisionCheck[1] + 1,
                np.inf,
            ]
            output["bins"] = 1
        else:
            output["divisions"] = sorted(set(divisionCheck))
            output["bins"] = len(set(divisionCheck)) - 3

    # Create a dataframe from the cuts and calculate frequencies
    df = cutDists(x, divisions=output["divisions"])
    if not numbers:
        df["Freq"] = df["Freq"] / df["Freq"].sum() * 100

    # Create the histogram plot
    plt.figure(figsize=(8, 6))

    sns.barplot(x="cuts", y="Freq", data=df)

    # Add axis labels, title, and formatting
    plt.xticks(rotation=90)
    plt.title("Quantile Trimmed Histogram", ha="center")
    plt.ylabel("Percentage" if not numbers else "Frequency")
    plt.xlabel("")

    # Remove top and right box lines
    sns.despine()

    # Display the plot
    plt.show()


# Example usage


# if __name__ == "__main__":
#     # width_bar()
#
#     x=  calc_widths(
#         "/home/bnt4me/virginia/repos/bedboss/scripts/stats/data/49666564c051c7ece0cacff6b74e1d90.bed.gz",
#     )
#     plotQTHist(x, bins= 30)
