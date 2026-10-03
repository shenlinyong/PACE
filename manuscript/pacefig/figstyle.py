"""Journal-style matplotlib defaults (Nature family: 89/183 mm widths, 5-7 pt text, vector output)."""
import matplotlib as mpl
import matplotlib.pyplot as plt

MM = 1 / 25.4
COL = dict(pace="#D55E00", abc="#0072B2", abc2="#56B4E9", dist="#8C8C8C", gbm="#009E73", pink="#CC79A7", gold="#E69F00",
           ink="#222222", mid="#6B6B6B", light="#D9D9D9", teal="#00796B", navy="#1F3A5F")
CASE = {"A": "#D55E00", "B_same": "#0072B2", "B_other": "#E69F00", "C": "#8C8C8C"}

def setup():
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Liberation Sans", "Arial", "DejaVu Sans"],
        "font.size": 6, "axes.labelsize": 6.5, "axes.titlesize": 6.5, "xtick.labelsize": 6, "ytick.labelsize": 6,
        "legend.fontsize": 5.5, "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5, "lines.linewidth": 1.0, "axes.spines.top": False,
        "axes.spines.right": False, "axes.edgecolor": COL["ink"], "axes.labelcolor": COL["ink"], "text.color": COL["ink"],
        "xtick.color": COL["ink"], "ytick.color": COL["ink"], "figure.dpi": 150, "savefig.dpi": 600,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none", "mathtext.default": "regular",
        "legend.frameon": False, "axes.grid": False,
    })

def panel(ax, letter, x=-0.18, y=1.06, size=8):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=size, fontweight="bold", va="bottom", ha="left")

def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(f"figures/{name}.{ext}", bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
