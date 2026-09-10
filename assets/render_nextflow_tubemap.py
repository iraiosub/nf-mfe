"""Render the nf-mfe workflow as a slide-ready tube map (SVG/PDF/PNG).

Run: python3 assets/render_nextflow_tubemap.py
Requires matplotlib. Topology follows main.nf and modules/local/*.nf;
this is a curated diagram, so update it when the workflow changes.
"""

import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "nf-mfe-matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch


OUT = Path(__file__).resolve().parent
INK, MUTED = "#243444", "#596A78"
GREEN, BLUE, PURPLE, ORANGE, REF = "#16866F", "#287BB5", "#8861AD", "#D28A23", "#567483"
PALE = "#A9BBB4"


def render():
    plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none", "pdf.fonttype": 42})
    fig = plt.figure(figsize=(16, 9), facecolor="white")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set(xlim=(0, 1920), ylim=(1080, 0), aspect="equal")
    ax.axis("off")

    def text(x, y, label, size=26, color=INK, weight="normal", **kwargs):
        return ax.text(x, y, label, fontsize=size * .6, color=color,
                       weight=weight, va="center", zorder=6, **kwargs)

    def line(points, color=GREEN, width=9, dashed=False, zorder=2):
        x, y = zip(*points)
        ax.plot(x, y, color=color, lw=width * .6, solid_capstyle="round",
                solid_joinstyle="round", linestyle=(0, (4, 4)) if dashed else "-",
                zorder=zorder)

    def station(x, y, color=GREEN, radius=11):
        ax.add_patch(Circle((x, y), radius, facecolor="white", edgecolor=color,
                            linewidth=3 * .6, zorder=5))

    def detail(x, y, title, rows, color=GREEN):
        line([(x, y), (x + 44, y)], color, 3)
        text(x, y + 32, title, 26, color, "bold")
        for i, row in enumerate(rows):
            text(x, y + 70 + 33 * i, row, 23, MUTED)

    def output(x, y, w, h, title, rows, color=GREEN, fill="#EFF8F4"):
        ax.add_patch(FancyBboxPatch((x, y), w, h,
            boxstyle="round,pad=0,rounding_size=14", facecolor=fill,
            edgecolor="#D9E4DF", linewidth=.8, zorder=3))
        text(x + 22, y + 34, title, 27, color, "bold")
        for i, row in enumerate(rows):
            text(x + 22, y + 74 + 32 * i, row, 23, MUTED)

    text(65, 65, "nf-mfe", 50, weight="bold")
    text(65, 120, "From chimeric interaction coordinates to RNA duplex energies and control comparisons", 29, MUTED)
    text(1855, 65, "NEXTFLOW WORKFLOW", 20, GREEN, "bold", ha="right")
    line([(65, 160), (1855, 160)], "#E0E6E4", 1)

    # Alternative input channels converge at the first process.
    line([(85, 360), (315, 360), (395, 440), (420, 440)], BLUE)
    line([(85, 535), (300, 535), (395, 440), (420, 440)], PURPLE)
    station(85, 360, BLUE, 9)
    station(85, 535, PURPLE, 9)
    text(75, 315, "Samplesheet TSV", 28, BLUE, "bold")
    text(75, 402, "Sample IDs + paths", 23, MUTED)
    text(75, 490, "File pattern", 28, PURPLE, "bold")
    text(75, 580, "Chimeric interaction tables", 23, MUTED)
    text(75, 225, "CHOOSE AN INPUT", 20, MUTED, "bold")

    # The genome reference is a required side input to extraction.
    line([(645, 285), (785, 285), (860, 360), (860, 440)], REF, 7)
    station(645, 285, REF, 9)
    text(630, 225, "Reference genome", 27, REF, "bold")
    text(630, 255, "FASTA + FAI index", 23, MUTED)

    # Shared per-chunk processing.
    line([(420, 440), (1160, 440)])
    for x in (420, 640, 860, 1080):
        station(x, 440)
    line([(420, 427), (420, 375)], PALE, 1.5)
    text(450, 350, "Split tables", 27, weight="bold", ha="center")
    text(640, 395, "Prepare BED", 27, weight="bold", ha="center")
    text(860, 487, "Extract sequences", 27, weight="bold", ha="center")
    text(1080, 395, "Add sequences", 27, weight="bold", ha="center")
    text(640, 360, "Left + right arms", 22, MUTED, ha="center")
    text(860, 521, "bedtools getfasta -s", 22, MUTED, ha="center")
    text(1080, 360, "lseq / rseq columns", 22, MUTED, ha="center")

    line([(420, 453), (420, 640), (385, 675)], PALE, 1.5)
    detail(360, 685, "Parallel chunks", ["10,000 rows by default", "Preserve each header", "Retain sample identity"])
    line([(860, 453), (860, 615), (795, 680)], PALE, 1.5)
    detail(730, 685, "Strand-aware sequences", ["BED-style coordinates", "0-based, half-open intervals", "Match both arms by read name"])

    # Exactly one MFE process runs. Both routes contain the observed fold.
    line([(1160, 440), (1280, 320), (1440, 320), (1560, 440)])
    line([(1160, 440), (1280, 560), (1440, 560), (1560, 440)], ORANGE)
    station(1340, 320)
    station(1340, 560, ORANGE)
    text(1340, 240, "Observed MFE", 28, GREEN, "bold", ha="center")
    text(1340, 278, "Default route", 23, MUTED, ha="center")
    text(1340, 395, "ViennaRNA", 24, MUTED, ha="center")
    text(1340, 429, "duplex folding", 24, MUTED, ha="center")
    text(1340, 610, "MFE + controls", 28, ORANGE, "bold", ha="center")
    text(1340, 646, "Shuffles + optional reversed arms", 22, MUTED, ha="center")
    detail(1150, 715, "Controls route enabled by", ["--shuffled_mfe or", "--flipped_arm_mfe", "Defaults: 100 shuffles · k-let = 2"], ORANGE)

    # Per-sample gather and published outputs.
    line([(1560, 440), (1640, 440)])
    station(1640, 440)
    text(1670, 487, "Merge chunks", 28, weight="bold")
    text(1670, 525, "Group by sample", 22, MUTED)
    text(1670, 558, "Keep one header", 22, MUTED)
    line([(1640, 440), (1740, 340), (1785, 340), (1785, 309)])
    output(1505, 185, 350, 124, "Sample MFE table", ["{sample_id}_mfe.tsv", "Sequences · MFE · structure"])

    # Plotting only follows concatenation when the controls route is selected.
    line([(1640, 453), (1640, 660), (1750, 770), (1785, 770), (1785, 830)], ORANGE, 7, True)
    station(1750, 770, ORANGE)
    text(1670, 635, "Plot summary", 28, ORANGE, "bold")
    text(1670, 674, "Controls route only", 22, MUTED)
    output(1545, 830, 310, 130, "Figures & statistics", ["Plots · PNG / PDF", "Summary statistics · TSV"], ORANGE, "#FCF6EC")

    text(75, 930, "Read left to right  →", 22, MUTED)
    text(730, 930, "Controls also publish per-chunk shuffle detail tables.", 22, MUTED)
    line([(65, 1000), (1855, 1000)], "#E0E6E4", 1)
    line([(75, 1040), (120, 1040)], GREEN, 7)
    station(98, 1040, GREEN, 7)
    text(140, 1040, "Workflow step", 22, MUTED)
    line([(390, 1040), (435, 1040)], ORANGE, 7)
    text(455, 1040, "Alternative MFE route", 22, MUTED)
    line([(785, 1040), (830, 1040)], ORANGE, 5, True)
    text(850, 1040, "Conditional output", 22, MUTED)
    text(1855, 1040, "Run reports · timeline · trace · execution DAG", 21, MUTED, ha="right")

    for ext in ("svg", "pdf", "png"):
        destination = OUT / f"nextflow_tubemap.{ext}"
        fig.savefig(destination, dpi=240, facecolor="white")
        print(destination)
    plt.close(fig)


if __name__ == "__main__":
    render()
