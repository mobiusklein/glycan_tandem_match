from ms_peak_picker.utils import draw_peaklist
from matplotlib import pyplot as plt
import csv


def write_annotated_peaklist(file_handle, peak_list):
    writer = csv.writer(file_handle)
    writer.writerow(["mz", "charge", "intensity", "annotation"])
    for peak in peak_list:
        for solution in peak.solution.name.split("|"):
            writer.writerow([peak.mz, peak.charge, peak.intensity, solution])


def annotate_spectrum(peaks, annotated_peaks, mz_region_start=120, mz_region_end=3000):
    ax = draw_peaklist(peaks, alpha=0.5)
    upper = max(ax.get_ylim())
    for solution in annotated_peaks:
        draw_peaklist(solution.fit.theoretical, alpha=0.8, ax=ax, color='red')
        y = max(p.intensity for p in solution.fit.experimental)
        y = min(y + 2000, upper * 0.8)
        label = "%s" % sorted(solution.solution.name.split("|"), key=len)[0]
        if solution.charge > 1:
            label += "$^{%d}$" % solution.charge
        ax.text(solution.mz, y, label, rotation=90, va='bottom', ha='center', fontsize=16)

    ax.set_xlim(mz_region_start, mz_region_end)
    max_xticks = 500
    xloc = plt.MaxNLocator(max_xticks)
    ax.xaxis.set_major_locator(xloc)
    for tl in ax.get_xticklabels():
        tl.set(rotation=90)
    fig = ax.get_figure()
    fig.set_figwidth(128)
    fig.set_figheight(64)
    return ax
