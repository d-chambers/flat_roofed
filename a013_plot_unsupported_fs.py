"""
Plot the unsupported factors of safety.
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import local



def init_fs_plot():
    """Create the figure and axis."""
    fig, ax = plt.subplots(2, 2, figsize=(7, 7), sharey=True)
    return fig, ax.flatten()


def plot_hist(ser, ax, title):
    """Plot the histograms."""
    bins = np.linspace(0, 5)
    ax.hist(ser.values, bins=bins)
    mean, median, std = ser.mean(), ser.median(), ser.std()
    gt_fs = ser >= local.DESIGN_FACTOR_OF_SAFETY
    fs_percent = int((gt_fs.sum() / len(ser)) * 100)
    new_title = f"{title} ({fs_percent}% passed)"
    ax.set_title(new_title)
    ax.axvline(local.DESIGN_FACTOR_OF_SAFETY, ls='--', color='red', alpha=0.4)


def plot(df, axes, fig, suptitle='Unsupported'):
    """First get """
    failed = np.isclose(df['buckling_limit'], 100)
    sub = df[~failed]

    solvable_percent = int(np.round((len(sub) / len(df)) * 100))
    overall_fs = local.get_overall_fs(sub)

    plot_hist(overall_fs, axes[0], 'Overall FS')
    plot_hist(sub['buckling_fs'], axes[1], 'Buckling FS')
    plot_hist(sub['sliding_fs'], axes[2], 'Sliding FS')
    plot_hist(sub['crushing_fs'], axes[3], 'Crushing FS')

    fig.supxlabel('Factor of Safety')
    fig.supylabel('Sample Count')
    fig.suptitle(f'{suptitle}\n({solvable_percent}% solvable)')


if __name__ == "__main__":
    df = pd.read_csv(local.unsupported_df_path)
    fig, axes = init_fs_plot()
    plot(df, axes, fig)
    plt.tight_layout()
    fig.savefig(local.unsupported_plot)
