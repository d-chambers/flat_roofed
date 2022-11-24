"""
Script to calculate barrier pillar width
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DDF = pd.Series(
    {0.5: 0.53, 0.6: 0.41, 0.7: 0.31, 0.8: 0.25, 0.9: 0.2, 1.0: 0.18, 1.1: 0.17, 1.2: 0.16}
)


def get_ldf(w, h):
    """Get the LDF factor based on assumptions made in the report."""
    ser = pd.Series(DDF)
    # get closest DDF value
    wth = w / h
    inds = np.searchsorted(ser.index, wth, side='left')
    inds[inds == len(ser)] = len(ser) - 1
    # Note: FF is 1 so it is not featured here
    return 1 - ser.iloc[inds].values


def get_barrier_pillar_strength(w, h=15, UCS=120, l=300):
    """Get the strength of the pillars in MPa."""
    const = 0.65  # science is done in m
    lbr = 1
    area = w * l
    circ = 2 * w + 2 * l
    effective_width = w + ((4 * area) / (circ) - w) * lbr
    ldf = get_ldf(effective_width, h)
    return const * UCS * ldf * (effective_width ** 0.3) / (h ** 0.59)


def get_abutment_stress(w, beta=21, h=250, rho=2700, l=300):
    """
    Get the stress acting on the pillars, in MPa
    """
    a = w * h + h * h * np.tan(np.deg2rad(beta))
    weight = a * rho * l * 9.81
    stress = weight / (w * l)
    return stress / 1_000_000


def plot_barrier_pillar_width(w=np.arange(1, 100, 0.5), fos=1.8):
    """Make a plot of various pillar widths. """
    stress = get_abutment_stress(w)
    strength = get_barrier_pillar_strength(w)
    sf = strength / stress
    fig, ax = plt.subplots(1, 1, figsize=(5,4))
    ax.plot(w, sf)
    ax.set_xlabel('Pillar Width (m)')
    ax.set_ylabel('Factor of Safety')
    meets_req = np.argmax(sf >= fos)
    selected_width = w[meets_req]
    ax.axvline(selected_width, color='black', ls='--')
    return ax


if __name__ == "__main__":
    plot_barrier_pillar_width()


