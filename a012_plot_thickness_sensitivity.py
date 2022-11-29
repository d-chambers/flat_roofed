"""
Plot the sensitivity to thickness for different failure mechanisms.
"""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

import local
from voussoir import DiederichBeam


def test_with_max_values():
    """
    Ensure stability is achieved with max values from scenario. IE, 5m thick beam.
    """
    out = []
    thicknesses = np.arange(0.5, 5, 0.01)
    for thickness in thicknesses:
        sample = {i: v[1] for i, v in local.diedrich_inputs.items()}
        sample['thickness'] = thickness
        # first calculate pressure on beam
        pressure = local.get_total_pressure(sample)
        vb = DiederichBeam(
            thickness=thickness,
            span=sample['span'],
            rockmass_ucs=sample['ucs'] * sample['ucs_coef'],
            stiffness=sample['stiffness'],
            joint_stiffness=sample['joint_stiffness'],
            joint_spacing=sample['joint_spacing'],
            joint_friction_angle=sample['joint_friction_angle'],
            density=0,
            pressure=pressure,
        )
        ser = vb.solve()
        sub = ser[[x for x in ser.index if x.endswith('_fs')]]
        sub['thickness'] = thickness
        out.append(sub)
    df = pd.DataFrame(out).set_index('thickness')
    return df


def plot_sensitivity(df):
    """Plot the sensitivity to thickness"""
    fig, ax = plt.subplots(1, 1)
    for label, ser in df.items():
        ax.plot(ser.index, ser.values, label=label.replace('_fs', ''))
    ax.legend()
    ax.set_xlabel('Beam Thickness (m)')
    ax.set_ylabel('Factor of Safety')
    ax.set_ylim([0, 5])
    ax.axhline(2.0, ls='--', color='red')
    return fig, ax


if __name__ == "__main__":
    df = test_with_max_values()
    fig, ax = plot_sensitivity(df)
    fig.savefig(local.unsupported_sensitivity_path)

