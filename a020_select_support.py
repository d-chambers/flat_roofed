"""
Select support by varying over reasonable bolt spacing/depths.
"""
import numpy as np

from voussoir import AbousleimanBeam

import pandas as pd

import local


def get_pressure(kwargs):
    """Get the pressure on the beam from overburden."""

    def get_sandstone_pressure(sandstone_height, thickness, sandstone_density, **kwargs):
        """Return the pressure on the beam exerted by sandstone above it."""
        height = sandstone_height - thickness
        return sandstone_density * height * local.GRAVITATIONAL_CONSTANT

    def get_shale_pressure(shale_height, shale_density, **kwargs):
        """Return pressure from shale on the beam."""
        return shale_density * shale_height * local.GRAVITATIONAL_CONSTANT

    return get_shale_pressure(**kwargs) + get_sandstone_pressure(**kwargs)


def sample(dist_dict, bolt_length):
    """Sample the stability using the Abousleim Beam."""
    out = []

    for _ in range(1000):
        sample = local.sample(dist_dict)
        # first calculate pressure on beam
        pressure = get_pressure(sample)
        vb = AbousleimanBeam(
            thickness=sample['thickness'],
            span=sample['span'],
            ucs=sample['ucs'] * sample['ucs_coef'],
            stiffness=sample['stiffness'],
            joint_stiffness=sample['joint_stiffness'],
            joint_spacing=sample['joint_spacing'],
            joint_friction_angle=sample['joint_friction_angle'],
            density=sample['sandstone_density'],
            pressure=pressure,
            bolt_length=bolt_length,
            horizontal_joint_spacing=sample['horizontal_joint_spacing'],
            brittleness_factor=sample['brittleness_factor'],
        )
        out.append(vb.solve())
    df = pd.DataFrame(out)
    return df


if __name__ == "__main__":
    # get the distribution dictionary
    dist_dict = local.create_uniform_distributions(local.diedrich_inputs)
    bolt_lengths = np.arange(1, 5.25, 0.25)
    passed_rates = []
    for bolt_length in bolt_lengths:
        df = sample(dist_dict, bolt_length=bolt_length)
        passed_rates.append(local.get_percent_passed(df))
        print(f"finished {bolt_length}")

    import matplotlib.pyplot as plt
    plt.plot(bolt_lengths, passed_rates)
    plt.show()


