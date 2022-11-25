"""
Calculate the naive solution without accounting for bedding planes
"""

from voussoir import DiederichBeam

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


if __name__ == "__main__":
    # get the distribution dictionary
    out = []
    dist_dict = local.create_uniform_distributions(local.diedrich_inputs)
    for _ in range(10_000):
        sample = local.sample(dist_dict)
        # first calculate pressure on beam
        pressure = get_pressure(sample)
        vb = DiederichBeam(
            thickness=sample['thickness'],
            span=sample['span'],
            ucs=sample['ucs'] * sample['ucs_coef'],
            stiffness=sample['stiffness'],
            joint_stiffness=sample['joint_stiffness'],
            joint_spacing=sample['joint_spacing'],
            joint_friction_angle=sample['joint_friction_angle'],
            density=sample['sandstone_density'],
            pressure=pressure,
        )
        out.append(vb.solve())
    df = pd.DataFrame(out)
    df.to_csv(local.unsupported_df_path, index=False)
