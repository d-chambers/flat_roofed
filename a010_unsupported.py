"""
Calculate the naive solution without accounting for bedding planes
"""

from voussoir import DiederichBeam

import pandas as pd

import local



if __name__ == "__main__":
    # get the distribution dictionary
    out = []
    dist_dict = local.create_uniform_distributions(local.diedrich_inputs)
    for _ in range(10_000):
        sample = local.sample(dist_dict)
        # first calculate pressure on beam
        pressure = local.get_pressure(sample)
        vb = DiederichBeam(
            thickness=sample['horizontal_joint_spacing'],
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
