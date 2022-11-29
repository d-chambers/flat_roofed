"""
Select support by varying over reasonable bolt spacing/depths.
"""
import numpy as np

from voussoir import AbousleimanBeam

import pandas as pd

import local


def sample(dist_dict, bolt_length, iterations=1000):
    """
    Sample the stability using the Abousleim Beam.

    Need to set density to 0 and use the entire column weight
    """
    out = []

    for _ in range(iterations):
        sample = local.sample(dist_dict)
        # first calculate pressure on beam
        vb = AbousleimanBeam(
            thickness=sample['thickness'],
            span=sample['span'],  # sample['span'],
            rockmass_ucs=sample['ucs'] * sample['ucs_coef'],
            stiffness=sample['stiffness'],
            joint_stiffness=sample['joint_stiffness'],
            joint_spacing=sample['joint_spacing'],
            joint_friction_angle=sample['joint_friction_angle'],
            density=0,
            pressure=local.get_total_pressure(sample),
            bolt_length=bolt_length,
            horizontal_joint_spacing=sample['horizontal_joint_spacing'],
            brittleness_factor=sample['brittleness_factor'],
        )
        breakpoint()
        out.append(vb.solve())
    return pd.DataFrame(out)


if __name__ == "__main__":
    # get the distribution dictionary
    dist_dict = local.create_uniform_distributions(local.diedrich_inputs)
    bolt_lengths = np.arange(1, 5.25, 0.25)
    dfs = []
    for bolt_length in [3]:  # bolt_lengths:  # bolt_lengths:  # bolt_lengths:
        breakpoint()
        df = sample(dist_dict, bolt_length=bolt_length, iterations=1)
        df['bolt_length'] = bolt_length
        dfs.append(df)
        print(f"finished {bolt_length}")
    results = pd.concat(dfs, axis=0).reset_index(drop=True)
    breakpoint()
    results.to_csv(local.bolt_length_df_path)
