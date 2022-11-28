"""
Script to vary the surcharge
"""
import itertools

import numpy as np

from voussoir import AbousleimanBeam

import pandas as pd

import local


def sample(dist_dict, bolt_length, surcharge, iterations=1000):
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
            ucs=sample['ucs'] * sample['ucs_coef'],
            stiffness=sample['stiffness'],
            joint_stiffness=sample['joint_stiffness'],
            joint_spacing=sample['joint_spacing'],
            joint_friction_angle=sample['joint_friction_angle'],
            density=0,
            pressure=local.get_total_pressure(sample) * surcharge / 100.,
            bolt_length=bolt_length,
            horizontal_joint_spacing=sample['horizontal_joint_spacing'],
            brittleness_factor=sample['brittleness_factor'],
        )
        out.append(vb.solve())
    return pd.DataFrame(out)


if __name__ == "__main__":
    # get the distribution dictionary
    dist_dict = local.create_uniform_distributions(local.diedrich_inputs)
    bolt_lengths = [1.0, 2.0, 3.0, 4.0]
    # bolt_lengths = np.linspace(1, 4, num=100)
    # surcharge_percentages = range(10, 110, 10)
    surcharge_percentages = np.linspace(10, 100, num=15)
    dfs = []
    for bolt_length, surchage in itertools.product(bolt_lengths, surcharge_percentages):
        df = sample(dist_dict, bolt_length=bolt_length, surcharge=surchage, iterations=500)
        df['bolt_length'] = bolt_length
        df['surcharge_percent'] = surchage
        print(f"finished {bolt_length}, {surchage}")
        dfs.append(df)
    results = pd.concat(dfs, axis=0).reset_index(drop=True)
    results.to_csv(local.surchage_varied_path)
