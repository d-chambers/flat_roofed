"""
Local variables used across scripts.
"""
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# --- File IO
input_path = Path('inputs')
output_path = Path('outputs')
output_path.mkdir(exist_ok=True, parents=True)

unsupported_df_path = output_path / 'unsupported_fs_df.csv'
unsupported_plot = output_path / 'unsupported_fs.png'

# --- Constants

GRAVITATIONAL_CONSTANT = 9.81  # m/s^2
DESIGN_FACTOR_OF_SAFETY = 2  # unitless

diedrich_inputs = dict(
    thickness=(0.5, 1.0),
    span=(15, 15),
    ucs=(20e6, 40e6),
    stiffness=(10e9, 14e9),
    joint_stiffness=(1000e9, 2000e9),
    joint_spacing=(1, 3),
    joint_friction_angle=(22, 28),
    shale_density=(2450, 2450),
    sandstone_density=(2450, 2450),
    sandstone_height=(5, 5),
    shale_height=(10, 10),
    brittleness_factor=(0.6, 0.6),
    ucs_coef=(0.3, 0.5),
)


# --- Helper Functions


def create_uniform_distributions(input_dict, num=3000):
    """Create a uniform sampling of input values in dict."""
    out = {}
    for i, v in input_dict.items():
        assert len(v) == 2, 'you must provide min/max tuple as key'
        out[i] = np.linspace(v[0], v[1], num=num)
    return out


def sample(dist_dict):
    """Given a dict of arrays, randomly pull one sample for each key."""
    out = {i: np.random.choice(v, 1)[0] for i, v in dist_dict.items()}
    return out


def plot_fs(df):
    """Plot the factor of safety from a dataframe."""
    breakpoint()



