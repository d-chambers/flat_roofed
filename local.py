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

unsupported_sensitivity_path = output_path / 'thickness_sensitivity.png'

unsupported_df_path = output_path / 'unsupported_fs_df.csv'
unsupported_plot = output_path / 'unsupported_fs.png'

# --- Constants

GRAVITATIONAL_CONSTANT = 9.81  # m/s^2
DESIGN_FACTOR_OF_SAFETY = 2  # unitless

diedrich_inputs = dict(
    thickness=(5, 5),
    span=(15, 15),
    ucs=(20e6, 40e6),
    stiffness=(10e9, 14e9),
    joint_stiffness=(1000e9, 2000e9),
    joint_spacing=(1, 3),
    # joint_friction_angle=(30, 40),
    joint_friction_angle=(25, 25),

    shale_density=(2450, 2450),
    sandstone_density=(2450, 2450),
    sandstone_height=(5, 5),
    shale_height=(10, 10),
    brittleness_factor=(0.6, 0.8),
    ucs_coef=(0.3, 0.5),
    horizontal_joint_spacing=(0.5, 1),
)

def get_percent_passed(df):
    """Determine how many of the test cases meet the FS."""
    total = len(df)
    unsolvable = df['buckling_limit'] >= 100
    sub = df[~unsolvable]
    fs = get_overall_fs(sub)
    return (fs >= DESIGN_FACTOR_OF_SAFETY).sum() / total


def get_overall_fs(df):
    """Get the overall factor of safety."""
    fs_cols = [x for x in df.columns if x.endswith('_fs')]
    fs_df = df[fs_cols]
    overall_fs = fs_df.min(axis=1)
    return overall_fs

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


def get_pressure(kwargs):
    """Get the pressure on the beam from overburden."""

    def get_sandstone_pressure(sandstone_height, thickness, sandstone_density, **kwargs):
        """Return the pressure on the beam exerted by sandstone above it."""
        height = sandstone_height - thickness
        return sandstone_density * height * GRAVITATIONAL_CONSTANT

    def get_shale_pressure(shale_height, shale_density, **kwargs):
        """Return pressure from shale on the beam."""
        return shale_density * shale_height * GRAVITATIONAL_CONSTANT

    return get_shale_pressure(**kwargs) + get_sandstone_pressure(**kwargs)



