"""
Local variables used across scripts.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --- File IO
input_path = Path('inputs')
output_path = Path('outputs')
output_path.mkdir(exist_ok=True, parents=True)

unsupported_sensitivity_path = output_path / 'a012_thickness_sensitivity.png'

unsupported_df_path = output_path / 'a013_unsupported_fs_df.csv'
unsupported_plot = output_path / 'a013_unsupported_fs.png'

bolt_length_df_path = output_path / 'a020_support_selection.csv'

surchage_varied_path = output_path / 'a030_varied_surcharges.csv'

# --- Constants

GRAVITATIONAL_CONSTANT = 9.81  # m/s^2
DESIGN_FACTOR_OF_SAFETY = 2  # unitless

diedrich_inputs = dict(
    thickness=(5, 5),
    span=(15, 15),
    ucs=(20e6, 40e6),
    stiffness=(10e9, 14e9),
    joint_stiffness=(1500e9, 4000e9),
    joint_spacing=(1, 3),
    # joint_friction_angle=(30, 40),
    joint_friction_angle=(28, 45),

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


def aggregate_passed_rates(df):
    """Aggregate the various passed rates."""
    out = {}
    df_len = len(df)
    for col in [x for x in df.columns if x.endswith('_fs')]:
        out[col] = (df[col] >= DESIGN_FACTOR_OF_SAFETY).sum() / df_len
    ser = pd.Series(out)
    ser['fs'] = ser.min()
    return ser


def sample(dist_dict, max_only=False):
    """Given a dict of arrays, randomly pull one sample for each key."""
    if max_only:
        out = {i: v[-1] for i, v in dist_dict.items()}
    else:
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

    return (get_shale_pressure(**kwargs) + get_sandstone_pressure(**kwargs)) * (2/3)


def get_total_pressure(kwargs):
    """Get the pressure on the beam from overburden."""

    def get_sandstone_pressure(sandstone_height, sandstone_density, **kwargs):
        """Return the pressure on the beam exerted by sandstone above it."""
        return sandstone_density * sandstone_height * GRAVITATIONAL_CONSTANT

    def get_shale_pressure(shale_height, shale_density, **kwargs):
        """Return pressure from shale on the beam."""
        return shale_density * shale_height * GRAVITATIONAL_CONSTANT

    return get_shale_pressure(**kwargs) + get_sandstone_pressure(**kwargs)
