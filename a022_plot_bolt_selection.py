"""
Plot the bolt selection.
"""
import matplotlib.pyplot as plt
import pandas as pd

import local


def aggregate_passed_rates(df):
    """Aggregate the various passed rates."""
    out = {}
    df_len = len(df)
    for col in [x for x in df.columns if x.endswith('_fs')]:
        out[col] = (df[col] >= local.DESIGN_FACTOR_OF_SAFETY).sum() / df_len
    return pd.Series(out)




if __name__ == "__main__":
    df = pd.read_csv(local.bolt_length_df_path)

    passed = df.groupby('bolt_length').apply(aggregate_passed_rates)
    breakpoint()



