"""
Plot the bolt selection.
"""
import matplotlib.pyplot as plt
import pandas as pd

import local



if __name__ == "__main__":
    df = pd.read_csv(local.bolt_length_df_path)

    passed = df.groupby('bolt_length').apply(local.aggregate_passed_rates)
    breakpoint()



