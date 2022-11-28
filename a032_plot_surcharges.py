"""
Make a plot of surcharges percentage and bolt length.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import local



def plot(df):
    bl_gb = df.groupby('bolt_length')
    fig, axes = plt.subplots(len(bl_gb), 1)


    for ax, (bl, sub) in zip(axes, df.groupby('bolt_length')):
        surcharge_percents = np.sort(sub['surcharge_percent'].unique())
        results = []
        for surchage in surcharge_percents:
            results.append(sub[sub['surcharge_percent'] == surchage].pipe(local.aggregate_passed_rates))
        sur_df = pd.DataFrame(results)
        for name, ser in sur_df.items():
            if name == 'fs':
                ax.plot(surcharge_percents, ser, label=name, alpha=0.3, ls='--')
            else:
                ax.plot(surcharge_percents, ser, label=name)
        ax.legend()
        ax.set_title(f"Bolt Length: {bl:.02f}m")





if __name__ == "__main__":
    df = pd.read_csv(local.surchage_varied_path)
    fig, ax = plot(df)

