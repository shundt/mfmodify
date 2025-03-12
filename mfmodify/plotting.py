# IMPORT
import matplotlib.pyplot as plt
from .utils import (
    lst_df_from_gwf_long
)

# VARIABLES

# FUNCTIONS
def plot_compare_lstbud_gwfs(gwf0, gwf1):
    lst_df_long_orig = lst_df_from_gwf_long(gwf0)
    lst_df_long_new = lst_df_from_gwf_long(gwf1)
    comps = lst_df_long_orig.bcompname.unique()
    for plot_comp in comps:
        idf_orig = (
            lst_df_long_orig
            .query(f'bcompname == "{plot_comp}"')
            .loc[:, ['sp', 'rate']]
        )
        idf_new = (
            lst_df_long_new
            .query(f'bcompname == "{plot_comp}"')
            .loc[:, ['sp', 'rate']]
        )
        # get fig and axes
        fig, (ax0, ax1, ax2) = plt.subplots(3, 1, figsize=(7.5, 10))
        # plot side by side
        ax0.plot(idf_orig.sp, idf_orig.iloc[:, 1], label='gwf0')
        ax0.plot(idf_new.sp, idf_new.iloc[:, 1], label='gwf1')
        ax0.set_title('values')
        ax0.set_xlabel('stress period')
        ax0.set_ylabel(f'{plot_comp} in volume per time')
        ax0.legend()
        # plot residual
        idf_resid = idf_orig.iloc[:, 1] - idf_new.iloc[:, 1]
        ax1.plot(idf_orig.sp, idf_resid)
        ax1.set_title('difference')
        ax1.set_xlabel('stress period')
        ax1.set_ylabel(f'{plot_comp} residual in volume per time')
        fig.suptitle(plot_comp, fontsize=16)
        fig.tight_layout()
        # plot percent residual
        idf_pct_resid = 100 * (idf_resid / idf_orig.iloc[:, 1]).fillna(0)
        ax2.plot(idf_orig.sp, idf_pct_resid, color='black')
        ax2.set_title('percent difference')
        ax2.set_xlabel('stress period')
        ax2.set_ylabel(f'{plot_comp} residual in percent')
        fig.suptitle(plot_comp, fontsize=16)
        fig.tight_layout()