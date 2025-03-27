# IMPORT
import os
import math
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import pandas as pd
import flopy
from flopy.utils.postprocessing import get_water_table
from flopy.mf6.utils.output_util import HeadFile
from .utils import (
    lst_df_from_gwf_long,
    get_final_total_lst
)

# VARIABLES

# FUNCTIONS
# def plot_compare_lstbud_gwfs(gwf0, gwf1, names=['gwf0', 'gwf1']):
#     # THIS IS BROKEN!!
#     lst_df_long_orig = lst_df_from_gwf_long(gwf0)
#     lst_df_long_new = lst_df_from_gwf_long(gwf1)
#     comps = lst_df_long_orig.bcompname.unique()
#     for plot_comp in comps:
#         idf_orig = (
#             lst_df_long_orig
#             .query(f'bcompname == "{plot_comp}"')
#             .loc[:, ['sp', 'rate']]
#         )
#         idf_new = (
#             lst_df_long_new
#             .query(f'bcompname == "{plot_comp}"')
#             .loc[:, ['sp', 'rate']]
#         )
#         # get fig and axes
#         fig, (ax0, ax1, ax2) = plt.subplots(3, 1, figsize=(7.5, 10))
#         # plot side by side
#         ax0.plot(idf_orig.sp, idf_orig.iloc[:, 1], label=names[0])
#         ax0.plot(idf_new.sp, idf_new.iloc[:, 1], label=names[1])
#         ax0.set_title('values')
#         ax0.set_xlabel('stress period')
#         ax0.set_ylabel(f'{plot_comp} in volume per time')
#         ax0.legend()
#         # plot residual
#         idf_resid = idf_orig.iloc[:, 1] - idf_new.iloc[:, 1]
#         ax1.plot(idf_orig.sp, idf_resid)
#         ax1.set_title('difference')
#         ax1.set_xlabel('stress period')
#         ax1.set_ylabel(f'{plot_comp} residual in volume per time')
#         fig.suptitle(plot_comp, fontsize=16)
#         fig.tight_layout()
#         # plot percent residual
#         idf_pct_resid = 100 * (idf_resid / idf_orig.iloc[:, 1]).fillna(0)
#         ax2.plot(idf_orig.sp, idf_pct_resid, color='black')
#         ax2.set_title('percent difference')
#         ax2.set_xlabel('stress period')
#         ax2.set_ylabel(f'{plot_comp} residual in percent')
#         fig.suptitle(plot_comp, fontsize=16)
#         fig.tight_layout()

def plot_compare_obs_sim(sim_path0, sim_path1, modelname=None, names=['sim0', 'sim1'], **kwargs):
    sim0 = flopy.mf6.MFSimulation.load(sim_ws=sim_path0, verbosity_level=0)
    # budget obs files
    if modelname is None:
        gwf0 = sim0.get_model()
    else:
        gwf0 = sim0.get_model(modelname)
    # get dictionary of output filenames
    obs_fileout_dict = {
        x.package_name: list(x.continuous.data.keys()) for x in gwf0.obs
    }
    # loop over and load
    pack_df_dict = {}
    for ipack, ifiles in obs_fileout_dict.items():
        ipack_df_list = []
        for ifile in ifiles:
            for vers, idir in zip(names, [sim_path0, sim_path1]):
                idf = (
                    pd
                    .read_csv(os.path.join(idir, ifile))
                    .melt(id_vars = ['time'], var_name='bound', value_name='boundobs')
                    .assign(package = ipack)
                    .assign(file = ifile)
                    .assign(version = vers)
                )
                ipack_df_list.append(idf)
        pack_df_dict[ipack] = pd.concat(ipack_df_list)
    # plot
    for ipack, ipack_df in pack_df_dict.items():
        nbounds = ipack_df.bound.nunique()
        nrows = min(nbounds, 3)
        ncols = int(np.ceil(nbounds / nrows))
        fig, axes_grid = plt.subplots(nrows, ncols, figsize=(10, 7.5))
        if isinstance(axes_grid, np.ndarray):
            axes = axes_grid.flatten()
        else:
            axes = list([axes_grid])
        i=0
        for ibound, idf in ipack_df.groupby('bound'):
            ax = axes[i]
            idf_piv = idf.pivot(index='time', columns='version', values='boundobs')
            idf_piv.plot(ax=ax, title=ibound, **kwargs)
            i += 1
        fig.tight_layout()
    return fig


def make_xs_line_along_dis_grid(gwf_dis, xy, line_len, along='row'):
    # get grid rotation angle
    angrot_deg = gwf_dis.get_package('dis').angrot.data
    # get angle in radians and shift quadrant if along columns
    if along == 'row':
        angrot_rad = math.radians(angrot_deg)
    elif along == 'column':
        angrot_rad = math.radians(angrot_deg + 90)
    else:
        raise(ValueError('Must be along either row or column'))
    # get line cooridates
    x0 = xy[0] + line_len * math.cos(angrot_rad)
    x1 = xy[0] + line_len * math.cos(angrot_rad + math.pi)
    y0 = xy[1] + line_len * math.sin(angrot_rad)
    y1 = xy[1] + line_len * math.sin(angrot_rad + math.pi)
    xs_line = [(x0, y0), (x1, y1)]
    return xs_line

# def plot_xs_across_pt(gwf, xs_line, xs_ylim=None, nlevs=10, zoom_rel_line=1.5, toplayer=0):
#     # TODO: THIS NEEDS SOME WORK. Doesn't choose good contour levels and more
#     # get data
#     kstpkper = gwf.output.head().get_kstpkper()[-1]
#     # get heads
#     hds = gwf.output.head().get_data(kstpkper=kstpkper)
#     # get water table at final step
#     wt = get_water_table(hds)
#     # create figure and gridspec
#     fig = plt.figure(figsize=(7.5, 10))
#     gs = gridspec.GridSpec(7, 1)
#     # create top and bottom axes
#     ax_map = fig.add_subplot(gs[:4, 0])
#     ax_xs = fig.add_subplot(gs[4:, 0])
    
#     # plot cross-section 
#     # get cross-section plotter
#     xs = flopy.plot.PlotCrossSection(
#         model=gwf,
#         line={'line': xs_line},
#         ax=ax_xs,
#         geographic_coords=True
#     )
#     # grid
#     lc = xs.plot_grid(zorder=-1)
#     # head fill
#     pc = xs.plot_array(hds, head=hds, alpha=0.5, masked_values=[1e30])
#     # head contours
#     # get levels
#     pc_values = pc.get_array().data
#     min_val = min(pc_values)
#     max_val = max(pc_values)
#     val_range = max_val - min_val
#     cont_int = int(np.round(val_range / nlevs))
#     cont_levels = np.arange(
#         np.round(min_val) - 5*cont_int, 
#         np.round(max_val) + 5*cont_int, 
#         cont_int
#     )
#     # plot
#     ctr = xs.contour_array(hds, levels=cont_levels, colors="b", masked_values=[1e30])
#     # water table surface
#     surf = xs.plot_surface(wt, masked_values=[1e30], color="blue", lw=2)
#     # head contour labels
#     labels = xs.ax.clabel(ctr, inline=0.25, fontsize=8, inline_spacing=0)

#     # plot location of line on grid and contours
#     # get a zoom window for the map
#     center_x = (xs_line[0][0] + xs_line[1][0]) / 2
#     center_y = (xs_line[0][1] + xs_line[1][1]) / 2
#     line_length = np.sqrt((xs_line[0][0] - xs_line[1][0])**2 + (xs_line[0][1] - xs_line[1][1])**2)
#     map_xlim = [center_x - zoom_rel_line/2*line_length, center_x + zoom_rel_line/2*line_length]
#     map_ylim = [center_y - zoom_rel_line/2*line_length, center_y + zoom_rel_line/2*line_length]
#     # get map plotter
#     pmv = flopy.plot.PlotMapView(gwf, ax=ax_map, layer=toplayer)
#     # add grid
#     lc = pmv.plot_grid(lw=0.5)
#     # head fill
#     pc = pmv.plot_array(hds, alpha=0.5, masked_values=[1e30])
#     # get levels
#     pc_values = pc.get_array().data
#     min_val = min(pc_values)
#     max_val = max(pc_values)
#     val_range = max_val - min_val
#     cont_int = int(np.round(val_range / nlevs))
#     cont_levels = np.arange(
#         np.round(min_val) - 5*cont_int, 
#         np.round(max_val) + 5*cont_int, 
#         cont_int
#     )
#     # add water table contours
#     ctr = pmv.contour_array(wt, levels=cont_levels, linewidths=0.75, cmap='winter_r')
#     # add line
#     ax_map.plot(
#         [xs_line[0][0], xs_line[1][0]], 
#         [xs_line[0][1], xs_line[1][1]], 
#         lw=0.75,
#         color='red'
#     )
#     # set zoom
#     ax_map.set_xlim(map_xlim)
#     ax_map.set_ylim(map_ylim)
#     # set equal scale
#     ax_map.set_aspect('equal', adjustable='box')
#     ax_xs.set_ylim(xs_ylim)
#     fig.tight_layout()
#     return fig

def plot_interpolated_ws(gwf, xs_line, ax, color='black', label='', iper=-1):
    # get data
    kstpkper = gwf.output.head().get_kstpkper()[iper]
    # get heads
    hds = gwf.output.head().get_data(kstpkper=kstpkper)
    # get water table at final step
    wt = get_water_table(hds)
    # get cross-section plotter
    xs = flopy.plot.PlotCrossSection(
        model=gwf,
        line={'line': xs_line},
        ax=ax,
        geographic_coords=True
    )
    # plot water table
    surf = xs.plot_surface(wt, masked_values=[1e30], color="black", lw=0.0)
    # get centers of lines for interpolated surface plot
    int_surf_xs = []
    int_surf_ys = []
    for lines in surf:
        line = lines[0]
        xs, ys = line.get_data()
        int_surf_xs.append(xs.mean())
        int_surf_ys.append(ys.mean())
    # plot interpolated surface
    int_surf = ax.plot(int_surf_xs, int_surf_ys, color=color, lw=1, label=label)
    return int_surf

def plot_interpolated_heads(gwf, xs_line, ax, layer=0, grid=True, idx=-1, **kwargs):
    # get head file path
    sim_ws = gwf.simulation.sim_path
    hd_name = gwf.oc.head_filerecord.get_data()[0][0]
    hd_file = os.path.join(sim_ws, hd_name)
    # get heads
    lay_hds = (
        HeadFile(hd_file, precision='double')
        .get_data(idx=idx)[layer, :, :]
    )
    # get cross-section plotter
    xs = flopy.plot.PlotCrossSection(
        model=gwf,
        line={'line': xs_line},
        ax=ax,
        geographic_coords=True
    )
    # plot grid
    if grid:
        lc = xs.plot_grid(linewidth=0.5, alpha=0.5)
    # plot heads
    surf = xs.plot_surface(lay_hds, masked_values=[1e30], color="black", lw=0.0)
    # get centers of lines for interpolated surface plot
    int_surf_xs = []
    int_surf_ys = []
    for lines in surf:
        line = lines[0]
        xs, ys = line.get_data()
        int_surf_xs.append(xs.mean())
        int_surf_ys.append(ys.mean())
    # make sure they are sorted by x value
    xs = np.array(int_surf_xs)
    ys = np.array(int_surf_ys)
    # Get the indices that would sort the x array
    sorted_indices = np.argsort(xs)
    # Sort x and reorder y based on the sorted indices
    sorted_x = xs[sorted_indices]
    sorted_y = ys[sorted_indices]
    # plot interpolated surface
    int_surf = ax.plot(sorted_x, sorted_y, **kwargs)
    return int_surf

def plot_compare_final_vols(gwf0, gwf1, names=['gwf1', 'gwf2']):
    # get total volume dataframes
    lst_tot_df0 = get_final_total_lst(gwf0)
    lst_tot_df1 = get_final_total_lst(gwf1)
    # inputs
    ins0 = lst_tot_df0.loc['in']
    ins1 = lst_tot_df1.loc['in'].set_index('btype')
    ins = (
        ins0
        .join(ins1, on='btype', rsuffix='1')
        .rename(columns={'total_volume': names[0], 'total_volume1': names[1]})
        .set_index('btype')
    )
    # outputs
    outs0 = lst_tot_df0.loc['out']
    outs1 = lst_tot_df1.loc['out'].set_index('btype')
    outs = (
        outs0
        .join(outs1, on='btype', rsuffix='1')
        .rename(columns={'total_volume': names[0], 'total_volume1': names[1]})
        .set_index('btype')
    )
    
    # get figure info
    fig, (ax_in, ax_out) = plt.subplots(2, 1, figsize=(7.5, 10))
    # plot both
    ins.plot(kind='bar', ax=ax_in)
    outs.plot(kind='bar', ax=ax_out)
    # format
    ax_in.set_xlabel('')
    ax_out.set_xlabel('')
    ax_in.set_ylabel('total volume')
    ax_out.set_ylabel('total volume')
    ax_in.set_title('IN')
    ax_out.set_title('OUT')
    fig.suptitle('Comparison of listing final cumulative volumes', fontsize=15)
    fig.tight_layout()
    return fig

def remove_axes(ax, x=True, y=True):
    if x:
        ax.xaxis.set_ticks([])
        ax.set_xlabel('')
    if y:
        ax.yaxis.set_ticks([])
        ax.set_ylabel('')
    return ax