# Packages
# Folder
from folder import *
import glob

# Handle data
import numpy as np
import pandas as pd

# Visualise
import matplotlib.pyplot as plt

# ---------------------------------------------------

def statistic_df_generation(sims_df):

    sims_copy_df = sims_df.copy(deep=True)
    sims_copy_df['mean'] = sims_df.mean(axis=1)
    sims_copy_df['sd'] = sims_df.std(axis=1)
    sims_copy_df['cv'] = sims_df.std(axis=1)/sims_df.mean(axis=1)
    sims_copy_df['min'] = sims_df.min(axis=1)
    sims_copy_df['max'] = sims_df.max(axis=1)

    return sims_copy_df

def plot_paras_each(kind, para_name, ind_for_zoom, depth_calculation=True):

    # Main path
    main_dir = fr"{bathy_path}\\bathy_analyses"

    # Get data
    if para_name == 'width':
        err_sd_sims = pd.read_csv(fr"{main_dir}\err_sd_width_{kind}.csv")
        err_meter_sims = pd.read_csv(fr"{main_dir}\err_meter_width_{kind}.csv")
        para_sims = pd.read_csv(fr"{main_dir}\para_width_{kind}.csv")

        ori_depth_bed_wse_df = pd.read_csv(fr"{main_dir}\ori_depth_bed_wse.csv")
        ori_para_df = pd.read_csv(fr"{main_dir}\original_para_widths_mean_0.25km.csv")

        if depth_calculation:
            depth_sims = pd.read_csv(fr"{main_dir}\depth_width_{kind}.csv")
        else:
            pass


    elif para_name == 'slope':
        err_sd_sims = pd.read_csv(fr"{main_dir}\err_sd_slope_{kind}.csv")
        err_meter_sims = pd.read_csv(fr"{main_dir}\err_meter_slope_{kind}.csv")
        para_sims = pd.read_csv(fr"{main_dir}\para_slope_{kind}.csv")

        ori_depth_bed_wse_df = pd.read_csv(fr"{main_dir}\ori_depth_bed_wse.csv")
        ori_para_df = pd.read_csv(fr"{main_dir}\original_para_slope_mean_2.0km.csv")

        if depth_calculation:
            depth_sims = pd.read_csv(fr"{main_dir}\depth_slope_{kind}.csv")
        else:
            pass

    else:
        err_sd_sims = pd.read_csv(fr"{main_dir}\err_sd_flow_{kind}.csv")
        err_meter_sims = pd.read_csv(fr"{main_dir}\err_meter_flow_{kind}.csv")
        para_sims = pd.read_csv(fr"{main_dir}\para_flow_{kind}.csv")

        ori_depth_bed_wse_df = pd.read_csv(fr"{main_dir}\ori_depth_bed_wse.csv")
        ori_para_df = pd.read_csv(fr"{main_dir}\original_para_flow.csv")

        if depth_calculation:
            depth_sims = pd.read_csv(fr"{main_dir}\depth_flow_{kind}.csv")
        else:
            pass

    # Data for visualisation
    coord_path = glob.glob(fr"{main_dir}\coord*")
    coord_df = pd.read_csv(coord_path[0])
    err_sd_sims_copy_df = statistic_df_generation(err_sd_sims)
    err_meter_sims_copy_df = statistic_df_generation(err_meter_sims)
    para_sims_copy_df = statistic_df_generation(para_sims)
    if depth_calculation:
        depth_sims_copy_df = statistic_df_generation(depth_sims)
    else:
        pass

    # Fontsize
    fontsize = 12

    # Visualise figure 1 for error in sd ---------
    fig1, ax1 = plt.subplots(figsize=(15, 6))

    # Multiple lines
    for i in range(1, 51, 1):
        ax1.plot(coord_df[0:680].distance, err_sd_sims[0:680][f'sim{i}'],
                 color='grey', linewidth=.4, alpha=.6, label='Simulated errors in sd' if i == 1 else "")

    # Max line
    ax1.plot(coord_df[0:680].distance, err_sd_sims_copy_df[0:680]['max'],
            color='black', linewidth=1, label='Max & Min of errors in sd')
    # Min line
    ax1.plot(coord_df[0:680].distance, err_sd_sims_copy_df[0:680]['min'],
            color='black', linewidth=1)
    # Mean line
    ax1.plot(coord_df[0:680].distance, err_sd_sims_copy_df[0:680]['mean'],
            color='black', linewidth=1, linestyle='--', label='Mean of errors in sd')

    ax1.legend(frameon=True, fontsize=fontsize, bbox_to_anchor=(.65, 1.2), facecolor='grey', framealpha=.2)
    ax1.set_xlabel('Distance (m)', fontsize=fontsize, labelpad=fontsize + 5)
    ax1.set_ylabel('WIDTH - Error in sd (m)', fontsize=fontsize, labelpad=fontsize + 2)

    # Save fig1
    fig1.savefig(
        fr"{main_dir}\\error_in_sd.png",
        bbox_inches='tight', dpi=400
    )

    # Visualise figure 2 for error in meter ---------
    fig2, ax2 = plt.subplots(figsize=(15, 6))

    # Multiple lines
    for i in range(1, 51, 1):
        ax2.plot(coord_df[0:680].distance, err_meter_sims[0:680][f'sim{i}'],
                color='grey', linewidth=.4, alpha=.6, label='Simulated errors in meter' if i == 1 else "")

    # Max line
    ax2.plot(coord_df[0:680].distance, err_meter_sims_copy_df[0:680]['max'],
            color='black', linewidth=1, label='Max & Min of errors in meter')
    # Min line
    ax2.plot(coord_df[0:680].distance, err_meter_sims_copy_df[0:680]['min'],
            color='black', linewidth=1)
    # Mean line
    ax2.plot(coord_df[0:680].distance, err_meter_sims_copy_df[0:680]['mean'],
            color='black', linewidth=1, linestyle='--', label='Mean of errors in meter')

    ax2.legend(frameon=True, fontsize=fontsize, bbox_to_anchor=(.7, 1.2), facecolor='grey', framealpha=.2)
    ax2.set_xlabel('Distance (m)', fontsize=fontsize, labelpad=fontsize + 5)
    ax2.set_ylabel('WIDTH - Error in meter (m)', fontsize=fontsize, labelpad=fontsize + 2)

    # Save fig2
    fig2.savefig(
        fr"{main_dir}\\error_in_meter.png",
        bbox_inches='tight', dpi=400
    )

    # Visualise figure 3 for para -----
    fig3, ax3 = plt.subplots(figsize=(15, 6))
    plot_name = para_name.upper()

    # Multiple lines
    for i in range(1, 51, 1):
        ax3.plot(coord_df[0:680].distance, para_sims[0:680][f'sim{i}'],
                color='grey', linewidth=.4, alpha=.6, label=f'Simulated {plot_name}' if i == 1 else "")

    # Max line
    ax3.plot(coord_df[0:680].distance, para_sims_copy_df[0:680]['max'],
            color='black', linewidth=1, label=f'Max & Min of {plot_name}')
    # Min line
    ax3.plot(coord_df[0:680].distance, para_sims_copy_df[0:680]['min'],
            color='black', linewidth=1)
    # Original line
    ax3.plot(coord_df[0:680].distance, ori_para_df[0:680]['z'],
            color='deepskyblue', linewidth=2.5, label=f'Original {plot_name}')
    # Mean line
    ax3.plot(coord_df[0:680].distance, para_sims_copy_df[0:680]['mean'],
            color='green', linewidth=1, linestyle='--', label=f'Mean of {plot_name}')

    ax3.legend(frameon=True, fontsize=fontsize, bbox_to_anchor=(.65, 1.27), facecolor='grey', framealpha=.2)
    ax3.set_xlabel('Distance (m)', fontsize=fontsize, labelpad=fontsize + 5)
    ax3.set_ylabel(f'{plot_name} (m)', fontsize=fontsize, labelpad=fontsize + 2)

    # Save fig3
    fig3.savefig(
        fr"{main_dir}\\{plot_name}_simulated_lines.png",
        bbox_inches='tight', dpi=400
    )

    # Visualise figure 4 for para - zooming-in version -------------
    fig4, ax4 = plt.subplots(figsize=(15, 6))

    # Original line
    ax4.plot(coord_df[ind_for_zoom[0]:ind_for_zoom[1]].distance,
             ori_para_df[ind_for_zoom[0]:ind_for_zoom[1]]['z'],
             color='deepskyblue', linewidth=2.5, label=f'Original {plot_name}')
    # Mean line
    ax4.plot(coord_df[ind_for_zoom[0]:ind_for_zoom[1]].distance,
             para_sims_copy_df[ind_for_zoom[0]:ind_for_zoom[1]]['mean'],
             color='green', linewidth=1, linestyle='--', label=f'Mean of {plot_name}')

    ax4.legend(frameon=True, fontsize=fontsize, bbox_to_anchor=(.65, 1.15), facecolor='grey', framealpha=.2)
    ax4.set_xlabel('Distance (m)', fontsize=fontsize, labelpad=fontsize + 5)
    ax4.set_ylabel(f'{plot_name} (m)', fontsize=fontsize, labelpad=fontsize + 2)

    # Save fig4
    fig4.savefig(
        fr"{main_dir}\\{plot_name}_simulated_lines_ZOOM.png",
        bbox_inches='tight', dpi=400
    )

    # Visualise figure 5 for para - AREA --------------
    fig5, ax5 = plt.subplots(figsize=(15, 6))

    # Calculate upper and lower
    mean_arr_para = para_sims_copy_df[0:680]['mean']
    sd_arr_para = para_sims_copy_df[0:680]['sd']
    upper_para = mean_arr_para + sd_arr_para
    lower_para = mean_arr_para - sd_arr_para

    # Set up plots
    ax5.fill_between(coord_df[0:680].distance, lower_para, upper_para, color='deepskyblue',
                    alpha=.5, edgecolor='green', linewidth=1.5, label='Standard deviation area')
    # ax5.errorbar(coord_df[0:680].distance, mean_arr_para, sd_arr_para,
    #             capsize=4, capthick=2, color='deepskyblue', zorder=0)
    # ax5.scatter(coord_df[0:680].distance, mean_arr_para,
    #            facecolor='green', s=5, zorder=1)
    ax5.plot(coord_df[0:680].distance, mean_arr_para,
            color='green', zorder=1, label=f'Mean of {plot_name}')

    ax5.legend(frameon=True, fontsize=fontsize, bbox_to_anchor=(.65, 1.15), facecolor='grey', framealpha=.2)
    ax5.set_xlabel('Distance (m)', fontsize=fontsize, labelpad=fontsize + 5)
    ax5.set_ylabel(f'{plot_name} (m)', fontsize=fontsize, labelpad=fontsize + 2)

    # Save fig5
    fig5.savefig(
        fr"{main_dir}\\{plot_name}_simulated_lines_AREA.png",
        bbox_inches='tight', dpi=400
    )

    if depth_calculation:
        # Visualise figure 6 for depth -------------------
        fig6, ax6 = plt.subplots(figsize=(15, 6))

        # Depth -----------------------------------------------------------------
        for i in range(1, 51, 1):
            ax6.plot(coord_df[0:680].distance, depth_sims[0:680][f'sim{i}'],
                    color='grey', linewidth=.4, alpha=.4, label="Simulated depths" if i == 1 else "")
        # Min
        ax6.plot(coord_df[0:680].distance, depth_sims_copy_df[0:680]['min'],
                color='black', linewidth=1, label='Min & max of depths')
        # Max
        ax6.plot(coord_df[0:680].distance, depth_sims_copy_df[0:680]['max'],
                color='black', linewidth=1, zorder=51)
        # Original
        ax6.plot(coord_df[0:680].distance, ori_depth_bed_wse_df[0:680]['depth'],
                color='deepskyblue', linewidth=2.5, label='Original depth')
        # Mean
        ax6.plot(coord_df[0:680].distance, depth_sims_copy_df[0:680]['mean'],
                color='green', linewidth=1.5, linestyle='--', label='Mean of depths')

        ax6.legend(frameon=True, fontsize=fontsize, bbox_to_anchor=(.65, 1.25), facecolor='grey', framealpha=.2)
        ax6.set_xlabel('Distance (m)', fontsize=fontsize, labelpad=fontsize + 5)
        ax6.set_ylabel('DEPTHS (m)', fontsize=fontsize, labelpad=fontsize + 2)

        # Save fig6
        fig6.savefig(
            fr"{main_dir}\\depths_simulated_lines.png",
            bbox_inches='tight', dpi=400
        )

        # Visualise figure 7 for depth - zooming-in version -------
        fig7, ax7 = plt.subplots(figsize=(15, 6))

        # Original
        ax7.plot(coord_df[ind_for_zoom[0]:ind_for_zoom[1]].distance,
                 ori_depth_bed_wse_df[ind_for_zoom[0]:ind_for_zoom[1]]['depth'],
                 color='deepskyblue', linewidth=2.5, label='Original depth')
        # Mean
        ax7.plot(coord_df[ind_for_zoom[0]:ind_for_zoom[1]].distance,
                 depth_sims_copy_df[ind_for_zoom[0]:ind_for_zoom[1]]['mean'],
                 color='green', linewidth=1.5, linestyle='--', label='Mean of depths')

        ax7.legend(frameon=True, fontsize=fontsize, bbox_to_anchor=(.65, 1.15), facecolor='grey', framealpha=.2)
        ax7.set_xlabel('Distance (m)', fontsize=fontsize, labelpad=fontsize + 5)
        ax7.set_ylabel('DEPTHS (m)', fontsize=fontsize, labelpad=fontsize + 2)

        # Save fig7
        fig7.savefig(
            fr"{main_dir}\\depths_simulated_lines_ZOOM.png",
            bbox_inches='tight', dpi=400
        )

        # Visualise figure 8 for depth - AREA
        fig8, ax8 = plt.subplots(figsize=(15, 6))

        # Calculate upper and lower
        mean_arr_depth = depth_sims_copy_df[0:680]['mean']
        sd_arr_depth = depth_sims_copy_df[0:680]['sd']
        upper_depth = mean_arr_depth + sd_arr_depth
        lower_depth = mean_arr_depth - sd_arr_depth

        # Set up plots
        ax8.fill_between(coord_df[0:680].distance, lower_depth, upper_depth, color='deepskyblue',
                        alpha=.5, edgecolor='green', linewidth=1.5, label='Standard deviation area')
        # ax8.errorbar(coord_df[0:680].distance, mean_arr_depth, sd_arr_depth,
        #             capsize=4, capthick=2, color='deepskyblue', zorder=0)
        # ax8.scatter(coord_df[0:680].distance, mean_arr_depth,
        #            facecolor='green', s=5, zorder=1)
        ax8.plot(coord_df[0:680].distance, mean_arr_depth,
                color='green', zorder=1, label='Mean of depths')

        ax8.legend(frameon=True, fontsize=fontsize, bbox_to_anchor=(.65, 1.15), facecolor='grey', framealpha=.2)
        ax8.set_xlabel('Distance (m)', fontsize=fontsize, labelpad=fontsize + 5)
        ax8.set_ylabel('DEPTH (m)', fontsize=fontsize, labelpad=fontsize + 2)

        # Save fig8
        fig8.savefig(
            fr"{main_dir}\\depths_simulated_lines_AREA.png",
            bbox_inches='tight', dpi=400
        )

    else:
        pass

def plot_paras_all(kind, ind_for_zoom):

    # Main path
    main_dir = fr"{bathy_path}\\bathy_analyses"

    # Get data
    ori_depth_bed_wse_df = pd.read_csv(fr"{main_dir}\ori_depth_bed_wse.csv")
    depth_sims = pd.read_csv(fr"{main_dir}\depth_combination_{kind}.csv")
    coord_df = pd.read_csv(fr"{main_dir}\coord_combination_{kind}.csv")

    # Calculate statistical simulations
    depth_sims_copy_df = statistic_df_generation(depth_sims)

    # Fontsize
    fontsize = 12

    # Visualise paras
    plot_paras_each(kind, 'width', [600, 680], False)
    plot_paras_each(kind, 'slope', [600, 680], False)
    plot_paras_each(kind, 'flow', [600, 680], False)


    # Visualise figure 1 for depth -------------------
    fig1, ax1 = plt.subplots(figsize=(15, 6))

    # Multiple lines
    for i in range(1, 51, 1):
        ax1.plot(coord_df[0:680].distance, depth_sims[0:680][f'sim{i}'],
                 color='grey', linewidth=.4, alpha=.4, label="Simulated depths" if i == 1 else "")
    # Min
    ax1.plot(coord_df[0:680].distance, depth_sims_copy_df[0:680]['min'],
             color='black', linewidth=1, label='Min & max of depths')
    # Max
    ax1.plot(coord_df[0:680].distance, depth_sims_copy_df[0:680]['max'],
             color='black', linewidth=1, zorder=51)
    # Original
    ax1.plot(coord_df[0:680].distance, ori_depth_bed_wse_df[0:680]['depth'],
             color='deepskyblue', linewidth=2.5, label='Original depth')
    # Mean
    ax1.plot(coord_df[0:680].distance, depth_sims_copy_df[0:680]['mean'],
             color='green', linewidth=1.5, linestyle='--', label='Mean of depths')

    ax1.legend(frameon=True, fontsize=fontsize, bbox_to_anchor=(.65, 1.25), facecolor='grey', framealpha=.2)
    ax1.set_xlabel('Distance (m)', fontsize=fontsize, labelpad=fontsize + 5)
    ax1.set_ylabel('DEPTHS (m)', fontsize=fontsize, labelpad=fontsize + 2)

    # Save fig1
    fig1.savefig(
        fr"{main_dir}\\depths_simulated_lines.png",
        bbox_inches='tight', dpi=400
    )

    # Visualise figure 2 for depth - zooming-in version -------
    fig2, ax2 = plt.subplots(figsize=(15, 6))

    # Original
    ax2.plot(coord_df[ind_for_zoom[0]:ind_for_zoom[1]].distance,
             ori_depth_bed_wse_df[ind_for_zoom[0]:ind_for_zoom[1]]['depth'],
             color='deepskyblue', linewidth=2.5, label='Original depth')
    # Mean
    ax2.plot(coord_df[ind_for_zoom[0]:ind_for_zoom[1]].distance,
             depth_sims_copy_df[ind_for_zoom[0]:ind_for_zoom[1]]['mean'],
             color='green', linewidth=1.5, linestyle='--', label='Mean of depths')

    ax2.legend(frameon=True, fontsize=fontsize, bbox_to_anchor=(.65, 1.15), facecolor='grey', framealpha=.2)
    ax2.set_xlabel('Distance (m)', fontsize=fontsize, labelpad=fontsize + 5)
    ax2.set_ylabel('DEPTHS (m)', fontsize=fontsize, labelpad=fontsize + 2)

    # Save fig7
    fig2.savefig(
        fr"{main_dir}\\depths_simulated_lines_ZOOM.png",
        bbox_inches='tight', dpi=400
    )

    # Visualise figure 3 for depth - AREA ---------------------
    fig3, ax3 = plt.subplots(figsize=(15, 6))

    # Calculate upper and lower
    mean_arr_depth = depth_sims_copy_df[0:680]['mean']
    sd_arr_depth = depth_sims_copy_df[0:680]['sd']
    upper_depth = mean_arr_depth + sd_arr_depth
    lower_depth = mean_arr_depth - sd_arr_depth

    # Set up plots
    ax3.fill_between(coord_df[0:680].distance, lower_depth, upper_depth, color='deepskyblue',
                     alpha=.5, edgecolor='green', linewidth=1.5, label='Standard deviation area')
    # ax3.errorbar(coord_df[0:680].distance, mean_arr_depth, sd_arr_depth,
    #             capsize=4, capthick=2, color='deepskyblue', zorder=0)
    # ax3.scatter(coord_df[0:680].distance, mean_arr_depth,
    #            facecolor='green', s=5, zorder=1)
    ax3.plot(coord_df[0:680].distance, mean_arr_depth,
             color='green', zorder=1, label='Mean of depths')

    ax3.legend(frameon=True, fontsize=fontsize, bbox_to_anchor=(.65, 1.15), facecolor='grey', framealpha=.2)
    ax3.set_xlabel('Distance (m)', fontsize=fontsize, labelpad=fontsize + 5)
    ax3.set_ylabel('DEPTH (m)', fontsize=fontsize, labelpad=fontsize + 2)

    # Save fig3
    fig3.savefig(
        fr"{main_dir}\\depths_simulated_lines_AREA.png",
        bbox_inches='tight', dpi=400
    )


