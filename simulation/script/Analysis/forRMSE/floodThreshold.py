# Prepare packages -----------------------------------------------------------------------------------------------------
from folder import *

# For data manipulation
import numpy as np                                                  # For data array manipulation

# For plotting
import matplotlib                                                   # For ticker control
import matplotlib.pyplot as plt                                     # For main/parent axis

# Packages for variation calculation
from statisticCalculation import calculation_dict                   # For calculating statistics
from impactCalculation import water_raster_nobackground, \
                              watermap_onepolygon_parallelism   # For re-generate and re-join polygons
# ----------------------------------------------------------------------------------------------------------------------


def flood_rate_comparison(
    flood_rate_range,
    flowdepth_csv_df,
    num_processes,
    resolution,
    building_path
):
    """
    @Definition:
                A function to calculate the change of areas/ buildings through different flood rates
    @References:
                https://stackoverflow.com/questions/6081008/dump-a-numpy-array-into-a-csv-file
                https://stackoverflow.com/questions/46614526/how-to-import-a-csv-file-into-a-data-array
    @Arguments:
                flood_rate_range (list):
                                        A list of 3 values includes lower, upper limits and step
                flowdepth_csv_df (pandas dataframe):
                                        Flowdepth dataframe which was written to csv file
                num_processes (int):
                                        A number of process for the parallelism
                resolution (int):
                                        Resolution value in meter
                building_path (string):
                                        Path of file containing building polygons
    @Returns:
                mean_building_arr, sd_building_arr, mean_area_arr, sd_area_arr (numpy array):
                                        Mean and sd arrays of areas/ buildings through different simulations through
                                        different flood rates
    """
    # Get flood rate
    flood_rates = np.arange(flood_rate_range[0], flood_rate_range[1], flood_rate_range[2])

    # Get empty mean and standard deviation lists
    # Building
    mean_building = []
    sd_building = []
    # Area
    mean_area = []
    sd_area = []

    # Collect lists of means and standard deviations
    for each_threshold in flood_rates:
        # Do impact calculation
        water_raster_nobackground(flowdepth_csv_df, each_threshold)
        watermap_onepolygon_parallelism(flowdepth_csv_df, each_threshold, num_processes)

        # Get dictionary of all statistics
        stat_dict = calculation_dict(
            flowdepth_csv_df, resolution,
            building_path,
            each_threshold,
            False
        )

        # Get means and standard deviations of areas/buildings
        # Building
        mean_building.append(stat_dict['building'].mean(axis=1)[0])
        sd_building.append(stat_dict['building'].std(axis=1)[0])
        # Area
        mean_area.append(stat_dict['area'].mean(axis=1)[0])
        sd_area.append(stat_dict['area'].std(axis=1)[0])

    # Convert to numpy array
    # Building
    mean_building_arr = np.array(mean_building).astype('float64')
    sd_building_arr = np.array(sd_building).astype('float64')
    # Area
    mean_area_arr = np.array(mean_area).astype('float64')
    sd_area_arr = np.array(sd_area).astype('float64')

    # Write into csv
    # Building
    np.savetxt(fr"{other_untransformation}\\mean_building_array.csv", mean_building_arr, delimiter=",")
    np.savetxt(fr"{other_untransformation}\\sd_building_array.csv", sd_building_arr, delimiter=",")
    # Area
    np.savetxt(fr"{other_untransformation}\\mean_area_array.csv", mean_area_arr, delimiter=",")
    np.savetxt(fr"{other_untransformation}\\sd_area_array.csv", sd_area_arr, delimiter=",")

    return mean_building_arr, sd_building_arr, mean_area_arr, sd_area_arr



def flood_rate_plotting(
    figsize,
    flood_rate_range,
    mean_arr, sd_arr, sd_rate,
    calculation_option,
    extract_name,
    color_list=['fuchsia', 'red', 'white', 'maroon']
):
    """
    @Definition:
                A function to plot probability density function plot
    @References:
                https://stackoverflow.com/questions/40292875/change-capstyle-for-errorbars-in-matplotlib
                https://stackoverflow.com/questions/61799259/plot-smoothing-matplotlib-and-seaborn
                https://stackoverflow.com/questions/7601334/how-to-set-the-line-width-of-error-bar-caps
                https://stackoverflow.com/questions/12957582/plot-yerr-xerr-as-shaded-region-rather-than-error-bars
                https://stackoverflow.com/questions/22481854/plot-mean-and-standard-deviation
    @Arguments:
                figsize (tuple):
                            A tuple of figsize in matplotlib subplot (width, height)
                flood_rate_range (list):
                            A list of 4 values includes lower, upper limits, step for them, and step for x axis
                mean_arr, sd_arr (numpy array):
                            Mean and sd arrays of areas/ buildings through different simulations through
                            different flood rates
                sd_rate (int):
                            Rate to get range of standard deviation. Ex: sd_rate = 2 means 2 * standard deviation
                calculation_option (string):
                            Statistical options includes mean, sd, cv, cell, area, and building
                extract_name (string):
                            Name of a specific output among all flood model outputs
                color_list (list):
                            A list of string represents 3 names of colors for <fill_between>, <errorbar> and <plot>,
                            and <scatter> functions
    @Returns:
                None
    """
    # Choose outputs
    if extract_name == 'out.max':
        plot_untransformation = wd_plot_untransformation
    else:
        plot_untransformation = wse_plot_untransformation

    # Set up axis
    fig, ax = plt.subplots(figsize=figsize)

    # Get flood rates
    flood_rates = np.arange(flood_rate_range[0], flood_rate_range[1], flood_rate_range[2])

    # Get upper and lower
    upper = mean_arr + sd_arr * sd_rate
    lower = mean_arr - sd_arr * sd_rate

    # Set up plots
    # Plot area of errors
    ax.fill_between(flood_rates, lower, upper, color=color_list[0], alpha=0.2, zorder=0)
    # Plot error bar
    ax.errorbar(flood_rates, mean_arr, sd_arr * sd_rate, linestyle='none', capsize=4, capthick=2, color=color_list[1],
                zorder=1)
    # Plot line
    # ax.plot(flood_rates, mean_arr, color=color_list[1], zorder=2)
    # Plot scatter
    ax.scatter(flood_rates, mean_arr, edgecolor=color_list[2], facecolor=color_list[3], s=50, zorder=3)

    # Fontsize
    fontsize = 17
    labelpad = 25

    # Design y labels
    if calculation_option == 'area':
        y_label = r'Areas (x100 $\mathrm{m}^2$)'
        ax.get_yaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x / 100))))
    else:
        y_label = "Number of buildings"

    # Adjust x and y labels
    ax.set_xlabel("Flood depth thresholds (m)", fontsize=fontsize, labelpad=labelpad)
    ax.set_ylabel(y_label, rotation=-270, fontsize=fontsize, labelpad=labelpad+5)

    # Set up ticks
    for item in (ax.get_xticklabels() + ax.get_yticklabels()):  # For x, y ticks' labels
        item.set_fontsize(fontsize - 3)
    ax.tick_params(direction='out', length=5, pad=labelpad - 17)

    # Set up x range
    xlabel_arr = np.array(np.round(flood_rates[::flood_rate_range[3]], 2), dtype='float')
    ax.set_xticks(xlabel_arr)

    # Save fig
    fig.savefig(
        fr"{plot_untransformation}\\threshold_{calculation_option}.png",
        bbox_inches='tight', dpi=330
    )
