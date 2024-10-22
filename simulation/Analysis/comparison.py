# Prepare packages -----------------------------------------------------------------------------------------------------
# Packages for data/arrays manipulation
import pandas as pd                                     # For reading data
import numpy as np                                      # For handling data (using clip function)


from statisticCalculation import calculation_dict       # For generating statistical dictionary
from scipy.stats import gaussian_kde, skew              # For kde calculation and skewness calculation
from mpl_toolkits.axes_grid1 import make_axes_locatable

# For plotting
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mc
import seaborn as sns
import colorsys

# ----------------------------------------------------------------------------------------------------------------------

# GET DATA #############################################################################################################

def get_datalist(
    list_filename,
    list_resolution,
    building_path,
    flood_rate
):
    """
    @Definition:
                A function to call all statistic data from csv files and store them into a list (a list of four
                dictionaries)
    @References:
                None.
    @Arguments:
                list_filename (list):
                                A list of four strings represent for four paths to csv files
                list_resolution (list):
                                A list of four resolutions
                building_path (string):
                                Path of file containing building polygons
                flood_rate (float):
                                A rate where the flood was defined

    @Returns:
                None.
    """


    # Create a list contains four dictionaries (mean, sd, cv, and cell)
    dictionary_data_list = []

    # Loop to get dataframe of each type of data
    for i in range(len(list_filename)):
        # Get data from csv file
        statdata_df = pd.read_csv(fr"{list_filename[i]}\\5_analysis\\wd\\bathy_csv\\all_simulations.csv")

        if len(list_resolution) != 1:
            # Get statistic data if there are many resolution values
            statdata = calculation_dict(
                statdata_df,
                list_resolution[i],
                building_path,
                flood_rate, 'out.max', -9999,
                False,
                list_filename[i],
                [
                    fr"{list_filename[i]}\\5_analysis\\wse\\bathy_csv",
                    fr"{list_filename[i]}\\5_analysis\\wd\\bathy_csv",
                    fr"{list_filename[i]}\\5_analysis\\elev\\bathy_csv"
                ]
            )

        else:
            # Get data if there is only one resolution
            statdata = calculation_dict(
                statdata_df,
                list_resolution[0],
                building_path,
                flood_rate, 'out.max', -9999,
                False,
                list_filename[i],
                [
                    fr"{list_filename[i]}\\5_analysis\\wse\\bathy_csv",
                    fr"{list_filename[i]}\\5_analysis\\wd\\bathy_csv",
                    fr"{list_filename[i]}\\5_analysis\\elev\\bathy_csv"
                ]
            )

        # Add into list
        dictionary_data_list.append(statdata)

    return dictionary_data_list


def statistic_dataframe(list_dataname, list_dataframe, calculation_option):
    """
    @Definition:
                A function to get the data for each statistics and store in a dataframe
    @References:
                None.
    @Arguments:
                list_dataname (list):
                                A list contains name of data
                list_dataframe (list):
                                A list contains dataframe according to name of data
                calculation_option (string):
                                "mean" means to calculate mean
                                "sd" means to calculate standard deviation
                                "cv" means to calculate coefficient of variation
                                "cell" means to calculate proportion of simulations of each cell being inundated
                                "area" means to calculate flooded areas
                                "building" means to calculate flooded buildings
                                "rmse" means to calculate RMSEs
    @Returns:
                comparison_dictionary (dictionary):
                                A dictionary contains all dataframes of type of data
    """
    # Get empty dictionary
    comparison_dictionary = dict()

    # Add keys and values into dictionary
    for i in range(len(list_dataname)):
        if calculation_option not in ["area", "building", "rmse", "mae"]:
            # Get dataframe of mean, sd, cv, and cell
            comparison_dictionary[f'{list_dataname[i]}'] = list_dataframe[i][calculation_option][
                calculation_option][list_dataframe[i][calculation_option][calculation_option] != -9999]
        else:
            # Get dataframe of area and building
            comparison_dictionary[f'{list_dataname[i]}'] = list_dataframe[i][calculation_option].iloc[0]

    return comparison_dictionary

def statistic_df_dict(list_dataname, list_dataframe):
    """
    @Definition:
                A function to get the data for each statistics and store in a dataframe
    @References:
                None.
    @Arguments:
                list_dataname (list):
                                A list contains name of data
                list_dataframe (list):
                                A list contains dataframe according to name of data
    @Returns:
                comparison_dictionary (dictionary):
                                A dictionary contains all dataframes of type of data
    """
    # Get an empty dictionary
    statistic_df_dictionary = dict()

    # Statistical selections
    stat_selection = ['mean', 'sd', 'cv', 'cell', 'area', 'building', 'rmse', 'mae']

    # Generate statistical selection and add into a dictionary
    for each_option in stat_selection:
        stat_selection_df = pd.DataFrame(statistic_dataframe(
            list_dataname,
            list_dataframe,
            each_option
        ))

        statistic_df_dictionary[f'{each_option}'] = stat_selection_df

    return statistic_df_dictionary

# END GET DATA #########################################################################################################


# PLOT COMPARISON ######################################################################################################

def lighten_color(color, amount=0.5):
    """
    @Definition:
                A function to lighten the color (by @IanHincks)
    @References:
                https://stackoverflow.com/questions/37765197/darken-or-lighten-a-color-in-matplotlib/49601444#49601444
                https://stackoverflow.com/questions/55656683/change-seaborn-boxplot-line-rainbow-color
    @Arguments:
                color (matplotlib color):
                                A tuple of 3 values of RGB color
                amount (float):
                                Level of brightness of color
    @Returns:
                A new tuple of 3 values of RGB color
    """
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])

def boxplotting(
    figsize,
    stat_df,
    y_label,
    calculation_option,
    violin=False,
    palette='husl'
):
    """
    @Definition:
                A function to plot boxplot (based on @JohanC)
    @References:
                https://stackoverflow.com/questions/55656683/change-seaborn-boxplot-line-rainbow-color
                https://stackoverflow.com/questions/36874697/how-to-edit-properties-of-whiskers-fliers-caps-etc-in-seaborn-boxplot/36893152#36893152
                https://stackoverflow.com/questions/37765197/darken-or-lighten-a-color-in-matplotlib/49601444#49601444

                https://stackoverflow.com/questions/332289/how-do-i-change-the-size-of-figures-drawn-with-matplotlib
                https://stackoverflow.com/questions/14908576/how-to-remove-frame-from-matplotlib-pyplot-figure-vs-matplotlib-figure-frame
                https://stackoverflow.com/questions/70954489/how-to-use-column-and-meanprops-in-seaborn-boxplot-with-showmeans-true

                https://cduvallet.github.io/posts/2018/03/boxplots-in-python
                https://seaborn.pydata.org/generated/seaborn.boxplot.html
                https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.boxplot.html#matplotlib.axes.Axes.boxplot
                https://matplotlib.org/stable/gallery/pyplots/boxplot_demo_pyplot.html#sphx-glr-gallery-pyplots-boxplot-demo-pyplot-py
                https://stackoverflow.com/questions/35131798/tweaking-seaborn-boxplot

                https://www.python-graph-gallery.com/33-control-colors-of-boxplot-seaborn
                https://www.python-graph-gallery.com/32-custom-boxplot-appearance-seaborn
                https://stackoverflow.com/questions/46735745/how-to-control-scientific-notation-in-matplotlib
    @Arguments:
                figsize (tuple):
                            A tuple of figsize in matplotlib subplot (width, height)
                stat_df (pandas dataframe):
                            A dataframe contains sets of simulations as columns. Ex: The columns are 'north
                            translation', 'east translation', ...
                y_label (string):
                            Y label represents for the statistical name of the boxplots
                calculation_option (string):
                            Statistical options includes mean, sd, cv, cell, area, and building
                palette (string/dictionary):
                            A string which is a name listed in palette color library of seaborn. Ex: 'husl' - also
                            the default
                            A dictionary represents the colors for each column in the pandas dataframe. Ex: {'north
                            translation': 'orange', 'east translation': 'dodgerblue', ...}
    @Returns:
                A new tuple of 3 values of RGB color
    """
    # Set up axis
    fig, ax = plt.subplots(figsize=figsize)

    if violin:

        # Boxplot
        sns.boxplot(
            data=stat_df,
            orient='h',  # Boxplots lie horizontally
            # showmeans=True,  # Turn on mean sign
            # meanprops=dict(marker='s', markersize=6),  # Use big dot to visualise mean sign
            flierprops=dict(marker='o', markersize=1.5),  # Visualise outliers
            width=0.08,  # Size/width of boxplots
            boxprops={'zorder': 2},
            palette=sns.husl_palette(4), saturation=1, ax=ax
        )

        sns.violinplot(data=stat_df, palette=palette, orient='h', width=1,
                       inner=None, linewidth=0, saturation=1)

    else:

        # Boxplot
        sns.boxplot(
            data=stat_df,
            orient='h',  # Boxplots lie horizontally
            # showmeans=True,  # Turn on mean sign
            # meanprops=dict(marker='s', markersize=6),  # Use big dot to visualise mean sign
            flierprops=dict(marker='o', markersize=1.5),  # Visualise outliers
            width=0.4,  # Size/width of boxplots
            palette=sns.husl_palette(4), saturation=1, ax=ax
        )

    # ax.plot([], [], '|', linewidth=1, color='black', label='median')
    # ax.plot([], [], 's', linewidth=1, color='black', label='mean')
    # ax.legend(fontsize='x-large', frameon=False, loc='lower right')

    # Colorise all lines of boxplots
    box_patches = [
        patch for patch in ax.patches if type(patch) == matplotlib.patches.PathPatch
    ]
    num_patches = len(box_patches)
    lines_per_boxplot = len(ax.lines) // num_patches
    for i, patch in enumerate(box_patches):
        # Set the linecolor on the patch to the facecolor, and set the facecolor to None
        col = lighten_color(patch.get_facecolor(), 1.5)
        patch.set_edgecolor(col)

        # Each box has associated Line2D objects (to make the whiskers, fliers, etc.)
        # Loop over them here, and use the same color as above
        for line in ax.lines[i * lines_per_boxplot: (i + 1) * lines_per_boxplot]:
            line.set_color(col)
            line.set_mfc(col)  # facecolor of fliers
            line.set_mec(col)  # edgecolor of fliers
            line.set_linewidth(0.7)

    # Design x labels
    if calculation_option == 'mean':
        x_label = "Mean (m)"
    elif calculation_option == 'sd':
        x_label = "Standard deviation (m)"
    elif calculation_option == 'cv':
        x_label = "Coefficient of variation (%)"
    elif calculation_option == 'cell':
        x_label = "Proportion (%)"
    elif calculation_option == 'area':
        # x_label = r'Areas (x100 $\mathrm{m}^2$)'
        x_label = r'Areas (x100 $\mathrm{m}^2$)'

        ax.axvline(x=2250800, color='blue', linestyle='--', linewidth=1) # CMR
        ax.axvline(x=2254300, color='black', linestyle='--', linewidth=1) # UF

        ax.plot([], [], '--', linewidth=1, color='blue', label='RMSE $_{CMR,\ no-error\ parameters}$ = 2250800 $\mathrm{m}^2$') # CMR
        ax.plot([], [], '--', linewidth=1, color='black', label='RMSE $_{UF,\ no-error\ parameters}$ = 2254300 $\mathrm{m}^2$') # UF
        ax.legend(fontsize=13, frameon=False, loc='upper right')

    elif calculation_option == 'building':
        x_label = "Number of buildings"
    elif calculation_option == 'rmse':
        x_label = "RMSEs"
        ax.axvline(x=0.538637, color='blue', linestyle='--', linewidth=1) # CMR
        ax.axvline(x=0.524923, color='black', linestyle='--', linewidth=1) # UF

        ax.plot([], [], '--', linewidth=1, color='blue', label='RMSE $_{CMR,\ no-error\ parameters}$ = 0.539') # CMR
        ax.plot([], [], '--', linewidth=1, color='black', label='RMSE $_{UF,\ no-error\ parameters}$ = 0.525') # UF
        ax.legend(fontsize=15, frameon=False, loc='upper right')

    elif calculation_option == 'mae':
        x_label = "MAEs"
        ax.axvline(x=0.353, color='blue', linestyle='--', linewidth=1) # CMR
        ax.axvline(x=0.361, color='black', linestyle='--', linewidth=1) # UF

        ax.plot([], [], '--', linewidth=1, color='blue', label='MAE $_{CMR,\ no-error\ parameters}$ = 0.353') # CMR
        ax.plot([], [], '--', linewidth=1, color='black', label='MAE $_{UF,\ no-error\ parameters}$ = 0.361') # UF
        ax.legend(fontsize=15, frameon=False, loc='upper right')

    else:
        x_label = calculation_option

    # Fontsize
    fontsize = 20
    labelpad = 25

    # Adjust x and y labels
    ax.set_xlabel(x_label, fontsize=fontsize, labelpad=labelpad)
    ax.set_ylabel(y_label, rotation=-270, fontsize=fontsize, labelpad=labelpad+5)

    # Control scientific notation
    # Refer here for more information
    # https://stackoverflow.com/questions/46735745/how-to-control-scientific-notation-in-matplotlib
    if calculation_option == 'area':
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x / 100))))
        # ax.get_xaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
    else:
        pass

    # Design size and style for ticks and labels
    for item in (ax.get_xticklabels() + ax.get_yticklabels()):  # For x, y ticks' labels
        item.set_fontsize(fontsize-3)
    ax.tick_params(direction='out', length=5, pad=labelpad-17)

    # Control the frame plot
    fig.set_size_inches(figsize)

def boxplotting2(
    figsize,
    stat_df,
    y_label,
    calculation_option,
    palette='husl'
):
    """
    @Definition:
                A function to plot boxplot (based on @JohanC)
    @References:
                https://stackoverflow.com/questions/55656683/change-seaborn-boxplot-line-rainbow-color
                https://stackoverflow.com/questions/36874697/how-to-edit-properties-of-whiskers-fliers-caps-etc-in-seaborn-boxplot/36893152#36893152
                https://stackoverflow.com/questions/37765197/darken-or-lighten-a-color-in-matplotlib/49601444#49601444

                https://stackoverflow.com/questions/332289/how-do-i-change-the-size-of-figures-drawn-with-matplotlib
                https://stackoverflow.com/questions/14908576/how-to-remove-frame-from-matplotlib-pyplot-figure-vs-matplotlib-figure-frame
                https://stackoverflow.com/questions/70954489/how-to-use-column-and-meanprops-in-seaborn-boxplot-with-showmeans-true

                https://cduvallet.github.io/posts/2018/03/boxplots-in-python
                https://seaborn.pydata.org/generated/seaborn.boxplot.html
                https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.boxplot.html#matplotlib.axes.Axes.boxplot
                https://matplotlib.org/stable/gallery/pyplots/boxplot_demo_pyplot.html#sphx-glr-gallery-pyplots-boxplot-demo-pyplot-py
                https://stackoverflow.com/questions/35131798/tweaking-seaborn-boxplot

                https://www.python-graph-gallery.com/33-control-colors-of-boxplot-seaborn
                https://www.python-graph-gallery.com/32-custom-boxplot-appearance-seaborn
                https://stackoverflow.com/questions/46735745/how-to-control-scientific-notation-in-matplotlib
    @Arguments:
                figsize (tuple):
                            A tuple of figsize in matplotlib subplot (width, height)
                stat_df (pandas dataframe):
                            A dataframe contains sets of simulations as columns. Ex: The columns are 'north
                            translation', 'east translation', ...
                y_label (string):
                            Y label represents for the statistical name of the boxplots
                calculation_option (string):
                            Statistical options includes mean, sd, cv, cell, area, and building
                palette (string/dictionary):
                            A string which is a name listed in palette color library of seaborn. Ex: 'husl' - also
                            the default
                            A dictionary represents the colors for each column in the pandas dataframe. Ex: {'north
                            translation': 'orange', 'east translation': 'dodgerblue', ...}
    @Returns:
                A new tuple of 3 values of RGB color
    """
    # Set up axis
    fig, ax = plt.subplots(figsize=figsize)
    fig.subplots_adjust(wspace=.03)

    # Divide the plot into two
    ax1 = ax
    divider = make_axes_locatable(ax)
    ax2 = divider.new_horizontal(size='10%', pad=0.2)
    fig.add_axes(ax2)

    # BOXPLOT
    # Boxplot 1
    sns.boxplot(
        data=stat_df,
        orient='h',  # Boxplots lie horizontally
        # showmeans=True,  # Turn on mean sign
        # meanprops=dict(marker='s', markersize=6),  # Use big dot to visualise mean sign
        flierprops=dict(marker='o', markersize=1.5),  # Visualise outliers
        width=0.4,  # Size/width of boxplots
        palette=sns.husl_palette(4), saturation=1, ax=ax1
    )

    # Boxplot 2
    sns.boxplot(
        data=stat_df,
        orient='h',  # Boxplots lie horizontally
        # showmeans=True,  # Turn on mean sign
        # meanprops=dict(marker='s', markersize=6),  # Use big dot to visualise mean sign
        flierprops=dict(marker='o', markersize=1.5),  # Visualise outliers
        width=0.4,  # Size/width of boxplots
        palette=sns.husl_palette(4), saturation=1, ax=ax2
    )

    # Zoom-in / limit the view to different portions of the data
    ax1.xaxis.set_ticks(np.arange(0, 35, step=5))
    ax1.set_xlim(left=-3, right=34) # left=-3, right=34

    ax2.xaxis.set_ticks(np.arange(245, 258, step=5)[1:])
    ax2.set_xlim(left=251, right=258)  # left=104


    ax1.margins(x=0)

    # Hide the spines between ax1 and ax2
    ax1.spines.right.set_visible(False)
    ax2.spines.left.set_visible(False)
    # ax1.xaxis.tick_top()
    # ax1.tick_params(labeltop=False) # don't put tick labels at the top

    d = .9  # proportion of vertical to horizontal extent of the slanted line
    kwargs = dict(marker=[(-1, -d), (1, d)], markersize=15,
                  linestyle="none", color='k', mec='k', mew=2, clip_on=False)
    ax1.plot([1, 1], [1, 0], transform=ax1.transAxes, **kwargs)
    ax2.plot([0, 0], [0, 1], transform=ax2.transAxes, **kwargs)


    # ax.plot([], [], '|', linewidth=1, color='black', label='median')
    # ax.plot([], [], 's', linewidth=1, color='black', label='mean')
    # ax.legend(fontsize='x-large', frameon=False, loc='lower right')

    # Colorise all lines of boxplots
    # AXIS 1
    box_patches = [
        patch for patch in ax1.patches if type(patch) == matplotlib.patches.PathPatch
    ]
    num_patches = len(box_patches)
    lines_per_boxplot = len(ax1.lines) // num_patches
    for i, patch in enumerate(box_patches):
        # Set the linecolor on the patch to the facecolor, and set the facecolor to None
        col = lighten_color(patch.get_facecolor(), 1.5)
        patch.set_edgecolor(col)

        # Each box has associated Line2D objects (to make the whiskers, fliers, etc.)
        # Loop over them here, and use the same color as above
        for line in ax1.lines[i * lines_per_boxplot: (i + 1) * lines_per_boxplot]:
            line.set_color(col)
            line.set_mfc(col)  # facecolor of fliers
            line.set_mec(col)  # edgecolor of fliers
            line.set_linewidth(0.7)

    # AXIS 2
    box_patches = [
        patch for patch in ax2.patches if type(patch) == matplotlib.patches.PathPatch
    ]
    num_patches = len(box_patches)
    lines_per_boxplot = len(ax2.lines) // num_patches
    for i, patch in enumerate(box_patches):
        # Set the linecolor on the patch to the facecolor, and set the facecolor to None
        col = lighten_color(patch.get_facecolor(), 1.5)
        patch.set_edgecolor(col)

        # Each box has associated Line2D objects (to make the whiskers, fliers, etc.)
        # Loop over them here, and use the same color as above
        for line in ax2.lines[i * lines_per_boxplot: (i + 1) * lines_per_boxplot]:
            line.set_color(col)
            line.set_mfc(col)  # facecolor of fliers
            line.set_mec(col)  # edgecolor of fliers
            line.set_linewidth(0.7)

    # Design x labels
    if calculation_option == 'mean':
        x_label = "Mean (m)"
    elif calculation_option == 'sd':
        x_label = "Standard deviation (m)"
    elif calculation_option == 'cv':
        x_label = "Coefficient of variation (%)"
    elif calculation_option == 'cell':
        x_label = "Proportion (%)"
    elif calculation_option == 'area':
        # x_label = r'Areas (x100 $\mathrm{m}^2$)'
        x_label = r'Areas ($\mathrm{m}^2$)'
    elif calculation_option == 'building':
        x_label = "Number of buildings"
    elif calculation_option == 'rmse':
        x_label = "RMSEs"
    else:
        x_label = calculation_option

    # Fontsize
    fontsize = 20
    labelpad = 25

    # Adjust x and y labels
    # ax1.set_xlabel(x_label, fontsize=fontsize, labelpad=labelpad)
    fig.text(.5, -.01, x_label, fontsize=fontsize, ha='center', va='center')
    ax1.set_ylabel(y_label, rotation=-270, fontsize=fontsize, labelpad=labelpad+5)

    # ax2.set_xlabel(x_label, fontsize=fontsize, labelpad=labelpad)
    # ax2.set_ylabel(y_label, rotation=-270, fontsize=fontsize, labelpad=labelpad+5)
    ax2.set_yticks([])

    # Control scientific notation
    # Refer here for more information
    # https://stackoverflow.com/questions/46735745/how-to-control-scientific-notation-in-matplotlib
    if calculation_option == 'area':
        # ax.get_xaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x / 100))))
        ax1.get_xaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x))))
    else:
        pass

    # Design size and style for ticks and labels
    for item in (ax1.get_xticklabels() + ax1.get_yticklabels()):  # For x, y ticks' labels
        item.set_fontsize(fontsize-3)
    ax1.tick_params(direction='out', length=5, pad=labelpad-17)

    if calculation_option == 'area':
        # ax.get_xaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x / 100))))
        ax2.get_xaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x))))
    else:
        pass

    # Design size and style for ticks and labels
    for item in (ax2.get_xticklabels()):  # For x, y ticks' labels
        item.set_fontsize(fontsize-3)
    ax2.tick_params(direction='out', length=5, pad=labelpad-17)


    # Control the frame plot
    fig.set_size_inches(figsize)


def kdeplotting(
    figsize,
    stat_df,
    x_axis_range,
    legend_title,
    legend_location,
    calculation_option,
    comparison_sign,
    palette='husl',
    mean_line=False
):
    """
    @Definition:
                A function to plot probability density function plot
    @References:
                https://stackoverflow.com/questions/37765197/darken-or-lighten-a-color-in-matplotlib/49601444#49601444

                https://stackoverflow.com/questions/4150171/how-to-create-a-density-plot-in-matplotlib/4152016#4152016
                https://stackoverflow.com/questions/29083429/python-scipy-kde-fit-scaling
                https://matplotlib.org/3.3.4/gallery/showcase/integral.html

                https://stackoverflow.com/questions/66903255/retrieve-values-from-scipy-gaussian-kde

                https://seaborn.pydata.org/generated/seaborn.set_style.html
                https://stackoverflow.com/questions/25540259/remove-or-adapt-border-of-frame-of-legend-using-matplotlib
                https://seaborn.pydata.org/generated/seaborn.kdeplot.html (multiple='stack')


    @Arguments:
                figsize (tuple):
                            A tuple of figsize in matplotlib subplot (width, height)
                stat_df (pandas dataframe):
                            A dataframe contains sets of simulations as columns. Ex: The columns are 'north
                            translation', 'east translation', ...
                x_axis_range (list)
                            A list of 5 values include lower and upper limit, step, bandwith, and ylimit
                legend_title (stirng):
                            Title of legend
                legend_location (list):
                            A list of two values represents for the location of legend box
                calculation_option (string):
                            Statistical options includes mean, sd, cv, cell, area, building, and rmse
                comparison_sign (r string):
                            A comparison sign will be added into the upper limit of x axis. Ex: for larger
                            comparison, add this <r'$\geq $'>
                palette (string/dictionary):
                            A string which is a name listed in palette color library of seaborn. Ex: 'husl' - also
                            the default
                            A dictionary represents the colors for each column in the pandas dataframe. Ex: {'north
                            translation': 'orange', 'east translation': 'dodgerblue', ...}
                mean_line (boolean):
                            If False, mean line for each distribution will not be appeared
                            If True, mean line for each distribution will be appeared
    @Returns:
                A new tuple of 3 values of RGB color
    """
    # Set up background
    sns.set_style('ticks')

    # Set up axis
    fig, ax = plt.subplots(figsize=figsize)

    # Fontsize
    fontsize = 18
    labelpad = 23

    for number_column in range(stat_df.shape[1]):
        # Plot kde
        sns.kdeplot(
            np.clip(stat_df[stat_df.columns[number_column]].dropna(), x_axis_range[0], x_axis_range[1]),
            fill=True, linewidth=1.25,
            clip=(x_axis_range[0], x_axis_range[1]),
            alpha=0.12,
            bw_adjust=x_axis_range[3],
            color=lighten_color(sns.color_palette(palette)[number_column], 1.5), # transformation (0.9),
                                                                                 # resolution (1.1)
            label=stat_df.columns[number_column],
            # multiple='stack', # to get white curve line
            ax=ax
        )

        if mean_line:
            # Calculate kde values by Gaussian
            kde_gauss = gaussian_kde(stat_df.iloc[:, number_column].dropna().to_numpy(), bw_method='scott')

            # Calculate mean and height
            mean_gauss = stat_df.iloc[:, number_column].dropna().to_numpy().mean()
            height_gauss = kde_gauss(mean_gauss)

            # Plot mean line
            ymean_line = np.linspace(0, height_gauss, 200)
            xmean_line = [mean_gauss]*len(ymean_line)
            ax.plot(xmean_line, ymean_line,
                    ls="-.", linewidth=1,
                    color=lighten_color(sns.color_palette(palette)[number_column], 1))

    # Legend
    leg = ax.legend(ncol=2 if stat_df.shape[1] > 3 else 1,
                    prop={'size': fontsize-4}, title=legend_title,
                    title_fontproperties={'size': fontsize-3, 'weight': 'bold'},
                    framealpha=1,
                    loc=legend_location)

    # Set ylim
    ax.set_ylim(top=None if x_axis_range[4]==0 else x_axis_range[4])

    # Ref: https://stackoverflow.com/questions/25540259/remove-or-adapt-border-of-frame-of-legend-using-matplotlib
    leg.get_frame().set_facecolor('white')
    leg.get_frame().set_linewidth(0.0)

    # Design labels
    if calculation_option == 'mean':
        x_label = "Mean (m)"
    elif calculation_option == 'sd':
        x_label = "Standard deviation (m)"
    elif calculation_option == 'cv':
        x_label = "Coefficient of variation (%)"
    else:
        x_label = "Proportion (%)"

    # Generate title and labels
    ax.set_xlabel(x_label, fontsize=fontsize-1, labelpad=labelpad)
    ax.set_ylabel("Probability density", rotation=-270, fontsize=fontsize-1, labelpad=labelpad+4)

    # Remove grid background lines (including x, y lines)
    # ax.grid(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

    # Set x label, for more information, please visit
    # https://stackoverflow.com/questions/26218704/matplotlib-histogram-with-collection-bin-for-high-values
    # For dtype: https://numpy.org/doc/stable/reference/arrays.dtypes.html

    # Set tick
    x_range = np.arange(x_axis_range[0], x_axis_range[1]+0.1, x_axis_range[2])
    xlabel_arr = np.array(np.round(x_range[:], 2), dtype='float')
    ax.set_xticks(xlabel_arr)

    # Set tick label
    # Refer here for information about signs of larger/less or equal to
    # https://matplotlib.org/stable/tutorials/text/mathtext.html
    # https://matplotlib.org/2.0.2/users/mathtext.html
    # https://stackoverflow.com/questions/18428823/python-matplotlib-less-than-or-equal-to-symbol-in-text
    # https://stackoverflow.com/questions/48646669/make-matplotlib-display-serif-font-in-math-mode
    if calculation_option != 'cell':
        xlabel_arr_cop = np.array(np.round(x_range[:], 2), dtype='str')
        xlabel_arr_cop[-1] = '{1}{0:.2f}'.format(xlabel_arr[-1], comparison_sign)
        ax.set_xticklabels(xlabel_arr_cop)
    else:
        xlabel_arr_cop = np.array(np.round(x_range[:], 2), dtype='str')
        xlabel_arr_cop[0] = '{1}{0:.2f}'.format(xlabel_arr[0], comparison_sign)
        ax.set_xticklabels(xlabel_arr_cop)

    # Design size and style for ticks and labels
    for item in (ax.get_xticklabels() + ax.get_yticklabels()):  # For x, y ticks' labels
        item.set_fontsize(fontsize-3)
    ax.tick_params(direction='out', length=5, pad=labelpad-17)


def comparison_calculation(stat_df, choice_calculation):
    """
    References:
                https://www.itl.nist.gov/div898/handbook/eda/section3/eda35b.htm#:~:text=Kurtosis%20is%20a%20measure%20of,would%20be%20the%20extreme%20case.
    """
    # List of results
    results = []

    # Get results of calculation
    for i in range(stat_df.shape[1]):
        text = "{:<24}: {:>5.4f}"
        q1 = np.nanpercentile(stat_df[stat_df.columns[i]], 25, method='hazen')
        q2 = np.nanpercentile(stat_df[stat_df.columns[i]], 50, method='hazen')
        q3 = np.nanpercentile(stat_df[stat_df.columns[i]], 75, method='hazen')

        if choice_calculation == 'galton skewness':
            galton_skewness = (q1 + q3 - (2 * q2))/(q3 - q1)
            # print(text.format(stat_df.columns[i], galton_skewness))
            results.append(galton_skewness)
        elif choice_calculation == 'quartile dev':
            quartile_dev = (q3 - q1)/2
            # print(text.format(stat_df.columns[i], quartile_dev))
            results.append(quartile_dev)
        elif choice_calculation == 'coefficient quartile dev':
            coef_quartile_dev = (q3 - q1)/(q3 + q1)
            # print(text.format(stat_df.columns[i], coef_quartile_dev))
            results.append(coef_quartile_dev)
        elif choice_calculation == 'quartile median':
            quartile_median = (q3 - q1)/q2
            # print(text.format(stat_df.columns[i], quartile_median))
            results.append(quartile_median)
        elif choice_calculation == 'fisher-pearson skewness':
            # print(text.format(stat_df.columns[i], skew(stat_df[stat_df.columns[i]].to_numpy(), nan_policy='omit')))
            results.append(skew(stat_df[stat_df.columns[i]].to_numpy(), nan_policy='omit'))
        elif choice_calculation == 'mean':
            # print(text.format(stat_df.columns[i], np.nanmean(stat_df[stat_df.columns[i]].to_numpy())))
            results.append(np.nanmean(stat_df[stat_df.columns[i]].to_numpy()))
        elif choice_calculation == 'sd':
            # print(text.format(stat_df.columns[i], np.nanstd(stat_df[stat_df.columns[i]].to_numpy())))
            results.append(np.nanstd(stat_df[stat_df.columns[i]].to_numpy()))
        elif choice_calculation == 'range':
            results.append(
                (stat_df.max(axis=0) - stat_df.min(axis=0))
            )
        elif choice_calculation == 'median':
            results.append(
                stat_df.median(axis=0)
            )
        else:
            # print(text.format(stat_df.columns[i], np.nanstd(stat_df[stat_df.columns[i]].to_numpy())/np.nanmean(stat_df[stat_df.columns[i]].to_numpy())))
            results.append(np.nanstd(stat_df[stat_df.columns[i]].to_numpy())/np.nanmean(stat_df[stat_df.columns[i]].to_numpy()))

    return results






