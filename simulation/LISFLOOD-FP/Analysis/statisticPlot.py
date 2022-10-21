# Prepare packages -----------------------------------------------------------------------------------------------------
from folder import *

# For data manipulation
import numpy as np                                                  # For data array manipulation

# For plotting
import matplotlib                                                   # For importing ticker
import matplotlib.pyplot as plt                                     # For main/parent axis
from matplotlib.patches import Polygon                              # For plotting rainbow polygon underline curve
from mpl_toolkits.axes_grid1.inset_locator import TransformedBbox, \
                                                  BboxPatch         # For creating the box for inset plot
from matplotlib_scalebar.scalebar import ScaleBar                   # For creating scale bar
import contextily as ctx                                            # For basemap
import seaborn as sns                                               # For plotting density

# For color
from colorDevelopment import *


# For raster manipulation
import xarray as xr                                                 # For raster manipulation
# ----------------------------------------------------------------------------------------------------------------------

### VARIATION MAP PLOTTING #############################################################################################
# References
"""
# Get data from plot
https://stackoverflow.com/questions/43565892/python-seaborn-distplot-y-value-corresponding-to-a-given-x-value

# Get map
https://server.arcgisonline.com/arcgis/rest/services
https://contextily.readthedocs.io/en/latest/intro_guide.html

# How to adjust location of box inside plot
https://matplotlib.org/stable/api/legend_api.html

# How to specify order of matplotlib plot layers
https://stackoverflow.com/questions/37246941/specifying-the-order-of-matplotlib-layers
https://matplotlib.org/3.1.1/gallery/misc/zorder_demo.html
https://stackoverflow.com/questions/44361658/zorder-specification-in-matplotlib-patch-collections

# How to make the contour line disappear
https://stackoverflow.com/questions/8263769/hide-contour-linestroke-on-pyplot-contourf-to-get-only-fills

# How to get rid of 1e6
https://www.codegrepper.com/code-examples/python/how+to+get+rid+of+1e6+in+matplotlib+plots

# How to add arrow to the map
https://stackoverflow.com/questions/58088841/how-to-add-a-north-arrow-on-a-geopandas-map
https://stackoverflow.com/questions/27598976/matplotlib-unknown-property-headwidth-and-head-width
https://stackoverflow.com/questions/44655006/can-i-change-the-arrowprops-properties-of-an-annotation-in-an-animation
https://stackoverflow.com/questions/30649329/python-anotation-arrow-direction
https://matplotlib.org/stable/gallery/lines_bars_and_markers/broken_barh.html#sphx-glr-gallery-lines-bars-and-markers-broken-barh-py
https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.annotate.html

# How to add ScaleBar into the map
https://stackoverflow.com/questions/70256031/how-do-you-display-the-scale-in-meters-the-north-arrow-and-the-axes-in-latitude
https://stackoverflow.com/questions/39786714/how-to-insert-scale-bar-in-a-map-in-matplotlib

# How to control the smooth of colors in xarray
https://docs.xarray.dev/en/stable/user-guide/plotting.html

# How to change the proportion of colors in map
https://stackoverflow.com/questions/68379452/how-to-create-a-n-color-cmap-in-plt-py
https://stackoverflow.com/questions/14777066/matplotlib-discrete-colorbar
https://stackoverflow.com/questions/72651896/how-to-fill-intervals-under-kde-curve-with-different-colors

# How to add inset plots in parent plot
https://stackoverflow.com/questions/45076945/matplotlib-mark-inset-with-different-edges-for-axes
https://matplotlib.org/stable/gallery/subplots_axes_and_figures/zoom_inset_axes.html#sphx-glr-gallery-subplots-axes-and-figures-zoom-inset-axes-py
https://matplotlib.org/stable/gallery/axes_grid1/inset_locator_demo.html
https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html#matplotlib.patches.Patch.set_fill
https://matplotlib.org/stable/gallery/subplots_axes_and_figures/axes_zoom_effect.html#sphx-glr-gallery-subplots-axes-and-figures-axes-zoom-effect-py
https://pythonmatplotlibtips.blogspot.com/2017/12/draw-axes-in-axes-using-zoomed-inset-axes.html
https://stackoverflow.com/questions/17458580/embedding-small-plots-inside-subplots-in-matplotlib
https://towardsdatascience.com/mastering-inset-axes-in-matplotlib-458d2fdfd0c0
https://scipython.com/blog/inset-plots-in-matplotlib/
https://stackoverflow.com/questions/24477220/use-subplots-to-zoom-into-timeseries-or-how-i-can-draw-lines-outside-of-axis-bor
https://regenerativetoday.com/some-tricks-to-make-matplotlib-visualization-even-better/

# How to hide axis text in matplotlib
https://stackoverflow.com/questions/2176424/hiding-axis-text-in-matplotlib-plots

# How to fill gradient colors horizontally under the curve
https://stackoverflow.com/questions/22081361/pyplot-vertical-gradient-fill-under-curve

# How to fill gradient colors vertically under the curve (RAINBOW)
https://stackoverflow.com/questions/18215276/how-to-fill-rainbow-color-under-a-curve-in-python-matplotlib

# How to plot kernel density line (not whole plot, just like a drawing line) by seaborn, scikit, and scipy
https://stackoverflow.com/questions/68396403/kernel-density-estimation-using-scipys-gaussian-kde-and-sklearns-kerneldensity
https://stackoverflow.com/questions/4150171/how-to-create-a-density-plot-in-matplotlib
https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gaussian_kde.html
https://matplotlib.org/stable/gallery/showcase/integral.html#sphx-glr-gallery-showcase-integral-py
https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Polygon.html#matplotlib.patches.Polygon

# How to control the smooth of kde density line
https://stackoverflow.com/questions/68396403/kernel-density-estimation-using-scipys-gaussian-kde-and-sklearns-kerneldensity
https://stackoverflow.com/questions/32900854/how-to-smooth-a-line-using-gaussian-kde-kernel-in-python-setting-a-bandwidth

# How to clip data
https://numpy.org/doc/stable/reference/generated/numpy.clip.html
https://stackoverflow.com/questions/32492518/numpy-clip-function-takes-twice-as-long-as-expected

# How to move left axis to right axis
https://stackoverflow.com/questions/13369888/matplotlib-y-axis-label-on-right-side

# How to colorise the edge of plot or frame of plot
https://stackoverflow.com/questions/7778954/elegantly-changing-the-color-of-a-plot-frame-in-matplotlib
https://towardsdatascience.com/handling-plot-axis-spines-in-python-f143b8554da2

# How to use **kwargs in matplotlib
https://stackoverflow.com/questions/47997772/how-to-use-kwargs-in-matplotlib-axes-axes-arrow-python-2-7

# How to create hillshade or terrain shade
https://www.youtube.com/watch?v=5dDZeEXws9Q
https://blog.datawrapper.de/shaded-relief-with-gdal-python/
https://www.geophysique.be/2014/02/25/shaded-relief-map-in-python/
https://www.l3harrisgeospatial.com/docs/topographicshading.html
"""

def polygon(axes, x1, y1, x2, y2, c):
    """
    @Definition:
                A function to create polygon
    @Arguments:
                axes (matplotlib axis):
                            Axis of matplotlib subplot
                x1, y1, x2, y2 (float):
                            Coordinates of head (x1, y1) and tail (x2, y2) in matplotlib plot
                c (each color in cmap):
                            Each color in cmap
    @Returns:
                None
    """
    # Create polygon
    polygon = Polygon(
        [(x1, y1), (x2, y2), (x2, 0), (x1, 0)], color=c
    )

    # Add polygon into patch
    axes.add_patch(polygon)


def rainbow_fill_sns(axes, X, Y, top, cmap):
    """
    @Definition:
                A function to create polygon
    @Arguments:
                axes (matplotlib axis):
                            Axis of matplotlib subplot
                x1, y1, x2, y2 (float):
                            Coordinates of head (x1, y1) and tail (x2, y2) in matplotlib plot
                c (each color in cmap):
                            Each color in cmap
    @Returns:
                None
    """
    # Get the size of the data (how many values)
    N = float(X.size)

    # Plot each vertical color
    for n, (x, y) in enumerate(zip(X, Y)):
        # Make the bottom x axis start at 0 of y axis
        axes.set_ylim(top=top, bottom=0)

        # Each color for each value
        color = cmap(n/N)

        # Use polygon
        if n+1 == N: continue
        polygon(axes, x, y, X[n+1], Y[n+1], color)

def mark_inset_noconnectedline(parent_axes, inset_axes, zorder=None, **kwargs):
    """
    @Definition:
                A function to create inset plot represent a zooming place in map
    @Arguments:
                parent_axes (matplotlib axis):
                            Axis of matplotlib subplot (the original axis)
                inset_axes (matplotlib axis):
                            Axis of matplotlib inset plot
                zorder (int):
                            Oridinal number of plot layer
                **kwargs:
                            Any other arguments
    @Returns:
                None
    """
    # Get the coordinates of the box based on the parent axes
    rect = TransformedBbox(inset_axes.viewLim, parent_axes.transData)

    # Create the box
    pp = BboxPatch(rect, fill=False, zorder=zorder, **kwargs)

    # Add to parent axes
    parent_axes.add_patch(pp)

def mapping(
    # general arguments
    figsize,
    name_statistics,
    building_data,
    terrain_data,
    # Map arguments
    map_level_range,
    cmap_terrain,
    cmap_flood,
    # Inset arguments
    inset_box_coords,
    background_inset_box_coords,
    clip_range,
    bw_adjust,
    density_rate,
    density_y_limit,
    pdf_xaxis_label,
    x_tick_range_pdf,
    comparison_sign,
    # Zoom arguments
    zoomed_coordinates=None,
    zoom=False
):
    """
    @Definition:
                A function to create a map with density
    @Arguments:
                figsize (tuple):
                            A tuple of figsize in matplotlib subplot (width, height)
                name_statistics (string):
                            Names of statistics include:
                                - Mean (mean)
                                - Standard deviation (sd)
                                - Coefficient of variation (cv)
                                - Proportion of simulations of each cell being inundated (cell)
                building_data (shapefile read by geopandas):
                            Polygons of building data
                terrain_data (raster read by xarray):
                            Raster of terrain shading

                map_level (list):
                            A list of three values for smoothing map level - upper and lower limit and step.
                            The step would help to smooth the color
                cmap_terrain, cmap_flood (cmap with format of matplotlib):
                            Selective cmap from matplotlib

                inset_box_coords (list):
                            A list of four values [position of parent axis x, position of parent axis y,
                                                   width of inset, height of inset]
                background_inset_box_coords (list)
                            A list of four values [position of parent axis x, position of parent axis y,
                                                   width of inset background, height of inset background]
                            This is for background that covers the inset plot
                clip_range (list):
                            A list with two values of upper and lower limits for clipping
                bw_adjust (float):
                            Value to adjust the smoothiness of the density line/plot
                density_rate (float):
                            Value to scale y axis of density plot
                density_y_limit (float):
                            Value to limit the y axis of density plot
                pdf_xaxis_label (string):
                            Label of x axis of probability density function plot
                x_tick_range_pdf (list):
                            A list of three values - upper and lower limits and step.
                            Ex: Normally the mean range will be around 0.1 to 6. Hence, the x range should be [0.1,
                            6.1, 1], in that the step equals to 1
                comparison_sign (r string):
                            A comparison sign will be added into the upper limit of x axis. Ex: for larger
                            comparison, add this <r'$\geq $'>.

                zoomed_coordinates (list)
                            A list of four values - xmin, xmax, ymin, ymax
                zoom (boolean):
                            If True, the zoom plot will be added into the parent axis
                            If False, the zoom plot will not be added
    @Returns:
                None
    """
    # Get data
    # Raster
    raster = xr.open_dataset(fr"{raster_bathy}\\{name_statistics}_flowdepth.nc")

    # Copy for density
    z_val = raster.Band1.values.flatten()
    dens = z_val[~np.isnan(z_val)]

    # Set up axes
    fig, parent_axis = plt.subplots(figsize=figsize)

    # (LAYER 0) Plot terrain shade
    terrain_data.plot(ax=parent_axis, cmap=cmap_terrain, add_colorbar=False, zorder=0)

    # (LAYER 2) Plot flood map -----------------------------------------------------------------
    # Plot map
    map_level = np.arange(map_level_range[0], map_level_range[1], map_level_range[2])
    raster.Band1.plot(ax=parent_axis, levels=map_level, cmap=cmap_flood, add_colorbar=False, alpha=0.8, zorder=2)

    # Trim map
    parent_axis.set_ylim(bottom=5471400)

    # Remove 1e6
    parent_axis.ticklabel_format(style='plain')

    # x, y label titles
    parent_axis.set_xlabel("NZTM, east (m)", fontsize=18, labelpad=17)
    parent_axis.set_ylabel("NZTM, north (m)", rotation=-270, fontsize=18, labelpad=17)

    # Design size and style for ticks and labels
    for item in (parent_axis.get_xticklabels() + parent_axis.get_yticklabels()):  # For x, y ticks' labels
        item.set_fontsize(10)
    parent_axis.tick_params(direction='out', length=8, pad=7)

    # Add scale bar
    parent_axis.add_artist(ScaleBar(
        1,
        font_properties={'weight': 'bold', 'size': 20},
        pad=0.5,
        length_fraction=0.1,
        box_color=None,
        box_alpha=0,
        color='black',
        scale_formatter=lambda value, unit: f'{value} {unit}',
        location="lower left"
    ))

    # (LAYER 4) Add arrow
    x_arrow, y_arrow, arrow_length = 0.04, 0.06, 0.08
    parent_axis.annotate(
        'N', fontweight='bold', color='black',
        xy=(x_arrow, y_arrow),
        xytext=(x_arrow, y_arrow+arrow_length),
        arrowprops=dict(arrowstyle='<-, head_width=0.5', facecolor='black', edgecolor='black', linewidth=4, mutation_scale=12),
        ha='center', va='center', fontsize=20,
        xycoords=parent_axis.transAxes, zorder=4
    )

    # Plot inset plots for density --------------------------------------------------------------
    # (LAYER 4) Create inset axis
    axins = parent_axis.inset_axes(inset_box_coords, zorder=5)
    parent_axis.add_patch(plt.Rectangle((background_inset_box_coords[0], background_inset_box_coords[1]),
                                        background_inset_box_coords[2], background_inset_box_coords[3],
                                        fc="w", alpha=0.7,
                                        transform=parent_axis.transAxes, zorder=4))

    # Create density
    density_plot = sns.kdeplot(
        np.clip(dens, clip_range[0], clip_range[1])*density_rate,
        bw_method='silverman',
        bw_adjust=bw_adjust,
        clip=(clip_range[0]*density_rate, clip_range[1]*density_rate)
    )

    # Get data from density
    x_dens = density_plot.lines[0].get_data()[0]/density_rate
    y_dens = density_plot.lines[0].get_data()[1]

    # Color shaded area below the density curve
    axins.plot(x_dens, y_dens, linewidth=0.01)
    rainbow_fill_sns(axins, x_dens, y_dens, density_y_limit, cmap_flood)

    # Generate title and labels
    axins.set_title("Probability density function", pad=18, fontsize=13, fontweight='bold')
    axins.set_xlabel(pdf_xaxis_label, fontsize=10.5, labelpad=9)
    axins.set_ylabel("Probability density", rotation=-270, fontsize=10.5, labelpad=12)

    # Set tick values
    x_range = np.arange(x_tick_range_pdf[0], x_tick_range_pdf[1], x_tick_range_pdf[2])
    xlabel_arr = np.array(np.round(x_range[:], 2), dtype='float')
    axins.set_xticks(xlabel_arr)

    # Set comparison sign
    xlabel_arr_cop = np.array(np.round(x_range[:], 2), dtype='str')
    xlabel_arr_cop[-1] = '{1}{0:.1f}'.format(xlabel_arr[-1], comparison_sign)
    axins.set_xticklabels(xlabel_arr_cop)

    # Design size and style for ticks and labels
    for item in (axins.get_xticklabels() + axins.get_yticklabels()):  # For x, y ticks' labels
        item.set_fontsize(9)
    axins.tick_params(direction='out', length=5, pad=6)

    # Remove grid background lines (including x, y lines)
    axins.spines['top'].set_visible(False)
    axins.spines['right'].set_visible(False)
    axins.spines['bottom'].set_visible(False)
    axins.spines['left'].set_visible(False)

    # Zooming -------------------------------------------------------------------------------------
    if zoom:
        # Get axis
        axins_zoom = parent_axis.inset_axes([0.70, 0.57, 0.3, 0.4], zorder=5)

        # Plot hillshade
        terrain_data.plot(ax=axins_zoom, cmap=cmap_terrain, add_colorbar=False, zorder=0)

        # Plot flood map for zooming out
        raster.Band1.plot(ax=axins_zoom, levels=map_level, cmap=cmap_flood, add_colorbar=False, alpha=0.7, zorder=2)

        # Set up coordinates to be zoomed in parent axis
        axins_zoom.set_xlim(zoomed_coordinates[0], zoomed_coordinates[1])
        axins_zoom.set_ylim(zoomed_coordinates[2], zoomed_coordinates[3])

        # Remove all labels and ticks of zoom axis
        axins_zoom.set_xticks([])
        axins_zoom.set_yticks([])
        axins_zoom.set_xlabel('')
        axins_zoom.set_ylabel('')

        # Remove black edge of plot background
        for spine in axins_zoom.spines.values():
            spine.set_edgecolor('blue')
            spine.set_linewidth(3)

        # Plot zoom
        mark_inset_noconnectedline(
            parent_axis,
            axins_zoom, zorder=4,
            fc="none", ec="blue", linewidth=2
        )

        # Plot building for zooming plot
        building_data.plot(
            ax=axins_zoom, color='w', edgecolor='0.2', alpha=1, zorder=1
        )

        # Plot basemap for zooming out
        ctx.add_basemap(
            ax=axins_zoom, crs=2193,
            source='https://server.arcgisonline.com/ArcGIS/rest/services/Canvas/World_Dark_Gray_Base/MapServer/tile/{z}/{y}/{x}',
            alpha=.3,
            zoom=16,
            zorder=0
        )

    else:
        pass


    # Plot basemap, building and save plot ---------------------------------------------------------
    # Save basemap
    ctx.add_basemap(
        ax=parent_axis, crs=2193,
        source='https://server.arcgisonline.com/ArcGIS/rest/services/Canvas/World_Dark_Gray_Base/MapServer/tile/{z}/{y}/{x}',
        alpha=0.5, zorder=1
    )

    # Plot building
    building_data.plot(
        ax=parent_axis, color='w', edgecolor='0.2', linewidth=0.05, alpha=0.8, zorder=3
    )

    # Save fig
    fig.savefig(
        fr"{plot_bathy}\\{name_statistics}.png",
        bbox_inches='tight', dpi=330
    )


def map_plotting(
    statistics,
    building_data,
    terrain_data
):
    """
    @Definition:
                A function to execute mapping
    @Arguments:
                statistics (list):
                            A list of string includes the name of each statistics ['mean', 'sd', 'cv', 'cell']
                building_data (shapefile read by geopandas):
                            Polygons of building data
                terrain_data (raster read by xarray):
                            Raster of terrain shading
    @Returns:
                None
    """
    # GENERAL --------------------------------------------------------------------
    # General arguments
    figsize = (25, 14)

    # Inset arguments
    inset_box_coords = [0.08, 0.74, 0.19, 0.169]
    background_inset_box_coords = [0.02, 0.68, 0.28, 0.29]


    # EACH STATISTIC -------------------------------------------------------------
    for each_stat in statistics:
        # Mean
        if each_stat == 'mean':
            # General arguments
            name_statistics = 'mean'

            # Map arguments
            map_level_range = [0.1, 5.1, 0.001]
            cmap_terrain = plt.get_cmap('gist_gray')
            # cmap_flood = plt.get_cmap('ocean')
            cmap_flood = get_gradient_cmap(hex_list110[::-1])

            # Inset arguments
            clip_range = [0.1, 5.1]
            bw_adjust = 1
            density_rate = 1
            density_y_limit = 0.75
            pdf_xaxis_label = "Means (m)"
            x_tick_range_pdf = [0.1, 6.1, 1]
            comparison_sign = r'$\geq $'

            # Zoom
            zoomed_coordinates = None
            zoom = False

        # STANDARD DEVIATION
        elif each_stat == 'sd':
            # General arguments
            name_statistics = 'sd'

            # Map arguments
            map_level_range = [0, 1.1, 0.001]
            cmap_terrain = plt.get_cmap('gist_gray')
            # cmap_flood = plt.get_cmap('turbo')
            cmap_flood = get_gradient_cmap(hex_list125)

            # Inset arguments
            clip_range = [0, 1]
            bw_adjust = 2
            density_rate = 10
            density_y_limit = 0.35
            pdf_xaxis_label = "Standard deviations (m)"
            x_tick_range_pdf = [0, 1.1, 0.2]
            comparison_sign = r'$\geq $'

            # Zoom
            zoomed_coordinates = None
            zoom = False

        # COEFFICIENT OF VARIATION
        elif each_stat == 'cv':
            # General arguments
            name_statistics = 'cv'

            # Map arguments
            map_level_range = [0, 200.1, 0.1]
            cmap_terrain = plt.get_cmap('gist_gray')
            # cmap_flood = plt.get_cmap('gnuplot')
            cmap_flood = get_gradient_cmap(hex_list108)

            # Inset arguments
            clip_range = [0, 200]
            bw_adjust = 1
            density_rate = 1
            density_y_limit = 0.035
            pdf_xaxis_label = "Coefficients of variations (m)"
            x_tick_range_pdf = [0, 201, 50]
            comparison_sign = r'$\geq $'

            # Zoom
            zoomed_coordinates = None
            zoom = False

        # PROPORTION OF SIMULATIONS OF EACH CELL BEING INUNDATED
        else:
            # General arguments
            name_statistics = 'cell'

            # Map arguments
            map_level_range = [0, 100, 5]
            cmap_terrain = plt.get_cmap('gist_gray')
            # cmap_flood = plt.get_cmap('turbo')
            cmap_flood = get_gradient_cmap(hex_list124)

            # Inset arguments
            clip_range = [0, 100]
            bw_adjust = 3
            density_rate = 1
            density_y_limit = 0.045
            pdf_xaxis_label = "Proportion (%)"
            x_tick_range_pdf = [0, 102, 20]
            comparison_sign = r''

            # Zoom
            zoomed_coordinates = [1771500, 1772100, 5472950, 5473600]
            zoom = True

        mapping(
            # General arguments
            figsize,
            name_statistics,
            building_data,
            terrain_data,
            # Map arguments
            map_level_range,
            cmap_terrain,
            cmap_flood,
            # Inset arguments
            inset_box_coords,
            background_inset_box_coords,
            clip_range,
            bw_adjust,
            density_rate,
            density_y_limit,
            pdf_xaxis_label,
            x_tick_range_pdf,
            comparison_sign,
            # Zoom arguments
            zoomed_coordinates,
            zoom
        )

# END VARIATION MAP PLOTTING ###########################################################################################


# IMPACT PLOTTING ######################################################################################################
def area_building_plotting(
    figsize,
    df,
    text_box_location,
    option
):
    """
    @Definition:
                A function to plot distribution of areas
    @References:
                https://stackoverflow.com/questions/69524514/how-to-modify-the-kernel-density-estimate-line-in-a-sns-histplot
                https://seaborn.pydata.org/generated/seaborn.distplot.html

                https://stackoverflow.com/questions/65400669/how-to-generate-two-separate-y-axes-for-a-histogram-on-the-same-figure-in-seabor
                https://stackoverflow.com/questions/26752464/how-do-i-align-gridlines-for-two-y-axis-scales-using-matplotlib
                https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.hist.html#matplotlib.axes.Axes.hist
                https://seaborn.pydata.org/generated/seaborn.kdeplot.html

                https://seaborn.pydata.org/generated/seaborn.histplot.html#seaborn.histplot
                https://seaborn.pydata.org/generated/seaborn.distplot.html
                https://stackoverflow.com/questions/27671748/how-to-print-y-axis-label-horizontally-in-a-matplotlib-pylab-chart
                https://stackoverflow.com/questions/24391892/printing-subscript-in-python

                https://stackoverflow.com/questions/45037386/trouble-aligning-ticks-for-matplotlib-twinx-axes (best
                answer for align two axis)
                https://stackoverflow.com/questions/12608788/changing-the-tick-frequency-on-x-or-y-axis-in-matplotlib
    @Arguments:
                figsize (tuple):
                            A tuple of figsize in matplotlib subplot (width, height)
                df (pandas dataframe):
                            Dataframe of simulations' areas
                text_box_location (list):
                            A list of text box coordinates [x, y]
                option (string):
                            "area" for area plotting
                            "building" for building plotting
    @Returns:
                None
    """
    # Set up axes
    fig, parent_axis = plt.subplots(figsize=figsize)

    # Get area's values under array format
    values = df.iloc[0].to_numpy()

    # Number of simulations
    num_bin = len(values)

    # Get another y axis - for 'Density'
    axis_density = parent_axis.twinx()

    # Plotting
    # AREA ---------------------------------------------------------
    if option == 'area':
        # Frequency plot
        sns.histplot(values / 100, bins=num_bin, stat='count',
                     legend=False,
                     edgecolor='darkgreen',
                     facecolor='springgreen',
                     ax=parent_axis)

        # Density plot
        sns.kdeplot(values / 100, ax=axis_density, color="r", linewidth=2)

        # X label
        x_label = f'Areas (x100 '
        x_label += r'$\mathrm{m}^2$)'

    # BUILDING ----------------------------------------------------
    else:
        # Frequency plot
        sns.histplot(values, bins=num_bin, stat='count',
                     legend=False,
                     edgecolor='navy',
                     facecolor='deepskyblue',
                     ax=parent_axis)

        # Density plot
        sns.kdeplot(values, ax=axis_density, color="r", linewidth=2)

        # Set y label for 'Frequency'
        parent_axis.set_ylabel("Number of simulations", fontsize=20, labelpad=38)

        # X label
        x_label = f'Number of buildings'

    # Set y label for 'Frequency'
    parent_axis.set_ylabel("Number of simulations", fontsize=20, labelpad=38)

    # Set y label for 'Density'
    axis_density.set_ylabel("Probability density", rotation=270, fontsize=20, labelpad=38)

    # Design ticks and x figures for 'Density'
    axis_density.tick_params(direction='out', length=8, pad=10)
    for item in (axis_density.get_xticklabels() + axis_density.get_yticklabels()):  # For x, y ticks' labels
        item.set_fontsize(15)

    # Approximately align with the first y axis. Please visit here for more information
    # https://stackoverflow.com/questions/45037386/trouble-aligning-ticks-for-matplotlib-twinx-axes
    # https://stackoverflow.com/questions/12608788/changing-the-tick-frequency-on-x-or-y-axis-in-matplotlib (for
    # changing interval in ylim)
    start, end = parent_axis.get_ylim()
    parent_axis.yaxis.set_ticks(np.arange(start, end, 1))

    len_axis1 = parent_axis.get_ylim()

    # Get lims of second axis - Density
    len_axis2 = axis_density.get_ylim()

    # Develop a function to get the general ticks for both x and y
    f = lambda x: len_axis2[0] + (x - len_axis1[0]) / (len_axis1[1] - len_axis1[0]) * (len_axis2[1] - len_axis2[0])

    # Get number of ticks
    num_ticks = f(parent_axis.get_yticks())

    # Set the ticks for second y axis
    axis_density.yaxis.set_major_locator(matplotlib.ticker.FixedLocator(num_ticks))

    # Set up for x axis to make sure it appears with integers
    parent_axis.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))

    # Ticks, title, and y label
    parent_axis.tick_params(direction='out', length=8, pad=10)
    parent_axis.set_xlabel(x_label, fontsize=20, labelpad=38)

    for item in (parent_axis.get_xticklabels() + parent_axis.get_yticklabels()):  # For x, y ticks' labels
        item.set_fontsize(15)

    # Create text string
    textstr = "mean = {0:.3f}".format(df.mean(axis=1)[0])
    textstr += "\nstdev = {0:.3f}".format(df.std(axis=1)[0])

    # place a text box in upper left in axes coords
    # Refer here for more information
    # https://stackoverflow.com/questions/50869424/increase-line-separation-in-matplotlib-annotation-text/50888491
    parent_axis.text(text_box_location[0], text_box_location[1], textstr, transform=parent_axis.transAxes, fontsize=13,
                     linespacing=1.5,
                     fontweight='normal',
                     fontstyle='italic',
                     horizontalalignment='left',
                     verticalalignment='top')

# END IMPACT PLOTTING ##################################################################################################


