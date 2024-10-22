# Prepare packages -----------------------------------------------------------------------------------------------------
from folder import *                                        # For paths of sub-folders

# For copy files
import shutil

# Packages for directory manipulation
import os                                                   # For path manipulation

# Packages for data manipulation
import numpy as np                                          # For data array manipulation
import pandas as pd                                         # For reading data from text or csv files

# Packages for spatial dataframe
import geopandas as gpd                                     # For spatial dataframe manipulation

# Packages for spatial data manipulation
import rioxarray as rxr                                     # For reading spatial data

# Packages for geometry manipulation
from shapely.geometry import LineString, MultiPoint         # For LineString, MultiPoint creation

# Packages for time manipulation
import datetime                                             # For converting numbers into datetime format

# ----------------------------------------------------------------------------------------------------------------------


# TIDE & FLOW DATA #####################################################################################################
def get_tide_data(tide_path):
    """
    @Definition:
                A function to get data from tide path (csv file)
    @References:
                https://stackoverflow.com/questions/13703720/converting-between-datetime-timestamp-and-datetime64
                https://stackoverflow.com/questions/16158795/cant-convert-dates-to-datetime64
                https://www.w3resource.com/pandas/series/series-dt-tz_convert.php#:~:text=The%20tz_convert()%20function%20is,one%20time%20zone%20to%20another.&text=Time%20zone%20for%20time.,and%20remove%20the%20timezone%20information
                https://stackoverflow.com/questions/26089670/unable-to-apply-methods-on-timestamps-using-series-built-ins
                https://docs.python.org/3/library/datetime.html#strftime-and-strptime-behavior
    @Arguments:
                coordinates_func (list):
                                        A list of xmin, ymin, xmax, ymax
                addition (int):
                                        A number used to extend the padding
    @Returns:
                (list):
                                        A list of x min, x max, y min and y max
    """
    # Get tide from file
    tide = pd.read_csv(
        tide_path,
        skiprows=np.arange(0, 9, 1),
        names=['DateChar', 'Level'],
        header=None
    )

    # Convert time data
    tide['DateTime'] = pd.to_datetime(tide['DateChar'], format="%Y-%m-%dT%H:%M:%SZ", utc=False)
    tide['DateTime'] = tide['DateTime'].dt.tz_localize(None) # Remove timezone

    return tide

def get_flow_data(flow_path):
    """
    @Definition:
                A function to get data from flow path (csv file)
    @References:
                https://stackoverflow.com/questions/13703720/converting-between-datetime-timestamp-and-datetime64
                https://stackoverflow.com/questions/16158795/cant-convert-dates-to-datetime64
                https://www.w3resource.com/pandas/series/series-dt-tz_convert.php#:~:text=The%20tz_convert()%20function%20is,one%20time%20zone%20to%20another.&text=Time%20zone%20for%20time.,and%20remove%20the%20timezone%20information
                https://stackoverflow.com/questions/26089670/unable-to-apply-methods-on-timestamps-using-series-built-ins
                https://docs.python.org/3/library/datetime.html#strftime-and-strptime-behavior
    @Arguments:
                tide_path (string):
                                        A string of directory path and tide file
    @Returns:
                tide (pandas dataframe):
                                        A dataframe of DateChar, Level, DateTime
    """


    # Get flow from file
    flow = pd.read_csv(
        flow_path,
        names=['DateChar', 'Flow', 'Quality'],
        skiprows=[0],
        header=None
    )

    # Convert time data
    flow['DateTime'] = pd.to_datetime(flow['DateChar'], format="%m/%d/%Y %H:%M:%S", utc=False)
    flow['DateTime'] = flow['DateTime'].dt.tz_localize(None)

    return flow

def tide_flow_data(tide_path, flow_path, date_start, date_end):
    """
    @Definition:
                A function to get data from tide path (csv file)
    @References:
                https://www.geeksforgeeks.org/how-to-add-time-onto-a-datetime-object-in-python/
                https://stackoverflow.com/questions/69283228/how-to-remove-hours-minutes-seconds-and-utc-offset-from-pandas-date-column-i
                https://stackoverflow.com/questions/40992976/python-convert-datetime-column-into-seconds
                https://gis.stackexchange.com/questions/367228/using-shapely-interpolate-to-evenly-re-sample-points-on-a-linestring-geodatafram

                https://stackoverflow.com/questions/42826388/using-time-zone-in-pandas-to-datetime
                https://stackoverflow.com/questions/22800079/converting-time-zone-pandas-dataframe
                https://stackoverflow.com/questions/25015711/time-data-does-not-match-format

                https://stackoverflow.com/questions/13703720/converting-between-datetime-timestamp-and-datetime64
                https://stackoverflow.com/questions/16158795/cant-convert-dates-to-datetime64
                https://www.w3resource.com/pandas/series/series-dt-tz_convert.php#:~:text=The%20tz_convert()%20function%20is,one%20time%20zone%20to%20another.&text=Time%20zone%20for%20time.,and%20remove%20the%20timezone%20information
                https://stackoverflow.com/questions/26089670/unable-to-apply-methods-on-timestamps-using-series-built-ins
                https://docs.python.org/3/library/datetime.html#strftime-and-strptime-behavior

    @Arguments:
                tide_path (string):
                                        A string of directory path and tide file
                flow_path (string):
                                        A string of directory path and flow file
                date_start (string):
                                        A string of starting date. Ex: <2005-01-05 00:00:00>
                date_end (string):
                                        A string of ending date. Ex: <2005-01-07 00:00:00>
    @Returns:
                tide_flow_data (pandas dataframe):
                                        A dataframe of DateTime, Flow, Level
    """

    # Get file values
    tide_data = get_tide_data(tide_path)
    flow_data = get_flow_data(flow_path)

    # Convert timezone of flow data
    time_change = datetime.timedelta(hours=12)
    flow_data['DateTime'] = flow_data['DateTime'] - time_change

    # Convert level values of tide data
    tide_data['Level'] = tide_data['Level'] - 0.11

    # Align flow and tide data
    tide_flow_data = pd.merge(tide_data, flow_data, how='outer', on='DateTime')[['DateTime', 'Flow', 'Level']]
    tide_flow_data = tide_flow_data.sort_values(by='DateTime').reset_index(drop=True)
    tide_flow_data['Level'] = tide_flow_data['Level'].interpolate(method='cubic')

    # Merge tide and flow
    tide_flow_data = tide_flow_data.loc[(tide_flow_data.DateTime >= date_start)
                                        & (tide_flow_data.DateTime <= date_end)]

    # Write out csv
    tide_flow_data.to_csv(fr"{other_data}\\tide_flow_data.csv", index=False)

# END TIDE & FLOW DATA #################################################################################################

# BDY FILE #############################################################################################################

def bdy_generation(
    resolution_func, domain_name
):
    """
    @Definition:
                A function to generate BDY file
    @References:
                https://stackoverflow.com/questions/51865367/cannot-convert-the-series-to-class-int
    @Arguments:
                resolution_func:
                                        Resolution of raster
                boundary_no_padding (list):
                                        A list contains xmin, ymin, xmax, ymax
    @Returns:
                None.
    """
    # Get tide flow data
    tide_flow_df = pd.read_csv(fr"{other_data}\\tide_flow_data.csv", sep=',')
    tide_flow_df['DateTime'] = pd.to_datetime(tide_flow_df['DateTime'], format="%Y-%m-%d %H:%M:%S", utc=False)

    # Convert to seconds
    tide_flow_df['seconds'] = pd.to_numeric(tide_flow_df.DateTime - np.min(tide_flow_df.DateTime))/(10**9)
    tide_flow_df['seconds'] = tide_flow_df['seconds'].astype('int')

    # Extract only flow data
    flow_df = tide_flow_df[['Flow', 'seconds']].copy(deep=True)
    flow_df['Flow'] = round(flow_df['Flow']/resolution_func, 4)
    flow_df['seconds'] = flow_df['seconds'].astype(int)

    # Extract only tide data
    tide_df = tide_flow_df[['Level', 'seconds']].copy(deep=True)
    tide_df['Level'] = round(tide_df['Level'], 3)
    tide_df['seconds'] = tide_df['seconds'].astype(int)

    # Write BDY file includes discharge and tide level
    with open(fr"{FPpara_path}\\bathymetry.bdy", "w") as discharge_tide:
        # Write discharge
        discharge_note = "LISFLOOD-FP setup\n"
        discharge_reference = discharge_note + domain_name + "\n"
        discharge_column_name = discharge_reference + '{0:<20}seconds\n'.format(tide_flow_df.Flow.shape[0])
        discharge_tide.write(discharge_column_name)
        for discharge_line in range(tide_flow_df.Flow.shape[0]):
            data_discharge = flow_df.iloc[discharge_line]
            text_river_discharge = '{0[0]:<20}{0[1]:.0f}\n'.format(data_discharge)
            discharge_tide.write(text_river_discharge)

        # Write tide
        tide_reference = 'Tide\n'
        tide_column_name = tide_reference + '{0:<20}seconds\n'.format(tide_flow_df.Level.shape[0])
        discharge_tide.write(tide_column_name)
        for tide_line in range(tide_flow_df.Level.shape[0]):
            data_tide = tide_df.iloc[tide_line]
            text_tide = '{0[0]:<20}{0[1]:.0f}\n'.format(data_tide)
            discharge_tide.write(text_tide)
# END BDY FILE #########################################################################################################


# SAMPLED TIDAL POINTS FOR BOUNDARY ####################################################################################
def tideboundary_points(
    boundary_no_padding,
    automatic_linestring=True
):
    """
    @Definition:
                A function to sample points on LineString. These points are tidal boundary
    @References:
                https://gis.stackexchange.com/questions/309251/how-to-get-equidistant-points-from-a-linestring-geographical-coordinates
    @Arguments:
                boundary_no_padding (list):
                                        A list contains xmin, ymin, xmax, ymax
    @Returns:
                None.
    """
    if automatic_linestring:
        # Get extent information of original no-padding raster
        no_padding_xmin, no_padding_ymin, no_padding_xmax, no_padding_ymax = boundary_no_padding

        # Create tidal LineString
        tideboundary_line = LineString([
            (no_padding_xmin + 3100 + 50, no_padding_ymax),
            (no_padding_xmin, no_padding_ymax - 2170 - 50)
        ])

        # Sample points on LineString
        tide_sampledis = 100
        tide_num_points = 40
        # Create a list comprehension to collect all coordinates of points
        tidepoints = [
            tideboundary_line.interpolate((1/(tideboundary_line.length/tide_sampledis)) * i, normalized=True) for i in range(tide_num_points)
        ]
        # Gather all points into MultiPoint geometry
        tidepoints_geo = MultiPoint(tidepoints)

    else:
        # Get tidal LineString
        path_linestring = fr"{other_data}\\linestring.shp"
        linestring_file = gpd.read_file(path_linestring)
        tideboundary_line = linestring_file['geometry'][0]

        # Sample points on LineString
        tide_sampledis = 6.5
        tide_num_points = 8
        # Create a list comprehension to collect all coordinates of points
        tidepoints = [
            tideboundary_line.interpolate((1 / (tideboundary_line.length / tide_sampledis)) * i, normalized=True) for i
            in range(tide_num_points)
        ]
        # Gather all points into MultiPoint geometry
        tidepoints_geo = MultiPoint(tidepoints)


    # List all points in MultiPoint
    tideboundary_points = list(tidepoints_geo.geoms)
    tideboundary_data = {'geometry': tideboundary_points}
    tidepts = gpd.GeoDataFrame(data=tideboundary_data, crs='EPSG:2193')

    # Write out tidal boundary points
    tidepts.to_file(
        fr"{other_data}\\tideboundary_points.geojson", driver="GeoJSON"
    )

# END SAMPLED TIDAL POINTS FOR BOUNDARY ################################################################################


# BCI FILE #############################################################################################################

def bci_generation(
    river_points, domain_name,
    ran_trans_i
):
    """
    @Definition:
                A function to generate BCI file
    @References:
                None.
    @Arguments:
                river_points (tuple):
                                        A tuple contains x and y coordinates
                domain_name (string):
                                        Name of domain. Ex: "Waikanae"
                ran_trans_i (array):
                                        A 3D array contains values of angle, x, and y coordinates of points in tiles
                                        (iterating variable generated from multiprocessing)
    @Returns:
                None.
    """
    # Get values
    number_simulation = ran_trans_i

    # RIVER POINT ----------------------------------------------------
    # Create array of point coordinates
    point_source = np.array([river_points]).astype('float64')

    # Construct river point injection
    injection = ['P', point_source[0][0], point_source[0][1], 'QVAR', domain_name]

    # TIDE BOUNDARY ---------------------------------------------------
    # Transformed folder creation
    pathlib.Path(fr"{FPpara_path}\\bathy_{number_simulation}").mkdir(parents=True, exist_ok=True)
    bci_path = fr"{FPpara_path}\\bathy_{number_simulation}"
    # Transform tidal boundary points
    tide_path_in = fr"{other_data}\\tideboundary_points.geojson"
    tide_path_out = fr"{bci_path}\\tideboundary_points_{number_simulation}.geojson"
    shutil.copy2(tide_path_in, tide_path_out)
    # Get dataframe of tidal boundary points
    tidepts = gpd.read_file(
        tide_path_out
    )

    # Write BCI file - injection point and tide
    with open(fr"{bci_path}\\bathy_{number_simulation}.bci", "w") as boundary:
        injection_text = '{0[0]:<5}{0[1]:<20}{0[2]:<20}{0[3]:<7}{0[4]:<5}\n'.format(injection)
        boundary.write(injection_text)
        for line in range(tidepts.shape[0]):
            tide_boundary = ['P', tidepts.iloc[line][0].x, tidepts.iloc[line][0].y,
                             'HVAR', 'Tide']
            tide_text = '{0[0]:<5}{0[1]:<20}{0[2]:<20}{0[3]:<7}{0[4]:<5}\n'.format(tide_boundary)
            boundary.write(tide_text)

# END BCI FILE #########################################################################################################


# PAR FILE #############################################################################################################

def par_generation(
    number_simulation
):
    """
    @Definition:
                A function to generate PAR file
    @References:
                None.
    @Arguments:
                number_simulation (string):
                                            A string to identify the order of simulation (should be angle, x, y)
    @Returns:
                None.
    """
    # Get directories
    # OUTPUT
    flood_output = fr"{FPoutput_path}\\bathy_{number_simulation}"
    pathlib.Path(flood_output).mkdir(parents=True, exist_ok=True)
    # BDY
    flood_bdy = fr"{FPpara_path}\\bathymetry.bdy"
    # BCI
    flood_bci = fr"{FPpara_path}\\bathy_{number_simulation}\\bathy_{number_simulation}.bci"
    # DEM
    flood_dem = fr"{dem_asc_path}\\generated_dem_bathy_{number_simulation}.asc"
    # MANNING'S N
    flood_n = fr"{n_asc_path}\\generated_n_bathy_{number_simulation}.asc"
    # STARTDEPTH
    flood_startdepth = fr"{startdepth_asc_path}\\startdepth_{number_simulation}.asc"

    # Create parameters list
    parameters_list = [
        ('resroot', 'out'),
        ('dirroot', flood_output),
        ('saveint', 3600),
        ('massint', 100),
        ('sim_time', 172800),
        ('initial_tstep', 2),
        ('bcifile', flood_bci),
        ('bdyfile', flood_bdy),
        ('DEMFile', flood_dem),
        ('manningfile', flood_n),
        ('startfile', flood_startdepth)
    ]

    # Write into array
    parameters_array = np.array(parameters_list)

    # Write PAR file
    with open(
        fr"{FPpara_path}\\bathy_{number_simulation}\\bathy_{number_simulation}.par", "w"
    ) as parameters:
        for each_parameter in range(parameters_array.shape[0]):
            data_parameter = parameters_array[each_parameter]
            text_parameter = '{0[0]:<20}{0[1]}\n'.format(data_parameter)
            parameters.write(text_parameter)
        parameters.write('acceleration\ndrain_nodata\n\n')
# END PAR FILE #########################################################################################################


# RUN FLOOD MODEL ######################################################################################################
def run_LISFLOOD(
    lisflood_software_path,
    number_simulation
):
    """
    @Definition:
                A function to run LISFLOOD
    @References:
                None.
    @Arguments:
                lisflood_software_path (string):
                                            Directory where LISFLOOD-FP software is stored
                number_simulation (string):
                                            A string to identify the order of simulation (should be angle, x, y)
    @Returns:
                None.
    """
    # Get directory of LISFLOOD-FP software function
    os.chdir(lisflood_software_path)

    # Get directory of par file
    lisflood_par_path = fr"{FPpara_path}\bathy_{number_simulation}\bathy_{number_simulation}.par"
    os.system(f"lisflood_v8_1_0.exe -v {lisflood_par_path}")

# END RUN FLOOD MODEL ##################################################################################################


# SIMULATE FLOOD EVENT #################################################################################################

def flood_simulation(
    domain_name,
    river_points,
    lisflood_software_path,
    ran_trans_i
):
    """
    @Definition:
                A function to generate BCI file
    @References:
                None.
    @Arguments:
                domain_name (string):
                                        Name of domain. Ex: "Waikanae"
                river_points (tuple):
                                        A tuple contains x and y coordinates
                lisflood_software_path (string):
                                        Directory where LISFLOOD-FP software is stored
                ran_trans_i (array):
                                        A 3D array contains values of angle, x, and y coordinates of points in tiles
                                        (iterating variable generated from multiprocessing)
    @Returns:
                None.
    """
    # Get values
    number_simulation = ran_trans_i

    # Create BCI file
    bci_generation(
        river_points, domain_name,
        ran_trans_i
    )

    # Create PAR file
    par_generation(
        number_simulation
    )

    # Run LISFLOOD-FP
    run_LISFLOOD(
        lisflood_software_path,
        number_simulation
    )
# END SIMULATE FLOOD PARAMETER FILES ###################################################################################



