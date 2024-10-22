# Prepare packages -----------------------------------------------------------------------------------------------------
from folder import *                                        # For paths of sub-folders

# Package for paths manipulation
import pathlib                                              # For manipulating the directory path

# Packages for manipulating tabular data
import pandas as pd                                         # For manipulating tabular data (e.g. dataframe)
import numpy as np                                          # For manipulating arrays/matrices

# Packages for manipulating spatial data (e.g. raster)
import rioxarray as rxr                                     # For manipulating spatial data under xarray array format
import json                                                 # For creating json text data format
from geofabrics import processor                            # For executing the conversion of LiDAR into a DEM raster

# Packages for manipulating geometries
import shapely.geometry                                     # For creating a polygon object
import geopandas as gpd                                     # For manipulating shape files

from valueChange import *
# ----------------------------------------------------------------------------------------------------------------------


# LIDAR-DERIVED DEMS ###################################################################################################
def dem_raster(resolution_func, chunk_size_func, processor_func,
               number_simulation, padding_func, lidar_dataset_name):
    """
    @Definition:
                A function to create a raster file from a las/laz file
    @References:
                https://github.com/rosepearson/GeoFabrics
    @Arguments:
                resolution_func (int or float):
                                            Resolution value in meter
                chunk_size_func (int):
                                            Size value of chunk
                processor_func (int):
                                            Number of processor
                number_simulation (string):
                                            A string to identify the order of simulation (should be angle, x, y)
                padding_func (list):
                                            A list of x min, x max, y min, and y max
                lidar_dataset_name (string):
                                            LiDAR name
    @Returns:
                None.
    """

    # Name of output
    output_name = fr"generated_dem_bathy_{number_simulation}.nc"

    # Get paths
    basepath = pathlib.Path(fr"{bathy_path}")
    result_folder = fr"bathy_{number_simulation}"
    bathy_dir = basepath / result_folder
    bathy_dir.mkdir(exist_ok=True)


    # Set crs
    h_crs = 2193
    v_crs = 7839

    # Bounding box
    x0 = padding_func[0]
    y0 = padding_func[1]
    x1 = padding_func[2]
    y1 = padding_func[3]
    catchment = shapely.geometry.Polygon([(x0, y0), (x1, y0), (x1, y1), (x0, y1)])
    catchment = gpd.GeoSeries([catchment])
    catchment = catchment.set_crs('epsg:2193')
    catchment.to_file(basepath / result_folder / "catchment_boundary.geojson", driver="GeoJSON")

    # Design JSON instructions --------------------------------------------
    # Design JSON instructions
    instruction_json = {}

    # DEM
    instruction_json['dem'] = {
        # outputs
        "output": {
            "crs": {
                "horizontal": h_crs,
                "vertical": v_crs
            },
            "grid_params": {
                "resolution": resolution_func
            }
        },

        # processing
        "processing": {
            "chunk_size": chunk_size_func,
            "number_of_cores": processor_func
        },

        # apis
        "apis": {
            "open_topography": {
                f"{lidar_dataset_name}": {
                    "crs": {
                        "horizontal": h_crs,
                        "vertical": v_crs
                    }
                }
            }
        },

        # data paths
        "data_paths": {
            "local_cache": str(basepath),
            "subfolder": result_folder,
            "catchment_boundary": "catchment_boundary.geojson",
            "raw_dem": "raw_dem.nc",

            # For land
            "land": "land.geojson",

            # For bathymetry
            "bathymetry_contours": "bathymetry_contours.geojson",
            "river_bathymetry": ["river_bathymetry.geojson", "fan_bathymetry.geojson",
                                 "open_drain_elevation_5m_width.geojson", "closed_drain_elevation_5m_width.geojson"],
            "river_polygons": ["river_polygon.geojson", "fan_polygon.geojson",
                               "open_drain_polygon_5m_width.geojson", "closed_drain_polygon_5m_width.geojson"],
            # Result
            "result_dem": fr"{dem_nc_path}\\{output_name}"
        },

        # general
        "general": {
            # For lidar
            "set_dem_shoreline": True,
            "drop_offshore_lidar": False,
            "lidar_classification_to_keep": [2, 9],
            "interpolation_method": "nearest",
            "lidar_interpolation_method": "idw",

            # For bathymetry
            "bathymetry_contours_z_label": "valdco",
            "bathymetry_points_type": ['rivers', 'rivers', 'drains', 'drains'],
            "bathymetry_points_z_label": ['bed_elevation_Rupp_and_Smart', "depths", "elevation", "elevation"]
        }
    }

    # Save the instructions
    with open(basepath / result_folder / "instructions.json", "w") as instruction:
        json.dump(instruction_json, instruction)

    runner = processor.RawLidarDemGenerator(instruction_json['dem'])
    runner.run()
    runner = processor.HydrologicDemGenerator(instruction_json['dem'])
    runner.run()
# END LIDAR-DERIVED DEMS ###############################################################################################


# ROUGHNESS ############################################################################################################
def roughness_raster(resolution_func, chunk_size_func, processor_func,
                     number_simulation, padding_func, lidar_dataset_name):
    """
    @Definition:
                A function to create a raster file from a las/laz file
    @References:
                https://github.com/rosepearson/GeoFabrics
    @Arguments:
                resolution_func (int or float):
                                            Resolution value in meter
                chunk_size_func (int):
                                            Size value of chunk
                processor_func (int):
                                            Number of processor
                number_simulation (string):
                                            A string to identify the order of simulation (should be angle, x, y)
                padding_func (list):
                                            A list of x min, x max, y min, and y max
                lidar_dataset_name (string):
                                            LiDAR name
    @Returns:
                None.
    """

    # Name of output
    output_roughness = fr"generated_roughness_bathy_{number_simulation}.nc"
    output_dem = fr"generated_dem_bathy_{number_simulation}.nc"

    basepath = pathlib.Path(fr"{bathy_path}")
    result_folder = f'bathy_roughness_{number_simulation}'

    # Set crs
    h_crs = 2193
    v_crs = 7839

    # Bounding box
    x0 = padding_func[0]
    y0 = padding_func[1]
    x1 = padding_func[2]
    y1 = padding_func[3]
    catchment = shapely.geometry.Polygon([(x0, y0), (x1, y0), (x1, y1), (x0, y1)])
    catchment = gpd.GeoSeries([catchment])
    catchment = catchment.set_crs(h_crs)
    catchment.to_file(basepath / result_folder / "catchment_boundary.geojson", crs="EPSG:2193", driver="GeoJSON")

    # Roughness
    instruction_roughness = {
        # roughness - output
        "output": {
            "crs": {
                "horizontal": h_crs,
                "vertical": v_crs
            },
            "grid_params": {
                "resolution": resolution_func
            }
        },

        # roughness - apis
        "apis": {
            "open_topography": {
                f"{lidar_dataset_name}": {
                    "crs": {
                        "horizontal": 2193,
                        "vertical": 7839
                    }
                }
            }
        },

        # roughness - processing
        "processing": {
            "chunk_size": chunk_size_func,
            "number_of_cores": processor_func
        },

        # roughness - data_paths
        "data_paths": {
            "local_cache": str(basepath),
            "subfolder": result_folder,
            "catchment_boundary": "catchment_boundary.geojson",
            "land": "land.geojson",
            "bathymetry_contours": "bathymetry_contours.geojson",
            "result_dem": fr"{dem_nc_path}\\{output_dem}",
            "result_geofabric": fr"{roughness_nc_path}\\{output_roughness}"
        },

        # roughness - general
        "general": {
            "set_dem_shoreline": True,
            "drop_offshore_lidar": False,
            "lidar_classifications_to_keep": [1, 2, 4, 9],
            "interpolation_method": "nearest",
            "lidar_interpolation_method": "idw"
        },

        # roughness - rivers
        "drains": {
            "width": 5
        }
    }

    # Save the instructions
    with open(basepath / result_folder / "instructions.json", "w") as instruction:
        json.dump(instruction_roughness, instruction)

    # Create ROUGHNESS raster
    runner = processor.RoughnessLengthGenerator(instruction_roughness)
    runner.run()

def zo_to_n(number_simulation, H):
    """
    @Definition:
                A function to create a raster file of Manning's n from roughness length
    @References:
                https://doi.org/10.1080/15715124.2017.1411923
    @Arguments:
                number_simulation (string):
                                            A string to identify the order of simulation (should be angle, x, y)
                H (int):
                                            Value of depth
    @Returns:
                A raster of Manning's n
    """
    # Extract roughness length
    full_roughness_nc = rxr.open_rasterio(
        fr"{roughness_nc_path}\\generated_roughness_bathy_{number_simulation}.nc"
    )
    zo = full_roughness_nc.zo

    # Convert roughness length (zo) to Manning's n
    manning_n = (0.41 * (H ** (1 / 6)) * ((H / zo) - 1)) / (np.sqrt(9.80665) * (1 + (H / zo) * (np.log(H / zo) - 1)))

    # Calibrate manning's n
    new_manning_n = manning_n * 1.7

    # Write out Manning's n
    new_manning_n.rio.to_raster(
        fr"{n_nc_path}\\generated_n_bathy_{number_simulation}.nc"
    )

# END ROUGHNESS ######################################################################################################


# STARTDEPTH ###########################################################################################################
def startdepth_generation(
    dem_onlychangepadding_path,
    number_simulation
):
    """
    @Definition:
                A function to generate startdepth file
    @References:
                None.
    @Arguments:
                dem_onlychangepadding_path (string):
                                Directory of dem that has no padding
                number_simulation (string):
                                A string to identify the order of simulation (should be angle, x, y)
    @Return:
                None.
    """
    # Get tide flow data
    tide_flow_df = pd.read_csv(fr"{other_data}\\tide_flow_data.csv", sep=',')
    tide_flow_df['DateTime'] = pd.to_datetime(tide_flow_df['DateTime'], format="%Y-%m-%d %H:%M:%S", utc=False)

    # Read dem
    dem_ori = rxr.open_rasterio(fr"{dem_onlychangepadding_path}")
    dem = dem_ori.z.copy(deep=True)

    # Select values for startdepth
    startdepth = dem.where(dem.values <= tide_flow_df.Level.iloc[0], other=np.nan)
    startdepth = startdepth * -1

    # Fill values np.nan as 0
    startdepth.fillna(9999)
    startdepth = startdepth.where(startdepth.values >= 0, other=0)
    startdepth = startdepth.where(startdepth.values != 9999, other=np.nan)

    # Create startdepth raster (GeoTiff)
    startdepth.rio.to_raster(fr"{startdepth_tiff_path}\\startdepth_{number_simulation}.tif", cache=False)

    # For checking original if in case
    startdepth.rio.to_raster(fr"{startdepth_nc_path}\\startdepth_{number_simulation}.nc")

# END STARTDEPTH #######################################################################################################



# SIMULATE RASTER GENERATION FUNCTION ##################################################################################
def raster_generation(
        resolution_func,
        chunk_size_func,
        processor_func,
        padding_func,
        lidar_dataset_name,
        ran_trans_i
):
    """
    @Definition:
                A function to create a raster file from a las/laz file
    @References:
                https://github.com/rosepearson/GeoFabrics
    @Arguments:
                resolution_func (int or float):
                                            Resolution value in meter
                chunk_size_func (int):
                                            Size value of chunk
                processor_func (int):
                                            Number of processor
                number_simulation (string):
                                            A string to identify the order of simulation (should be angle, x, y)
                padding_func (list):
                                            A list of x min, x max, y min, and y max
                lidar_dataset_name (string):
                                            LiDAR name
                ran_trans_i (array):
                                            A 3D array contains values of angle, x, and y coordinates of points in tiles
                                            (iterating variable generated from multiprocessing)
    @Returns:
                None
    """
    # Get values
    number_simulation = ran_trans_i

    # Create DEM
    dem_raster(
        resolution_func,
        chunk_size_func,
        processor_func,
        number_simulation,
        padding_func,
        lidar_dataset_name
    )

    # Create MANNING'S N
    # Roughness
    roughness_raster(
        resolution_func,
        chunk_size_func,
        processor_func,
        number_simulation,
        padding_func,
        lidar_dataset_name
    )
    # Manning's n
    zo_to_n(
        number_simulation,
        1
    )

    # Create STARTDEPTH
    startdepth_generation(
        fr"{dem_nc_path}\\generated_dem_bathy_{number_simulation}.nc",
        number_simulation
    )

    # Convert netCDF to GeoTIFF
    nc_to_tiff(
        number_simulation
    )

    # Change value of sea
    value_sea_change(
        number_simulation,
        polygon=True # need to adjust manually
    )

    # Convert GeoTiff into ASCII
    convert_to_asc(number_simulation)

# END SIMULATE RASTER GENERATION FUNCTION ##############################################################################


