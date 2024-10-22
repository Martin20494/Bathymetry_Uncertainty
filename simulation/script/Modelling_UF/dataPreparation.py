# Prepare packages -------------------------------------------------------------------------
# Folder packages
from folder import *

# Packages for path creation
import pathlib
import shutil

# Packages for downloading lidar
import geoapis
import geoapis.lidar

# Packages for converting lidar point cloud las/laz file into DEM raster netCDF file
import numpy as np
import shapely.geometry
import geopandas
import json
from geofabrics import processor

# Packages for manipulate raster
import rioxarray as rxr
from osgeo import gdal

import multiprocessing
# ------------------------------------------------------------------------------------------

# NESTED PARALLELISM --------------------------------------------------------------------------------------------------
"""
@Definition:
            A class to perform nested parallelism
@References:
            https://stackoverflow.com/questions/50937362/multiprocessing-on-python-3-jupyter
            https://stackoverflow.com/questions/66420735/simple-multiprocessing-pool-hangs-in-jupyter-notebook
            https://stackoverflow.com/questions/23641475/multiprocessing-working-in-python-but-not-in-ipython/23641560#23641560
            https://stackoverflow.com/questions/5442910/how-to-use-multiprocessing-pool-map-with-multiple-arguments
            https://stackoverflow.com/questions/20886565/using-multiprocessing-process-with-a-maximum-number-of-simultaneous-processes

            https://stackoverflow.com/questions/6974695/python-process-pool-non-daemonic
            https://www.geeksforgeeks.org/partial-functions-python/
"""


class NoDaemonProcess(multiprocessing.Process):
    @property
    def daemon(self):
        return False

    @daemon.setter
    def daemon(self, value):
        pass


class NoDaemonContext(type(multiprocessing.get_context())):
    Process = NoDaemonProcess


# We sub-class multiprocessing.pool.Pool instead of multiprocessing.Pool
# because the latter is only a wrapper function, not a proper class.
class MyPool(multiprocessing.pool.Pool):
    def __init__(self, *args, **kwargs):
        kwargs['context'] = NoDaemonContext()
        super(MyPool, self).__init__(*args, **kwargs)


# END NESTED PARALLELISM ----------------------------------------------------------------------------------------------


def download_lidar(boundary_coordinates_func, name_dataset, simulation=False):
    """
    @Definition:
                A function to download las/laz files from Open Topography
    @References:
                https://github.com/niwa/geoapis/wiki/Basic-Usage
    @Arguments:
                boundary_coordinates_func (list):
                                                    A list of xmin, ymin, xmax, ymax of rectangular boundary
                name_dataset (string):
                                                    Name of the dataset
                simulation (boolean):
                                                    False if downloading data for preparation
                                                    True if downloading data for simulation
    @Returns:
                None.
    """
    # Get crs
    h_crs = 2193

    # Get xmin, ymin, xmax, ymax of rectangular boundary
    x0_func = boundary_coordinates_func[0]
    y0_func = boundary_coordinates_func[1]
    x1_func = boundary_coordinates_func[2]
    y1_func = boundary_coordinates_func[3]

    # Assign those coordinates into shapely geometry polygon and store under geopandas format
    lidar_bound_coordinates = shapely.geometry.Polygon(
        [(x0_func, y0_func), (x1_func, y0_func), (x1_func, y1_func), (x0_func, y1_func)])
    lidar_bound_coordinates = geopandas.GeoSeries([lidar_bound_coordinates])
    lidar_bound_coordinates = lidar_bound_coordinates.set_crs(h_crs)

    # Create test_catchment.zip
    if simulation:
        test_path = pathlib.Path(f"{bathy_path}\\test_catchment")
    else:
        test_path = pathlib.Path(f"{original_lidar_path}\\test_catchment")
    lidar_bound_coordinates.to_file(test_path)
    shutil.make_archive(base_name=test_path, format='zip', root_dir=test_path)
    shutil.rmtree(test_path)
    lidar_bound_path = pathlib.Path(str(test_path) + ".zip")

    # Create polygon
    lidar_polygon = geopandas.read_file(lidar_bound_path)
    lidar_polygon.to_crs(h_crs)

    # Download lidar
    if simulation:
        lidar_fetcher = geoapis.lidar.OpenTopography(cache_path=fr"{bathy_path}",
                                                     search_polygon=lidar_polygon,
                                                     verbose=True)
    else:
        lidar_fetcher = geoapis.lidar.OpenTopography(cache_path=fr"{original_lidar_path}",
                                                     search_polygon=lidar_polygon,
                                                     verbose=True)
    lidar_fetcher.run(name_dataset)


def distance_calculation(point_1, point_2):
    """
    @Definition:
                A function to calculate the distance of two points using Euclidean method
    @Reference:
                https://orion.math.iastate.edu/dept/links/formulas/form2.pdf
    @Arguments:
                point_1 (list):
                                    First point
                point_2 (list):
                                    Second point
    @Returns:
                distance (float):
                                    Euclidean distance between two points
    """
    # Calculating the euclidean distance between two points
    return np.sqrt(np.power((point_1[0] - point_2[0]), 2) + np.power((point_1[1] - point_2[1]), 2))

def get_divisible_number(changed_number, divisible_number=16):
    """
    @Definition:
                A function to find the padding number that is divisible by 16
    @Reference:
                https://stackoverflow.com/questions/8002217/how-do-you-check-whether-a-number-is-divisible-by-another-number-python
    @Arguments:
                changed_number (float):
                                    Number that needs changing into the number that can be divisible by 16
                divisible_number (float):
                                    Default is 16
    @Return:
                changed_number (float):
                                    A new number that can divisible by 16
    """
    # Get the number that is divisible by 16
    while changed_number % divisible_number != 0:
        changed_number = changed_number + 1

    return changed_number

def padding_combination(coordinates_func, addition):
    """
    @Definition:
                A function to create the coordinates of padding/grid size
    @References:
                None.
    @Arguments:
                coordinates_func (list):
                                        A list of xmin, ymin, xmax, ymax
                addition (int):
                                        A number used to extend the padding
    @Returns:
                (list):
                                        A list of x min, x max, y min and y max
    """
    # Get xmin, ymin, xmax, ymax
    x_min = coordinates_func[0]
    y_min = coordinates_func[1]
    x_max = coordinates_func[2]
    y_max = coordinates_func[3]

    # Calculate the center of padding
    center_x_func = ((x_min + x_max) / 2)
    center_y_func = ((y_min + y_max) / 2)
    center_padding_func = np.array([center_x_func, center_y_func])

    # Calculate the radius
    # Choose randomly a corner coordinates of rectangle point cloud to calculate the distance with center point
    radius = distance_calculation(np.array([x_min, y_min]), center_padding_func)

    # Get the coordinates of middle points of 4 sides of the rectangle point cloud
    middle_point_1 = np.array([x_min, center_y_func])
    middle_point_2 = np.array([center_x_func, y_max])
    middle_point_3 = np.array([x_max, center_y_func])
    middle_point_4 = np.array([center_x_func, y_min])

    # Calculate the distances between these points and the center point
    distance_center_1 = distance_calculation(middle_point_1, center_padding_func)
    distance_center_2 = distance_calculation(middle_point_2, center_padding_func)
    distance_center_3 = distance_calculation(middle_point_3, center_padding_func)
    distance_center_4 = distance_calculation(middle_point_4, center_padding_func)

    # Calculate the distances between the above middle points and the middle points of 4 sides of the padding
    distance_1 = radius - distance_center_1
    distance_2 = radius - distance_center_2
    distance_3 = radius - distance_center_3
    distance_4 = radius - distance_center_4

    # Convert into numbers which will be divisible by 16
    distance_padding_1 = get_divisible_number(np.ceil(distance_1))
    distance_padding_2 = get_divisible_number(np.ceil(distance_2))
    distance_padding_3 = get_divisible_number(np.ceil(distance_3))
    distance_padding_4 = get_divisible_number(np.ceil(distance_4))

    # Calculate the coordinates of the middle points of 4 sides of the padding
    padding_middle_point_1 = np.array((middle_point_1[0] - distance_padding_1, center_y_func))
    padding_middle_point_2 = np.array((center_x_func, middle_point_2[1] + distance_padding_2))
    padding_middle_point_3 = np.array((middle_point_3[0] + distance_padding_3, center_y_func))
    padding_middle_point_4 = np.array((center_x_func, middle_point_4[1] - distance_padding_4))

    # Calculate the coordinates of 4 corners of the padding
    padding_x_min = np.min(np.array([padding_middle_point_1[0],
                                     padding_middle_point_2[0],
                                     padding_middle_point_3[0],
                                     padding_middle_point_4[0]])) - addition

    padding_x_max = np.max(np.array([padding_middle_point_1[0],
                                     padding_middle_point_2[0],
                                     padding_middle_point_3[0],
                                     padding_middle_point_4[0]])) + addition

    padding_y_min = np.min(np.array([padding_middle_point_1[1],
                                     padding_middle_point_2[1],
                                     padding_middle_point_3[1],
                                     padding_middle_point_4[1]])) - addition

    padding_y_max = np.max(np.array([padding_middle_point_1[1],
                                     padding_middle_point_2[1],
                                     padding_middle_point_3[1],
                                     padding_middle_point_4[1]])) + addition

    return [padding_x_min, padding_y_min, padding_x_max, padding_y_max]


def dem_raster_reference(resolution_func,
                         chunk_size_func, processor_func,
                         padding_func, lidar_dataset_name, element_name):
    """
    @Definition:
                A function to create a reference raster file from a las/laz file used for center calculation and
                padding reference
    @References:
                https://github.com/rosepearson/GeoFabrics
    @Arguments:
                resolution_func (int or float):
                                            Resolution value in meter
                chunk_size_func (int):
                                            Size value of chunk
                processor_func (int):
                                            Number of processor
                padding_func (list):
                                            A list of x min, x max, y min and y max
                lidar_dataset_name (string):
                                            LiDAR name
                element_name (string):
                                            Name of element folder
    @Returns:
                None.
    """

    # Set crs
    h_crs = 2193
    v_crs = 7839

    # Paths
    basepath = pathlib.Path(original_lidar_path)
    result_folder = element_name
    data_dir = basepath / result_folder
    data_dir.mkdir(parents=True, exist_ok=True)

    # Bounding box/ Catchment boundary
    x0 = padding_func[0]
    y0 = padding_func[1]
    x1 = padding_func[2]
    y1 = padding_func[3]
    catchment = shapely.geometry.Polygon([(x0, y0), (x1, y0), (x1, y1), (x0, y1)])
    catchment = geopandas.GeoSeries([catchment])
    catchment = catchment.set_crs(h_crs)
    catchment.to_file(basepath / result_folder / "catchment_boundary.geojson", crs="EPSG:2193", driver="GeoJSON")

    # Design JSON instructions
    instruction_json = {}

    # DRAINS
    instruction_json['drains'] = {
        # drain - output
        "output": {
            "crs": {
                "horizontal": h_crs,
                "vertical": v_crs
            },
            "grid_params": {
                "resolution": 1
            }
        },

        # drain - processing
        "processing": {
            "chunk_size": 1000,
            "number_of_cores": 4
        },

        # drain - data_paths
        "data_paths": {
            "local_cache": str(basepath),
            "subfolder": result_folder,
            "catchment_boundary": "catchment_boundary.geojson"
        },

        # drain - apis
        "apis": {
            "open_topography": {
                "Wellington_2013": {
                    "crs": {
                        "horizontal": h_crs,
                        "vertical": v_crs
                    }
                }
            },
            "linz": {
                "key": "857413f41ce446ed8961e2f1e960a24b",
                "land": {
                    "layers": [51559]
                }
            }
        },

        # drain - general
        "general": {
            "lidar_classifications_to_keep": [2, 9]
        },

        # drain - drains
        "drains": {
            "width": 5
        }
    }

    # RIVERS
    instruction_json['rivers'] = {
        # river - output
        "output": {
            "crs": {
                "horizontal": h_crs,
                "vertical": v_crs
            },
            "grid_params": {
                "resolution": 1
            }
        },

        # river - processing
        "processing": {
            "chunk_size": 1000,
            "number_of_cores": 4
        },

        # river - data_paths
        "data_paths": {
            "local_cache": str(basepath),
            "subfolder": result_folder,
        },

        # river - apis
        "apis": {
            "open_topography": {
                "Wellington_2013": {
                    "crs": {
                        "horizontal": h_crs,
                        "vertical": v_crs
                    }
                }
            },
            "linz": {
                "key": "857413f41ce446ed8961e2f1e960a24b",
                "land": {
                    "layers": [51559]
                },
                "bathymetry_contours": {
                    "layers": [50554]
                }
            }
        },

        # river - general
        "general": {
            "set_dem_shoreline": True,
            "bathymetry_contours_z_label": "valdco",
            "drop_offshore_lidar": False,
            "lidar_classifications_to_keep": [2, 9],
            "interpolate_missing_values": True
        },

        # river - rivers
        "rivers": {
            "osm_id": 132793862,
            "veg_lidar_classifications_to_keep": [2, 3, 4, 5, 9],
            "max_channel_width": 120,
            "min_channel_width": 10,
            "max_bank_height": 2,
            "rec_alignment_tolerance": 65,
            "width_centre_smoothing": 10,
            "channel_area_threshold": 19000000,
            "channel_rec_id": 9253579,
            "cross_section_spacing": 10,
            "min_bank_height": 0.75,
            "rec_file": str(basepath / "rec2_3.geojson"),
            "flow_file": str(basepath / "flow_and_friction.csv.gz")
        }
    }


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
            },
            "linz": {
                "key": "857413f41ce446ed8961e2f1e960a24b",
                "land": {
                    "layers": [51559]
                },
                "bathymetry_contours": {
                    "layers": [50554]
                }
            }
        },

        # data paths
        "data_paths": {
            "local_cache": str(basepath),
            "subfolder": result_folder,
            "catchment_boundary": "catchment_boundary.geojson",
            "raw_dem": "raw_dem.nc",

            # For bathymetry
            "river_bathymetry": ["river_bathymetry.geojson", "fan_bathymetry.geojson",
                                 "open_drain_elevation_5m_width.geojson", "closed_drain_elevation_5m_width.geojson"],
            "river_polygons": ["river_polygon.geojson", "fan_polygon.geojson",
                               "open_drain_polygon_5m_width.geojson", "closed_drain_polygon_5m_width.geojson"],
            # Result
            "result_dem": f"{element_name}.nc"
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
            "bathymetry_points_z_label": ['bed_elevation_Neal_et_al', "depths", "elevation", "elevation"]
        }
    }

    # Save the instructions
    with open(basepath/ result_folder / "instructions.json", "w") as instruction:
        json.dump(instruction_json, instruction)

    # Create DEM raster
    runner = processor.RiverBathymetryGenerator(instruction_json['rivers'])
    runner.run()
    runner = processor.DrainBathymetryGenerator(instruction_json['drains'])
    runner.run()
    runner = processor.RawLidarDemGenerator(instruction_json['dem'])
    runner.run()
    runner = processor.HydrologicDemGenerator(instruction_json['dem'])
    runner.run()


def terrain_shading(
    altitude,
    azimuth
):
    """
    @Definition:
                A function to create terrain shading
    @References:
                https://www.youtube.com/watch?v=5dDZeEXws9Q
                https://blog.datawrapper.de/shaded-relief-with-gdal-python/
                https://www.geophysique.be/2014/02/25/shaded-relief-map-in-python/
                https://www.l3harrisgeospatial.com/docs/topographicshading.html
    @Arguments:
                azimuth and altitude (float):
                            Values to customise the terrain shading. Ex: altitude = 45, azimuth = 355
    @Returns:
                None.
    """
    # Convert nc to tiff
    terrain = rxr.open_rasterio(fr"{original_lidar_path}\\shading\\shading.nc")
    terrain.z.rio.to_raster(fr"{original_lidar_path}\\shading\\shading.tiff")

    # Create shading
    terrain_tiff = gdal.Open(fr"{original_lidar_path}\\shading\\shading.tiff")
    terrain_shade = gdal.DEMProcessing(
        fr"{original_lidar_path}\\shading\\terrain_shading.tiff",
        terrain_tiff,
        "hillshade",
        computeEdges=True,
        altitude=altitude, azimuth=azimuth
    )

    # Close the file
    terrain_shade = None
    terrain_tiff = None


