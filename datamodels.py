import pathlib
from typing import Dict, List, Optional, Union

from pydantic import BaseModel


class Equi7GridTiles(BaseModel):
    grid_path: pathlib.Path
    land_path: pathlib.Path
    prj_path: pathlib.Path
    tile_extent: float


class Equi7Grid(BaseModel):
    tiles: Dict[str, Equi7GridTiles]


class SensorProduct(BaseModel):
    snippet: Union[Dict, str]
    bands: Optional[List[str]]
    resolution: float
    edge_size: int
    start_date: Optional[str]
    end_date: Optional[str]


class Sensors(BaseModel):
    sensors: Dict[str, SensorProduct]


def get_equi7_tiles(
    tile_names: Optional[List[str]] = [
        "AF",
        "AN",
        "AS",
        "EU",
        "NA",
        "OC",
        "SA",
    ],
    tile_extent: float = 500 * 128,
) -> Equi7Grid:
    # Create a dictionary of tiles
    tiles = dict()

    # Add the tiles to the dictionary
    if tile_names is not None:
        for tile_name in tile_names:
            tiles[tile_name] = {
                "grid_path": pathlib.Path(
                    f"./equi7grid/{tile_name}/{tile_name}_T6.geojson"
                ),
                "land_path": pathlib.Path(
                    f"./equi7grid/{tile_name}/LAND_{tile_name}.geojson"
                ),
                "prj_path": pathlib.Path(
                    f"./equi7grid/{tile_name}/WKT_{tile_name}.prj"
                ),
                "tile_extent": tile_extent,
            }

    return Equi7Grid(tiles=tiles)


def get_sensor_products(
    product: Optional[List[str]] = [
        "MODIS_061_MCD43A4",
        "WorldCLIM_v2",
        "MODIS_061_MOD11A1",
        "BNU_FGS_CCNL_v1",
        "fabdem",
        "geomorpho90",
        "s1gbm",
        "isric",
        "global_salinity",
        "ai0",
        "gwa",
        "ESA_WorldCover_v200",
        "OpenLandMap",
        "Oxford_MAP",
    ]
) -> Sensors:
    # Create a dictionary of sensors
    sensors = dict()

    # Add the sensors to the dictionary
    if "MODIS_061_MCD43A4" in product:
        sensors["MODIS_061_MCD43A4"] = {
            "snippet": "MODIS/061/MCD43A4",
            "bands": [
                "Nadir_Reflectance_Band1",
                "Nadir_Reflectance_Band2",
                "Nadir_Reflectance_Band3",
                "Nadir_Reflectance_Band4",
                "Nadir_Reflectance_Band5",
                "Nadir_Reflectance_Band6",
                "Nadir_Reflectance_Band7",
                "BRDF_Albedo_Band_Mandatory_Quality_Band1",
                "BRDF_Albedo_Band_Mandatory_Quality_Band2",
                "BRDF_Albedo_Band_Mandatory_Quality_Band3",
                "BRDF_Albedo_Band_Mandatory_Quality_Band4",
                "BRDF_Albedo_Band_Mandatory_Quality_Band5",
                "BRDF_Albedo_Band_Mandatory_Quality_Band6",
                "BRDF_Albedo_Band_Mandatory_Quality_Band7",
            ],
            "resolution": 500,
            "edge_size": 128,
            "start_date": "2000-02-24",
            "end_date": None,
        }

    if "WorldCLIM_v2" in product:
        sensors["WorldCLIM_v2"] = {
            "snippet": "WORLDCLIM/V1/MONTHLY",
            "bands": ["tavg", "tmin", "tmax", "prec"],
            "resolution": 1000,
            "edge_size": 64,
            "start_date": "1960-01-01",
            "end_date": "1990-01-01",
        }

    if "MODIS_061_MOD11A1" in product:
        sensors["MODIS_061_MOD11A1"] = {
            "snippet": "MODIS/061/MOD11A1",
            "bands": ["LST_Day_1km", "LST_Night_1km", "QC_Day", "QC_Night"],
            "resolution": 1000,
            "edge_size": 64,
            "start_date": "2000-02-24",
            "end_date": None,
        }

    if "BNU_FGS_CCNL_v1" in product:
        sensors["BNU_FGS_CCNL_v1"] = {
            "snippet": "BNU/FGS/CCNL/v1",
            "bands": ["b1"],
            "resolution": 1000,
            "edge_size": 64,
            "start_date": "1992-01-01",
            "end_date": "2014-01-01",
        }

    if "fabdem" in product:
        sensors["fabdem"] = {
            "snippet": "projects/sat-io/open-datasets/FABDEM",
            "bands": [],
            "resolution": 30,
            "edge_size": 2133,
            "start_date": None,
            "end_date": None,
        }

    if "geomorpho90" in product:
        sensors["geomorpho90"] = {
            "snippet": "projects/sat-io/open-datasets/Geomorpho90m/geom",
            "bands": [],
            "resolution": 250,
            "edge_size": 256,
            "start_date": None,
            "end_date": None,
        }

    if "s1gbm" in product:
        sensors["s1gbm"] = {
            "snippet": {
                "VH": "projects/sat-io/open-datasets/S1GBM/normalized_s1_backscatter_VH",
                "VV": "projects/sat-io/open-datasets/S1GBM/normalized_s1_backscatter_VV",
            },
            "bands": [],
            "resolution": 100,
            "edge_size": 2560,
            "start_date": None,
            "end_date": None,
        }

    if "isric" in product:
        sensors["isric"] = {
            "snippet": {
                key: f"projects/soilgrids-isric/{key}"
                for key in [
                    "bdod_mean",
                    "cec_mean",
                    "cvfo_mean",
                    "clay_mean",
                    "sand_mean",
                    "silt_mean",
                    "nitrogen_mean",
                    "phh2o_mean",
                    "soc_mean",
                    "ocd_mean",
                    "ocs_mean",
                ]
            },
            "bands": [],
            "resolution": 250,
            "edge_size": 256,
            "start_date": None,
            "end_date": None,
        }

    if "global_salinity" in product:
        sensors["global_salinity"] = {
            "snippet": "projects/sat-io/open-datasets/global_soil_salinity",
            "bands": [],
            "resolution": 250,
            "edge_size": 256,
            "start_date": None,
            "end_date": None,
        }

    if "ai0" in product:
        sensors["ai0"] = {
            "snippet": {
                key: f"projects/sat-io/open-datasets/global_ai/{key}"
                for key in ["global_ai_yearly", "global_ai_monthly"]
            },
            "bands": [],
            "resolution": 1000,
            "edge_size": 64,
            "start_date": None,
            "end_date": None,
        }

    if "gwa" in product:
        sensors["gwa"] = {
            "snippet": {
                key: f"projects/earthengine-legacy/assets/projects/sat-io/open-datasets/global_wind_atlas/{key}"
                for key in [
                    "air-density",
                    "capacity-factor",
                    "power-density",
                    "ruggedness-index",
                    "wind-speed",
                ]
            },
            "bands": [],
            "resolution": 1000,
            "edge_size": 64,
            "start_date": None,
            "end_date": None,
        }

    if "ESA_WorldCover_v200" in product:
        sensors["ESA_WorldCover_v200"] = {
            "snippet": "ESA/WorldCover/v200",
            "bands": ["Map"],
            "resolution": 10,
            "edge_size": 6400,
            "start_date": "2021-01-01",
            "end_date": "2022-01-01",
        }

    if "OpenLandMap" in product:
        sensors["OpenLandMap"] = {
            "snippet": "OpenLandMap/PNV/PNV_BIOME-TYPE_BIOME00K_C/v01",
            "bands": ["biome_type"],
            "resolution": 1000,
            "edge_size": 64,
            "start_date": "2001-01-01",
            "end_date": "2002-01-01",
        }

    if "Oxford_MAP" in product:
        sensors["Oxford_MAP"] = {
            "snippet": "projects/sat-io/open-datasets/oxford_map",
            "bands": ["map"],
            "resolution": 5000,
            "edge_size": 64,
            "start_date": "2001-01-01",
            "end_date": "2002-01-01",
        }

    return Sensors(sensors=sensors)
