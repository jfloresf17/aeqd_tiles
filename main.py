import pathlib

import geopandas as gpd

from utils import composite_collection, generate_grids

tile_name = "SA"  # From ["AF", "AN", "AS", "EU", "NA", "OC", "SA"]
tile_id = 292  # Aleatory shape number
dst_path = pathlib.Path("./output")

# Generate tile grid
generate_grids(tile_name=tile_name, tile_extent=500 * 128, dst_path=dst_path)

#  Load the generated grid
grid = gpd.read_file(f"./output/{tile_name}/{tile_name}_buffers.geojson")

# Filter by tile_id
grid = grid[grid["id"] == tile_id]

# Get the shapely geometry
geom = grid.geometry.values[0]

# Testing the composite collection for:

# - MODIS_061_MCD43A4
composite_collection(
    geom=geom,
    tile_name=tile_name,
    tile_id=tile_id,
    product="MODIS_061_MCD43A4",
    start_date="2014-01-01",
    end_date="2022-12-31",
    dst_folder=dst_path,
)

# - MODIS_061_MOD11A1
composite_collection(
    geom=geom,
    tile_name=tile_name,
    tile_id=tile_id,
    product="MODIS_061_MOD11A1",
    start_date="2014-01-01",
    end_date="2022-12-31",
    dst_folder=dst_path,
)
# - s1gbm
composite_collection(
    geom=geom,
    tile_name=tile_name,
    tile_id=tile_id,
    product="s1gbm",
    subproduct="VV",
    dst_folder=dst_path,
)

# - global_salinity
composite_collection(
    geom=geom,
    tile_name=tile_name,
    tile_id=tile_id,
    product="global_salinity",
    dst_folder=dst_path,
)
