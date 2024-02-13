import pathlib

import geopandas as gpd
import tqdm

from utils import generate_grids, get_datacube

tile_name = "SA"  # From ["AF", "AN", "AS", "EU", "NA", "OC", "SA"]
shape_n = 292  # Aleatory shape number

# Generate tile grid
generate_grids(tile_name=tile_name, tile_extent=500 * 128, dst_path="./output")

# Testing the datacube with single shape of the grid
grid = gpd.read_file(f"./output/{tile_name}/{tile_name}_buffers.geojson")
grid = grid[grid["id"] == shape_n]
shape = grid.geometry.values[0]

## Get the datacube
cube = get_datacube(
    buffer_shape=shape,
    tile_name=tile_name,
    product="MODIS_061_MCD43A4",
    start_date="2019-01-01",
    end_date="2019-03-01",
)

# Get dates
dates = cube["time"].dt.strftime("%Y_%m_%d").values.tolist()

# Download all images from the datacube (59 images in approximately 1 second)
for i in tqdm.tqdm(range(len(cube))):
    date = dates[i]
    folder_path = f"./output/{tile_name}/{shape_n}"
    pathlib.Path(folder_path).mkdir(parents=True, exist_ok=True)

    cube[i].rio.to_raster(
        f"{folder_path}/{date}.tif",
        driver="GTiff",
        dtype="float32",
        tiled=True,
        blockxsize=256,
        blockysize=256,
        compress="LZW",
    )


