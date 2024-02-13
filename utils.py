import pathlib
from typing import List

import ee
import geopandas as gpd
import numpy as np
import pyproj
import shapely
import tqdm
import xarray as xr

from datamodels import get_equi7_tiles, get_sensor_products

# Initialize Earth Engine with high volume
ee.Initialize(opt_url="https://earthengine-highvolume.googleapis.com")


def generate_grids(
    tile_name: str = "AF",
    tile_extent: float = 500 * 128,
    dst_path: str = "./output",
) -> List[List[float]]:
    """
    Generate a grid of points for a given Equi7Grid tile name.

    Args:
    tile_name (str): The Equi7Grid tile name.
    tile_extent (float): The extent of the grid cell in meters.
    dst_path (str): The destination path to save the grid.

    Returns:
    List[List[float]]: A list of lists with the coordinates of the
    grid cells.
    """

    # Get the tile grid
    equi7 = get_equi7_tiles()
    tiles = equi7.tiles[tile_name]
    gdf = gpd.read_file(tiles.grid_path)

    # Get the PROJ4-string
    with open(tiles.prj_path, "r") as file:
        crs = file.read()

    # Make a grid from the layer
    xmin, ymin, xmax, ymax = gdf.total_bounds

    # create the cells in a loop
    grid_cells = []
    for x0 in np.arange(xmin, xmax + tile_extent, tile_extent):
        for y0 in np.arange(ymin, ymax + tile_extent, tile_extent):
            # bounds
            x1 = x0 - tile_extent
            y1 = y0 + tile_extent

            # create the shapely geometry
            cell_box = shapely.geometry.box(x0, y0, x1, y1)

            # append the centroid to the list
            grid_cells.append(cell_box)

    # Create a GeoDataFrame
    cell = gpd.GeoDataFrame(geometry=grid_cells)

    # Intersect with the land tile layer
    gdf_land = gdf[gdf["COVERSLAND"] == 1]
    cell_land = cell[cell.intersects(gdf_land.unary_union)]

    # Intersect with the countries layer
    countries = gpd.read_file(tiles.land_path)
    cell_countries = cell_land[cell_land.intersects(countries.unary_union)]

    # Get centroids
    cell_points = cell_countries.copy()
    cell_points["geometry"] = cell_points["geometry"].centroid

    # Set projection
    cell_points.set_crs(crs, inplace=True)

    # Generate a buffer around each point
    buffer = cell_points.buffer(tile_extent / 2).envelope
    buffer_list = []
    for geom in tqdm.tqdm(buffer):
        minx, miny, maxx, maxy = geom.bounds
        polygon = [[minx, miny], [minx, maxy], [maxx, maxy], [maxx, miny]]

        # Convert to latlon
        transformer = pyproj.Transformer.from_crs(
            crs, "EPSG:4326", always_xy=True
        )
        latlon_coords = [
            list(transformer.transform(x[0], x[1])) for x in polygon
        ]

        # append to list
        buffer_list.append(latlon_coords)

    # Save buffer list as GeoJSON
    gdf_buffer = gpd.GeoDataFrame(
        data={"id": range(len(buffer_list))},
        geometry=gpd.GeoSeries(
            [shapely.Polygon(coords) for coords in buffer_list]
        ),
        crs="EPSG:4326",
    )

    dst_path = pathlib.Path(dst_path)
    folder_path = dst_path / tile_name
    folder_path.mkdir(parents=True, exist_ok=True)

    gdf_buffer.to_file(
        f"{folder_path}/{tile_name}_buffers.geojson", driver="GeoJSON"
    )


def get_datacube(
    buffer_shape: shapely.geometry.Polygon,
    tile_name: str = "AF",
    product: str = "MODIS_061_MCD43A4",
    start_date: str = "2019-01-15",
    end_date: str = "2019-12-31",
) -> xr.Dataset:
    """
    Get a datacube from Earth Engine for a given buffer shape.

    Args:
    buffer_shape (shapely.geometry.Polygon): The buffer shape.
    tile_name (str): The Equi7Grid tile name.
    product (str): The product name.
    start_date (str): The start date.
    end_date (str): The end date.

    Returns:
    xr.Dataset: The datacube for the given buffer polygon.
    """

    # Convert shapely.geometry.Polygon to ee.Geometry.Polygon
    buffer_shape = list(buffer_shape.exterior.coords.xy)
    buffer_coords = [[x, y] for x, y in zip(buffer_shape[0], buffer_shape[1])]
    polygon = ee.Geometry.Polygon(buffer_coords)

    # Get AOI
    aoi = polygon.getInfo()

    # Generate BBox
    west = aoi["coordinates"][0][0][0]
    south = aoi["coordinates"][0][0][1]
    east = aoi["coordinates"][0][2][0]
    north = aoi["coordinates"][0][2][1]

    BBox = ee.Geometry.BBox(west, south, east, north)

    # Get the MODIS collection
    sensor_metadata = get_sensor_products()
    snippet = sensor_metadata.sensors[product].snippet
    bands = sensor_metadata.sensors[product].bands
    resolution = sensor_metadata.sensors[product].resolution
    edge_size = sensor_metadata.sensors[product].edge_size

    collection_filter = (
        ee.ImageCollection(snippet)
        .filterBounds(BBox)
        .filterDate(start_date, end_date)
        .select(bands)
    )

    # Get the CRS
    equi7 = get_equi7_tiles()
    tiles = equi7.tiles[tile_name]
    with open(tiles.prj_path, "r") as file:
        crs = file.read()

    # Get datacube from Earth Engine
    cube = xr.open_dataset(
        collection_filter, engine="ee", geometry=BBox, scale=resolution, crs=crs
    )

    # Rename and transpose dims
    cube = (
        cube.rename(Y="y", X="x")
        .to_array("band")
        .transpose("time", "band", "y", "x")
    )

    # Delete attributes to the cube
    cube.attrs = dict()

    # Add attributes to the cube
    collection_id = collection_filter.get("system:id").getInfo()

    cube.attrs = dict(
        collection=collection_id,
        stac="https://earthengine-stac.storage.googleapis.com/catalog/catalog.json",
        crs=crs,
        resolution=resolution,
        edge_size=edge_size,
        central_lat=(north + south) / 2,
        central_lon=(east + west) / 2,
        time_start=start_date,
        time_end=end_date,
    )

    # Rename the cube
    cube.name = collection_id

    return cube
