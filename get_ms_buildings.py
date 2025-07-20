# %%
import pandas as pd
import geopandas as gpd
from shapely import geometry
import mercantile
from tqdm import tqdm
import os
import tempfile

# %%
cities_used = gpd.read_file(
    r"D:/Public Transport Study/Data/Boundaries/Boundaries.gpkg", layer="cities_used"
)
# Filter to Oslo, Norway
cities_used = cities_used[cities_used["Name"] == "Stockholm, Sweden"]
cities_used = cities_used.to_crs("EPSG:32633")
print(cities_used.head())

aoi_shape = cities_used.unary_union
minx, miny, maxx, maxy = aoi_shape.bounds
print(minx, miny, maxx, maxy)

cities_used.explore()

# %%
quad_keys = set()
for tile in list(mercantile.tiles(minx, miny, maxx, maxy, zooms=9)):
    quad_keys.add(mercantile.quadkey(tile))
quad_keys = list(quad_keys)
print(f"The input area spans {len(quad_keys)} tiles: {quad_keys}")

# %%
df = pd.read_csv(
    "https://minedbuildings.z5.web.core.windows.net/global-buildings/dataset-links.csv",
    dtype=str,
)
df.head()

# %%
idx = 0
combined_gdf = gpd.GeoDataFrame()
with tempfile.TemporaryDirectory() as tmpdir:
    # Download the GeoJSON files for each tile that intersects the input geometry
    tmp_fns = []
    for quad_key in tqdm(quad_keys):
        rows = df[df["QuadKey"] == quad_key]
        if rows.shape[0] == 1:
            url = rows.iloc[0]["Url"]

            df2 = pd.read_json(url, lines=True)
            df2["geometry"] = df2["geometry"].apply(geometry.shape)

            gdf = gpd.GeoDataFrame(df2, crs=4326)
            fn = os.path.join(tmpdir, f"{quad_key}.geojson")
            tmp_fns.append(fn)
            if not os.path.exists(fn):
                gdf.to_file(fn, driver="GeoJSON")
        elif rows.shape[0] > 1:
            print(f"Multiple rows found for QuadKey: {quad_key}, using the first one.")
            url = rows.iloc[0]["Url"]

            df2 = pd.read_json(url, lines=True)
            df2["geometry"] = df2["geometry"].apply(geometry.shape)

            gdf = gpd.GeoDataFrame(df2, crs=4326)
            fn = os.path.join(tmpdir, f"{quad_key}.geojson")
            tmp_fns.append(fn)
            if not os.path.exists(fn):
                gdf.to_file(fn, driver="GeoJSON")
        else:
            print(f"QuadKey not found in dataset: {quad_key}, skipping.")
            continue

    # Merge the GeoJSON files into a single file
    for fn in tmp_fns:
        gdf = gpd.read_file(fn)  # Read each file into a GeoDataFrame
        gdf = gdf[gdf.geometry.within(aoi_shape)]  # Filter geometries within the AOI
        gdf["id"] = range(idx, idx + len(gdf))  # Update 'id' based on idx
        idx += len(gdf)
        combined_gdf = pd.concat([combined_gdf, gdf], ignore_index=True)

# %%
# combined_gdf.explore()

# Write to building_footprints geopackage
combined_gdf.to_file(
    "./data/building_footprints.gpkg",
    layer="building_footprints_stockholm",
    driver="GPKG",
)

# %%
combined_gdf.size()
# %%
