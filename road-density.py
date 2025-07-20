import geopandas as gpd
import osmnx as ox
import pandas as pd
import fiona

# Configuration
city_name = "Oslo, Norway"
city_short_name = "oslo"
target_crs = "EPSG:27393"

print(f"Testing road density calculation for {city_name}")

# Load city boundary
print("Loading city boundary...")
cities_used = gpd.read_file(
    r"D:/Public Transport Study/Data/Boundaries/Boundaries.gpkg", layer="cities_used"
)
cities_used = cities_used.to_crs(target_crs)
selected_city_polygon = cities_used[cities_used["Name"] == city_name]

# Load population density data
print("Loading population density data...")
population_density_cities_used = gpd.read_file(
    r"./data/population_density.gpkg", layer="population_density_cities_used"
)
cities_used_boundaries = gpd.read_file(r"./data/urban_areas.gpkg", layer="cities_used")
selected_city_boundary = cities_used_boundaries[
    cities_used_boundaries["Name"] == city_name
]

# Filter population density to city boundary
if population_density_cities_used.crs != selected_city_boundary.crs:
    population_density_cities_used = population_density_cities_used.to_crs(
        selected_city_boundary.crs
    )

population_density_city = gpd.overlay(
    population_density_cities_used, selected_city_boundary, how="intersection"
)
print(f"Population hexagons: {len(population_density_city)}")

# Add hex_id if not exists
if "hex_id" not in population_density_city.columns:
    population_density_city["hex_id"] = range(len(population_density_city))

# Query road network
print("Querying road network from OSM...")
query_polygon_4326 = selected_city_polygon.to_crs("EPSG:4326")["geometry"].values[0]

road_types = {
    "highway": [
        "motorway",
        "trunk",
        "primary",
        "secondary",
        "tertiary",
        "unclassified",
        "residential",
        "motorway_link",
        "trunk_link",
        "primary_link",
        "secondary_link",
        "tertiary_link",
    ]
}

city_roads_4326 = ox.features_from_polygon(query_polygon_4326, tags=road_types)
city_roads_4326 = city_roads_4326[city_roads_4326.geometry.type == "LineString"]

print(f"Queried {len(city_roads_4326)} road segments")

# Reproject to target CRS
city_roads = city_roads_4326.to_crs(target_crs)
city_roads["length_m"] = city_roads.geometry.length

print(f"Total road length: {city_roads['length_m'].sum() / 1000:.2f} km")

# Calculate road density per hexagon
print("Calculating road density per hexagon...")

# Reproject roads to match population density CRS for overlay
target_overlay_crs = population_density_city.crs
city_roads_proj = city_roads.to_crs(target_overlay_crs)
city_roads_proj["road_id"] = range(len(city_roads_proj))

# Perform overlay
road_overlay_result = gpd.overlay(
    city_roads_proj[["road_id", "geometry"]],
    population_density_city,
    how="intersection",
    keep_geom_type=False,
)

# Calculate intersected lengths
road_overlay_result["intersected_length_m"] = road_overlay_result.geometry.length

# Aggregate per hexagon
road_stats_per_hex = road_overlay_result.groupby("hex_id").agg(
    total_road_length_m=("intersected_length_m", "sum"),
    road_segment_count=("road_id", pd.Series.nunique),
)

# Merge back to population data
population_density_city = population_density_city.merge(
    road_stats_per_hex, on="hex_id", how="left"
)

# Fill NaN values
population_density_city["total_road_length_m"] = population_density_city[
    "total_road_length_m"
].fillna(0)
population_density_city["road_segment_count"] = (
    population_density_city["road_segment_count"].fillna(0).astype(int)
)

# Calculate road density
hex_area_sqm = population_density_city.geometry.iloc[0].area
hex_area_sqkm = hex_area_sqm / 1_000_000

population_density_city["road_density_km_per_sqkm"] = (
    population_density_city["total_road_length_m"] / 1000
) / hex_area_sqkm

# Print results
print(f"\nRoad density results:")
print(
    f"Hexagons with roads: {len(population_density_city[population_density_city['road_segment_count'] > 0])}"
)
print(
    f"Total road length: {population_density_city['total_road_length_m'].sum() / 1000:.2f} km"
)
print(
    f"Average road density: {population_density_city['road_density_km_per_sqkm'].mean():.2f} km/sq km"
)
print(
    f"Max road density: {population_density_city['road_density_km_per_sqkm'].max():.2f} km/sq km"
)

# Save road data properly
print("Saving road data...")
roads_gpkg_path = "./data/road_networks_test.gpkg"
roads_layer_name = f"{city_short_name}_roads"

city_roads.to_file(roads_gpkg_path, layer=roads_layer_name, driver="GPKG")

# Save enhanced population data
pop_density_gpkg_path = "./data/population_density.gpkg"
pop_density_layer_name = f"population_density_with_buildings_roads_{city_short_name}"

population_density_city.to_file(
    pop_density_gpkg_path, layer=pop_density_layer_name, driver="GPKG"
)

print("Data saved successfully!")
print(
    f"Road density columns available: {[col for col in population_density_city.columns if 'road' in col.lower()]}"
)
