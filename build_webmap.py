from pathlib import Path
import sys
import json
import math
import logging
from typing import Optional, Tuple, List, Dict

import geopandas as gpd
from shapely.geometry import Point, Polygon, MultiPolygon
from shapely.ops import unary_union
import fiona
import folium
from folium import FeatureGroup
from folium.plugins import HeatMap
from shapely.wkb import dumps as wkb_dumps, loads as wkb_loads


def setup_logger(log_path: Path) -> logging.Logger:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("webmap_builder")
    logger.setLevel(logging.INFO)
    # Clear existing handlers when re-running
    logger.handlers.clear()
    fh = logging.FileHandler(log_path, encoding="utf-8")
    fh.setLevel(logging.INFO)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    fmt = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    fh.setFormatter(fmt)
    ch.setFormatter(fmt)
    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger


def list_layers(gpkg: Path) -> List[str]:
    try:
        return list(fiona.listlayers(str(gpkg)))
    except Exception:
        return []


def pick_layer(layers: List[str], include_keys: List[str]) -> Optional[str]:
    keys = [k.lower() for k in include_keys]
    for lyr in layers:
        name = lyr.lower()
        if all(k in name for k in keys):
            return lyr
    # fallback: any one key match
    for lyr in layers:
        name = lyr.lower()
        if any(k in name for k in keys):
            return lyr
    return None


def load_urban_area(data_dir: Path, logger: logging.Logger) -> Optional[gpd.GeoDataFrame]:
    # First try urban_areas.gpkg
    ua_gpkg = data_dir / "urban_areas.gpkg"
    if ua_gpkg.exists():
        layers = list_layers(ua_gpkg)
        cand = (
            pick_layer(layers, ["urban", "area"])
            or pick_layer(layers, ["urban", "boundary"])
            or (layers[0] if layers else None)
        )
        if cand:
            try:
                gdf = gpd.read_file(ua_gpkg, layer=cand)
                logger.info(f"Loaded urban area from {ua_gpkg.name} layer '{cand}' with {len(gdf)} features")
                return gdf
            except Exception as e:
                logger.warning(f"Failed to read layer {cand} from urban_areas.gpkg: {e}")
    # Fallback: union of urban_cells*.geojson
    cells_files = [
        data_dir / "urban_cells_oslo_growing_refined.geojson",
        data_dir / "urban_cells_oslo_growing.geojson",
    ]
    polys = []
    for p in cells_files:
        if p.exists():
            try:
                cells = gpd.read_file(p)
                if cells.crs is None:
                    # assume projected from notebook target; leave as-is
                    pass
                geom = unary_union(cells.geometry)
                polys.append(geom)
                logger.info(f"Loaded urban cells from {p.name} ({len(cells)} features)")
                break
            except Exception as e:
                logger.warning(f"Failed to read {p.name}: {e}")
    if polys:
        gdf = gpd.GeoDataFrame(geometry=[polys[0]], crs=cells.crs)
        return gdf
    logger.warning("Urban area not found; map will render without it.")
    return None


def load_population_cells(data_dir: Path, logger: logging.Logger) -> Optional[gpd.GeoDataFrame]:
    # Prefer refined geojson cells if present
    for name in ["urban_cells_oslo_growing_refined.geojson", "urban_cells_oslo_growing.geojson"]:
        fp = data_dir / name
        if fp.exists():
            try:
                gdf = gpd.read_file(fp)
                if "population" not in gdf.columns:
                    logger.warning(f"'{name}' missing 'population' column; density/tooltip may be limited.")
                return gdf
            except Exception as e:
                logger.warning(f"Failed to read {name}: {e}")
    # Try population_density.gpkg if exists
    pd_gpkg = data_dir / "population_density.gpkg"
    if pd_gpkg.exists():
        layers = list_layers(pd_gpkg)
        cand = pick_layer(layers, ["population"]) or (layers[0] if layers else None)
        if cand:
            try:
                return gpd.read_file(pd_gpkg, layer=cand)
            except Exception as e:
                logger.warning(f"Failed to read {pd_gpkg.name}:{cand}: {e}")
    logger.warning("Population cells not found.")
    return None


def load_rail(data_dir: Path, logger: logging.Logger) -> Tuple[Optional[gpd.GeoDataFrame], Optional[gpd.GeoDataFrame]]:
    rail_gpkg = data_dir / "rail_lines.gpkg"
    if not rail_gpkg.exists():
        logger.warning("rail_lines.gpkg not found.")
        return None, None
    layers = list_layers(rail_gpkg)
    lines_layer = pick_layer(layers, ["rail", "lines"]) or pick_layer(layers, ["rail", "edge"]) or None
    stops_layer = pick_layer(layers, ["rail", "station"]) or pick_layer(layers, ["stop"]) or None
    lines = None
    stops = None
    if lines_layer:
        try:
            lines = gpd.read_file(rail_gpkg, layer=lines_layer)
            logger.info(f"Loaded rail lines: {lines_layer} ({len(lines)})")
        except Exception as e:
            logger.warning(f"Failed to read lines {lines_layer}: {e}")
    if stops_layer:
        try:
            stops = gpd.read_file(rail_gpkg, layer=stops_layer)
            logger.info(f"Loaded rail stops: {stops_layer} ({len(stops)})")
        except Exception as e:
            logger.warning(f"Failed to read stops {stops_layer}: {e}")
    return lines, stops


def load_isochrones(data_dir: Path, logger: logging.Logger) -> Optional[gpd.GeoDataFrame]:
    iso_gpkg = data_dir / "station_isochrones.gpkg"
    if not iso_gpkg.exists():
        logger.warning("station_isochrones.gpkg not found.")
        return None
    layers = list_layers(iso_gpkg)
    cand = pick_layer(layers, ["isochrone"]) or pick_layer(layers, ["15", "min"]) or (layers[0] if layers else None)
    if cand:
        try:
            gdf = gpd.read_file(iso_gpkg, layer=cand)
            logger.info(f"Loaded isochrones: {cand} ({len(gdf)})")
            return gdf
        except Exception as e:
            logger.warning(f"Failed to read isochrone layer {cand}: {e}")
    return None


def to_wgs84(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    if gdf is None:
        return None
    if gdf.crs is None:
        return gdf
    try:
        return gdf.to_crs("EPSG:4326")
    except Exception:
        return gdf


def add_urban_area_layer(m: folium.Map, gdf: gpd.GeoDataFrame) -> FeatureGroup:
    group = FeatureGroup(name="Urban Area", show=True)
    folium.GeoJson(
        to_wgs84(gdf),
        name="Urban Area",
        style_function=lambda f: {"color": "#1f78b4", "weight": 2, "fillColor": "#1f78b4", "fillOpacity": 0.1},
        tooltip=folium.features.GeoJsonTooltip(fields=[c for c in gdf.columns if c != "geometry"][:5]),
    ).add_to(group)
    group.add_to(m)
    return group


def add_population_layers(
    m: folium.Map, cells: gpd.GeoDataFrame
) -> Tuple[Optional[FeatureGroup], Optional[FeatureGroup]]:
    cells = cells.copy()
    # Compute density per km^2 if possible
    try:
        if cells.crs and cells.crs.is_projected:
            area_km2 = cells.geometry.area / 1_000_000.0
        else:
            tmp = cells.to_crs(3857)
            area_km2 = tmp.geometry.area / 1_000_000.0
        if "population" in cells.columns:
            cells["density_km2"] = cells["population"] / area_km2.replace(0, math.nan)
        else:
            cells["density_km2"] = area_km2 * 0.0
    except Exception:
        cells["density_km2"] = 0
    # Choropleth-like GeoJson layer with tooltip
    choropleth_group = FeatureGroup(name="Population (cells)", show=False)
    cells_wgs = to_wgs84(cells)

    def style_fn(feat):
        d = feat["properties"].get("density_km2", 0) or 0
        # simple quantile-ish color ramp
        if d > 20000:
            col = "#800026"
        elif d > 10000:
            col = "#BD0026"
        elif d > 5000:
            col = "#E31A1C"
        elif d > 2500:
            col = "#FC4E2A"
        elif d > 1000:
            col = "#FD8D3C"
        elif d > 500:
            col = "#FEB24C"
        elif d > 250:
            col = "#FED976"
        else:
            col = "#FFEDA0"
        return {"color": "none", "fillColor": col, "fillOpacity": 0.6}

    folium.GeoJson(
        cells_wgs,
        name="Population Cells",
        style_function=style_fn,
        tooltip=folium.features.GeoJsonTooltip(
            fields=[c for c in ["population", "density_km2"] if c in cells_wgs.columns],
            aliases=["Population", "Density/km²"],
            localize=True,
        ),
    ).add_to(choropleth_group)
    choropleth_group.add_to(m)

    # HeatMap layer using centroids weighted by population
    heat_group = FeatureGroup(name="Population HeatMap", show=False)
    # Compute centroids in a metric CRS, then convert to WGS84 for mapping
    try:
        cent_metric = cells.to_crs(3857).geometry.centroid
        pts = gpd.GeoSeries(cent_metric, crs=3857).to_crs(4326)
    except Exception:
        pts = cells_wgs.geometry.centroid
    weights = cells_wgs["population"] if "population" in cells_wgs.columns else None
    heat_data = []
    for i, p in enumerate(pts):
        if p.is_empty:
            continue
        if weights is not None:
            w = float(weights.iloc[i]) if not math.isnan(float(weights.iloc[i])) else 0.0
        else:
            w = 1.0
        heat_data.append([p.y, p.x, max(0.0, w)])
    if heat_data:
        HeatMap(heat_data, radius=12, blur=15, min_opacity=0.3).add_to(heat_group)
        heat_group.add_to(m)
    return choropleth_group, heat_group


def add_rail_layers(
    m: folium.Map, lines: gpd.GeoDataFrame, stops: gpd.GeoDataFrame
) -> Tuple[Optional[FeatureGroup], Optional[FeatureGroup]]:
    color_by_type = {
        "subway": "#1f77b4",
        "metro": "#1f77b4",
        "light_rail": "#2ca02c",
        "tram": "#ff7f0e",
        "train": "#9467bd",
        "rail": "#8c564b",
        "monorail": "#e377c2",
    }

    def color_for(val: Optional[str]) -> str:
        if not val:
            return "#444444"
        v = str(val).lower()
        for key, col in color_by_type.items():
            if key in v:
                return col
        return "#444444"

    line_group = FeatureGroup(name="Rail Lines", show=True)
    lines_wgs = to_wgs84(lines) if (lines is not None and not lines.empty) else None
    if lines_wgs is not None:
        folium.GeoJson(
            lines_wgs,
            name="Rail Lines",
            style_function=lambda f: {"color": color_for(f["properties"].get("rail_type")), "weight": 3},
            tooltip=folium.features.GeoJsonTooltip(
                fields=[c for c in ["name", "rail_type"] if c in lines_wgs.columns],
                aliases=["Name", "Rail Type"],
            ),
        ).add_to(line_group)
        line_group.add_to(m)

    stop_group = FeatureGroup(name="Rail Stops", show=True)
    stops_wgs = to_wgs84(stops) if (stops is not None and not stops.empty) else None
    if stops_wgs is not None:
        for _, r in stops_wgs.iterrows():
            geom = r.geometry
            if geom is None or geom.is_empty:
                continue
            name = r.get("name", "Stop")
            rtype = r.get("rail_type", r.get("railway", ""))
            color = color_for(rtype)
            folium.CircleMarker(
                location=[geom.y, geom.x],
                radius=4,
                color=color,
                fill=True,
                fill_color=color,
                fill_opacity=0.9,
                tooltip=f"{name} ({rtype})",
            ).add_to(stop_group)
        stop_group.add_to(m)
    return line_group, stop_group


def add_isochrones_layer(m: folium.Map, iso: gpd.GeoDataFrame) -> FeatureGroup:
    group = FeatureGroup(name="15-min Walk Isochrones", show=True)
    folium.GeoJson(
        to_wgs84(iso),
        name="Isochrones",
        style_function=lambda f: {"color": "#3182bd", "weight": 1, "fillColor": "#3182bd", "fillOpacity": 0.25},
    ).add_to(group)
    group.add_to(m)
    return group


def compute_reachability_stats(
    cells: Optional[gpd.GeoDataFrame], iso: Optional[gpd.GeoDataFrame], logger: logging.Logger
) -> Tuple[Optional[float], Optional[int], Optional[int]]:
    if cells is None or iso is None or cells.empty or iso.empty:
        return None, None, None
    # Project to a metric CRS for safe area/centroid operations
    metric_crs = 3857
    cc = cells.to_crs(metric_crs)
    ii = iso.to_crs(metric_crs)

    # Robust union of isochrone geometries
    def robust_union(geoms):
        parts = []
        for g in geoms:
            if g is None:
                continue
            try:
                if g.is_empty:
                    continue
            except Exception:
                continue
            try:
                g2 = wkb_loads(wkb_dumps(g, output_dimension=2))
            except Exception:
                g2 = g
            gt = getattr(g2, "geom_type", "")
            if gt in ("Polygon", "MultiPolygon"):
                parts.append(g2)
            elif gt == "GeometryCollection" and hasattr(g2, "geoms"):
                for p in g2.geoms:
                    if p is not None and not p.is_empty and getattr(p, "geom_type", "") in ("Polygon", "MultiPolygon"):
                        parts.append(p)
        if not parts:
            return None
        items = parts
        while len(items) > 1:
            next_items = []
            for i in range(0, len(items), 256):
                chunk = items[i : i + 256]
                u = None
                for q in chunk:
                    try:
                        u = q if u is None else u.union(q)
                    except Exception:
                        try:
                            qf = q.buffer(0)
                            u = qf if u is None else u.union(qf)
                        except Exception:
                            continue
                if u is not None and not u.is_empty:
                    next_items.append(u)
            items = next_items if next_items else items[:1]
        return items[0]

    geom_iso = robust_union(ii.geometry)
    if geom_iso is None or geom_iso.is_empty:
        logger.warning("Isochrone geometry union failed or empty; skipping stats.")
        return None, None, None
    cc = cc.copy()
    cc["_centroid"] = cc.geometry.centroid
    within_mask = cc["_centroid"].within(geom_iso)
    in_cells = cc[within_mask]
    tot_pop = int(cc["population"].sum()) if "population" in cc.columns else None
    in_pop = int(in_cells["population"].sum()) if "population" in cc.columns else None
    pct = (in_pop / tot_pop * 100.0) if (in_pop is not None and tot_pop and tot_pop > 0) else None
    logger.info(f"Reachability: {in_pop} / {tot_pop} people within 15-min walk ({pct:.2f}% if available)")
    return pct, in_pop, tot_pop


def add_stats_panel(m: folium.Map, pct: Optional[float], in_pop: Optional[int], tot_pop: Optional[int]):
    val_line = (
        f"<div><strong>PT Reachability:</strong> {pct:.2f}% of population within 15-min walk</div>"
        if pct is not None
        else "<div><strong>PT Reachability:</strong> N/A</div>"
    )
    pop_line = (
        f"<div><strong>Population covered:</strong> {in_pop:,} / {tot_pop:,}</div>"
        if in_pop is not None and tot_pop is not None
        else ""
    )
    html = f"""
    <div style='position: absolute; bottom: 20px; right: 20px; z-index: 9999; background: rgba(255,255,255,0.9); padding: 10px 12px; border: 1px solid #aaa; border-radius: 4px; box-shadow: 0 1px 4px rgba(0,0,0,0.3); font-family: Arial, sans-serif; font-size: 12px;'>
        <div style='font-weight: 600; margin-bottom: 4px;'>Public Transport Reachability</div>
        {val_line}
        {pop_line}
        <div style='margin-top:6px; font-size:11px;'>Toggle layers via the control panel (top-right).</div>
    </div>
    """
    m.get_root().html.add_child(folium.Element(html))


def main():
    # Paths
    script_dir = Path(__file__).resolve().parent
    root_dir = script_dir.parent
    data_dir = script_dir / "data"
    report_dir = root_dir / "Report"
    report_dir.mkdir(parents=True, exist_ok=True)
    log_path = report_dir / "webmap_build.log"
    logger = setup_logger(log_path)

    logger.info("Starting webmap builder...")

    # Load datasets
    ua_gdf = load_urban_area(data_dir, logger)
    pop_cells = load_population_cells(data_dir, logger)
    rail_lines, rail_stops = load_rail(data_dir, logger)
    iso_gdf = load_isochrones(data_dir, logger)

    # Determine initial map center
    center = [59.9139, 10.7522]  # default Oslo center as fallback
    if ua_gdf is not None and not ua_gdf.empty:
        try:
            cen = to_wgs84(ua_gdf).geometry.unary_union.centroid
            center = [float(cen.y), float(cen.x)]
        except Exception:
            pass
    elif pop_cells is not None and not pop_cells.empty:
        try:
            tb = to_wgs84(pop_cells).total_bounds
            center = [float((tb[1] + tb[3]) / 2.0), float((tb[0] + tb[2]) / 2.0)]
        except Exception:
            pass

    # Build map
    m = folium.Map(location=center, zoom_start=11, tiles="CartoDB positron")
    folium.TileLayer("OpenStreetMap").add_to(m)
    # Add a terrain-like layer with attribution
    folium.TileLayer(
        tiles="https://stamen-tiles.a.ssl.fastly.net/terrain/{z}/{x}/{y}.jpg",
        name="Stamen Terrain",
        attr="Map tiles by Stamen Design, CC BY 3.0 — Map data © OpenStreetMap contributors",
    ).add_to(m)

    # Add layers
    if ua_gdf is not None and not ua_gdf.empty:
        add_urban_area_layer(m, ua_gdf)
    if pop_cells is not None and not pop_cells.empty:
        add_population_layers(m, pop_cells)
    if rail_lines is not None and not rail_lines.empty or rail_stops is not None and not rail_stops.empty:
        add_rail_layers(m, rail_lines, rail_stops)
    if iso_gdf is not None and not iso_gdf.empty:
        add_isochrones_layer(m, iso_gdf)

    # Stats panel
    pct, in_pop, tot_pop = compute_reachability_stats(pop_cells, iso_gdf, logger)
    add_stats_panel(m, pct, in_pop, tot_pop)

    folium.LayerControl(position="topright", collapsed=False).add_to(m)

    out_html = report_dir / "webmap.html"
    m.save(str(out_html))
    logger.info(f"Saved web map to: {out_html}")
    print(str(out_html))


if __name__ == "__main__":
    main()
