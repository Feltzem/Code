# Public Transport Study - Code Components

This directory contains the data processing, analysis, and modeling components for analyzing urban areas and their metro train systems. The project focuses on defining city boundaries based on urban characteristics and evaluating the quality and coverage of their metro networks.

## Project Overview

The analysis workflow consists of three main components:

### 1. Data Collection and Processing
- **Population Density**: Processing census and population data into standardized hexagonal grids
- **Building Footprints**: Retrieving and processing building data from Microsoft Buildings and OpenStreetMap
- **Road Networks**: Extracting and analyzing road density from OpenStreetMap
- **Metro Systems**: Collecting metro line and station data from OpenStreetMap

### 2. Urban Area Detection
- **Spatial Analysis**: Aggregating multiple urban indicators (population, buildings, roads) per hexagonal cell
- **Region Growing Algorithm**: Identifying contiguous urban areas using multi-criteria analysis
- **Boundary Definition**: Creating precise urban area boundaries for metro coverage analysis

### 3. Metro System Analysis
- **Network Coverage**: Calculating population coverage of metro systems
- **Accessibility Analysis**: Determining walking distances to stations
- **Quality Metrics**: Assessing metro network characteristics and performance

## Script and Notebook Descriptions

### Core Analysis Notebooks
- **`calc_urban_area.ipynb`** (Calculate Urban Area): Interactive notebook for comprehensive urban area analysis using population density, building footprints, and road networks. Implements a region-growing algorithm to identify contiguous urban areas with visualization and detailed reporting.
- **`model_census.ipynb`** (Model Census Data): Processing and modeling of population/census data, to eventually find Functional Urban Areas (FUA's) if more appropriate for PT coverage.
- **`get_metro_data.ipynb`** (Get Metro Data): Retrieval and processing of metro system data from OpenStreetMap, including lines, stations, and network topology.
- **`osm_experiments.ipynb`** (OSM Experiments): Experimental notebook for testing and refining OpenStreetMap data extraction methods and exploring different query approaches.

### Supporting Scripts
- **`road-density.py`**: Calculates road density per hexagonal grid cell for a specific city (Oslo). Queries OpenStreetMap for road segments, overlays them onto population hexagons, and computes road density metrics.
- **`get_ms_buildings.py`**: Retrieves building footprint data from Microsoft Buildings dataset, processes geometries, and prepares building data for urban analysis.

### Data Flow
1. **Input Data**: Raw census, building, and boundary data from external sources
2. **Processing**: Scripts convert data into standardized hexagonal grids and extract OSM features
3. **Analysis**: Notebooks perform spatial aggregation and urban area detection
4. **Output**: Processed datasets saved to `./data/` folder in GeoPackage format
5. **Visualization**: Data imported into ArcGIS Pro for cartographic analysis and presentation

### Output Data Structure
All processed data is saved to the `./data/` folder in GeoPackage (.gpkg) format for compatibility with both Python geospatial libraries and ArcGIS Pro:
- **Population density grids**: Hexagonal cells with population statistics
- **Building aggregations**: Building counts and areas per hexagon
- **Road networks**: Road segments with length and density metrics
- **Urban area boundaries**: Detected urban areas and individual cells
- **Metro networks**: Processed metro lines and stations

### Directory Structure
- **`cache/`**: Stores cached data from external APIs (OSM, Microsoft Buildings) to improve performance on repeated runs
- **`data/`**: Contains all processed datasets in GeoPackage format, ready for analysis in ArcGIS Pro
- **`Workspace.code-workspace`**: VS Code workspace configuration for the project

## Visualization and Mapping
Final visualization, cartographic analysis, and presentation mapping are performed in **ArcGIS Pro** using the processed datasets from the `./data/` folder. This separation allows for:
- Advanced cartographic styling and layout design
- Integration with additional spatial datasets
- Professional-quality map production for reports and presentations
- Leveraging ArcGIS Pro's specialized urban analysis tools
