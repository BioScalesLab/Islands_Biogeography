# Data Dictionary

## Overview
This directory contains datasets used in the Islands Biogeography analysis project. Data is organized by ecosystem type (marine/terrestrial) and diversity metric (taxonomic/functional).

## File Descriptions

### Marine Ecosystems

#### `model_data_coral_taxonomic_arch.csv`
**Description**: Coral species composition (taxonomic diversity) data across island archipelagos

**Variables**:
- `Archipelago_ID`: Unique identifier for archipelago
- `Island_ID`: Unique identifier for island
- `Distance`: Geographic distance between islands (km)
- `diff_isolation`: Quaternary isolation difference (time since last connection)
- `diff_past_isolation`: Past Quaternary isolation metric
- `diff_reef_area`: Difference in reef area between islands (km²)
- `diff_past_reef_area`: Historical reef area difference
- `beta_sor`: Sørensen beta-diversity index (species turnover)
- `Region`: Geographic region classification

**Data Type**: Community composition matrix
**Number of Observations**: Island-pair comparisons
**Source**: Reef survey data

---

#### `model_data_coral_functional_arch.csv`
**Description**: Coral trait-based (functional) diversity data across island archipelagos

**Variables**: (Same structure as taxonomic, with functional diversity metrics)
- `beta_func_sor`: Functional Sørensen diversity index
- Trait-based metrics of coral communities

**Data Type**: Trait-based community composition
**Number of Observations**: Island-pair comparisons

---

#### `model_data_fish_taxonomic_arch.csv`
**Description**: Fish species composition (taxonomic diversity) data

**Variables**:
- Same structural variables as coral data
- `beta_sor`: Fish species turnover
- Fish-specific community metrics

**Data Type**: Community composition matrix
**Number of Observations**: Island-pair comparisons

---

#### `model_data_fish_functional_arch_VERSION2.csv`
**Description**: Fish trait-based (functional) diversity data (latest version)

**Note**: VERSION2 represents the latest processed dataset with corrected trait values

**Variables**: Fish functional diversity metrics
**Data Type**: Trait-based community composition

---

### Terrestrial Ecosystems

#### `bird_taxonomic_only_native.csv`
**Description**: Native bird species composition data across islands

**Variables**:
- `Island_ID`: Island identifier
- `Species_ID`: Bird species identifier
- `Presence`: Binary presence/absence (0/1)
- `Distance`: Inter-island distance
- `Island_Area`: Area of island (km²)
- `Elevation`: Maximum elevation (m)
- `beta_sor`: Bird species turnover
- Geographic and ecological covariates

**Data Type**: Native species occurrence data
**Units**: Species counts and area in km², elevation in meters
**Filtering**: Native species only (non-introduced)

---

#### `bird_functional_only_native.csv`
**Description**: Native bird trait-based (functional) diversity data

**Variables**:
- Community-weighted mean trait values
- Trait diversity indices
- `beta_func_sor`: Functional diversity turnover
- Same geographic/ecological covariates as taxonomic data

**Data Type**: Trait-based diversity metrics
**Traits Included**: Body size, diet, foraging behavior, habitat preference

---

#### `plant_taxonomic_only_native_final.csv`
**Description**: Native plant species composition data (final processed version)

**Variables**:
- `Plot_ID`: Vegetation plot identifier
- `Species_ID`: Plant species identifier
- `Abundance`: Species abundance/cover (%)
- `Distance`: Inter-island distance
- `Island_Area`: Island area (km²)
- `Precipitation`: Mean annual precipitation (mm)
- `Temperature`: Mean annual temperature (°C)
- `beta_sor`: Plant species turnover

**Data Type**: Plant community composition
**Units**: Cover percentage (%), climate in standard units
**Filtering**: Native species only

---

#### `plant_functional_only_native_final.csv`
**Description**: Native plant trait-based (functional) diversity data (final version)

**Variables**:
- Community-weighted mean trait values
- Functional diversity indices
- `beta_func_sor`: Functional turnover
- Environmental covariates

**Data Type**: Trait-based diversity metrics
**Traits Included**: Leaf area, wood density, height, seed size, growth rate

---

## Data Quality Notes

### Data Standardization
- All numeric variables are in consistent units
- Geographic distances: kilometers
- Areas: square kilometers
- Temperature: degrees Celsius
- Precipitation: millimeters

### Missing Data
- Indicated as `NA` in CSV files
- Handled within R scripts using `na.omit()` or multiple imputation
- Documented in analysis scripts

### Scaling and Transformation
- Distance variables may be log-transformed in models
- Beta-diversity indices are typically bounded [0, 1]
- Predictor variables are centered and scaled in statistical models

## Data Access and Privacy

- All data are de-identified and aggregated at island-level
- No individual organism-level data included
- Data collection permits archived (not included in repository)

## Citation

If using any datasets from this project, please cite the original source and this repository:

```
Waechter, L. et al. (2026). Islands Biogeography dataset: Beta-diversity 
patterns in island archipelagos. BioScales Lab. 
https://github.com/BioScalesLab/Islands_Biogeography
```

## Data Maintenance

**Last Updated**: January 30, 2026
**Data Version**: Final (v1.0)
**Curator**: Luiza Waechter

For data updates or corrections, please contact: luizawaechter.s@gmail.com

---

**Note**: Raw data files are managed separately. This repository contains processed, 
analysis-ready datasets. Large datasets (>50MB) are managed with Git LFS.
