# Islands Biogeography

A comprehensive research project investigating beta-diversity patterns across island archipelagos using taxonomic and functional approaches for marine and terrestrial ecosystems.

## Overview

This repository contains statistical models and data analysis for studying biodiversity patterns in island systems. The project examines:

- **Marine ecosystems**: Coral and fish communities
- **Terrestrial ecosystems**: Plant and bird communities

Each ecosystem is analyzed using both **taxonomic** (species composition) and **functional** (trait-based) diversity metrics.

## Project Structure

```
Islands_Biogeography/
├── Data/                          # Raw and processed datasets
│   ├── Coral data (taxonomic & functional)
│   ├── Fish data (taxonomic & functional)
│   ├── Plant data (taxonomic & functional)
│   └── Bird data (taxonomic & functional)
├── RCode/                         # R analysis scripts
│   ├── Island_Betadiversity_models.R    # Main statistical models
│   ├── Spearman_correlations.R          # Correlation analyses
│   └── Density_Plots.R                  # Visualization scripts
├── docs/                          # Documentation
├── README.md                      # This file
├── CONTRIBUTING.md                # Contribution guidelines
└── Islands_Biogeography.Rproj     # RStudio project file
```

## Data Files

### Marine Ecosystems
- `model_data_coral_taxonomic_arch.csv` - Coral taxonomic diversity data
- `model_data_coral_functional_arch.csv` - Coral functional diversity data
- `model_data_fish_taxonomic_arch.csv` - Fish taxonomic diversity data
- `model_data_fish_functional_arch_VERSION2.csv` - Fish functional diversity data

### Terrestrial Ecosystems
- `bird_taxonomic_only_native.csv` - Bird taxonomic diversity data
- `bird_functional_only_native.csv` - Bird functional diversity data
- `plant_taxonomic_only_native_final.csv` - Plant taxonomic diversity data
- `plant_functional_only_native_final.csv` - Plant functional diversity data

## Requirements

### R Packages

Required packages (automatically loaded via `pacman`):
- `tidyverse` - Data manipulation and visualization
- `glmmTMB` - Generalized Linear Mixed Models
- `ggplot2` - Advanced plotting
- `brms` - Bayesian regression modeling
- `tidybayes` - Bayesian analysis tools
- `bayestestR` - Bayesian model testing
- `caret` - Machine learning preprocessing
- `xgboost` - Gradient boosting models
- `marginaleffects` - Model effect estimation
- And many more (see script headers for complete list)

### System Requirements
- R 4.0 or higher
- RStudio recommended

## Quick Start

1. **Clone the repository**
   ```bash
   git clone https://github.com/BioScalesLab/Islands_Biogeography.git
   cd Islands_Biogeography
   ```

2. **Open in RStudio**
   ```bash
   open Islands_Biogeography.Rproj
   ```

3. **Run the analysis scripts**
   ```r
   source("RCode/Island_Betadiversity_models.R")
   ```

## Analysis Pipeline

### 1. Data Loading and Preparation
- Import biodiversity and environmental data
- Format data for modeling (factors, scaling)
- Data quality checks

### 2. Statistical Modeling
- Generalized Linear Mixed Models (GLMM) with `glmmTMB`
- Bayesian regression models with `brms`
- Model validation and diagnostics

### 3. Correlation Analysis
- Spearman rank correlations between variables
- Multivariate analysis of diversity patterns

### 4. Visualization
- Density plots and distribution analysis
- Effect size plots
- Comparative visualizations across taxonomic groups

## Key Analyses

### Beta Diversity Models
The main analysis examines beta-diversity (species turnover) in relation to:
- Geographic distance between islands
- Temporal isolation (Quaternary period)
- Environmental factors (reef area, habitat heterogeneity)

Models are fitted separately for:
- Each ecosystem type (marine vs. terrestrial)
- Each diversity metric (taxonomic vs. functional)

## Configuration

- **Parallel Processing**: Enabled using all available cores
  ```r
  options(mc.cores = parallel::detectCores())
  ```

## Contributing

Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on how to contribute to this project.

## Authors

- **Luiza Waechter** - Lead researcher
- BioScales Lab

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this code or data in your research, please cite:

```
Waechter, L. et al. (2026). Islands Biogeography: Beta-diversity patterns 
in island archipelagos. [Unpublished manuscript]. BioScales Lab.
```

## Acknowledgments

- BioScales Lab research group
- Data contributors and field teams
- R development community

## Contact

For questions or inquiries, please contact:
- Luiza Waechter: luizawaechter.s@gmail.com

## Status

🔄 **Active Development** - Ongoing analysis and model refinement

---

**Last Updated**: January 30, 2026
