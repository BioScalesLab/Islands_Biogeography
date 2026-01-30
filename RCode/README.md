# RCode Directory

## Overview
This directory contains all R analysis scripts for the Islands Biogeography project.

## Scripts

### [Island_Betadiversity_models.R](Island_Betadiversity_models.R)
**Purpose**: Main statistical modeling pipeline for analyzing beta-diversity patterns

**Contents**:
1. **Data Loading**: Import taxonomy and functional diversity data for all taxa
2. **Data Preparation**: Scaling, transformations, and variable engineering
3. **Marine Models**:
   - Coral taxonomic diversity models
   - Coral functional diversity models
   - Fish taxonomic diversity models
   - Fish functional diversity models
4. **Terrestrial Models**:
   - Bird taxonomic diversity models
   - Bird functional diversity models
   - Plant taxonomic diversity models
   - Plant functional diversity models
5. **Model Types**:
   - GLMM with `glmmTMB`
   - Bayesian models with `brms`
   - Model diagnostics and validation

**Key Analysis Variables**:
- Response: Beta-diversity (Sørensen index)
- Predictors: Geographic distance, temporal isolation, environmental factors
- Random Effects: Archipelago groups

**Output**: Model summaries, predictions, effect plots

**Dependencies**: 
- `glmmTMB`, `brms`, `tidyverse`, `bayesplot`, `tidybayes`

---

### [Spearman_correlations.R](Spearman_correlations.R)
**Purpose**: Correlation analysis between beta-diversity and environmental variables

**Contents**:
1. **Data Preparation**: Format data for correlation analysis
2. **Spearman Rank Correlations**: Non-parametric correlation tests
3. **Correlation Matrices**: Create correlation heatmaps
4. **Significance Testing**: Calculate p-values and confidence intervals
5. **Visualization**: Correlation plots and network diagrams

**Analysis Variables**:
- Correlations between taxonomic and functional diversity
- Environmental predictors (distance, isolation, area, climate)
- Cross-ecosystem comparisons

**Output**: Correlation tables, heatmaps, p-value matrices

**Dependencies**: 
- `corrplot`, `Hmisc`, `tidyverse`, `ggplot2`

---

### [Density_Plots.R](Density_Plots.R)
**Purpose**: Distribution analysis and visualization of diversity metrics

**Contents**:
1. **Distribution Analysis**: Examine distribution shapes of beta-diversity
2. **Density Plots**: Create kernel density estimates
3. **Comparative Plots**: Compare distributions across taxa/ecosystems
4. **Statistical Tests**: Normality tests, distribution comparisons
5. **Visualization**: Publication-quality density plots with ggplot2

**Analysis Variables**:
- Beta-diversity distributions
- Environmental variable distributions
- By ecosystem type and diversity metric

**Output**: Density plots, distribution comparison figures

**Dependencies**: 
- `ggplot2`, `ggridges`, `tidyverse`, `patchwork`

---

## Running the Scripts

### Sequential Execution (Recommended)
```r
# Load main analysis (takes 2-4 hours for full Bayesian models)
source("RCode/Island_Betadiversity_models.R")

# Run supplementary analyses
source("RCode/Spearman_correlations.R")
source("RCode/Density_Plots.R")
```

### Individual Script Execution
Each script can be run independently, but ensure data is loaded first:
```r
# Load data manually if running scripts individually
setwd("Data/")
coral_tax <- read.csv("model_data_coral_taxonomic_arch.csv")
# ... load other data files
```

## Script Structure

Each analysis script follows this standard structure:

```
1. File Header
   - Script name
   - Author and date
   - Purpose
   - Input/output description

2. Library Loading
   - Required packages
   - Custom function definitions

3. Configuration
   - Set working directory
   - Options (parallel processing, etc.)

4. Data Loading & Preparation
   - Read CSV files
   - Data cleaning and transformation
   - Variable scaling/coding

5. Analysis Sections
   - Clearly marked with comments
   - Logical grouping by analysis type

6. Visualization
   - Plot creation
   - Export to files

7. Output & Reporting
   - Summary statistics
   - Model diagnostics
   - Publication-ready tables
```

## Code Style Guidelines

- **Naming**: Use snake_case for variables and functions
- **Comments**: Use `#` for inline comments, `# ----` for section markers
- **Line Length**: Maximum 100 characters
- **Spacing**: 2-space indentation

## Performance Notes

- **Bayesian Models**: Computationally intensive (2-4 hours for full run)
  - Uses parallel processing across all cores
  - Save model objects after fitting
  
- **Data Processing**: Vectorized operations where possible
- **Memory**: Large datasets may require adequate RAM (8GB+ recommended)

## Troubleshooting

### Common Issues

**Package Installation Errors**
```r
# Use pacman for automatic installation
if (!require("pacman")) install.packages("pacman")
pacman::p_load(package_name)
```

**Memory Issues**
```r
# Reduce parallel chains in brms models
options(mc.cores = 4)  # Instead of all cores
```

**Data Not Found**
```r
# Set working directory explicitly
setwd("/path/to/Data/")
# Or use relative paths from RStudio project
```

## Contributing

When adding new scripts:
1. Follow naming convention: `Descriptive_Analysis_Name.R`
2. Include comprehensive file header
3. Document all sections with comments
4. Add reference to this README.md
5. Ensure reproducibility

## Version History

- **v1.0** (Jan 30, 2026): Initial analysis pipeline
  - Core models implemented
  - Correlation and visualization analyses

---

**Last Updated**: January 30, 2026
**Maintainer**: Luiza Waechter
