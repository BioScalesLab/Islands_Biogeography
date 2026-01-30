# Quick Start Guide

## 5-Minute Setup

### Prerequisites
- R 4.0+
- RStudio (recommended)
- Git (for cloning)

### Installation

1. **Clone the repository**
```bash
git clone https://github.com/BioScalesLab/Islands_Biogeography.git
cd Islands_Biogeography
```

2. **Open in RStudio**
```bash
open Islands_Biogeography.Rproj
# or double-click the .Rproj file
```

3. **Install required packages**
```r
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse, glmmTMB, brms, ggplot2, 
  tidybayes, bayestestR, bayesplot, ggeffects
)
```

4. **Set working directory**
```r
# Automatically set by RStudio project, or:
setwd("Data/")
```

## Running Your First Analysis

### Option 1: Full Analysis (2-4 hours)
```r
source("RCode/Island_Betadiversity_models.R")
```

### Option 2: Quick Analysis (10 minutes)
```r
# Load just the data and run correlations
source("RCode/Spearman_correlations.R")
source("RCode/Density_Plots.R")
```

### Option 3: Specific Dataset
```r
# Load a single dataset and explore
library(tidyverse)
setwd("Data/")

# Example: Load coral data
corals <- read.csv("model_data_coral_taxonomic_arch.csv") %>%
  mutate_if(is.character, as.factor)

summary(corals)
head(corals)
```

## Project Structure

```
├── RCode/                      # Analysis scripts
│   └── Island_Betadiversity_models.R  # Main analysis
├── Data/                       # Data files (12 CSV files)
├── outputs/                    # Generated results
├── README.md                   # Full documentation
└── PROJECT_CONFIG.md           # Detailed configuration
```

## Key Files

| File | Purpose | Run Time |
|------|---------|----------|
| `Island_Betadiversity_models.R` | Statistical models (GLMM, Bayesian) | 2-4 hours |
| `Spearman_correlations.R` | Correlation analysis | 10 min |
| `Density_Plots.R` | Distribution visualization | 5 min |

## Common Tasks

### View Data Summary
```r
library(tidyverse)

# Load example dataset
corals <- read.csv("Data/model_data_coral_taxonomic_arch.csv")

# Quick summary
glimpse(corals)
summary(corals$beta_sor)
```

### Create a Simple Plot
```r
library(ggplot2)

corals %>%
  ggplot(aes(x = Distance, y = beta_sor)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal() +
  labs(
    x = "Geographic Distance (km)",
    y = "Beta Diversity (Sørensen)",
    title = "Coral Species Turnover vs. Distance"
  )
```

### Fit a Simple Model
```r
library(glmmTMB)

# Load data
corals <- read.csv("Data/model_data_coral_taxonomic_arch.csv")

# Scale predictors
corals_scaled <- corals %>%
  mutate(
    distance_sc = scale(Distance)[,1],
    isolation_sc = scale(diff_isolation)[,1]
  )

# Fit model
model <- glmmTMB(
  beta_sor ~ distance_sc + isolation_sc + (1|Archipelago_ID),
  family = beta_family(),
  data = corals_scaled
)

summary(model)
```

## Troubleshooting

### Issue: Package installation fails
```r
# Try installing dependencies first
install.packages(c("tidyverse", "glmmTMB"), 
                 repos = "https://cloud.r-project.org/")
```

### Issue: Working directory not found
```r
# Check current directory
getwd()

# Navigate to project folder explicitly
setwd("/path/to/Islands_Biogeography/Data/")

# Or use relative path if in RStudio project
setwd("Data/")
```

### Issue: Memory error on large models
```r
# Reduce computational demands
options(mc.cores = 2)  # Use fewer cores

# Run simpler models first
# Increase chains and iterations gradually
```

### Issue: Package conflicts or "masking" warnings
```r
# This is normal. Manage with explicit namespace:
ggplot2::ggplot(...)
dplyr::select(...)
```

## Getting Help

1. **Documentation**
   - [README.md](README.md) - Full project overview
   - [CONTRIBUTING.md](CONTRIBUTING.md) - Guidelines
   - [PROJECT_CONFIG.md](PROJECT_CONFIG.md) - Detailed info
   - [Data/README.md](Data/README.md) - Data dictionary

2. **RStudio Help**
   - `?function_name` - Function documentation
   - `help(package="pacman")` - Package help

3. **Contact**
   - Email: luizawaechter.s@gmail.com
   - GitHub Issues: [Report a bug](../../issues)

## Next Steps

- [ ] Clone and open the project
- [ ] Install required packages
- [ ] Explore data files
- [ ] Run a simple correlation analysis
- [ ] Read [PROJECT_CONFIG.md](PROJECT_CONFIG.md) for full details
- [ ] Review [CONTRIBUTING.md](CONTRIBUTING.md) if planning contributions

## Performance Tips

💡 **For faster results:**
- Start with `Spearman_correlations.R` (10 min)
- Use a subset of data for testing
- Reduce MCMC iterations in Bayesian models

⚡ **For parallel processing:**
```r
# Automatically uses all available cores:
options(mc.cores = parallel::detectCores())

# Or specify cores:
options(mc.cores = 4)
```

🖥️ **System requirements for full run:**
- CPU: 4+ cores recommended
- RAM: 8+ GB
- Storage: 5+ GB free

---

**Version**: 1.0  
**Last Updated**: January 30, 2026  
**Contact**: luizawaechter.s@gmail.com
