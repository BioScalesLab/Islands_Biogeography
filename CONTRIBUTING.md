# Contributing to Islands Biogeography

Thank you for your interest in contributing! Please follow these guidelines to maintain code quality and project organization.

## Code of Conduct

- Be respectful and professional
- Maintain scientific rigor and reproducibility
- Provide clear documentation for all changes

## How to Contribute

### 1. Fork and Clone
```bash
git clone https://github.com/BioScalesLab/Islands_Biogeography.git
cd Islands_Biogeography
```

### 2. Create a Feature Branch
```bash
git checkout -b feature/description-of-feature
# or
git checkout -b fix/description-of-bug
```

### 3. Make Your Changes

#### Code Style Guidelines

**R Code:**
- Use descriptive variable names (snake_case)
- Add comments for complex logic blocks
- Maximum line length: 100 characters
- Use `tidyverse` style conventions

**Example:**
```r
# Good
model_data_scaled <- raw_data %>%
  mutate(
    distance_scaled = scale(distance)[,1],
    area_scaled = scale(area)[,1]
  )

# Avoid
data<-raw_data; data$dist_sc=scale(data$Distance)
```

#### Documentation

- Add file headers explaining script purpose
- Document functions with parameter descriptions
- Update relevant README sections
- Include inline comments for non-obvious logic

**Script Header Template:**
```r
# ============================================================================
# Script: Analysis Description
# Author: Your Name
# Date: YYYY-MM-DD
# Purpose: Brief description of what this script does
# 
# Input: Data files used
# Output: Results/figures produced
# ============================================================================
```

### 4. Data Management

- **Large Files**: Use Git LFS for files > 50 MB
  ```bash
  git lfs track "*.RData"
  ```
- **Raw Data**: Never commit raw data without documentation
- **Derived Data**: Include data dictionaries in CSV headers
- Add `.gitignore` entries for local files (temporary outputs, cache, etc.)

### 5. Testing

Before committing:
- Run scripts and verify outputs
- Check for errors and warnings
- Test with clean workspace (`rm(list=ls())`)
- Ensure reproducibility

### 6. Commit Messages

Write clear, descriptive commit messages:

```
# Good
Add temperature variables to beta-diversity models
Implement Bayesian hierarchical model for coral taxonomy
Fix data scaling issue in fish functional analysis

# Avoid
Update script
Fix bug
Changes
```

Format:
- First line: 50 characters or less, concise summary
- Blank line
- Detailed description (if needed)

### 7. Push and Create Pull Request

```bash
git push origin feature/your-feature-name
```

Then create a Pull Request on GitHub with:
- Clear title and description
- Reference any related issues
- Explain your changes and their impact

## Project Organization Guidelines

### File Naming
- Scripts: `Descriptive_Analysis_Name.R`
- Data: `data_source_grouping_version.csv`
- Outputs: `figure_or_table_descriptive_name_date.png/pdf`

### Folder Structure
```
RCode/              # All analysis scripts
Data/               # All data files
docs/               # Additional documentation
outputs/            # Results, figures, tables
```

### Comments and Structure

```r
# ============================================================================
# SECTION 1: DATA LOADING AND PREPARATION
# ============================================================================

# Load libraries
library(tidyverse)

# Load data
data <- read.csv("data_file.csv")

# ============================================================================
# SECTION 2: STATISTICAL ANALYSIS
# ============================================================================

# Fit model
model <- glmmTMB(formula = response ~ predictor + (1|group), data = data)
```

## Issues and Bug Reports

### Reporting Issues

Include:
1. **Description**: What doesn't work?
2. **Steps to Reproduce**: How to replicate the issue
3. **Expected vs. Actual**: What should happen vs. what happens
4. **Environment**: R version, OS, package versions
5. **Error Message**: Full error output

**Example:**
```
Title: Coral model throws error with temperature variable

Description:
When running Island_Betadiversity_models.R with the temperature variable 
included, the script crashes at line 245.

Steps to reproduce:
1. Load model_data_coral_taxonomic_arch.csv
2. Run the "Corals Beta Taxonomic" section
3. Fit model with temperature predictor

Error:
Error in glmmTMB: dimension mismatch between fixed and random effects

Environment:
- R 4.3.0
- glmmTMB 1.1.8
- Ubuntu 20.04
```

## Documentation

- Update README.md for major changes
- Add docstrings for new functions
- Create/update analysis documentation in `docs/`
- Comment tricky or non-obvious code

## Review Process

1. **Automated Checks**: Code style and basic tests
2. **Peer Review**: At least one other team member reviews code
3. **Approval**: Project lead approves and merges

## Questions?

- Open an issue for discussion
- Contact: luizawaechter.s@gmail.com

---

**Last Updated**: January 30, 2026
