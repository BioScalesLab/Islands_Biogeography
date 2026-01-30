# Project Configuration

## Project Information

**Title**: Islands Biogeography: Beta-Diversity Patterns in Island Archipelagos

**Lead Researcher**: Luiza Waechter  
**Institution**: BioScales Lab  
**Start Date**: 2024  
**Current Status**: Active Development  
**Last Updated**: January 30, 2026

## Research Objectives

1. Quantify beta-diversity (species turnover) patterns across island systems
2. Compare taxonomic (species-based) and functional (trait-based) diversity metrics
3. Examine environmental and geographic drivers of biodiversity patterns
4. Develop predictive models for biodiversity in isolated island communities
5. Compare marine (coral, fish) vs. terrestrial (plant, bird) ecosystem patterns

## Study Scope

### Geographic Coverage
- Multiple island archipelagos globally
- Range of island sizes and isolation distances
- Varying environmental conditions

### Taxonomic Groups
- **Marine**: Corals (Scleractinia), Fish (multiple families)
- **Terrestrial**: Birds (native species), Plants (vascular plants)

### Diversity Metrics
- **Taxonomic**: Species composition, species richness, Sørensen index
- **Functional**: Trait-based diversity, functional richness, community-weighted means

## Data Summary

| Ecosystem | Taxonomic | Functional | Records |
|-----------|-----------|-----------|---------|
| Corals | ✓ | ✓ | 150-200 island pairs |
| Fish | ✓ | ✓ | 150-200 island pairs |
| Birds | ✓ | ✓ | 100-150 islands |
| Plants | ✓ | ✓ | 100-150 islands |

## Analysis Methods

### Statistical Models
- **Generalized Linear Mixed Models** (glmmTMB)
- **Bayesian Hierarchical Models** (brms)
- **Spearman Rank Correlations**
- **Distribution Analysis**

### Key Predictor Variables
- Geographic distance between islands
- Quaternary isolation (time since last connection)
- Island area and environmental factors
- Climate variables (temperature, precipitation)

### Response Variables
- Sørensen beta-diversity index
- Species turnover (Beta-SOR)
- Functional diversity metrics
- Community composition similarity

## Computing Environment

### Required Software
- **R**: Version 4.0 or higher
- **RStudio**: Latest version (recommended)
- **Operating System**: Linux, macOS, or Windows

### Recommended Hardware
- **CPU**: Multi-core processor (4+ cores)
- **RAM**: 8+ GB for full Bayesian models
- **Storage**: 10+ GB for raw data and model outputs

### Key R Packages
See `Island_Betadiversity_models.R` for complete package list

## Directory Organization

```
Islands_Biogeography/
├── RCode/                    # Analysis scripts
│   ├── Island_Betadiversity_models.R       # Main analysis
│   ├── Spearman_correlations.R             # Correlation analysis
│   ├── Density_Plots.R                     # Visualization
│   └── README.md                           # Script documentation
├── Data/                     # Data files
│   ├── Coral data (4 files)
│   ├── Fish data (2 files)
│   ├── Bird data (2 files)
│   ├── Plant data (2 files)
│   └── README.md                           # Data dictionary
├── outputs/                  # Analysis results (generated)
├── docs/                     # Additional documentation
├── README.md                 # Project overview
├── CONTRIBUTING.md           # Contribution guidelines
├── LICENSE                   # MIT License
├── .gitignore               # Git configuration
└── Islands_Biogeography.Rproj # RStudio project
```

## Project Timeline

- **Phase 1 (2024)**: Data collection and compilation
- **Phase 2 (2024-2025)**: Model development and testing
- **Phase 3 (2025-2026)**: Final analyses and manuscript preparation
- **Phase 4 (2026)**: Publication and public release

## Collaboration Guidelines

### Code Review
- All changes reviewed before merging
- Clear commit messages required
- Documentation updates mandatory

### Contribution Process
See [CONTRIBUTING.md](../CONTRIBUTING.md)

### Communication
- Email: luizawaechter.s@gmail.com
- GitHub Issues: For bugs and feature requests
- Discussions: For general questions

## Data Management

### Data Storage
- Original raw data: Archived separately
- Processed data: Version controlled in `Data/`
- Large files: Managed with Git LFS (>50 MB)

### Reproducibility
- All analyses documented in R scripts
- Seed values set for random number generation
- Model summaries and diagnostics exported
- Session info recorded for computational environment

## Quality Assurance

### Code Quality
- ✓ Consistent style and formatting
- ✓ Comprehensive comments
- ✓ Error handling
- ✓ Reproducible results

### Data Quality
- ✓ Validation checks
- ✓ Missing data documentation
- ✓ Outlier assessment
- ✓ Units standardization

### Model Validation
- ✓ Convergence diagnostics
- ✓ Residual analysis
- ✓ Cross-validation
- ✓ Sensitivity analysis

## Publication Strategy

### Planned Outputs
1. **Primary manuscript**: Beta-diversity comparative analysis
2. **Data paper**: Dataset description and methods
3. **Supplementary materials**: Extended analyses
4. **Code repository**: Open-source analysis pipeline

### Expected Impacts
- Advance understanding of island biogeography
- Demonstrate taxonomic-functional diversity relationships
- Provide tools for conservation planning
- Enable hypothesis testing in island systems

## Funding and Acknowledgments

**Research Support**: BioScales Lab  
**Data Contributors**: [Field teams and institutions]  
**Software**: R and community packages

## Related Publications

List of related papers (when available):
- [Publication 1]
- [Publication 2]

## License

This project is licensed under the MIT License. See [LICENSE](../LICENSE) for details.

## Contact

**Lead Researcher**: Luiza Waechter  
Email: luizawaechter.s@gmail.com  
GitHub: [@BioScalesLab](https://github.com/BioScalesLab)

---

**Last Updated**: January 30, 2026  
**Version**: 1.0
