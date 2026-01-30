# Frequently Asked Questions (FAQ)

## General Questions

### What is this project about?

**Islands Biogeography** is a research project investigating biodiversity patterns in island archipelagos. We examine both marine ecosystems (corals and fish) and terrestrial ecosystems (plants and birds), using taxonomic (species-based) and functional (trait-based) approaches.

**Key Questions**:
- How do island distance and isolation affect biodiversity?
- Do taxonomic and functional diversity patterns differ?
- What environmental factors drive species turnover?

### Who is behind this project?

**Lead Researcher**: Luiza Waechter  
**Affiliation**: BioScales Lab  
**Collaborators**: [Field teams and research institutions]

### How can I cite this project?

```bibtex
@software{waechter2026islands,
  author = {Waechter, Luiza},
  title = {Islands Biogeography: Beta-diversity patterns in island archipelagos},
  year = {2026},
  url = {https://github.com/BioScalesLab/Islands_Biogeography}
}
```

Or in text format:
```
Waechter, L. (2026). Islands Biogeography: Beta-diversity patterns in island 
archipelagos. BioScales Lab. Retrieved from 
https://github.com/BioScalesLab/Islands_Biogeography
```

### Is this project open source?

Yes! The code is licensed under the MIT License, making it free to use, modify, and distribute for both commercial and non-commercial purposes.

---

## Getting Started

### How do I install the required packages?

The easiest way is to use the `pacman` package:

```r
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, glmmTMB, brms, ggplot2, tidybayes, bayestestR)
```

This will automatically install any missing packages.

### How long does the full analysis take to run?

- **Full Bayesian models**: 2-4 hours (depending on hardware)
- **Correlation analysis**: ~10 minutes
- **Visualization**: ~5 minutes

You can run just the correlations and plots for a quick overview.

### What are the system requirements?

**Minimum:**
- R 4.0+
- 4 GB RAM
- 5 GB disk space

**Recommended:**
- R 4.2+
- 8+ GB RAM
- Multi-core processor
- 10+ GB disk space

### Can I run this on Mac/Windows/Linux?

Yes! The project is fully cross-platform. It runs on:
- ✅ macOS
- ✅ Windows
- ✅ Linux (including Ubuntu, CentOS, Fedora)

### How do I update the packages?

```r
# Update all packages
update.packages()

# Or update specific packages
pacman::p_update()
```

---

## Data Questions

### Where can I find the data?

The processed, analysis-ready data is in the `Data/` directory:
- 12 CSV files with taxonomic and functional diversity data
- Full data dictionary in [Data/README.md](Data/README.md)

### Can I access the raw data?

Raw field data is archived separately for data security and privacy reasons. Contact luizawaechter.s@gmail.com for access requests.

### What are the data sources?

**Marine**: 
- Coral and fish survey data from reef monitoring programs
- Multiple island archipelagos globally

**Terrestrial**:
- Bird and plant surveys from ecological monitoring
- Native species only (non-introduced)

### Can I use this data in my research?

Yes, with proper citation. See the citation section above.

### Are there size limits for the data files?

Large files (>50 MB) are managed with Git LFS (Git Large File Storage). They download automatically but require LFS to be installed.

---

## Analysis Questions

### What do the models do?

The models examine relationships between:
- **Response**: Beta-diversity (species turnover between islands)
- **Predictors**: Geographic distance, isolation time, island size, climate
- **Random Effects**: Archipelago grouping

### What does beta-diversity mean?

**Beta-diversity** measures how different species communities are between two locations. High beta-diversity means communities are very different; low beta-diversity means they're similar.

We use the **Sørensen index**, which ranges from 0 (identical) to 1 (completely different).

### Why do you compare taxonomic and functional diversity?

Different metrics reveal different patterns:
- **Taxonomic diversity**: Species presence/absence
- **Functional diversity**: Ecological traits and roles

Functional diversity can be maintained even if species composition changes.

### What's the difference between GLMM and Bayesian models?

- **GLMM** (glmmTMB): Frequentist approach, faster computation
- **Bayesian** (brms): Incorporates prior information, better uncertainty estimates

We run both for comparison.

---

## Technical Questions

### How do I increase computation speed?

```r
# Use more parallel processing cores
options(mc.cores = parallel::detectCores())

# Reduce MCMC iterations (less precise but faster)
# In brms models, reduce iter and chains parameters
```

### What if I get memory errors?

```r
# Reduce parallel processing
options(mc.cores = 2)

# Run models one at a time instead of all at once
# Use a subset of data for testing
```

### How do I debug errors in the code?

1. Check the error message carefully
2. Verify your working directory: `getwd()`
3. Confirm all data files are present: `list.files("Data/")`
4. Check package versions: `packageVersion("glmmTMB")`
5. Open an issue on GitHub if you can't resolve it

### Can I modify the code?

Absolutely! See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### How do I contribute my changes?

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Submit a pull request

See [CONTRIBUTING.md](CONTRIBUTING.md) for detailed instructions.

---

## Troubleshooting

### I'm getting "package not found" errors

Install the missing package:
```r
install.packages("package_name", repos = "https://cloud.r-project.org/")
```

Or use pacman:
```r
pacman::p_load(package_name)
```

### The data won't load

Check your working directory:
```r
getwd()

# Set to Data folder
setwd("Data/")

# List available files
list.files()
```

### Models won't converge

This is common with complex data. Try:
1. Reduce model complexity (fewer predictors)
2. Increase iterations: `iter = 4000, warmup = 2000`
3. Use priors to guide the model
4. Check for multicollinearity in predictors

### I get different results each time

This is normal for Bayesian models (random sampling). To reproduce exactly:

```r
set.seed(12345)
# Run model
```

### The plots don't look right

Check that ggplot2 is loaded:
```r
library(ggplot2)
```

And verify the data format:
```r
str(your_data)
```

---

## Collaboration & Contributing

### Can I contribute to this project?

Yes! We welcome contributions from the community. See [CONTRIBUTING.md](CONTRIBUTING.md).

### What types of contributions are welcome?

- Bug reports and fixes
- Code improvements
- Documentation enhancements
- New analyses or features
- Data corrections
- Visualization improvements

### How do I report a bug?

Open a GitHub issue with:
1. Description of the bug
2. Steps to reproduce
3. Expected vs. actual behavior
4. Error message (if applicable)
5. Your environment (R version, OS, package versions)

### What if I find a security issue?

Please report it to luizawaechter.s@gmail.com rather than opening a public issue. See [SECURITY.md](SECURITY.md).

### Can I use this code in my publication?

Yes! Please cite the repository and our work. See citation guidelines above.

---

## Advanced Questions

### How do I add new data to the analysis?

1. Add your CSV file to `Data/`
2. Update [Data/README.md](Data/README.md) with file description
3. Load it in your analysis script
4. Submit a pull request if you want to merge changes

### How do I create custom models?

1. Load your data
2. Scale variables if needed
3. Use `glmmTMB()` or `brms::brm()` to fit
4. Use `summary()`, `plot()`, `predict()` for results

See [Island_Betadiversity_models.R](RCode/Island_Betadiversity_models.R) for examples.

### How do I modify the plots?

The plots use `ggplot2`. You can:
1. Change themes: `+ theme_minimal()`
2. Adjust colors: `+ scale_color_manual()`
3. Add facets: `+ facet_wrap(~group)`
4. Customize labels: `+ labs(title = "Custom Title")`

---

## Contact & Support

### How do I contact the project lead?

**Email**: luizawaechter.s@gmail.com  
**Response Time**: 24-48 hours

### Where can I ask questions?

1. **GitHub Discussions**: For project questions
2. **GitHub Issues**: For bugs and feature requests
3. **Email**: For urgent or sensitive matters
4. **Documentation**: Check [README.md](README.md) and [QUICKSTART.md](QUICKSTART.md)

### How do I stay updated?

- **Watch** the GitHub repository for updates
- **Star** to show support
- **Subscribe** to GitHub notifications
- **Check** [CHANGELOG.md](CHANGELOG.md) for version updates

---

## Useful Links

- [Quick Start Guide](QUICKSTART.md) - Get running in 5 minutes
- [Contribution Guidelines](CONTRIBUTING.md) - How to contribute
- [Data Dictionary](Data/README.md) - Understand your data
- [Project Configuration](PROJECT_CONFIG.md) - Detailed project info
- [Code of Conduct](CODE_OF_CONDUCT.md) - Community guidelines
- [Security Policy](SECURITY.md) - Report vulnerabilities

---

**Last Updated**: January 30, 2026  
**Contact**: luizawaechter.s@gmail.com  
**Repository**: [BioScalesLab/Islands_Biogeography](https://github.com/BioScalesLab/Islands_Biogeography)

**Can't find what you're looking for?** Open an issue or email us!
