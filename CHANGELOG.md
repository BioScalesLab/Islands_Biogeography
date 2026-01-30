# Changelog

All notable changes to the Islands Biogeography project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-01-30

### Added

#### Documentation
- Comprehensive `README.md` with project overview and quick start guide
- `QUICKSTART.md` for new users (5-minute setup)
- `CONTRIBUTING.md` with code style guidelines and contribution workflow
- `PROJECT_CONFIG.md` with detailed project configuration and research objectives
- `Data/README.md` with complete data dictionary for all 12 datasets
- `RCode/README.md` with documentation for each analysis script
- `CHANGELOG.md` (this file) for version tracking

#### Project Infrastructure
- MIT License for open-source distribution
- Comprehensive `.gitignore` with R, IDE, and system-specific patterns
- Git LFS configuration for managing large data files
- `.gitattributes` for tracking large RData files with Git LFS

#### Code Organization
- Restructured repository with professional directory layout
- Organized documentation in English throughout
- Established code style guidelines (snake_case, line length limits)
- Added file headers and section markers in R scripts

### Changed
- Updated `.gitignore` from minimal to comprehensive (100+ patterns)
- Reorganized project documentation structure
- Improved file naming and documentation consistency
- Enhanced data file descriptions and usage instructions

### Fixed
- Resolved GitHub repository file size limit issue by implementing Git LFS
- Cleaned up git history to remove oversized data files
- Fixed encoding issues in data dictionary

## [Unreleased]

### Planned Features
- [ ] Add automated testing framework
- [ ] Implement continuous integration (GitHub Actions)
- [ ] Create Docker container for reproducible environment
- [ ] Add model comparison benchmarking
- [ ] Develop Shiny app for interactive exploration
- [ ] Create manuscript generation pipeline

### Upcoming Changes
- [ ] Expand Bayesian model diagnostics
- [ ] Add cross-validation framework
- [ ] Implement sensitivity analysis pipeline
- [ ] Create publication-ready figure templates

---

## Version Information

### Current Release: v1.0.0
- **Release Date**: January 30, 2026
- **Status**: Stable
- **R Requirement**: R ≥ 4.0.0
- **Key Dependencies**:
  - glmmTMB 1.1+
  - brms 2.17+
  - tidyverse 1.3+
  - bayestestR 0.12+

### Previous Versions
- **v0.x.x**: Beta testing and development phases

---

## Guidelines for Reporting Changes

### Types of Changes
- **Added**: New features, functions, or documentation
- **Changed**: Changes in existing functionality
- **Deprecated**: Features or functions no longer recommended
- **Removed**: Removed features or functions
- **Fixed**: Bug fixes or error corrections
- **Security**: Security vulnerability fixes

### Version Numbers

We follow Semantic Versioning (MAJOR.MINOR.PATCH):
- **MAJOR**: Significant changes or breaking changes to API
- **MINOR**: New features that are backwards compatible
- **PATCH**: Bug fixes and small improvements

Example: `1.2.3`
- Major: 1 (major version)
- Minor: 2 (new features)
- Patch: 3 (bug fix)

---

## Migration Guides

### Upgrading to v1.0.0 from v0.x

1. **Update code structure**
   - Ensure scripts follow new style guidelines
   - Update relative paths if directory structure changed

2. **Package updates**
   - Update all packages to latest compatible versions
   - Run `pacman::p_load()` to ensure all packages installed

3. **Data format**
   - No breaking changes to data format
   - All previous data files compatible

---

## Roadmap

### Q1 2026
- [x] Professional repository structure
- [ ] Extended model diagnostics
- [ ] Publication of v1.0

### Q2 2026
- [ ] Automated testing suite
- [ ] Docker environment setup
- [ ] GitHub Actions CI/CD

### Q3-Q4 2026
- [ ] Interactive Shiny app
- [ ] Manuscript submission
- [ ] Community engagement

---

## Contributors

### v1.0.0 Contributors
- **Luiza Waechter** (@LuizaWaechter) - Lead developer
- **BioScales Lab** - Research group

### Future Contributors
We welcome contributions! See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

---

## Support

- **Issues & Bug Reports**: [GitHub Issues](../../issues)
- **Discussions**: [GitHub Discussions](../../discussions)
- **Email**: luizawaechter.s@gmail.com
- **Documentation**: See [README.md](README.md)

---

## License

All changes and versions are licensed under the [MIT License](LICENSE).

---

**Last Updated**: January 30, 2026  
**Maintained by**: Luiza Waechter  
**Repository**: [BioScalesLab/Islands_Biogeography](https://github.com/BioScalesLab/Islands_Biogeography)
