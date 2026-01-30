# Security Policy

## Reporting Security Vulnerabilities

If you discover a security vulnerability in the Islands Biogeography project, please report it responsibly to:

**Email**: luizawaechter.s@gmail.com  
**Subject**: Security Vulnerability Report - Islands Biogeography

### Please Include:
1. Description of the vulnerability
2. Steps to reproduce (if applicable)
3. Potential impact
4. Suggested fix (if you have one)

### Response Timeline
- **Initial Response**: Within 48 hours
- **Investigation**: Within 1 week
- **Resolution**: Varies based on severity

We appreciate responsible disclosure and will work with you to address any issues.

## Supported Versions

### Current Release
| Version | Supported       | Release Date   | End of Support |
|---------|-----------------|----------------|----------------|
| 1.x     | ✅ Yes          | Jan 30, 2026   | Jan 30, 2027   |

### Previous Versions
Older versions may have security issues and are not actively supported. We recommend updating to the latest version.

## Security Best Practices

### For Users
1. Keep R and all packages updated
2. Verify data sources before use
3. Use strong authentication for GitHub
4. Enable two-factor authentication on GitHub account
5. Review code changes before running

### For Developers
1. Never commit sensitive information (API keys, passwords)
2. Use `.gitignore` to prevent accidental commits
3. Review code before pushing
4. Keep dependencies updated
5. Use environment variables for configuration

### For Data Management
1. Keep raw data separate from version control
2. Use Git LFS for large files
3. Document data sources and privacy considerations
4. Verify data integrity
5. Use HTTPS for all downloads

## Dependencies and Package Security

### Regular Updates
We monitor package updates for security issues and update dependencies regularly.

### Known Issues
If a critical security issue is discovered in dependencies:
1. We will investigate immediately
2. Release a patched version
3. Communicate via GitHub security advisories

### Package Management
```r
# Keep packages updated
update.packages()

# Or use pacman for controlled updates
pacman::p_update()
```

## Data Privacy

### Handling User Data
- Only use public, de-identified data
- No personal information stored in repository
- Data collection follows ethical guidelines
- Respect data contributor privacy

### Third-Party Services
- GitHub for code hosting
- GitHub Releases for distribution
- No external analytics or tracking in this repository

## Compliance

This project aims to comply with:
- ✅ Open science best practices
- ✅ FAIR data principles
- ✅ Research ethics guidelines
- ✅ GitHub security standards

## Vulnerability Disclosure

When a vulnerability is confirmed:

1. **Assessment**: Severity rating (Critical/High/Medium/Low)
2. **Notification**: Security advisories issued
3. **Resolution**: Fix developed and tested
4. **Release**: Patched version released
5. **Communication**: Users notified of update

## Security Contact

**Lead Security Researcher**: Luiza Waechter  
**Email**: luizawaechter.s@gmail.com  
**Response Window**: 24-48 hours

## Additional Resources

- [GitHub Security Advisories](https://github.com/BioScalesLab/Islands_Biogeography/security/advisories)
- [Dependabot Alerts](https://github.com/BioScalesLab/Islands_Biogeography/security/dependabot)
- [GitHub Security Best Practices](https://docs.github.com/en/code-security)

---

**Last Updated**: January 30, 2026  
**Policy Version**: 1.0
