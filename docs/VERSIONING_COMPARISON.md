# Detailed Versioning Comparison: setuptools-scm vs Static

## setuptools-scm Disadvantages

### 1. **Git Dependency**
- **Problem**: Requires `.git` directory to determine version
- **Impact**:
  - Breaks if someone downloads source tarball (not git clone)
  - Requires `.git_archival.txt` workaround for non-git distributions
  - Users who download "Download ZIP" from GitHub get version detection failures

**Example failure:**
```bash
# User downloads ZIP from GitHub
wget https://github.com/user/repo/archive/refs/heads/main.zip
unzip main.zip
cd repo-main
pip install .
# ERROR: No git repository found, can't determine version
```

### 2. **Build-Time Dependency**
- **Problem**: setuptools-scm must be installed *before* building
- **Impact**:
  - Adds to `build-system.requires` in pyproject.toml
  - Increases build time (need to download/install setuptools-scm first)
  - More dependencies to maintain and potentially fail

**Before pip can even start installing your package:**
```toml
[build-system]
requires = ["hatchling", "hatch-vcs"]  # Must download these first
build-backend = "hatchling.build"
```

### 3. **Shallow Clone Issues**
- **Problem**: Many CI/CD systems use shallow clones by default
- **Impact**:
  - GitHub Actions: `fetch-depth: 0` required to get tag history
  - GitLab CI: `GIT_DEPTH: 0` needed
  - Without tags, version detection fails or gives wrong version

**CI Configuration required:**
```yaml
# .github/workflows/test.yml
- uses: actions/checkout@v4
  with:
    fetch-depth: 0  # MUST include this or setuptools-scm fails
```

### 4. **Development Version Confusion**
- **Problem**: Between tags, you get versions like `2.0.0.dev42+g1234abc.d20241111`
- **Impact**:
  - Confusing for users/developers
  - Long version strings in logs
  - Hard to correlate to actual code state
  - `pip list` shows cryptic versions

**What users see:**
```bash
$ pip list
intronIC  2.0.0.dev42+g9412cb2.d20241111  # What does this mean?

$ intronIC --version
intronIC 2.0.0.dev42+g9412cb2.d20241111  # Confusing!
```

### 5. **Dirty Working Tree Versions**
- **Problem**: Uncommitted changes add `.dYYYYMMDD` suffix
- **Impact**:
  - Installing from development directory gives version like `2.0.0.dev5+gabc.d20241111`
  - Every install might have different version string
  - Hard to track what version someone is actually using

### 6. **Not Obvious in Source Code**
- **Problem**: Version isn't visible in any source file
- **Impact**:
  - Can't grep for version in code
  - Must understand git tagging to know version
  - New contributors confused about version management

**Comparison:**
```python
# Static: Clear in pyproject.toml
[project]
version = "2.0.0"  # ← Right there!

# setuptools-scm: Where's the version?
[project]
dynamic = ["version"]  # ← Not helpful
```

### 7. **Tag Naming Requirements**
- **Problem**: Must follow specific tag format conventions
- **Impact**:
  - Tags must be `v1.2.3` or similar (configurable but requires setup)
  - Inconsistent tagging breaks version detection
  - Team must understand and follow tagging rules

### 8. **Source Distribution Complexity**
- **Problem**: Source distributions (sdist) need special handling
- **Impact**:
  - Must include `.git_archival.txt` file
  - Must configure `hatch-vcs` or similar to embed version in sdist
  - Extra files and configuration

**Required files:**
```
.git_archival.txt  # For non-git installs
pyproject.toml     # With special configuration
_version.py        # Generated at build time
```

### 9. **Installation from PyPI Edge Case**
- **Problem**: If version isn't properly embedded in wheel/sdist
- **Impact**:
  - Users installing from PyPI might get version detection failures
  - Requires careful testing of distribution process

### 10. **Learning Curve**
- **Problem**: Contributors need to understand version management
- **Impact**:
  - Must know: "To release, tag with v2.0.0"
  - Must understand git tagging
  - Must know how to fix version detection issues
  - Not obvious for Python-only developers

---

## Static Versioning Disadvantages

For balance, here are static versioning disadvantages:

### 1. **Manual Updates Required**
- **Problem**: Must remember to update version in files
- **Impact**: Can forget to bump version before release
- **Mitigation**: Scripts, checklists, pre-commit hooks

### 2. **Multiple Files Risk**
- **Problem**: Version might be in multiple places
- **Impact**: Can get out of sync
- **Mitigation**: Single source in pyproject.toml, read elsewhere via importlib.metadata

### 3. **No Automatic Dev Versions**
- **Problem**: Between releases, all installs show same version
- **Impact**: Can't distinguish different commits
- **Mitigation**: Use `-dev` suffix, or don't worry about it

---

## Head-to-Head Comparison

| Aspect | setuptools-scm | Static |
|--------|----------------|--------|
| **Visibility** | ❌ Hidden in git tags | ✅ Clear in pyproject.toml |
| **Simplicity** | ❌ Complex setup | ✅ Simple |
| **Dependencies** | ❌ Requires hatch-vcs/setuptools-scm | ✅ None (built-in importlib) |
| **CI/CD** | ❌ Needs fetch-depth config | ✅ Works everywhere |
| **Download ZIP** | ❌ Fails without git_archival | ✅ Always works |
| **Dev versions** | ✅ Automatic | ❌ Manual |
| **Git tags** | ✅ Automatic from tags | ❌ Manual tagging + version bump |
| **Learning curve** | ❌ Higher | ✅ Lower |
| **Error-prone** | ❌ More edge cases | ❌ Forgot to bump |
| **PyPI-ready** | ✅ Yes (if configured correctly) | ✅ Yes (always) |

---

## Real-World Example: setuptools-scm Failure

Here's what happens when things go wrong:

```bash
# Scenario: CI/CD with shallow clone
git clone --depth=1 https://github.com/user/intronic.git
cd intronic
pip install .

# Error:
LookupError: setuptools-scm was unable to detect version for intronic.

The following have been tried:
• git
  Traceback (most recent call last):
    ...
  LookupError: git repository has no tags

Make sure you're either building from a fully intact git repository
or PyPI tarballs. Most other sources (such as GitHub's tarballs, a
git checkout without the .git directory) do not contain the necessary
metadata and will not work.
```

**User reaction:** 😡 "Why can't I just install this?"

---

## When setuptools-scm Makes Sense

Despite disadvantages, setuptools-scm IS appropriate for:

1. **Large teams** where manual version bumps are often forgotten
2. **Frequent releases** where automation saves significant time
3. **Projects with many contributors** where consistency is hard to enforce
4. **When you want build metadata** in version strings for debugging
5. **Projects where dev versions matter** (e.g., pre-release testing)

---

## Recommendation for intronIC

### Use **Static Versioning** Because:

1. **Small team** (mostly you) - manual versioning is fine
2. **Infrequent releases** - not making releases daily
3. **Scientific software** - users prefer stable, clear versions
4. **Easy for contributors** - obvious how versioning works
5. **No surprises** - always works, everywhere
6. **Simple is better** - no build-time dependencies or edge cases

### Can Migrate Later If:
- Team grows significantly
- Release frequency increases
- Automated versioning becomes valuable
- You get tired of manual version bumps

**Migration is easy:**
```bash
# Switch takes ~10 minutes
# 1. Add hatch-vcs to pyproject.toml
# 2. Change version to dynamic
# 3. Add .git_archival.txt
# 4. Tag current state
# Done!
```

---

## Hybrid Approach (Best of Both?)

**Compromise:** Static version + git tag validation

```python
# In CI, check version matches tag
import tomli
from pathlib import Path

with open("pyproject.toml", "rb") as f:
    version = tomli.load(f)["project"]["version"]

git_tag = os.environ["GITHUB_REF"].split("/")[-1]  # refs/tags/v2.0.0

if git_tag != f"v{version}":
    raise ValueError(f"Tag {git_tag} doesn't match version {version}")
```

**Benefits:**
- Version in source (static)
- Enforced consistency (automated check)
- No setuptools-scm dependency
- Clear error if you forget to update

---

## Summary

**setuptools-scm is not bad**, but it's **overkill for most projects**.

**Main disadvantages:**
1. Git dependency (breaks for source downloads)
2. CI/CD complexity (shallow clones)
3. Build-time dependencies
4. Version string confusion
5. Learning curve

**For intronIC:** Static versioning is **simpler, more reliable, and sufficient**.

You can always switch to setuptools-scm later if the project grows to the point where automation is valuable.
