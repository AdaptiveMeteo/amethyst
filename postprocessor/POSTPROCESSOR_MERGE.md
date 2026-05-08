# Postprocessor Merge Documentation

**Status**: In progress | **Branch**: devel | **Current Head**: 9929c84 | **Last Updated**: October 28, 2024

## Overview

The postprocessor component underwent a major merge integrating two previously separate implementations. The merge involved:

1. **Consolidating validation logic** from two independent postprocessing pipelines
2. **Creating a unified atmospheric utilities library** (`utilities/atmos/atmos_tools.py`)
3. **Refactoring file naming conventions** (removing "mirto" prefixes)
4. **Standardizing dependencies** across validation workflows

## Merge Timeline

| Commit | Message | Key Changes |
|--------|---------|-------------|
| f38e8b5 | "Merged two versions of the postprocessing folder. Still to be tested" | Initial merge combining two implementations |
| 9929c84 | "FIRST EDIT of the merged postprocessor: removed mirto from filenames and changed calls to atmos library" | Refactoring & atmos library integration |

## Scope of Changes

### Files Deleted
- `postprocessor/validation/mirto_plot_validation_profiles.py` (623 lines, consolidated into other files)

### Files Renamed (mirto → generic naming)
- `mirto_validation_skewT_map.py` → `plot_validation_skewT_map.py`
- `mirto_wrf_validation_skewT_map.py` → `plot_wrf_validation_skewT_map.py`

### Files Modified
| File | Changes |
|------|---------|
| `postprocessor/validation/mirto_plot_validation_profiles_single.py` | Updated imports to use atmos_tools library |
| `postprocessor/validation/plot_validation_profiles.py` | Refactored ~161 lines; now uses unified atmos utilities |
| `postprocessor/validation/validate_retrievals.py` | 10 lines changed; imports from `utilities.atmos.atmos_tools` |
| `postprocessor/validation/validate_wrf.py` | 10 lines changed; imports from `utilities.atmos.atmos_tools` |
| `postprocessor/validation/validate_wrf_single.py` | 12 lines changed; imports from `utilities.atmos.atmos_tools` |

### Files Added
- **`utilities/atmos/atmos_tools.py`** (352 lines) — New unified atmospheric utilities library

## New Atmospheric Utilities Library

**Location**: [utilities/atmos/atmos_tools.py](../../utilities/atmos/atmos_tools.py)

This library consolidates atmospheric calculations previously scattered across validation scripts:

### Core Functions

| Function | Purpose | Input | Output |
|----------|---------|-------|--------|
| `check_units()` | Normalize temperature/pressure units | TEMP, PRES, optional units | Standardized arrays (K, hPa) |
| `compute_temp()` | Calculate temperature from potential temperature | THETA, PRES, optional units | Temperature (K) |
| `compute_rh()` | Calculate relative humidity from vapor & saturation | TEMP, PRES, QVAPOR | RH (%) |
| `compute_dewp()` | Calculate dew point from vapor & temperature | TEMP, PRES, QVAPOR | Dew point (K) |
| `interpolate_pressure_grid()` | Interpolate profiles to standard pressure levels | Profiles, current pressures, target pressures | Interpolated profiles |

**Units Supported**:
- **Temperature**: Kelvin (K), Celsius (C)
- **Pressure**: hPa, millibars (mb), Pascals (Pa)

### Example Usage

```python
from utilities.atmos.atmos_tools import compute_temp, compute_rh, check_units

# Normalize units
temp_k, pres_hpa = check_units(temp_c, pres_pa, units=['Pa', 'C'])

# Compute derived variables
temperature = compute_temp(theta, pres_hpa)
relative_humidity = compute_rh(temperature, pres_hpa, qvapor)
```

## Validation Workflow Changes

### Before Merge
- Multiple validation scripts with duplicated atmospheric calculations
- Mixed naming conventions (mirto_* prefix for MIRTO-specific code)
- Scattered pressure interpolation logic
- Inconsistent unit handling across scripts

### After Merge
- Single source of truth for atmospheric calculations (`atmos_tools.py`)
- Standardized naming conventions (generic `plot_*` and `validate_*` functions)
- Centralized pressure grid interpolation
- Consistent unit checking across all validation workflows

### Updated Import Pattern

**Old**:
```python
# Scattered functions across mirto_* modules
from mirto_plot_tools import compute_rh, compute_dewp
```

**New**:
```python
# Unified atmos library
from utilities.atmos.atmos_tools import compute_rh, compute_dewp, check_units
```

## Known Issues & Testing Status

⚠️ **Status**: "Still to be tested" (per merge commit message)

### Outstanding Validation Tasks

1. **End-to-end validation workflow**
   - [ ] Run `validate_retrievals.py` against test dataset
   - [ ] Run `validate_wrf.py` against WRF output
   - [ ] Run `validate_wrf_single.py` for single-profile validation
   - [ ] Verify output plots generate correctly

2. **Atmospheric calculations**
   - [ ] Compare RH calculations against canonical reference (e.g., NCAR formulae)
   - [ ] Test dew point calculations across temperature extremes (-50°C to +50°C)
   - [ ] Verify pressure interpolation accuracy (particularly near extrema)
   - [ ] Cross-check against radiosondes (existing sondes in [validation_sondes/](validation_sondes/))

3. **Backward compatibility**
   - [ ] Verify renamed `plot_*` functions produce identical output to `mirto_*` versions
   - [ ] Check file I/O compatibility with existing datasets
   - [ ] Ensure all instrument channels are still supported

4. **Legacy code**
   - [ ] Deprecate remaining `mirto_*` files (still present: `mirto_plot_tools.py`, `mirto_plot_tr.py`, `mirto_stats.py`, etc.)
   - [ ] Migrate or remove deprecated code paths

## Remaining "mirto" Files

The following files still use the "mirto" naming convention and should be reviewed for consolidation or deprecation:

```
postprocessor/validation/
├── mirto_plot_tools.py          (18 KB) - plot utilities, likely candidates for atmos_tools
├── mirto_plot_tr.py             (22 KB) - transmission/radiance plots
├── mirto_plot_validation_profiles_single.py (29 KB) - refactored but name unchanged
├── mirto_stats.py               (17 KB) - statistics functions
├── mirto_validation_rawinsonde.py (4.8 KB) - sonde comparison
└── mirto_wrf_xSection.py        (24 KB) - WRF cross-section plots
```

**Action Items**:
- [ ] Rename `mirto_plot_validation_profiles_single.py` → `plot_validation_profiles_single.py`
- [ ] Extract plot utilities from `mirto_plot_tools.py` into `atmos_tools.py` or new `plot_tools.py`
- [ ] Migrate `mirto_stats.py` functions into statistical utilities module
- [ ] Consolidate radiance/transmission plotting into `plot_tr.py`

## Development Workflow for Future Postprocessor Changes

### When Adding New Validation Features

1. **Add atmospheric calculations to `utilities/atmos/atmos_tools.py`** if they compute thermodynamic variables
2. **Import from atmos_tools in validation scripts**, don't duplicate logic
3. **Use generic naming** (e.g., `plot_*` or `validate_*`), avoid tool-specific prefixes
4. **Test against reference data** (see Outstanding Validation Tasks above)

### When Modifying Existing Validation

1. **Search for both old and new function names** — legacy code may still exist
2. **Update imports to use `utilities.atmos.atmos_tools`**
3. **Check for deprecated `mirto_*` modules** that might conflict
4. **Run validation suite** against test dataset before committing

## Testing Data & Reference Files

| Category | Location | Size | Purpose |
|----------|----------|------|---------|
| Radiosondes | [validation_sondes/](validation_sondes/) | ~1 MB | Reference for validation comparisons |
| Station coordinates | [stations.txt](stations.txt) | ~805 KB | Geographic reference |
| Region definitions | [regions_dictionary.py](regions_dictionary.py) | 46 KB | Domain mask definitions |
| Wyoming Sondes | [wyoming_sonde_downloader.py](wyoming_sonde_downloader.py) | 16 KB | Script to fetch reference data |

## Next Steps

### Short-term (Immediate)
1. **Complete testing** against all validation workflows
2. **Fix any failures** in atmospheric calculations
3. **Consolidate remaining mirto_* files** (refactor or deprecate)
4. **Update this documentation** with test results

### Medium-term (Sprint)
1. **Create validation test suite** (unit tests for atmos_tools functions)
2. **Document instrument-specific validation** (CRIS vs IASI differences)
3. **Automate validation workflows** (CI/CD integration)

### Long-term (Roadmap)
1. **Performance optimization** for large-scale validation
2. **Real-time diagnostic plotting** for operational monitoring
3. **Integration with postprocessor diagnostic module** ([diagnostic/](diagnostic/))

## References

- **Main merge**: [f38e8b5](https://github.com/Adaptive-Meteo/amethyst/commit/f38e8b5)
- **Post-merge refactoring**: [9929c84](https://github.com/Adaptive-Meteo/amethyst/commit/9929c84)
- **Atmospheric utilities**: [utilities/atmos/atmos_tools.py](../../utilities/atmos/atmos_tools.py)
- **Validation workflows**: [postprocessor/validation/](.)

---

**Maintainer**: Paolo Scaccia (<paolo.scaccia@adaptivemeteo.com>)  
**Last Reviewed**: October 28, 2024
