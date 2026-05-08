# Instrument Configuration Agent

**Specialization**: Instrument setup, channel selection, and multi-instrument support  
**Location**: [ancillary/instrument/](../instrument/) + [amethyst_config.py](../../amethyst_config.py)  
**Supported Instruments**: CRIS (default), IASI (configurable)  
**Status**: CRIS fully supported; IASI recently added (commit f903961)

## Purpose

This agent specializes in configuring AMETHYST for different satellite instruments. Instrument configuration touches every stage of the pipeline (preprocessing, processing, postprocessing) and determines which spectral channels and state variables are retrieved.

## Quick Start: Switching Instruments

To enable a new instrument or switch between CRIS and IASI:

```python
# In amethyst_config.py, modify common_vars:
common_vars = {
    "instrument": "iasi",  # Change from "cris" to "iasi"
    # ... rest of config
}
```

**Then**:
1. Verify channel list exists: [ancillary/instrument/](../instrument/)
2. Update state vector if needed (see State Vector Configuration below)
3. Recompile forward model if using new spectral ranges: `cd processor/ossfm/v3 && make`

## Project Structure for Instruments

```
ancillary/instrument/
├── CRIS                         # Primary instrument (HDF5 format)
│   ├── [channel data files]
├── IASI                         # Secondary instrument (NetCDF format)
│   ├── iasi_chList.dat         # Full IASI channel list
│   ├── iasi_chList_tr.dat      # Thermal (IR) channels only
│   └── iasi_selected_channels.py # Python module with channel selection logic
├── [other instruments/]
```

## Key Configuration Files

### 1. Central Config: amethyst_config.py

**Location**: [amethyst_config.py](../../amethyst_config.py)  
**Role**: Single source of truth for all instrument settings

**Instrument-Specific Settings**:
```python
common_vars = {
    "instrument": "cris",           # Active instrument
    "channels_list": "/path/to/channels.txt",  # Channel file for active instrument
    "selected_channels": [list of channel indices],  # Subset for retrieval
    "selected_state_vector_variables": [... state codes ...],  # What to retrieve
}

processor_vars = {
    "forward_model_version": 3,     # Must match compiled OSS module
    # Instrument-specific forward model settings
}
```

### 2. CRIS Configuration

**Location**: [ancillary/instrument/CRIS/](../instrument/CRIS/)  
**Format**: HDF5 or binary channel descriptors  
**Status**: Production-ready  
**Channel Range**: Typical ~2000 channels (full resolution)

**Configuration Pattern**:
```python
# In amethyst_config.py
common_vars = {
    "instrument": "cris",
    "channels_list": "ancillary/instrument/CRIS/cris_channels.txt",
    "selected_channels": list(range(0, 100)),  # Select first 100 channels
}
```

### 3. IASI Configuration (Recently Added)

**Location**: [ancillary/instrument/IASI/](../instrument/IASI/)  
**Format**: Text-based channel lists (`.dat` files) + Python selection module  
**Status**: Recently added (commit f903961); under validation  
**Channel Range**: 8461 channels (full resolution)

**Files**:
- `iasi_chList.dat` — Full IASI channel list
- `iasi_chList_tr.dat` — Thermal IR channels only (recommended for atmospheric retrieval)
- `iasi_selected_channels.py` — Python module for intelligent channel selection

**Configuration Pattern**:
```python
# In amethyst_config.py
common_vars = {
    "instrument": "iasi",
    "channels_list": "ancillary/instrument/IASI/iasi_chList_tr.dat",
    "selected_channels": list(range(0, 50)),  # Select first 50 thermal IR channels
}
```

## State Vector Configuration

The state vector defines what atmospheric variables are retrieved. Configuration is instrument-independent but must match selected channels.

### State Vector Codes

| Code | Variable | Abbreviation | Units | Typical Range | Notes |
|------|----------|--------------|-------|----------------|-------|
| -2 | Surface Temperature | T_skin | K | 250–320 K | Land/ocean surface |
| -1 | Surface Emissivity | ε | dimensionless | 0.7–1.0 | Wavelength-dependent |
| 0 | Temperature Profile | T(p) | K | 200–320 K | Atmosphere |
| 1 | Water Vapor | q (VMCR) | ppmv | 0–20,000 ppmv | Major uncertainty source |
| 2 | CO₂ | CO₂ (VMCR) | ppmv | 300–450 ppmv | Often fixed |
| 3 | Ozone | O₃ (VMCR) | ppmv | 0–10 ppmv | Retrieved if channels available |

### Example State Vector Configuration

```python
common_vars = {
    "selected_state_vector_variables": [
        -2,  # Surface temperature
        -1,  # Emissivity (at select channels)
        0,   # Temperature profile
        1,   # Water vapor
        # 2,  # CO2 (commented out = fixed, not retrieved)
        # 3,  # O3 (commented out = fixed)
    ],
}
```

### Creating Custom State Vectors

**Rules**:
1. Must start with -2 (T_skin) or -1 (emissivity) for radiative transfer to work
2. Must include 0 (temperature profile) for any atmospheric retrieval
3. Order matters: output jacobian dimensions follow state vector order
4. Adding variable → retrieve more DOF but slower inversion & more noise

**Process**:
1. Edit `selected_state_vector_variables` in [amethyst_config.py](../../amethyst_config.py)
2. Recompile forward model to ensure jacobian shape propagates: `cd processor/ossfm/v3 && make`
3. Test on single observation via [processor/main/amethyst_code_main.py](../../processor/main/amethyst_code_main.py)

## Channel Selection Strategies

### Full Resolution (Default)
**Use Case**: Research, when computational budget is unlimited  
**Channels**: 2000+ (CRIS) or 8461 (IASI)  
**CPU**: Slowest  
**Noise**: Highest (redundancy)

```python
common_vars = {
    "selected_channels": list(range(0, n_total_channels)),
}
```

### Optimal Spectral Sampling (OSS, Recommended)
**Use Case**: Operational, fast inversion with minimal noise  
**Channels**: 300–500 (pre-selected for information content)  
**CPU**: Fast  
**Noise**: Minimal (optimized)  
**Implementation**: Already built into forward model (v3); just select OSS-subset channels

```python
common_vars = {
    "selected_channels": [indices of OSS-selected channels],  # Provided by forward model
}
```

### Thermal IR Only
**Use Case**: Atmospheric profile retrieval (temperature, humidity)  
**Channels**: 1000–2000 (CRIS thermal range) or 4000–8461 (IASI thermal range)  
**Advantage**: Avoid solar contribution issues; cleaner inversion

```python
# Example for IASI thermal
common_vars = {
    "channels_list": "ancillary/instrument/IASI/iasi_chList_tr.dat",  # tr = thermal
}
```

### By Wavenumber Range
**Use Case**: Targeted retrieval (e.g., CO₂ at 15 μm = 667 cm⁻¹)  
**Channels**: 10–100 (user-specified)  
**Implementation**: Write custom Python script to select by wavenumber

```python
import numpy as np

# Example: Select CRIS channels between 600–700 cm-1
wavenumbers = np.loadtxt("ancillary/instrument/CRIS/wavenumbers.txt")
ch_indices = np.where((wavenumbers > 600) & (wavenumbers < 700))[0]
common_vars = {"selected_channels": ch_indices.tolist()}
```

## Observation Error Configuration

Observation error (measurement uncertainty) must match selected channels.

### Current Setup

**Location**: [utilities/create_cris_obs_err.py](../../utilities/create_cris_obs_err.py)  
**Status**: CRIS-only; IASI version needed  
**Purpose**: Generates observation error covariance matrix from instrument specs

**Issue**: If you select different channels, observation error matrix dimension mismatch!

### When to Update Observation Error

| Scenario | Action |
|----------|--------|
| **Switch instruments (CRIS → IASI)** | Run IASI version of obs_err generator |
| **Change selected_channels** | Regenerate covariance matrix for new channel subset |
| **Update instrument specs** | Modify [utilities/create_cris_obs_err.py](../../utilities/create_cris_obs_err.py) with new noise specs |

### Creating Observation Error for New Instrument

**Example: IASI observation error**

```bash
# Run existing script for CRIS
python utilities/create_cris_obs_err.py

# For IASI, create new script: utilities/create_iasi_obs_err.py
# Copy create_cris_obs_err.py and modify:
# 1. Load IASI channel specs (noise, wavelength, etc.)
# 2. Generate covariance for IASI channels
# 3. Save to ancillary/instrument/IASI/obs_error_cov.nc
```

## Adding a New Instrument

**Step-by-step process**:

### 1. Create Instrument Directory

```bash
mkdir ancillary/instrument/[NEW_INSTRUMENT]/
```

### 2. Prepare Channel Information

**Collect**:
- Channel numbers (1 to N)
- Wavenumbers or wavelengths (cm⁻¹ or μm)
- Noise characteristics (NEDT in K, for brightness temperature channels)
- Spectral response functions (if available)

**Create channel list file**:
```
# NEW_INSTRUMENT channel list
ch_number,wavenumber_cm-1,nedt_k,wavelength_um
1,400.0,0.5,25.0
2,400.5,0.5,24.98
...
```

### 3. Add to amethyst_config.py

```python
common_vars = {
    "instrument": "new_instrument",
    "channels_list": "ancillary/instrument/NEW_INSTRUMENT/channels.txt",
    "selected_channels": [list of indices],
}
```

### 4. Create Observation Error File

```bash
# 1. Create utilities/create_new_instr_obs_err.py (based on CRIS version)
# 2. Run to generate obs_error_cov.nc
# 3. Verify dimensions match selected_channels
python utilities/create_new_instr_obs_err.py
```

### 5. Update Forward Model (if needed)

If spectral range, channel structure, or OSS coefficients differ:
```bash
cd processor/ossfm/v3
# Edit oss_ir.F90 if forward model logic changes
make clean && make  # Recompile
```

### 6. Test Integration

```python
# Test script: test_new_instrument.py
from preprocessor.amethyst_preprocessor import launch_preprocessing
from processor.amethyst_processor import launch_processing
from postprocessor.validation.validate_retrievals import ValidationDataset

# 1. Preprocess test data
launch_preprocessing(config_instrument="new_instrument")

# 2. Run retrieval
launch_processing(config_instrument="new_instrument")

# 3. Validate output
# (See postprocessor/validation/ for examples)
```

## Common Configuration Mistakes

⚠️ **Critical errors** that cause crashes:

1. **Channel count mismatch**
   ```
   Error: Jacobian shape (N_channels, M_state) does not match observation error (K, K)
   ```
   **Fix**: Verify `selected_channels` dimension matches observation error covariance

2. **Instrument not found**
   ```
   Error: FileNotFoundError: ancillary/instrument/[INSTRUMENT]/channels.txt
   ```
   **Fix**: Create instrument directory and channel list file

3. **Forward model version mismatch**
   ```
   Error: ImportError: ossir module version does not match fm_version config (expected 3, got 2)
   ```
   **Fix**: Recompile forward model: `cd processor/ossfm/v3 && make`

4. **State vector channel mismatch**
   ```
   Error: Emissivity channel -1 selected but only 100 channels available
   ```
   **Fix**: Ensure `selected_state_vector_variables` matches channel count (e.g., skip emissivity if <1000 channels)

## Validation Against Reference Data

### Radiosondes (for temperature/humidity)
**Location**: [postprocessor/validation/validation_sondes/](../../postprocessor/validation/validation_sondes/)  
**Use**: Compare retrieved profiles against ground truth  
**Script**: [postprocessor/validation/validate_retrievals.py](../../postprocessor/validation/validate_retrievals.py)

### WRF Model Output (for spatial comparison)
**Location**: [postprocessor/validation/](../../postprocessor/validation/)  
**Use**: Validate spatial patterns and fields  
**Scripts**: `validate_wrf.py`, `validate_wrf_single.py`

## Performance by Instrument

| Instrument | Channels | Typical Selection | Inversion Time | Notes |
|------------|----------|-------------------|-----------------|-------|
| **CRIS** | ~2000 | 300–500 (OSS) | 5–10 sec | Production-ready; most optimized |
| **IASI** | 8461 | 300–500 (OSS) | 10–20 sec | Recently added; may need tuning |

## Related Documentation

- **Central config**: [amethyst_config.py](../../amethyst_config.py)
- **Preprocessing**: [preprocessor/amethyst_preprocessor.py](../../preprocessor/amethyst_preprocessor.py)
- **Forward model**: [processor/ossfm/FORWARD_MODEL_AGENT.md](../ossfm/FORWARD_MODEL_AGENT.md)
- **Observation error**: [utilities/create_cris_obs_err.py](../../utilities/create_cris_obs_err.py)
- **Validation**: [postprocessor/POSTPROCESSOR_MERGE.md](../../postprocessor/POSTPROCESSOR_MERGE.md)

---

**Last Updated**: May 8, 2026  
**Maintained by**: AMETHYST Development Team
