# Forward Model Development Agent

**Specialization**: Optimal Spectral Sampling (OSS) Forward Model  
**Location**: [processor/ossfm/](../ossfm/)  
**Language**: Python + Fortran F90  
**Compilation**: f2py + gfortran

## Purpose

This agent specializes in developing, optimizing, and debugging the Optimal Spectral Sampling (OSS) forward model—the computationally intensive radiative transfer core of AMETHYST. The forward model is the bottleneck for inversion performance (~1000 radiative transfer calls per satellite observation).

## Forward Model Architecture

```
processor/ossfm/
├── v2/                          # Previous version (legacy)
├── v3/                          # Current version (active)
│   ├── oss_ir.F90              # Fortran 90 OSS radiative transfer
│   ├── ossFM.py                # Python wrapper for F90 module
│   ├── Makefile                # Compilation via f2py
│   └── Makefile.bkp            # Backup
└── doc/                         # Documentation
```

### Current Version (v3)

**Language**: Fortran 90 + Python wrapper  
**Status**: Active (set in [amethyst_config.py](../../amethyst_config.py) with `fm_version: 3`)  
**Purpose**: Fast radiative transfer using pre-computed OSS coefficients

### Legacy Version (v2)

**Status**: Do not use; retained for reference  
**Why**: API compatibility; switch requires `fm_version: 2` in config + OSS recompilation

## Key Files & Responsibilities

### Core Fortran Module: oss_ir.F90

**Location**: [processor/ossfm/v3/oss_ir.F90](../../processor/ossfm/v3/oss_ir.F90)  
**Lines of Code**: ~100-500 (estimate, typical for OSS)  
**Purpose**: Implements fast radiative transfer via optimal spectral sampling

**Key Subroutines** (typical OSS patterns):
- `oss_forward()` — Main forward radiative transfer operator
- `oss_jacobian()` — Compute jacobians (sensitivity of radiances to state variables)
- `oss_init()` — Load OSS coefficients (precomputed tables)
- `oss_od()` — Optical depth calculations

**Physical Calculations**:
- Absorber optical depths (function of pressure, temperature, concentration)
- Atmospheric transmittance (layer-by-layer integration)
- Surface/clouds contributions
- Jacobians for retrieved state variables

### Python Wrapper: ossFM.py

**Location**: [processor/ossfm/v3/ossFM.py](../../processor/ossfm/v3/ossFM.py)  
**Purpose**: f2py-generated interface to Fortran module  
**Compiled Output**: `ossir.cpython-311-x86_64-linux-gnu.so` (or similar, depending on Python version)

**Typical Methods**:
```python
ossir.oss_forward(...)      # Radiance forward model
ossir.oss_jacobian(...)     # Jacobian computation
```

### Makefile: Compilation Instructions

**Location**: [processor/ossfm/v3/Makefile](../../processor/ossfm/v3/Makefile)  
**Purpose**: Compile Fortran code to Python-callable shared object

**Key Build Variables**:
- `f2py3` — f2py command (Python 3 variant)
- `gfortran` — GNU Fortran compiler (must support F90)
- `MKL_FLAGS` — Intel MKL BLAS/LAPACK link flags (for performance)

**Build Process**:
```bash
cd processor/ossfm/v3
make clean
make              # Generates ossir.cpython-*.so
```

**Troubleshooting Build**:
1. **f2py not found**: `pip install numpy` (f2py comes with NumPy)
2. **gfortran not found**: `apt-get install gfortran` (Linux) or `brew install gcc` (macOS)
3. **MKL not found**: Install Intel MKL or link against system BLAS (modify Makefile `MKL_FLAGS`)

## Common Development Tasks

### Debugging Radiative Transfer Issues

**Symptom**: Retrieved profiles unrealistic, large residuals  
**Root Causes** (in priority order):

1. **OSS Coefficient Mismatch**
   - Check that ancillary OSS coefficients match forward model version
   - Location: [ancillary/forward_model/](../../ancillary/forward_model/)
   - Compare coefficient file version against `fm_version` in config

2. **Numerical Precision**
   - Verify Fortran code uses appropriate REAL precision (typically REAL*8 / double precision)
   - Check for underflow/overflow in optical depth calculations
   - Add debug output in Fortran code (via `PRINT *` statements)

3. **Input Array Dimensions**
   - Jacobian computation: shape mismatch between state vector and output jacobian
   - Trace through [processor/main/ForwardModel.py](../../processor/main/ForwardModel.py) for shape transformations

### Optimizing Forward Model Performance

**Bottleneck Location**: Radiative transfer loop (typical: 1000 calls/observation, ~10 sec/observation on 2019 hardware)

**Optimization Strategies**:
1. **Vectorize Fortran loops** if not already done (f2py passes NumPy arrays efficiently)
2. **Use BLAS/LAPACK** for matrix operations via MKL (already compiled in Makefile)
3. **Cache OSS coefficients** (already done at module load time)
4. **Profile with Fortran profiler**: `gfortran -pg` + `gprof`

### Updating OSS Coefficients

**When**: When instrument channels or spectral range changes  
**How**:
1. Regenerate OSS coefficients using (typically external) spectroscopic database
2. Place coefficients in ancillary format under [ancillary/forward_model/](../../ancillary/forward_model/)
3. Update [amethyst_config.py](../../amethyst_config.py) path reference
4. **Recompile** Fortran module if coefficient loading code changed: `cd processor/ossfm/v3 && make`

### Adding New Jacobian Variables

**Example**: Adding jacobian w.r.t. a new absorber (e.g., SO₂)

**Steps**:
1. **Fortran implementation**:
   - Add subroutine `oss_jacobian_so2()` in `oss_ir.F90`
   - Integrate with main `oss_jacobian()` routine
   
2. **Update state vector** in [amethyst_config.py](../../amethyst_config.py)
   - Add SO₂ code (e.g., 4) to `selected_state_vector_variables`
   
3. **Recompile**:
   ```bash
   cd processor/ossfm/v3 && make clean && make
   ```
   
4. **Test**:
   - Run single observation through [processor/main/amethyst_code_main.py](../../processor/main/amethyst_code_main.py)
   - Verify jacobian output dimensions

## Integration Points

### Where Forward Model is Called

| Caller | Purpose | File |
|--------|---------|------|
| **1DVar Inversion Loop** | Compute radiance/jacobian at each iteration | [processor/main/amethyst_code_main.py](../../processor/main/amethyst_code_main.py) |
| **Observation Error Class** | Precompute jacobians for error covariance | [processor/dobjects/ObservationError.py](../../processor/dobjects/ObservationError.py) |
| **ForwardModel Class** | High-level radiative transfer interface | [processor/main/ForwardModel.py](../../processor/main/ForwardModel.py) |

### Import Pattern

```python
# In processor/main/ForwardModel.py (typical usage)
from processor.ossfm.v3 import ossir    # Compiled Fortran module

# Usage
radiance, jacobian = ossir.oss_forward(...)
```

**Gotcha**: If `fm_version` is set to 2 in config, import will fail unless v2 module is compiled. Always sync `fm_version` and compiled module version.

## Testing Forward Model Changes

### Unit Test Pattern

```python
# Test OSS forward model independently
import numpy as np
from processor.ossfm.v3 import ossir

# Input state vector
temperature = np.array([250, 260, 270, 280, 290, 300, 310, 320])
pressure = np.array([100, 150, 200, 300, 500, 700, 850, 1000])
h2o_vmr = np.array([1e-6, 2e-6, 5e-6, 1e-5, 2e-5, 5e-5, 1e-4, 1e-3])

# Call forward model
radiance = ossir.oss_forward(temperature, pressure, h2o_vmr, ...)

# Validate output
assert radiance.shape == (expected_channels,), f"Shape mismatch: {radiance.shape}"
assert np.all(radiance > 0), "Negative radiances detected (nonphysical)"
```

### Integration Test Pattern

1. **Use validation observations** from [postprocessor/validation/](../../postprocessor/validation/)
2. **Compare jacobians** against numerical finite-difference approximations
3. **Check jacobian signs** (e.g., dRad/dTemp should be positive for most IR channels)

## Troubleshooting Guide

| Issue | Diagnosis | Solution |
|-------|-----------|----------|
| ImportError: `ossir` module not found | Forward model not compiled | `cd processor/ossfm/v3 && make` |
| TypeError: Dimension mismatch in `oss_jacobian()` | State vector shape incorrect | Check `selected_state_vector_variables` in config |
| Negative radiances | Numerical instability in Fortran | Check for underflow; verify OSS coefficients |
| 10x slower than expected | Not using optimized BLAS | Recompile with MKL: check Makefile `MKL_FLAGS` |
| Jacobian NaN values | Division by zero or invalid state | Check for zero pressure or negative temperature in input |

## Performance Profiling

### Profile Fortran Code

```bash
cd processor/ossfm/v3
gfortran -pg -O2 oss_ir.F90 -o oss_profile
./oss_profile < input.txt
gprof oss_profile gmon.out
```

### Profile Python Wrapper

```python
import cProfile
from processor.main.amethyst_code_main import core

profiler = cProfile.Profile()
profiler.enable()
core().invert()
profiler.disable()
profiler.print_stats(sort='cumulative')
```

## Related Documentation

- **Main inversion loop**: [processor/main/amethyst_code_main.py](../../processor/main/amethyst_code_main.py)
- **Forward model wrapper**: [processor/main/ForwardModel.py](../../processor/main/ForwardModel.py)
- **Observation error setup**: [processor/dobjects/ObservationError.py](../../processor/dobjects/ObservationError.py)
- **Ancillary coefficients**: [ancillary/forward_model/](../../ancillary/forward_model/)
- **Build environment**: [misc/set_env.sh](../../misc/set_env.sh)

---

**Last Updated**: May 8, 2026  
**Maintained by**: AMETHYST Development Team
