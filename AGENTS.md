# AMETHYST Agent Instructions

**AMETHYST** is a satellite atmospheric retrieval system that inverts hyperspectral infrared observations (CRIS/IASI instruments) into atmospheric profiles using a 1DVar variational algorithm with Optimal Spectral Sampling (OSS) radiative transfer.

## Project Architecture

The system is organized as a four-stage pipeline:

1. **Preprocessor** → Generates First Guess (FG), A Priori, and Field-of-View data from WRF models
2. **Processor** → Parallel 1DVar inversion engine (main computation)
3. **Postprocessor** → Validation diagnostics and mapping (minimal)
4. **Ancillary** → Static data: climatologies, instrument specs, OSS coefficients
5. **Utilities** → IASI reader library, geometry helpers, standalone tools

## Key Entry Points

| File | Purpose |
|------|---------|
| [processor/amethyst_processor.py](processor/amethyst_processor.py) | Main orchestrator—spawns parallel worker processes |
| [processor/main/amethyst_code_main.py](processor/main/amethyst_code_main.py) | Core 1DVar inversion logic (`core` class, `invert()` method) |
| [preprocessor/amethyst_preprocessor.py](preprocessor/amethyst_preprocessor.py) | Data preparation entry point |
| [amethyst_config.py](amethyst_config.py) | Central configuration—all paths, instrument settings, state vector config |

## Project Conventions

- **Main modules**: `amethyst_[component].py`
- **Data classes**: CamelCase in [processor/dobjects/](processor/dobjects/) (e.g., `Solar`, `Hitran`, `ObservationError`)
- **Functions**: `lowercase_with_underscores`
- **Configuration**: Module-level dicts (`common_vars`, `preprocessor_vars`, `processor_vars`, `postprocessor_vars`)
- **State vector codes**: -2=Surface Temperature, -1=Emissivity, 0=Temperature, 1=H₂O, 2=CO₂, 3=O₃
- **File formats**: NetCDF4 (internal), HDF5 (CRIS input)
- **Parallelization**: `multiprocessing.Process` with queue communication

## Environment Setup & Build

```bash
# 1. Set up conda environment
conda env create -f misc/environment.yml

# 2. Configure Python environment and paths
cd misc && ./set_env.sh

# 3. Compile Fortran forward model (requires gfortran, f2py, MKL)
cd processor/ossfm/v3 && make  # Creates ossir.cpython-*.so

# 4. Run processor (example)
python processor/amethyst_processor.py [arguments]
```

**Critical**: The `set_env.sh` script must be sourced to configure `PYTHONPATH` before running Python code.

## State of Development (devel branch, latest: 9929c84)

Recent work on the `devel` branch focused on **merging two postprocessor versions** and migrating to the **atmos library**:
- Removed `mirto` from filenames
- Updated function calls to use `atmos` library instead of legacy code
- Postprocessor merge still under testing

## Dependencies

- **Python**: 3.11
- **Scientific**: numpy 1.26.4, scipy 1.13.1, pandas 2.2.2, netCDF4 1.6.2
- **Instruments**: piasi_reader 0.9.7 (custom package in [utilities/piasi_reader-0.9.7](utilities/piasi_reader-0.9.7/))
- **Forward model**: Fortran 90 (compiled via f2py), requires gfortran and MKL
- **Data**: HDF5 1.10.6, NetCDF 4.8.1

## Critical Pitfalls & Gotchas

⚠️ **Must address before code changes:**

1. **Hard-coded paths** → Config expects `/mnt/satellite/amethyst_test_data/` (reconfigure if different)
2. **Forward model version mismatch** → `fm_version` setting must match compiled OSS module version
3. **Silent worker process crashes** → Always check logs; multiprocessing failures are hard to debug
4. **PYTHONPATH setup required** → Missing `set_env.sh` → import failures
5. **Ancillary file dependencies** → Missing `.nc` files cause immediate crashes (126MB+ data required)
6. **Channel list compatibility** → Must match instrument channels to `selected_state_vector_variables`
7. **WRF input format strict** → FG generator expects specific NetCDF structure
8. **Deprecated code everywhere** → `/deprecated/` directories are legacy; do not extend or rely on
9. **Configuration reload** → Changes to [amethyst_config.py](amethyst_config.py) require Python restart
10. **Observation error matrix dimensions** → Must exactly match selected channels

## Common Development Tasks

### Adding a new atmospheric variable

1. Define state vector code in [amethyst_config.py](amethyst_config.py) (e.g., 4=SO₂)
2. Add to `selected_state_vector_variables` dict
3. Implement retrieval logic in [processor/main/amethyst_code_main.py](processor/main/amethyst_code_main.py)
4. Update ancillary files if needed

### Debugging inversion failures

1. Check `processor/` logs for worker process crashes
2. Verify [amethyst_config.py](amethyst_config.py) paths exist
3. Confirm ancillary files match instrument/channel selection
4. Use `-v` flag in processor for verbose output

### Supporting a new instrument

1. Add channel list to [ancillary/instrument/](ancillary/instrument/)
2. Update `common_vars["instrument"]` in [amethyst_config.py](amethyst_config.py)
3. If radiative transfer differs, regenerate OSS coefficients or modify [processor/ossfm/](processor/ossfm/)
4. Update observation error matrix in [utilities/create_cris_obs_err.py](utilities/create_cris_obs_err.py)

## Testing & Validation

- Unit tests are sparse; see [preprocessor/fg_generator/deprecated/wrf2firstguess/tests/](preprocessor/fg_generator/deprecated/wrf2firstguess/tests/) for examples (legacy)
- Validation patterns in [postprocessor/validation/](postprocessor/validation/) are the reference for output checking
- Manual profiling with `cProfile` is available via flags in [processor/amethyst_processor.py](processor/amethyst_processor.py)

## Performance Considerations

- **Parallelization via multiprocessing** (not threading) because Fortran forward model is not thread-safe
- **Bottleneck**: Forward model radiative transfer (~1000 calls per observation)
- **OSS library** pre-computes optimal spectral sampling to accelerate radiative transfer
- Monitor queue communication overhead in [processor/amethyst_processor.py](processor/amethyst_processor.py) for large datasets

## Architecture Decision Points

- **Why pipeline architecture?** Decoupling enables independent scaling: preprocessor runs once, processor parallelizes per observation
- **Why Python + Fortran?** Python orchestrates I/O & data; Fortran handles numerically intensive radiative transfer
- **Why central config dict?** Single source of truth for paths/params, simplifies parameter passing through nested objects
- **Why multiprocessing?** Fortran forward model designed for isolated processes, not threads

## Specialized Agents

For domain-specific development, refer to these specialized agent guides:

| Agent | Location | Purpose |
|-------|----------|---------|
| **Forward Model Agent** | [processor/ossfm/FORWARD_MODEL_AGENT.md](processor/ossfm/FORWARD_MODEL_AGENT.md) | OSS radiative transfer development, Fortran/Python optimization, debugging radiative transfer issues |
| **Instrument Config Agent** | [ancillary/instrument/INSTRUMENT_CONFIG_AGENT.md](ancillary/instrument/INSTRUMENT_CONFIG_AGENT.md) | Instrument setup, channel selection, adding new instruments (CRIS, IASI, custom) |
| **Postprocessor Merge Docs** | [postprocessor/POSTPROCESSOR_MERGE.md](postprocessor/POSTPROCESSOR_MERGE.md) | Understanding the postprocessor consolidation, atmospheric utilities library (atmos_tools), validation workflow changes |

---

**Questions?** Refer to the [amethyst_config.py](amethyst_config.py) for all configuration options and [misc/](misc/) directory for environment setup details.
