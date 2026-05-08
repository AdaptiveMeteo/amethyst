# AMETHYST

AMETHYST is the evolution of the proprietary Mirto retrieval system for hyperspectral infrared observations. Mirto used an Optimal Spectral Sampling (OSS) forward model; AMETHYST 1.0 uses RTTOV for radiative transfer, and AMETHYST 2.0 is planned to transition to Helios-RT.

## Overview

The AMETHYST pipeline is organized into four main components:

- `preprocessor/`: prepares input fields such as first guess, a priori state, and field-of-view data from model output
- `processor/`: executes the core inversion engine, including the 1DVar algorithm and radiative transfer (RTTOV in AMETHYST 1.0, with Helios-RT planned for AMETHYST 2.0)
- `postprocessor/`: validates retrievals, produces diagnostic plots, and compares results with reference data
- `ancillary/`: contains static data, instrument definitions, climatologies, and forward model coefficients
- `utilities/`: helper tools, IASI reader support, atmospheric utilities, and observation error generation

## Key Files

- `amethyst_config.py`: central configuration file for instrument, state vector, and processing settings
- `processor/amethyst_processor.py`: main processor orchestrator, handles multiprocessing and inversion workflow
- `processor/main/amethyst_code_main.py`: core retrieval logic and 1DVar inversion
- `processor/main/ForwardModel.py`: forward model interface used by the inversion core
- `processor/ossfm/v3/`: legacy OSS forward model implementation retained for reference; AMETHYST 1.0 uses RTTOV and AMETHYST 2.0 will use Helios-RT
- `postprocessor/POSTPROCESSOR_MERGE.md`: documentation for the merged postprocessor and atmospheric utilities
- `processor/ossfm/FORWARD_MODEL_AGENT.md`: forward model development guide
- `ancillary/instrument/INSTRUMENT_CONFIG_AGENT.md`: instrument configuration guide

## System Description

AMETHYST ingests satellite radiances and ancillary model data, then performs atmospheric retrievals using a variational method. The project combines:

- physics-based radiative transfer via RTTOV in AMETHYST 1.0, with a planned transition to Helios-RT in AMETHYST 2.0
- inversion of temperature, moisture, and trace gas state variables
- instrument-dependent channel selection and observation error handling
- multiprocessing to scale retrievals across many profiles

The system is designed to support multiple instruments, with CRIS as the default and IASI support included.

## Setup

1. Create the conda environment:

   ```bash
   conda env create -f misc/environment.yml
   ```

2. Configure the environment:

   ```bash
   cd misc
   ./set_env.sh
   ```

3. Compile the OSS forward model:

   ```bash
   cd processor/ossfm/v3
   make
   ```

## Usage

Typical workflow:

1. Configure instrument and retrieval settings in `amethyst_config.py`
2. Run the preprocessor to generate FG/apriori/FOV inputs
3. Run the processor to perform retrievals
4. Run postprocessor validation as needed

## Getting Started

Example command sequence:

```bash
cd /home/paoloa/repository/amethyst
conda env create -f misc/environment.yml
cd misc && ./set_env.sh
cd processor/ossfm/v3 && make
cd /home/paoloa/repository/amethyst
python preprocessor/amethyst_preprocessor.py
python processor/amethyst_processor.py
python postprocessor/validation/validate_retrievals.py
```

Adjust the last three commands based on your configuration and available input data.

## Notes

- The `processor/ossfm/v3` module must be compiled before running the processor.
- The `misc/set_env.sh` script must be sourced to set up `PYTHONPATH` and environment variables.
- Many validation and diagnostic scripts still use legacy naming and are under consolidation.

## Useful Documentation

- `AGENTS.md`: agent-oriented guide to the repository and development conventions
- `postprocessor/POSTPROCESSOR_MERGE.md`: postprocessor merge and validation documentation
- `processor/ossfm/FORWARD_MODEL_AGENT.md`: forward model development documentation
- `ancillary/instrument/INSTRUMENT_CONFIG_AGENT.md`: instrument configuration documentation

## License

The repository appears to be licensed under GPL-style terms as indicated by the source files.
