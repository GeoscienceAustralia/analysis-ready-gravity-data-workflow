# Analysis Ready Gravity Data Workflow

This repository introduces an open-source tool for computing regional gravimetric geoids using gravity observations. Our primary goal is to create a platform for analysis-ready gravity data, where algorithms and code for working with gravity data are openly shared and improved. The initial release features a tile-wise least-squares collocation (LSC) method based on gravity anomaly observations.

## Overview

This self-contained codebase provides three major steps for geoid calculation:
1. **Topographic corrections to gravity anomalies**
2. **Calculation of covariance functions**
3. **Performing least squares collocation (LSC)**

The code is designed as a starting point for research and development in geoid calculation using a variety of gravity data, such as terrestrial gravity anomalies, gravity anomalies from satellite altimetry, airborne gravity anomalies, and gravity gradients.

## Workflow

The workflow for generating the geoid from gravity observations can be found in `RunDemo.m`. This script is highly flexible and can be adapted to different regions and processing requirements. Designed to be run on a personal computer, `RunDemo.m` allows users to specify an area of interest and directories where the final plots will be saved.

In `RunDemo.m`, all necessary specifications are defined. The script first calls the functions needed to apply topographic corrections to the gravity observations within the area of interest. It then performs the geoid calculation by dividing the area into tiles, which reduces the matrix sizes for easier inversion.

At the end of the process, `RunDemo.m` combines all the tiles and outputs the geoid in TIFF format for the specified region.

## Remove-Predict-Restore

In summary, for computing a geoid from gravity anomalies, the process can be conceptualized as a sequence of "remove-predict-restore" operations, where the Global Gravity Model (GGM) and topographic effects are removed, a geoid is predicted (here with LSC), and then the effects are restored to obtain the final geoid model. The `functions` folder provides all the MATLAB functions to perform these three steps for geoid calculations.

## Usage

The code can be run using the `Run*.m` scripts.

## Citations

If you use this code in your research, please contact [neda.darbeheshti@ga.gov.au](mailto:neda.darbeheshti@ga.gov.au) for citation information.



