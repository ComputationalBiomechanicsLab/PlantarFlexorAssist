# PlantarFlexorAssist
This repository contains the simulation results and optimization setups associated with the study:

"Quantifying the Need for Assistance in Patients with Plantar Flexor Muscle Weakness to Improve Gait Outcomes." by Amini et al..

The files are provided to facilitate reproducibility, enable additional analyses, and serve as starting points for future SCONE optimization studies.

Here is a revised repository structure description reflecting the new organization:

---

# Repository Structure

## 1. SCONE Setup Files

This folder contains all files required to reproduce the SCONE optimization studies presented in the manuscript.

The optimization setups include the baseline gait controllers as well as the assistance controller used for the assisted-gait simulations. The assistance controller was implemented in Lua and integrated with the original gait controller, operating alongside the baseline control framework during optimization and simulation.

## 2. Visualization

This folder contains the MATLAB scripts used to process, analyze, and visualize the SCONE simulation results presented in the manuscript.

The folder includes:

* **SCONE Results** – A compressed archive containing the SCONE simulation results required by the visualization scripts.
* MATLAB scripts for extracting gait data, performing post-processing, and generating the figures and analyses presented in the manuscript.

> **Important:** **Before running any MATLAB visualization script, you must first extract (unzip) the `SCONE Results` archive. The scripts are configured to access the simulation results from the extracted folder, and they will not run correctly if the archive remains compressed.**

