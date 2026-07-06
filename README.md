# PlantarFlexorAssist
This repository contains the simulation results and optimization setups associated with the study:

"Quantifying the Need for Assistance in Patients with Plantar Flexor Muscle Weakness to Improve Gait Outcomes." by Amini et al..

The files are provided to facilitate reproducibility, enable additional analyses, and serve as starting points for future SCONE optimization studies.

## Repository Structure

### 1. SCONE 

The SCONE folder contains two main components:

#### Results

This directory contains simulation results for unimpaired and plantarflexor-weakness (PF) gait models.

The results of the **unimpaired** and **PF** models were originally obtained from the study:

*"Predicting gait adaptations due to ankle plantarflexor muscle weakness and contracture using physics-based musculoskeletal simulations."*

These models were used as baseline models in the present study.

In addition, two subject-specific models were generated in this work:

* **95th-percentile male model**
* **5th-percentile female model**

These models were developed to investigate the effects of plantarflexor muscle weakness across different anthropometric characteristics and to generate the corresponding impaired gait simulations.

#### Optimization Setups

This directory contains all files required to reproduce the SCONE optimization studies presented in the manuscript.

For the assisted-gait simulations, an additional assistance controller was integrated into the original gait controller. The assistance controller was implemented in **Lua** and operated alongside the baseline gait control framework during optimization and simulation.

### 2. Visualization Folder

The Visualization folder contains MATLAB scripts used to process and extract data from SCONE simulation results.

These scripts were developed to generate the figures and analyses presented in the manuscript and may also be used for further post-processing of the simulation outputs.

