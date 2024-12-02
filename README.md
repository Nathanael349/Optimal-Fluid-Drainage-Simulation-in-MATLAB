# Fluid Mechanics Optimization: MATLAB Simulation for Drainage Efficiency
This repository contains a MATLAB-based simulation to analyze and optimize the drainage performance of a simple fluidic system. The project aims to determine the optimal discharge pipe length that minimizes drain time while maintaining efficient water travel distance.

## Overview

The project models fluid dynamics using engineering principles, including:
- **Bernoulli's Equation**
- **Frictional Losses** (via Haaland and Darcy-Weisbach equations)
- **Flow Regime Transitions** (laminar, turbulent, and transitional flow)

The simulation evaluates system behavior under varying conditions to identify an optimal design balancing drainage speed and flow distance.

## Features

- **Dynamic Flow Simulation**: Accounts for time-varying pressure, velocity, and Reynolds number.
- **Loss Calculations**: Includes frictional and minor losses in the system.
- **Optimization Analysis**: Determines the ideal pipe length and evaluates performance using a balanced ratio of drain time and travel distance.

## Usage

### Prerequisites
- MATLAB (tested with R2022a, but should work on other versions with minor adjustments)

### Running the Simulation
1. Clone this repository:
    ```bash
    git clone https://github.com/<your-username>/<repo-name>.git
    cd <repo-name>
    ```
2. Open the MATLAB script (`fluid_drainage_simulation.m`) in MATLAB.
3. Run the script to execute the simulation.

### Outputs
- Optimal discharge length
- Drain time for the selected configuration
- Graphical results comparing theoretical and experimental data

## Assumptions

1. Pipe material is cast iron with a roughness of 0.0005 m.
2. Flow regime transitions between laminar and turbulent based on the Reynolds number.
3. Quasi-steady-state analysis is applied, assuming constant velocity and height at small time intervals.
4. Losses are dominated by frictional and minor losses; entrance and exit losses are neglected.

## Results

- The optimal pipe length: **0.5 meters**
- Corresponding drain time: **4 minutes and 35 seconds**
- The system achieves a balanced trade-off between drainage speed and flow distance.
