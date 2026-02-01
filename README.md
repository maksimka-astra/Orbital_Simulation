# Orbital_Simulation

This repository contains the Python code for an IB Mathematics Internal Assessment.  
It simulates the motion of Earth around the Sun to study the accuracy and stability of numerical integration methods.

## Overview

The simulation implements two numerical methods:

- Euler method (first-order)
- Fourth-order Runge–Kutta method (RK4)

It compares the methods in terms of:
- orbital trajectories
- relative and absolute radius errors

The total simulation time is fixed at 10 years, and the time step Δt is specified by the user.

## Requirements

- **Python 3.x**
- **NumPy** (for numerical calculations)
- **Matplotlib** (for plotting figures)
- **csv** (standard Python library for CSV output)

These libraries must be installed to run the code successfully.

## Input

The user is prompted to enter:
- time step Δt (in years)

The number of iterations is computed automatically to maintain the fixed simulation duration.

## Output (CSV Files)

The code generates two CSV files:

- `euler.csv` — results from the Euler method  
- `rk4.csv` — results from the RK4 method

Each file includes:
- time (years)
- x-position (AU)
- y-position (AU)
- relative radius error
- absolute radius error

These files were used to create the figures included in the report.

## Visualisation

The script produces plots for:
- orbital trajectories
- relative radius error over time
- absolute radius error over time

## Notes

This code is for educational and analytical purposes and supports the mathematical analysis in the Internal Assessment.
