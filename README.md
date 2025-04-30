# BabyBCal Simulation Analysis

Analysis code for the BabyBCal Beam Test Simulation at FTBF.

## Description

The `histos.cpp` file peforms the primary analysis, creates histograms for fitting and plotting, and writes the histograms to a `Histos.root` file.
The remaining `.cpp` files produce fits and various for various scans and require the necessary histograms be stored in the `Histos.root` file.

## Getting Started

### Dependencies

* The analysis scripts are run in a root terminal in the eic-shell.

### Executing program

#### Perform the primary analysis.
* Create a `sims/` directory that contains the `MCSimulation.edm4hep.root` file.
* Edit `histos.cpp` line 506 to include the name of your `MCSimulation.edm4hep.root` file.
* Edit `histos.cpp` line 517 to include the appropriate parameters of your simulation in the `createHistos()` function call.
* Run `histos.cpp` from a root terminal.
```
root -b -l -q histos.cpp
```
#### Perform the fits and plots for a paricular scan, example: `energyScanPlot.cpp`.
* Create a `FTBF_data/` directory that contains the `Data.root` file with histograms.
* Edit `energyScanPlot.cpp` lines 384-386 to include the correct histogram names from your `Data.root` file.
* Edit `energyScanPlot.cpp` lines 393-405 to include the correct histogram names from your `Histos.root` file.
* Run `energyScanPlot.cpp` from a root terminal.
```
root -b -l -q energyScanPlot.cpp
```
