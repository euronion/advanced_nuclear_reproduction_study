<!--
SPDX-FileCopyrightText: 2022 Johannes Hampp <johannes.hampp@zeu.uni-giessen.de>

SPDX-License-Identifier: MIT
-->

Original model and model files from Lei Duan on Feb 16, 2022, https://github.com/LDuan3008/Advanced_nuclear_2021/.

Own implementation aiming to reproduce the same results and explore further scenarios based on same methodology.

## Installation

1. Install Anaconda/Miniconda/ ... for environment management using `conda` or `mamba` commands. `mamba` recommended, e.g. by installing [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge)
2. Setup environment
```
conda env create -f envs/environment.yaml
# or
mamba env create -f envs/environment.yaml
```
3. Setup gurobi license, see [documentation](https://www.gurobi.com/documentation/).
4. Download the folder "Model and input data" from the [original paper repository](https://github.com/LDuan3008/Advanced_nuclear_2021/) and place it into the main folder.
    The folder structure should look like this:
    ```
    /duan_advanced_nuclear_2021$ tree -d -L 2
    .
    ├── .reuse
    ├── config
    ├── data
    │   ├── network_templates
    │   └── technologydata
    ├── envs
    ├── Model and input data
    │   └── Input_Data
    ├── resources
    ├── results
    │   └── validation
    ├── scenarios
    │   └── validation
    └── workflow
        ├── notebooks
        ├── rules
        └── scripts
    ```

## Running

Every action is supposed to happen in the installed activated environment.
Activate the environment by typing in command line (with `conda` or `mamba` setup):
```
conda activate advanced-nuclear
# or
mamba activate advanced-nuclear
```

### Run a single scenario

Use `snakemake` and specify scenario file(s) to be solved and be generated:

```
snakemake -call results/template_<template>/years_<years>/country_<country>/technologydata_<technologydata>/co2limit_<co2limit>/results.nc
```

or for varying PGP costs:
```
snakemake -call results/template_<template>/years_<years>/country_<country>/technologydata_<technologydata>/co2limit_<co2limit>/pgpfactor_<pgpfactor>/results.nc
```

where `<...>` indicate wildcards to substitute:

| wildcard | values / meaning | examples |
| ---      | ---              | ---      |
| `<template>` | network template to use. Defines which technologies are included and how they are connected. Any file name in `data/network_templates/` | advanced_nuclear |
| `<years>` | Weather and demand year(s) time-series to use for VRES feed-in. | "2019", "1980-2019" |
| `<country>` | Two letter country code to select demand and VRES availability time-series on. | "US", "JP" |
| `<technologydata>` | Technology data (efficiencies, cost assumptions, discount rate, CO2 reference emissions) to use. Any file name in `data/technologydata/` | "EIA" |
| `<co2limit>` | CO2 emission limit to apply. Calculated based on CO2 reference emissions. Between 1.0 (=100% emissions allowed) and 0.0 (=No emissions allowed). | "1.0","0.0","0.05" |
| `pgpfactor>` | (Optional) Factor to multiple PGP technology CAPEX with | "1.0", "0.1", "0.3" |

e.g.

```
snakemake -call results/template_advanced_nuclear_baseload/years_2019/country_US/technologydata_AdvancedNuclearEIA/co2limit_0/network.nc
```

### Running scenarios

By default the file `scenarios/default.csv` is used.
Each row entry in this file represents a scenario which will be run.
Wildcard values are set in their respective columns.
For possible wildcard values, see table above in section "Single scenario".

To calculate all scenario contained therein:

```
snakemake -call --restart-times 3 all_scenarios
```

> Note:
> The `--restart-times 3` flag is used for `snakemake` as this repository relies on `jupyter notebook` for some `snakemake` rules.
> As a lot of `jupyter notebook servers` are started, some of them conflict on trying to reserve the same `port`, causing the
> `snakemake` `rule` to fail. `--restart-times 3` tries to re-run all failed rules `3` times before considering them `failed` and
> aborting the workflow.


### Reproduce / explore scenarios similar to original paper figures

The workflow is currently setup to create scenarios and figures for two types of plot from the original paper:

* Figure 2 (`system_cost` plot)
* Supplementary Materials Figure 18 (`pgp_sensitivity` plot)

The three plots (figure 2 left side, right side, suppl. figure 18) can be recreated using the following command:

```
snakemake -call --restart-times 3 plot_validation_figures
```

The results will be plotted and stored in the folder `results/validation`.
The scenarios used for the validation are defined in `scenarios/validation/*.csv`.

### Exploring additional scenarios

Currently the following scenarios similar to those from the original paper by Duan et al. are supported:

* System costs (Figure 2 from original paper). Run:
    ```
    snakemake -call --restart-times 3 plot_all_system_costs
    ```
    * Scenarios calculating system costs with different `co2limit` for: network template / year / country / technology data
    * Scenarios defined in file: `scenarios/system_costs.csv`
    * Scenarios defining original figure from Figure 2:
    
        | network template | year | CC | technology assumptions | CO2 | 
        | ---              | ---  | -- | ---                    | --- | 
        | advanced_nuclear | 2019 | US | AdvancedNuclearEIA     | 0   |
        | advanced_nuclear | 2019 | CN | AdvancedNuclearEIA     | 0   |
        | advanced_nuclear | 2019 | DE | AdvancedNuclearEIA     | 0   |
        | advanced_nuclear | 2019 | ZA | AdvancedNuclearEIA     | 0   |
        | advanced_nuclear | 2019 | AU | AdvancedNuclearEIA     | 0   |
        | advanced_nuclear | 2019 | BR | AdvancedNuclearEIA     | 0   |
        | advanced_nuclear | 2019 | US | AdvancedNuclearEIA4000 | 0   |
        | advanced_nuclear | 2019 | CN | AdvancedNuclearEIA4000 | 0   |
        | advanced_nuclear | 2019 | DE | AdvancedNuclearEIA4000 | 0   |
        | advanced_nuclear | 2019 | ZA | AdvancedNuclearEIA4000 | 0   |
        | advanced_nuclear | 2019 | AU | AdvancedNuclearEIA4000 | 0   |
        | advanced_nuclear | 2019 | BR | AdvancedNuclearEIA4000 | 0   |
    
    * Plots for each scenario are called `pgp_sensitivity.{png,pdf}` and are created in the respective `results/` subdirectories

* PGP scenarios (Supplemantary Figure 18 from original paper). Run:
    ```
    snakemake -call --restart-times 3 plot_all_pgp_sensitivities
    ```
    * Scenarios with different PGP costs based on: network template / year / country / technology data / co2limit
    * Scenarios defined in file: `scenarios/pgp_sensitivities.csv`
    * Scenarios defining original figure from supplementary materials:
    
        | network template     | year | CC | technology assumptions     | CO2 | 
        | ---                  | ---  | -- | ---                        | --- | 
        | advanced_nuclear_PGP | 2019 | US | AdvancedNuclearEIA4000_PGP | 0   |
        | advanced_nuclear_PGP | 2019 | CN | AdvancedNuclearEIA4000_PGP | 0   |
        | advanced_nuclear_PGP | 2019 | DE | AdvancedNuclearEIA4000_PGP | 0   |
        | advanced_nuclear_PGP | 2019 | ZA | AdvancedNuclearEIA4000_PGP | 0   |
        | advanced_nuclear_PGP | 2019 | AU | AdvancedNuclearEIA4000_PGP | 0   |
        | advanced_nuclear_PGP | 2019 | BR | AdvancedNuclearEIA4000_PGP | 0   |
    
    * Plots for each scenario are called `pgp_sensitivity.{png,pdf}` and are created in the respective `results/` subdirectories

## Network templates (`data/network_templates` folder)

Different types of network templates used as basis to construct PyPSA network on which the optimisation is conducted.

* `advanced_nuclear`: Nuclear with nuclear heat source, TES and generator seperated
* `conventional_nuclear`: Nuclear with en-block nuclear, i.e. electricity is direct output of nuclear technology

Multiple suffixes allow for networks with additional technologies:
* `_baseload`: The nuclear technology is forced to run in baseload operation at 100% availability. Without all nuclear technologies can be freely dispatched
* `_PGP`: Allow for PGP (H2 electrolysis, storage and H2E)
* `_load-shifting`: Allow for load shifting in the network
* `_resisitveheater`: Allow for a resistive heater to feed surplus electricity from grid into a TES storage included with the advanced nuclear configuration

## Technology data (`data/technologydata` folder)

Different technology data assumptions (cost, efficiencies, specific emissions, discount rate).
Assumptions are such that they aim to emulate scenarios from the Duan et al. paper.
The original model uses cost input as determined as "full hourly cost", while our model
uses CAPEX+FOM+VOM+Marginal cost+lifetime as input.
The cost numbers in the files in the `technologydata` should recreate the cost assumptions used in the original model.

In general they are determined as follows:

* `var_cost` from the original model files are considered as `fuel cost`, `VOM`` are set to zero
* ``CAPEX` is calculated based on CRF (discount rate, lifetime) from `fixed_cost` and includes `FOM`; `FOM` are set to zero.

**Scenarios:**

* *ConventionalNuclearEIA.csv* (original paper)
    * Representing input data from original paper "U_Case_ConventionalNuclearEIA.csv"
    * Technologies:
        * load shedding
        * Natural Gas
        * Natural Gas with CCS
        * Wind
        * PV
        * Storage (battery)
        * Regular nuclear (~6300 USD/kW_e)

* *AdvancedNuclearEIA.csv* (original paper)
    * Representing input data from original paper "U_Case_AdvancedNuclearEIA.csv"
    * Nuclear seperated into heat source + TES + generator
    * EIA nuclear cost of ~6300 USD/kW_e as baseline
    * generator and TES cost taken from EIA SAM and heat source cost adjusted to match EIA total cost
    * Scenario for calculating Figure 2 (a-f) from original paper
    
* *AdvancedNuclearEIA4000* (original paper)
    * Representing input data from original paper "U_Case_AdvancedNuclearEIA4000.csv"
    * Nuclear seperated into heat source + TES + generator
    * Nuclear cost of ~ 4000 USD/kW_e as baseline
    * generator and TES cost taken from EIA SAM and heat source cost adjusted to match reduced total cost of ~4000 USD/kW_e
    * Scenario for calculating Figure 2 (g-l) from original paper
    
* *AdvancedNuclearEIA4000_PGP.csv* (original paper)
    * Representing input data from original paper "U_Case_AdvancedNuclearEIA4000_PGP.csv"
    * Based on *AdvancedNuclearEIA.csv* with reduced cost for advanced nuclear setup (cf. *AdvancedNuclearEIA4000*)
    * Additional options for storing electricity as hydrogen
    * Power to gas (P2G): Electrolysis
    * Gas storage: Cavern
    * Gas to power (G2P): Fuel cell
    * Scenario for calculating Supplementary Figure 18 from original paper
    
* *ConventionalNuclearEIA_flexible-NPP.csv* (new)
    * Based on *ConventionalNuclearEIA*
    * NPPs are not considered baseload, but model economics allow them to operate flexibly if advantageous
    * Nuclear cost disaggregated into investment, FOM and variable cost
    * Cost disaggregation based on "MASTER_cost_assumptions.xlsx" file from original paper repository

* *AdvancedNuclearEIA_IRENA-PV.csv* (new)
    * Based on *AdvancedNuclearEIA*
    * Lower PV cost assumed (based on IRENA: RES cost 2020 report)
    
* *AdvancedNuclearEIA_IRENA-wind.csv* (new)
    * Based on *AdvancedNuclearEIA*
    * Lower wind cost assumed (based on IRENA: RES cost 2020 report)

* *AdvancedNuclearEIA_IRENA-RES.csv* (new)
    * Based on *AdvancedNuclearEIA*
    * Lower RES (PV, wind) cost assumed (based on IRENA: RES cost 2020 report)

* *AdvancedNuclearEIA4000_PGP_IEA-H2.csv* (new)
    * Based on *AdvancedNuclearEIA4000_PGP*
    * Lower PGP cost based on IEA "Future of Hydrogen" "Today" numbers
        * Water electrolysis for today
        * Compressed hydrogen storage for long-term storage and Gas-to-power
    
* *AdvancedNuclearEIA4000_PGP_IRENA-RES_IEA-H2.csv* (new)
    * Based on *AdvancedNuclearEIA4000_PGP*
    * Lower RES (PV, wind) cost assumed (based on IRENA: RES power generation cost 2020 report)
    * Lower PGP cost based on IEA "Future of Hydrogen" "Today" numbers
        * Water electrolysis for today
        * Compressed hydrogen storage for long-term storage and Gas-to-power
    
* *AdvancedNuclearEIA_PGP.csv* (new)
    * Based on *AdvancedNuclearEIA4000_PGP.csv*
    * PGP option based on *AdvancedNuclearEIA4000_PGP*
    * Higher advanced nuclear cost corresponding EIA cost at ~6300 USD/kW_e (cf. *AdvancedNuclearEIA*)

* *AdvancedNuclearEIA_PGP_IRENA-RES_IEA-H2.csv* (new)
    * Based on *AdvancedNuclearEIA4000_PGP.csv*
    * PGP option based on lower cost assumptions from *AdvancedNuclearEIA4000_PGP_IRENA-RES_EIA-H2*
    * Higher advanced nuclear cost corresponding EIA cost at ~6300 USD/kW_e (cf. *AdvancedNuclearEIA*)

> TODO: Describe additional scenarios added

## Caveats

* natural gas with CCS implemented in a quick and dirty way by increasing the gas to electricity efficiency. Doesn't affect costs with which the technology is included in the model, only gas consumption (and thus CO2 constraint)

## Configuration

The file `config/config.default.yaml` allows for some configuration of the workflow.
Current options:

* a flag to reproduce behaviour specific to original Duan et al paper.
* Resolution (x-axis steps) used for PGP and system cost plots
* Colors used for plots

## LICENSE and COPYRIGHT

(c) 2022 Johannes Hampp (johannes.hampp@zeu.uni-giessen.de)

Source code is licensed under MIT license.
CC0-1.0 license on insignificant files.
CC-BY-4.0 on all input data files.

For licensing information refer to the SPDX identifiers of the individual files or if not present
the information in `.reuse/h5dep`.
