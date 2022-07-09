# SPDX-FileCopyrightText: 2022 Johannes Hampp <johannes.hampp@zeu.uni-giessen.de>
#
# SPDX-License-Identifier: MIT

# -*- coding: utf-8 -*-
import pypsa

n = pypsa.Network(snakemake.input["network"])
n.consistency_check()

solver_name = "gurobi"
solver_options = {
    "method": 2,  # barrier
    "crossover": 0,
    "BarConvTol": 1.0e-5,
    "AggFill": 0,
    "PreDual": 0,
    "GURO_PAR_BARDENSETHRESH": 200,
    "threads": snakemake.threads,
}

formulation = "kirchhoff"

# NB no extra_functionality necessary since battery power has no cost
status, termination_condition = n.lopf(
    pyomo=False,
    solver_name=solver_name,
    solver_options=solver_options,
    formulation=formulation,
)
print(status, termination_condition)

n.export_to_netcdf(snakemake.output["network"])
