# SPDX-FileCopyrightText: 2022 Johannes Hampp <johannes.hampp@zeu.uni-giessen.de>
#
# SPDX-License-Identifier: MIT

import numpy as np

pgp_sensitivities = Paramspace(
    pd.read_csv("scenarios/pgp_sensitivities.csv", dtype=str), param_sep="_"
)


rule modify_pgp_cost:
    input:
        network=f"resources/attach_co2limit/{pgp_sensitivities.wildcard_pattern}/network.nc",
    output:
        network=f"resources/modify_pgp_cost/{pgp_sensitivities.wildcard_pattern}/pgpfactor_{{pgpfactor}}/network.nc",
    notebook:
        "../notebooks/modify_pgp_cost.py.ipynb"


rule solve_pgp_sensitivity_network:
    input:
        network=f"resources/modify_pgp_cost/{pgp_sensitivities.wildcard_pattern}/pgpfactor_{{pgpfactor}}/network.nc",
    output:
        network=f"results/{pgp_sensitivities.wildcard_pattern}/pgpfactor_{{pgpfactor}}/network.nc",
    script:
        "../scripts/solve_network.py"


rule plot_pgp_sensitivity:
    input:
        results=expand(
            "results/{wcp}/pgpfactor_{pgpfactor}/results.csv",
            pgpfactor=config["pgp_sensitivities"]["pgp_factors"],
            wcp=pgp_sensitivities.wildcard_pattern,
        ),
    output:
        figure=expand(
            "results/{wcp}/pgp_sensitivitity.{ext}",
            wcp=pgp_sensitivities.wildcard_pattern,
            ext=["pdf", "png"],
        ),
    notebook:
        "../notebooks/plot_pgp_sensitivity.py.ipynb"


# Rule to create all plots for all scenarios defined via 'pgp_sensitivitites' Paramspace
# and combine them for the same countries
aux_ps_pgp = Paramspace(
    pgp_sensitivities.dataframe.drop(columns="country").drop_duplicates(),
    filename_params="*",
    param_sep="_",
)


rule plot_all_pgp_sensitivities:
    input:
        # Create a figure for all scenarios defined using paramspace "pgp_sensitivities"
        expand(
            "results/{params}/pgp_sensitivitity.{ext}",
            ext=["pdf"],
            params=pgp_sensitivities.instance_patterns,
        ),
    output:
        expand(
            "results/pgp_sensitivitity_{params}.pdf",
            params=aux_ps_pgp.instance_patterns,
        ),
    notebook:
        "../notebooks/combine_pdf_plots.py.ipynb"
