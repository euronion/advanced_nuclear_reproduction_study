# SPDX-FileCopyrightText: 2022 Johannes Hampp <johannes.hampp@zeu.uni-giessen.de>
#
# SPDX-License-Identifier: MIT

# Load definition of scenarios to run and plot
system_costs = Paramspace(
    pd.read_csv("scenarios/system_costs.csv", dtype=str), param_sep="_"
)


rule plot_system_costs:
    message:
        "Plotting co2limit vs. system cost / technology cost shares."
    input:
        # Create all networks with configured "co2limit" resolution related to this scenario plot
        networks=expand(
            "results/{wcp}/co2limit_{co2limit}/network.nc",
            co2limit=config["system_costs"]["co2limits"],
            wcp=system_costs.wildcard_pattern,
        ),
    output:
        figure=[
            f"results/{system_costs.wildcard_pattern}/system_costs.png",
            f"results/{system_costs.wildcard_pattern}/system_costs.pdf",
        ],
    log:
        f"log/{system_costs.wildcard_pattern}/plot_system_costs.log",
    notebook:
        "../notebooks/plot_system_costs.py.ipynb"


aux_ps = Paramspace(
    system_costs.dataframe.drop(columns="country").drop_duplicates(),
    filename_params="*",
    param_sep="_",
)


# Rule to create all plots for all scenarios defined via 'system_costs' Paramspace
rule plot_all_system_costs:
    input:
        # Create a figure for all scenarios defined using paramspace "system_costs"
        expand(
            "results/{params}/system_costs.{ext}",
            ext=["pdf"],
            params=system_costs.instance_patterns,
        ),
    output:
        expand(
            "results/system_costs_{params}.pdf",
            params=aux_ps.instance_patterns,
        ),
    notebook:
        "../notebooks/combine_pdf_plots.py.ipynb"
