# SPDX-FileCopyrightText: 2022 Johannes Hampp <johannes.hampp@zeu.uni-giessen.de>
#
# SPDX-License-Identifier: MIT

# Load definition of scenarios to run and plot
fig2 = Paramspace(
    pd.read_csv("scenarios/validation/fig2.csv", dtype=str), param_sep="_"
)

sfig18 = Paramspace(
    pd.read_csv("scenarios/validation/sfig18.csv", dtype=str), param_sep="_"
)

aux_ps_fig2 = Paramspace(
    fig2.dataframe.drop(columns="country").drop_duplicates(),
    filename_params="*",
    param_sep="_",
)

# Rule to create all plots for figure 2 and then merge all country-specific plots per scenario
rule plot_fig2:
    input:
        expand(
            "results/{params}/system_costs.{ext}",
            ext=["pdf"],
            params=fig2.instance_patterns,
        ),
    output:
        expand(
            "results/validation/fig2_{params}.pdf",
            params=aux_ps_fig2.instance_patterns,
        ),
    notebook:
        "../notebooks/combine_pdf_plots.py.ipynb"

# Rule to create all plots for supplementary figure 18 and then merge all country-specific plots per scenario
aux_ps_sfig18 = Paramspace(
    sfig18.dataframe.drop(columns="country").drop_duplicates(),
    filename_params="*",
    param_sep="_",
)

rule plot_sfig18:
    input:
        expand(
            "results/{params}/pgp_sensitivitity.{ext}",
            ext=["pdf"],
            params=sfig18.instance_patterns,
        ),
    output:
        expand(
            "results/validation/sfig18_{params}.pdf",
            params=aux_ps_sfig18.instance_patterns,
        ),
    notebook:
        "../notebooks/combine_pdf_plots.py.ipynb"

rule plot_validation_figures:
    input:
        rules.plot_fig2.output,
        rules.plot_sfig18.output,