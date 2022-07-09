# SPDX-FileCopyrightText: 2022 Johannes Hampp <johannes.hampp@zeu.uni-giessen.de>
#
# SPDX-License-Identifier: MIT

# -*- coding: utf-8 -*-
import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pypsa

logger = logging.getLogger(__name__)

n = pypsa.Network(snakemake.input.network)

df = pd.read_csv(snakemake.input.demand, skiprows=5)
index = pd.to_datetime(df[["year", "month", "day", "hour"]])
# Change of convention: an hour is represented by its beginning
index = index - pd.Timedelta("1h")
demand = pd.Series(df.demand.values, index)

years = snakemake.wildcards["years"].split("-")
if int(years[-1]) == 2019 and snakemake.config["emulate_original_paper"]:
    # This is only used if the year is 2019
    # Take the 2018 demand data, create a copy with a shifted time index
    # and append it to the demand time-series.
    logger.warning(
        "'emulate_original_paper' enabled. "
        "Using demand data labelled '2018' also for year 2019."
    )

    demand_2018 = demand[-8760:]
    demand_2018.index = demand_2018.index + pd.Timedelta("8760h")
    demand = pd.concat([demand, demand_2018])


demand = demand[n.snapshots]

n.add("Load", "load", bus="electricity", p_set=demand)

n.export_to_netcdf(snakemake.output.network)
