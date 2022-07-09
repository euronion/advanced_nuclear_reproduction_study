# SPDX-FileCopyrightText: 2022 Johannes Hampp <johannes.hampp@zeu.uni-giessen.de>
#
# SPDX-License-Identifier: MIT

# -*- coding: utf-8 -*-
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pypsa

n = pypsa.Network(snakemake.input.network)


for carrier in ["wind", "solar"]:
    df = pd.read_csv(snakemake.input[carrier], skiprows=3)
    index = pd.to_datetime(df[["year", "month", "day", "hour"]])
    # Change of convention: an hour is represented by its beginning
    index = index - pd.Timedelta("1h")
    ds = pd.Series(df[carrier[0] + "_cfs"].values, index)
    ds = ds[n.snapshots]

    n.add(
        "Generator",
        carrier,
        bus="electricity",
        p_nom_extendable=True,
        p_max_pu=ds,
    )


n.export_to_netcdf(snakemake.output.network)
